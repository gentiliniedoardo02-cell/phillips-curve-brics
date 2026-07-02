library(tidyverse)
library(lmtest)
library(sandwich)
library(zoo)
library(patchwork)
library(tseries)
library(strucchange)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
brics_raw <- read_csv("brics_panel.csv")

cat("Rows:", nrow(brics_raw), "\n")
cat("Countries:", unique(brics_raw$country), "\n")
cat("Period:", format(min(brics_raw$date), "%Y-%m"), "—",
    format(max(brics_raw$date), "%Y-%m"), "\n")

brics <- brics_raw %>%
  group_by(country) %>%
  arrange(date) %>%
  mutate(
    inflation_lag  = lag(inflation, 1),
    d_log_exch = c(NA, diff(log(exchange_rate))) * 100,
    exports_growth = c(NA, diff(log(exports))) * 100,
    imports_growth = c(NA, diff(log(imports))) * 100,
    period = ifelse(date < as.Date("2018-07-01"), "Pre-2018", "Post-2018"),
    post_trade = ifelse(date >= as.Date("2018-07-01"), 1, 0),
    unemployment_c = unemployment - mean(unemployment, na.rm = TRUE),
    lambda= unemployment_c * post_trade
  ) %>%
  filter(!is.na(inflation_lag), !is.na(d_log_exch),
         !is.na(exports_growth), !is.na(imports_growth)) %>%
  ungroup()

COUNTRIES <- c("Brazil", "Russia", "India", "China", "South_Africa")
BREAK_DATE <- as.Date("2018-07-01")

country_colors <- c(
  "Brazil" = "blue",
  "Russia" = "orange",
  "India" = "red",
  "China" = "lightblue",
  "South_Africa" = "grey50"
)

library(urca)

# unit root tests 
# ADF + pp + KPSS
for (p in COUNTRIES) {
  cat("\n", p, "\n")
  df_p <- brics %>% filter(country == p)
  #stationarity on inflation
  x_inf <- df_p$inflation
  adf_inf <- ur.df(x_inf, type = "drift", 
                   lags = 4, selectlags = "AIC")
  pp_inf <- pp.test(x_inf, alternative = "stationary")
  kpss_inf <- kpss.test(x_inf, null = "Level")
  
  cat("INFLATION\n ADF stat:", 
      round(adf_inf@teststat[1], 3),
      "- PP p:", round(pp_inf$p.value, 4),
      "- KPSS p:", round(kpss_inf$p.value, 4), "\n")
  
  #stationarity on unemployment
  x_un <- df_p$unemployment
  adf_un <- ur.df(x_un, type = "drift",
                  lags = 4, selectlags = "AIC")
  pp_un <- pp.test(x_un, alternative = "stationary")
  kpss_un <- kpss.test(x_un, null = "Level")
  
  cat("UNEMPLOYMENT\n ADF stat:",
      round(adf_un@teststat[1], 3),
      "- PP p:", round(pp_un$p.value, 4),
      "- KPSS p:", round(kpss_un$p.value, 4), "\n")
}

# Figure 1: inflation
ggplot(brics, aes(x = date, y = inflation, color = country)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = BREAK_DATE, linetype = "dashed",
             color = "gray", linewidth = 0.5) +
  annotate("text", x = as.Date("2018-10-01"), y = Inf,
           label = "Trade war", vjust = 2, hjust = 0, size = 3, color = "gray") +
  scale_color_manual(values = country_colors) +
  labs(title = "BRICS inflation 2010–2023",
       x = NULL, y = "Inflation (%)", color = NULL,
       caption = "Figure 1: Monthly inflation, dashed = July 2018") +
  theme_minimal()+
  theme(legend.position = "bottom")

# Figure 2: unemployment
ggplot(brics, aes(x = date, y = unemployment, color = country)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = BREAK_DATE, linetype = "dashed",
             color = "gray", linewidth = 0.5) +
  scale_color_manual(values = country_colors) +
  labs(title = "BRICS unemployment 2010–2023",
       x = NULL, y = "Unemployment (%)", color = NULL,
       caption = "Figure 2: Monthly unemployment rate") +
  theme_minimal()+
  theme(legend.position = "bottom")

# Figure 3: exchange rate
ggplot(brics, aes(x = date, y = d_log_exch, color = country)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = BREAK_DATE, linetype = "dashed",
             color = "gray", linewidth = 0.6) +
  scale_color_manual(values = country_colors) +
  labs(title = "BRICS exchange rate variation (%)",
       x = NULL, y = "Δlog exchange rate (%)", color = NULL,
       caption = "Figure 3: Proxy for global trade shock") +
  theme_minimal()+
  theme(legend.position = "bottom")

# Figure 4: Phillips Curve scatter
ggplot(brics, aes(x = unemployment, y = inflation, color = country)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = country_colors) +
  facet_wrap(~ country, scales = "free", ncol = 3) +
  labs(title = "Phillips Curve BRICS 2010–2023",
       x = "Unemployment", y = "Inflation",
       caption = "Figure 4: OLS regression line, full sample") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# M1 and M2 estimation
models_m1 <- list()
models_m2 <- list()
aic_tab <- data.frame(country = character(),
                        r2_m1 = numeric(), r2_m2 = numeric(),
                        aic_m1 = numeric(), aic_m2 = numeric())

for (p in COUNTRIES) {
  df_p <- brics %>% filter(country == p)
  m1 <- lm(inflation ~ unemployment, data = df_p)
  m2 <- lm(inflation ~ inflation_lag + unemployment + d_log_exch, data = df_p)
  
  models_m1[[p]] <- m1
  models_m2[[p]] <- m2
  
  cat("\n", p, "\n")
  cat("M1 R2:", round(summary(m1)$r.squared, 3),
      " M2 R2:", round(summary(m2)$r.squared, 3), "\n")
  cat("AIC M1:", round(AIC(m1), 1),
      " M2:", round(AIC(m2), 1), "\n")
  
  aic_tab <- rbind(aic_tab, data.frame(
    country = p,
    r2_m1  = round(summary(m1)$r.squared, 3),
    r2_m2 = round(summary(m2)$r.squared, 3),
    aic_m1  = round(AIC(m1), 1),
    aic_m2  = round(AIC(m2), 1)
  ))
}

# Figure 5 - 6: AIC  and R2 comparison 
aic_long <- aic_tab %>%
  dplyr::select(country, aic_m1, aic_m2) %>%
  pivot_longer(cols = c(aic_m1, aic_m2),
               names_to = "modello",
               values_to = "aic") %>%
  mutate(modello = ifelse(modello == "aic_m1", "M1 baseline", "M2 open economy"))

ggplot(aic_long, aes(x = country, y = aic, fill = modello)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("M1 baseline" = "gray", "M2 open economy" = "lightblue")) +
  labs(title = "AIC: M1 vs M2",
       subtitle = "Lower = better",
       x = NULL, y = "AIC", fill = "Modello",
       caption = "Figure 5") +
  theme(legend.position = "bottom")

r2_long <- aic_tab %>%
  dplyr::select(country, r2_m1, r2_m2) %>%
  pivot_longer(-country, names_to = "model", values_to = "r2") %>%
  mutate(model = ifelse(model == "r2_m1", "M1", "M2"))

ggplot(r2_long, aes(x = country, y = r2, fill = model)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("M1" = "gray", "M2" = "lightblue")) +
  labs(title = "R2: M1 vs M2",
       x = NULL, y = "R2", fill = NULL,
       caption = "Figure 6") +
  theme(legend.position = "bottom")

# Misspecification tests on M2
for (p in COUNTRIES) {
  cat("\n", p, "\n")
  m2 <- models_m2[[p]]
  
  bp <- bptest(m2)
  white <- bptest(m2, ~ fitted(m2) + I(fitted(m2)^2))
  cat("Breusch-Pagan p =", round(bp$p.value, 4),
      " White p =", round(white$p.value, 4), "\n")
  dw  <- dwtest(m2)
  bg4 <- bgtest(m2, order = 4)
  cat("Durbin-Watson p =", round(dw$p.value, 4),
      " BG(4) p =", round(bg4$p.value, 4), "\n")
  reset <- resettest(m2, power = 2:3, type = "fitted")
  cat("RESET p =", round(reset$p.value, 4), "\n")
}

# Figure 7: M2 residuals
resid_df <- do.call(rbind, lapply(COUNTRIES, function(p) {
  data.frame(
    country     = p,
    observation = seq_along(residuals(models_m2[[p]])),
    residual    = residuals(models_m2[[p]])
  )
}))

ggplot(resid_df, aes(x = observation, y = residual)) +
  geom_line(color = "lightblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  facet_wrap(~ country, scales = "free_y", ncol = 3) +
  labs(title = "M2 residuals",
       x = "Observation", y = "Residual",
       caption = "Figure 7") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

# Figure 8: ACF residuals
acf_plots <- lapply(COUNTRIES, function(p) {
  resid <- residuals(models_m2[[p]])
  acf_out <- acf(resid, lag.max = 24, plot = FALSE)
  acf_df <- data.frame(lag = as.numeric(acf_out$lag[-1]),
                        acf = as.numeric(acf_out$acf[-1]))
  ci <- qnorm(0.975) / sqrt(length(resid))
  ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_col(fill = "lightblue", width = 0.5) +
    geom_hline(yintercept = c(ci, -ci), linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, color = "black") +
    labs(title = p, x = "Lag", y = "ACF") +
    theme_minimal()
})
wrap_plots(acf_plots, ncol = 3) +
  plot_annotation(caption = "Figure 8: ACF of M2 residuals")

# Figure 8b: PACF residuals M2
pacf_plots <- lapply(COUNTRIES, function(p) {
  resid <- residuals(models_m2[[p]])
  pacf_out <- pacf(resid, lag.max = 24, plot = FALSE)
  pacf_df  <- data.frame(lag = as.numeric(pacf_out$lag),
                         pacf = as.numeric(pacf_out$acf))
  ci <- qnorm(0.975) / sqrt(length(resid))
  ggplot(pacf_df, aes(x = lag, y = pacf)) +
    geom_col(fill = "lightblue", width = 0.5) +
    geom_hline(yintercept = c(ci, -ci), 
               linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, color = "black") +
    labs(title = p, x = "Lag", y = "PACF") +
    theme_minimal()
})
wrap_plots(pacf_plots, ncol = 3) +
  plot_annotation(caption = "Figure 8b: PACF of M2 residuals")
# HAC correction — Newey-West (lag = 4)
for (p in COUNTRIES) {
  cat("\n", p, "\n")
  print(coeftest(models_m2[[p]], vcov = NeweyWest(models_m2[[p]], lag = 4)))
}
# Structural break tests: Chow - QA 
# Chow test
cat("\nChow test (break = July 2018)\n")
for (p in COUNTRIES) {
  df_p <- brics %>% filter(country == p)
  bp   <- which(df_p$date == BREAK_DATE)
  chow <- sctest(inflation ~ inflation_lag + unemployment + d_log_exch,
                 data = df_p, type = "Chow", point = bp)
  cat(p, "F =", round(chow$statistic, 3), " p =", round(chow$p.value, 4), "\n")
}
# Quandt-Andrews test
cat("\nQuandt-Andrews test\n")
for (p in COUNTRIES) {
  df_p <- brics %>% filter(country == p)
  qa   <- sctest(inflation ~ inflation_lag + unemployment + d_log_exch,
                 data = df_p, type = "supF")
  cat(p, "supF =", round(qa$statistic, 3), " p =", round(qa$p.value, 4), "\n")
}

# M3 — interaction model with structural break
# lambda = change in Phillips Curve slope post-2018 (u+lambda)
models_m3 <- list()
for (p in COUNTRIES) {
  df_p <- brics %>% filter(country == p)
  m3 <- lm(inflation ~ inflation_lag + unemployment_c + d_log_exch +
             post_trade + lambda, data = df_p)
  models_m3[[p]] <- m3
  
  cat("\n", p, "\n")
  print(summary(m3))
  cat("HAC (lag = 4):\n")
  print(coeftest(m3, vcov = NeweyWest(m3, lag = 4)))
}
library(car)
vif(models_m3[["Brazil"]])

# Table 2: lambda summary
lambda_tab <- do.call(rbind, lapply(COUNTRIES, function(p) {
  m3  <- models_m3[[p]]
  hac <- coeftest(m3, vcov = NeweyWest(m3, lag = 4))
  data.frame(
    country    = p,
    kappa      = round(coef(m3)["unemployment_c"], 3),      #controlla errore 
    lambda     = round(hac["lambda", "Estimate"], 3),
    slope_post = round(coef(m3)["unemployment_c"] + coef(m3)["lambda"], 3),  
    p_value    = round(hac["lambda", "Pr(>|t|)"], 4)
  )
}))

cat("\nTable 2: change in Phillips slope post-2018 (M3, HAC)\n")
print(lambda_tab)

# Figure 9: lambda with HAC confidence intervals
lambda_se <- data.frame(
  country = COUNTRIES,
  lambda  = lambda_tab$lambda,
  se      = sapply(COUNTRIES, function(p) {
    coeftest(models_m3[[p]],
             vcov = NeweyWest(models_m3[[p]], lag = 4))["lambda", "Std. Error"]
  })
)

ggplot(lambda_se, aes(x = reorder(country, lambda), y = lambda)) +
  geom_col(fill = "lightblue", width = 0.6) +
  geom_errorbar(aes(ymin = lambda - 1.96 * se,
                    ymax = lambda + 1.96 * se),
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  coord_flip() +
  labs(title = "Change in Phillips Curve slope post-2018",
       subtitle = "Lambda from M3 — HAC 95% CI",
       x = NULL, y = "Lambda",
       caption = "Figure 9: CI not crossing zero = significant change") 

# Subsample M2 pre/post 
kappa_pre  <- c()
kappa_post <- c()

for (p in COUNTRIES) {
  df_pre  <- brics %>% filter(country == p, date < BREAK_DATE)
  df_post <- brics %>% filter(country == p, date >= BREAK_DATE)
  kappa_pre[p]  <- round(coef(lm(inflation ~ inflation_lag + unemployment +
                                   d_log_exch, data = df_pre))["unemployment"], 3)
  kappa_post[p] <- round(coef(lm(inflation ~ inflation_lag + unemployment +
                                   d_log_exch, data = df_post))["unemployment"], 3)
}
# Figure 10: Phillips Curve pre vs post
ggplot(brics, aes(x = unemployment, y = inflation, color = period)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Pre-2018" = "blue", "Post-2018" = "red")) +
  facet_wrap(~ country, scales = "free", ncol = 3) +
  labs(title = "Phillips Curve pre and post 2018 trade war",
       subtitle = "Subsamples — illustrative only",
       x = "Unemployment", y = "Inflation", color = NULL,
       caption = "Figure 10") +
  theme_minimal()+
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

sessionInfo()