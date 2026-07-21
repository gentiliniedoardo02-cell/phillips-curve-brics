# Did the Phillips Curve change shape in BRICS after the 2018 trade war?

Group project, Macroeconometrics I - ISEG Lisbon

## The question

We wanted to check whether a specific shock, the 2018 US-China trade war, actually moved the needle on the Phillips Curve relationship in the five BRICS economies, or whether it's just noise like everything else.

Monthly data, January 2010 to December 2023, five countries, 840 observations total.

## How we approached it

Built it up in three steps rather than jumping straight to the complicated model:

- **M1** - the textbook Phillips curve, inflation on unemployment, nothing else. Mostly here as a baseline to show how little it explains on its own (R² under 0.06 everywhere, basically zero for China).
- **M2** - added lagged inflation (as a stand-in for adaptive expectations, following Friedman-Phelps) and the change in the exchange rate, since currency depreciation is basically the same transmission mechanism as a tariff imported goods get more expensive either way.
- **M3** - the actual test: an interaction term between unemployment and a post-July-2018 dummy, to see if the slope itself changed after the trade war started.

Every model gets run per country, then corrected with Newey-West HAC standard errors once we'd checked for heteroskedasticity and autocorrelation. We also ran Chow and Quandt-Andrews tests to check for a clean structural break at the trade war date, and unit root tests up front to make sure everything was stationary.

One methodological choice worth flagging: we used the raw unemployment rate rather than an unemployment gap (which would need a NAIRU estimate) for BRICS labour markets.
## What we found

No clean structural break anywhere: the Chow test comes back insignificant for all five countries, so the trade war didn't cause an abrupt shift on a single date.

But the interaction model (M3) tells a more interesting story once you let the slope change gradually rather than looking for a hard break: **India and China both show a significant slope change at the 10% level**, in opposite directions.

- **India**: slope steepens (λ = −0.437, p = 0.061) — counterintuitive at first, but plausible if trade diversion from China pushed extra demand into Indian manufacturing
- **China**: slope flips sign, from negative to positive (λ = +0.348, p = 0.098), consistent with a stagflationary shock: exports got hit by tariffs while supply-chain disruption pushed costs up at the same time
- **Brazil, Russia, South Africa**: no significant change either way

So the effect wasn't uniform across the bloc, but it shows up specifically in the two economies most directly exposed to the US-China trade relationship, and it works through different channels in each.

## Files

- `Phillip_curve_brics.R` - full pipeline: data construction, unit root tests, M1/M2/M3 estimation, misspecification tests, HAC correction, Chow/Quandt-Andrews, all plots
- `brics_panel.csv` — monthly panel, 2010–2023
- `Phillips_Curve_BRICS_paper.pdf` - full writeup with literature review and country-by-country discussion


Edoardo Gentilini — [LinkedIn](https://linkedin.com/in/edoardo-gentilini) — gentilini.edoardo02@gmail.com
