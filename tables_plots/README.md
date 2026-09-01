Alphabetically:

bib_table.tex outputs a LaTeX table via xtable package presenting the cumulative lag coefficients for each bin of the binary and days models

binary_plot.pdf outputs a plot showing the three cumulative lag coefficients of the binary model

binary_vcovs.pdf outputs a plot showing the three cumulative lag coefficients of the binary model with four different variance-covariance specifications: County-Year two-way clustered, and Conley spatial errors (100km, 200km, and 400km cutoffs)

binary_vcovs_dk.pdf outputs a plot adds Driscoll-Kraay variance-covariance specification to the existing four specifications from the binary_vcovs.pdf plot

cont_0001_asinh_plot.pdf outputs a plot of the 24 exposure coefficients part of the event study (lead and lag endpoints, 11 lags, 10 leads, and 1 simultaneous)

cumcoeff_summary.tex outputs a LaTeX table via xtable package presenting the cumulative lag coefficients for the three primary models (continuous, binary, and days specifications)

endpoint_plot.pdf outputs a plot showing the individual coefficients of the continuous model with and without binned endpoints

parameter_plot.pdf outouts a plot showing the continuous model with five different exponential decay parameters chosen

primary_results.tex outputs a LaTeX table via fixest::etable package presenting the results for each of the three models; since all of the coefficients are presented, the table contains over 360 rows so I construct alternatives

summary_stats.tex outputs a LaTeX table via xtable package presenting summary statistics for each of the variables used in the regressions
