#for data cleaning and manipulation
library(tidyverse)
#for panel fixed effects models
library(fixest)
#for hypothesis testing
library(marginaleffects)
#for sourcing data
library(here)
#load panel data
panel_all <- readRDS(here("source_data", "panel_all.rds"))

#calculate binned endpoints following Schmidheiny & Seigloch (2023)
#the lag endpoint calculates the cumulative sum of the variable starting at lag 12 and moving forward
lag_endpoint_bin <- function(x) {lag(cumsum(as.numeric(x)), 12)}
#the lead endpoint calculates the cumulative sum of the variable starting at lead 12 and moving backward
lead_endpoint_bin <- function(x) {lead(rev(cumsum(rev(as.numeric(x)))), 12)}
panel_all <- panel_all %>% 
  #run functions for each variable for each city
  group_by(CITY) %>%
  mutate(
    across(starts_with("asinh_"), lag_endpoint_bin, .names = "lag_{.col}_end"),
    across(starts_with("asinh_"), lead_endpoint_bin, .names = "lead_{.col}_end"),
    across(starts_with("cont_"), lag_endpoint_bin, .names = "lag_{.col}_end"),
    across(starts_with("cont_"), lead_endpoint_bin, .names = "lead_{.col}_end"),
    across(starts_with("binary_"), lag_endpoint_bin, .names = "lag_{.col}_end"),
    across(starts_with("binary_"), lead_endpoint_bin, .names = "lead_{.col}_end"),
    across(starts_with("days_"), lag_endpoint_bin, .names = "lag_{.col}_end"),
    across(starts_with("days_"), lead_endpoint_bin, .names = "lead_{.col}_end"),
  ) %>%
  ungroup()

#define functions for varying robust inference specifications
countyyear_cluster_vcov <- function(model) {
  vcov(model, vcov = ~COUNTY + year)
}
#the Conley spatial errors follow Conley (1999) and incorporate a distance cutoff
#I choose cutoffs of 100km, 200km, and 400km in this analysis
conley100_vcov <- function(model) {
  vcov(model, vcov = vcov_conley(model, cutoff = 100))
}
conley200_vcov <- function(model) {
  vcov(model, vcov = vcov_conley(model, cutoff = 200))
}
conley400_vcov <- function(model) {
  vcov(model, vcov = vcov_conley(model, cutoff = 400))
}

#define functions to create strings for hypothesis testing (joint nullity for cumulative terms)
#where horizon = 'l' or 'f', trans = 'asinh' or the empty string, and parameter is the appropriate theta
cont_terms <- function(horizon, trans, parameter) {
  horizon_range <- if (horizon == "l") 1:11 else 2:11
  paste0("`", horizon, "(", trans, "_cont_", parameter, ", ", horizon_range, ")`")
}
cont_hyp_string <- function(horizon, trans, parameter) {
  paste(paste(cont_terms(horizon, trans, parameter), collapse = " + "), "= 0")
}

#define function for cumulative coefficients of continuous models
#where 'horizons' can be a vector (c("l", "f")) or individual terms
cont_coeffs <- function(model, trans, parameter, vcov, horizons) {
  for (horizon in horizons) {
    print(
      hypotheses(
        model, 
        cont_hyp_string(horizon, trans, parameter), 
        vcov))
    }}

#define function to collect outputs from continuous exposure models 
cont_output <- function(model, trans, parameter, vcov, label) {
  #extract coefficients of the model to a data frame using a function to simplify plotting
  coeffs_df <- as.data.frame(summary(model, vcov = vcov)$coeftable)
  coeffs_df$name <- rownames(coeffs_df)
  coeffs_df %>%
    #filter only coefficients of the continuous exposure (treatment) variable, including lags, leads, and binned endpoints
    filter(grepl(pattern = "cont", name)) %>%
    #label each coefficient with a number (-12 to 12) in order to graph on a numeric x-axis
    mutate(horizon = 
             case_when(
               name == paste0(trans, "_cont_", parameter) ~ 0,
               #lag terms to 1 to 11 (right side of graph)
               grepl("l\\(", name) ~ as.numeric(sub(".*, (\\d+)\\)", "\\1", name)),
               #lead terms to -11 to -2 (left side of graph)
               grepl("f\\(", name) ~ -as.numeric(sub(".*, (\\d+)\\)", "\\1", name)),
               #binned endpoints to -12 and 12
               name == paste0("lead_", trans, "_cont_", parameter, "_end") ~ -12,
               name == paste0("lag_", trans, "_cont_", parameter, "_end") ~ 12),
           Model = label) %>%
    rename(SE = 'Std. Error')
}

#run regression for continuous exposure model, excluding first lead and including endpoint bins
cont_0001_asinh_model <- feols(price_log ~ 
                                 l(asinh_cont_0001, 1:11) + asinh_cont_0001 + f(asinh_cont_0001, 2:11) + lag_asinh_cont_0001_end + lead_asinh_cont_0001_end + temp + ppt | CITY + month_year,
                               data = panel_all, panel.id = ~CITY + month_year)
#print cumulative lag and cumulative lead coefficient estimates for the continuous exposure model
cont_coeffs(cont_0001_asinh_model, "asinh", "0001", countyyear_cluster_vcov(cont_0001_asinh_model), c("l", "f"))
#create data frame presenting individual coefficients
cont_0001_asinh_output <- cont_output(cont_0001_asinh_model, "asinh", "0001", countyyear_cluster_vcov(cont_0001_asinh_model), label = "Model with Parameter = 0.0001")

#check robustness with other specifications:
#run model with no endpoint binning and 12 lags/leads
cont_0001_asinh_model2 <- feols(price_log ~
                                  l(asinh_cont_0001, 1:12) + asinh_cont_0001 + f(asinh_cont_0001, 1:12) +
                                  temp + ppt | CITY + month_year, data = panel_all, panel.id = ~CITY + month_year)
#test joint nullity of cumulative lag coefficients
hypotheses(cont_0001_asinh_model2,
           paste(paste(paste0("`l(asinh_cont_0001, ", 1:12, ")`"), collapse = " + "), "= 0"), 
           vcov = countyyear_cluster_vcov(cont_0001_asinh_model2))
#test joint nullity of cumulative lead coefficients
hypotheses(cont_0001_asinh_model2,
           paste(paste(paste0("`f(asinh_cont_0001, ", 1:12, ")`"), collapse = " + "), "= 0"), 
           vcov = countyyear_cluster_vcov(cont_0001_asinh_model2))
endpoint_outputs <- bind_rows(cont_output(cont_0001_asinh_model, "asinh", "0001", countyyear_cluster_vcov(cont_0001_asinh_model), "with Endpoint Binning"),
                              cont_output(cont_0001_asinh_model2, "asinh", "0001", countyyear_cluster_vcov(cont_0001_asinh_model2), "w/o Endpoint Binning"))

#run regressions with various parameters and plot results for comparison
cont_0004_asinh_model <- feols(price_log ~ l(asinh_cont_0004, 1:11) + asinh_cont_0004 + f(asinh_cont_0004, 2:11) + lag_asinh_cont_0004_end + lead_asinh_cont_0004_end +
                                 temp + ppt | CITY + month_year, data = panel_all, panel.id = ~CITY + month_year)
cont_0002_asinh_model <- feols(price_log ~ l(asinh_cont_0002, 1:11) + asinh_cont_0002 + f(asinh_cont_0002, 2:11) + lag_asinh_cont_0002_end + lead_asinh_cont_0002_end +
                                 temp + ppt | CITY + month_year, data = panel_all, panel.id = ~CITY + month_year)
cont_00005_asinh_model <- feols(price_log ~ l(asinh_cont_00005, 1:11) + asinh_cont_00005 + f(asinh_cont_00005, 2:11) + lag_asinh_cont_00005_end + lead_asinh_cont_00005_end +
                                  temp + ppt | CITY + month_year, data = panel_all, panel.id = ~CITY + month_year)
cont_000025_asinh_model <- feols(price_log ~ l(asinh_cont_000025, 1:11) + asinh_cont_000025 + f(asinh_cont_000025, 2:11) + lag_asinh_cont_000025_end + lead_asinh_cont_000025_end +
                                   temp + ppt | CITY + month_year, data = panel_all, panel.id = ~CITY + month_year)
#find the cumulative lag coefficients for each of the models
cont_coeffs(cont_0004_asinh_model, "asinh", "0004", countyyear_cluster_vcov(cont_0004_asinh_model), "l")
cont_coeffs(cont_0002_asinh_model, "asinh", "0002", countyyear_cluster_vcov(cont_0002_asinh_model), "l")
cont_coeffs(cont_00005_asinh_model, "asinh", "00005", countyyear_cluster_vcov(cont_00005_asinh_model), "l")
cont_coeffs(cont_000025_asinh_model, "asinh", "000025", countyyear_cluster_vcov(cont_000025_asinh_model), "l")

parameter_outputs <- bind_rows(
  cont_output(cont_0004_asinh_model, "asinh", "0004", countyyear_cluster_vcov(cont_0004_asinh_model), "Theta = 0.0004"),
  cont_output(cont_0002_asinh_model, "asinh", "0002", countyyear_cluster_vcov(cont_0002_asinh_model), "Theta = 0.0002"),
  cont_output(cont_0001_asinh_model, "asinh", "0001", countyyear_cluster_vcov(cont_0001_asinh_model), "Theta = 0.0001"),
  cont_output(cont_00005_asinh_model, "asinh", "00005", countyyear_cluster_vcov(cont_00005_asinh_model), "Theta = 0.00005"),
  cont_output(cont_000025_asinh_model, "asinh", "000025", countyyear_cluster_vcov(cont_000025_asinh_model), "Theta = 0.000025")
)

#define distances vector for use in loops
distances <- c("0_10", "10_50", "50_100")
#where horizon = 'l' or 'f', bin = 'binary' or 'days', and distance is '0_10', '10_50', or '50_100'
bin_terms <- function(horizon, bin, distance) {
  horizon_range <- if (horizon == "l") 1:11 else 2:11
  paste0("`", horizon, "(", bin, "_", distance, ", ", horizon_range, ")`")
}
bin_hyp_string <- function(horizon, bin, distance) {
  paste(paste(bin_terms(horizon, bin, distance), collapse = " + "), "= 0")
}

#prints hypotheses outputs (cumulative coefficient estimates)
#where 'distances' can be a vector (c("0_10", "10_50", "50_100")) or each term individually and horizon can be 'l', 'f', or both
bin_coeffs <- function(model, bin, vcov, horizons, distances)
  for (distance in distances) {
    for (horizon in horizons) {
      print(
        hypotheses(
          model,
          bin_hyp_string(horizon, bin, distance),
          vcov))
    }}

#put cumulative coefficient estimates in a data frame
bin_output <- function(model, bin, vcov, horizons, distances, label) {
  coeffs <- list()
  #adds each of the three distances to the data frame
  for (horizon in horizons) {
    for (distance in distances) {
      coeffs[[paste0(horizon, "_", distance)]] <- hypotheses(model, bin_hyp_string(horizon, bin, distance), vcov) %>%
        mutate(horizon = if (horizon == "l") "Lag" else "Lead",
               distance_bin = distance)
    }
  }
  bind_rows(coeffs) %>%
    #use the middle of the bin for easier visualization when plotted
    mutate(distance_mean = case_when(
      distance_bin == "0_10" ~ 5,
      distance_bin == "10_50" ~ 30,
      distance_bin == "50_100" ~ 75),
      Model = label)
}

#run regression for binary bin exposure model
binary_model <- feols(price_log ~
                        l(binary_0_10, 1:11) + binary_0_10 + f(binary_0_10, 2:11) + lag_binary_0_10_end + lead_binary_0_10_end +
                        l(binary_10_50, 1:11) + binary_10_50 + f(binary_10_50, 2:11) + lag_binary_10_50_end + lead_binary_10_50_end +
                        l(binary_50_100, 1:11) + binary_50_100 + f(binary_50_100, 2:11) + lag_binary_50_100_end + lead_binary_50_100_end +
                        temp + ppt | CITY + month_year, data = panel_all, panel.id = ~CITY + month_year)
#prints the cumulative lag coefficient of the binary model for each distance bin
bin_coeffs(binary_model, "binary", countyyear_cluster_vcov(binary_model), "l", c("0_10", "10_50", "50_100"))
#create data frame presenting cumulative coefficients
binary_output <- bin_output(binary_model, "binary", countyyear_cluster_vcov(binary_model), c("l","f"), c("0_10", "10_50", "50_100"), "Binary Bins Model")

#as a robustness check, I run regressions and plot results of the binary model with varying variance-covariance matrices
vcov_outputs <- bind_rows(
  bin_output(binary_model, "binary", countyyear_cluster_vcov(binary_model), "l", c("0_10", "10_50", "50_100"), "Two-Way: County-Year"),
  bin_output(binary_model, "binary", conley100_vcov(binary_model), "l", c("0_10", "10_50", "50_100"), "Spatial: Conley (100km Cutoff)"),
  bin_output(binary_model, "binary", conley200_vcov(binary_model), "l", c("0_10", "10_50", "50_100"), "Spatial: Conley (200km Cutoff)"),
  bin_output(binary_model, "binary", conley400_vcov(binary_model), "l", c("0_10", "10_50", "50_100"), "Spatial: Conley (400km Cutoff)")
)
vcov_outputs_dk <- vcov_outputs %>% bind_rows(
  bin_output(binary_model, "binary", vcov_DK(binary_model, lag = 12), "l", c("0_10", "10_50", "50_100"), "Mix: Driscoll-Kraay")
)

#run regression for days bin exposure model
days_model <- feols(price_log ~
                      l(days_0_10, 1:11) + days_0_10 + f(days_0_10, 2:11) + lag_days_0_10_end + lead_days_0_10_end +
                      l(days_10_50, 1:11) + days_10_50 + f(days_10_50, 2:11) + lag_days_10_50_end + lead_days_10_50_end +
                      l(days_50_100, 1:11) + days_50_100 + f(days_50_100, 2:11) + lag_days_50_100_end + lead_days_50_100_end +
                      temp + ppt | CITY + month_year, data = panel_all, panel.id = ~CITY + month_year)
#prints the cumulative lag coefficient of the days model for each distance bin
bin_coeffs(days_model, "days", countyyear_cluster_vcov(days_model), "l", c("0_10", "10_50", "50_100"))

#plots the cumulative lag coefficient estimates for days bin model
days_output <- bin_output(days_model, "days", countyyear_cluster_vcov(days_model), c("l","f"), c("0_10", "10_50", "50_100"), "Days Bins Model")

saveRDS(cont_0001_asinh_output, here("source_data", "cont_0001_asinh_output.rds"))
saveRDS(cont_0001_asinh_output2, here("source_data", "cont_0001_asinh_output2.rds"))
saveRDS(endpoint_outputs, here("source_data", "endpoint_outputs.rds"))
saveRDS(parameter_outputs, here("source_data", "parameter_outputs.rds"))
saveRDS(binary_output, here("source_data", "binary_output.rds"))
saveRDS(vcov_outputs, here("source_data", "vcov_outputs.rds"))
saveRDS(vcov_outputs_dk, here("source_data", "vcov_outputs_dk.rds"))
saveRDS(days_output, here("source_data", "days_output.rds"))
