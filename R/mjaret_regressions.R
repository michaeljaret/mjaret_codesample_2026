#for data cleaning and manipulation
library(tidyverse)
#for panel fixed effects models
library(fixest)
#for hypothesis testing
library(marginaleffects)
#for sourcing data
library(here)
#for mapping functions across alternative models
library(purrr)
#load panel data
panel_all <- readRDS(here("intermediate_data", "panel_all.rds"))

#calculate binned endpoints following Schmidheiny & Seigloch (2023)
#the lag endpoint calculates the cumulative sum of the variable starting at lag 12 and moving forward
lag_endpoint <- function(x, lag) {lag(cumsum(as.numeric(x)), lag)}
#the lead endpoint calculates the cumulative sum of the variable starting at lead 12 and moving backward
lead_endpoint <- function(x, lead) {lead(rev(cumsum(rev(as.numeric(x)))), lead)}
panel_all <- panel_all %>% 
  arrange(CITY, month_year) %>%
  #run functions for each variable for each city
  group_by(CITY) %>%
  mutate(
    across(starts_with("asinh_"), ~ lag_endpoint(.x, 12), .names = "lag_{.col}_end_12"),
    across(starts_with("asinh_"), ~ lead_endpoint(.x, 12), .names = "lead_{.col}_end_12"),
    across(starts_with("cont_"), ~ lag_endpoint(.x, 12), .names = "lag_{.col}_end_12"),
    across(starts_with("cont_"), ~ lead_endpoint(.x, 12), .names = "lead_{.col}_end_12"),
    lag_binary_end_4 = lag_endpoint(binary_0_10 + binary_10_50 + binary_50_100, 4),
    lead_binary_end_2 = lag_endpoint(binary_0_10 + binary_10_50 + binary_50_100, 2)) %>%
  ungroup()
  
panel_all <- panel_all %>% panel(panel.id = ~CITY + month_year) #set panel ID for fixest regressions
  
countyyear_cluster_vcov <- function(model) {
  vcov(model, vcov = ~COUNTY + year)}
conley_vcov <- function(model, cutoff) {
  vcov(model, vcov = vcov_conley(model, lat = "lat", lon = "lon", cutoff = cutoff))}

#define function to create event study model for continuous variable
cont_model <- function(dep, trans, parameter, lags, leads, endpoints, controls, fe) {
  trans[[1]] <- if_else(trans[1] == "asinh", "asinh_", trans[1])
  variable <- paste0(trans[1], "cont_", trans[2], "_", parameter) #creates variable name
  lag_terms <- paste0("l(", variable, ", ", lags, ")")
  lead_terms <- paste0("f(", variable, ", ", leads, ")") 
  #combines all non-control terms into an equation; includes endpoints if argument is true
  main_terms <- if_else(endpoints == TRUE, 
                        paste(c(lag_terms, variable, lead_terms, paste0("lag_", variable, "_end_", max(lags) + 1), paste0("lead_", variable, "_end_", max(leads) + 1)), collapse = " + "),
                        paste(c(lag_terms, variable, lead_terms), collapse = " + "))
  controls_terms <- paste(controls, collapse = " + ")
  fixed_effects <- paste(fe, collapse = " + ")
  equation <- paste0(dep, " ~ ", main_terms, " + ", controls_terms, " | ", fixed_effects) #creates fixest-friendly equation string
  return(as.formula(equation)) #turns equation string into formula
  }

#define functions to create strings for hypothesis testing (joint nullity for cumulative terms)
#where horizon = 'l' or 'f', trans = permutations of c('asinh'/'', 'log'/'raw'), and the decay parameter (i.e. '0001')
cont_model_terms <- function(horizon, trans, parameter, lags = NULL, leads = NULL) {
  trans[[1]] <- if_else(trans[1] == "asinh", "asinh_", trans[1])
  lead_terms <- paste0("`f(", trans[1], "cont_", trans[2], "_", parameter, ", ", leads, ")`")
  lag_terms <- paste0("`l(", trans[1], "cont_", trans[2], "_", parameter, ", ", lags, ")`")
  if (horizon == "l") {
    return(lag_terms)
  } else if (horizon == "f") {
    return(lead_terms)
  } else {
    return(NULL)}
  }

cont_hyp_string <- function(horizon, trans, parameter, lags = NULL, leads = NULL) {
  paste(paste(cont_model_terms(horizon, trans, parameter, lags = lags, leads = leads), collapse = " + "), "= 0")
  }

#define function to print cumulative coefficients of continuous models to console
#where 'horizons' can be a vector (c("l", "f")) or individual terms
cont_coeffs <- function(model, trans, parameter, vcov, horizons, lags = NULL, leads = NULL) {
  for (horizon in horizons) {
    print(
      hypotheses(
        model, 
        cont_hyp_string(horizon, trans, parameter, lags = lags, leads = leads), 
        vcov))}
  }

#define function to collect outputs from continuous exposure models 
cont_output <- function(model, trans, parameter, vcov, label) {
  #extract coefficients of the model to a data frame using a function to simplify plotting
  trans[[1]] <- if_else(trans[1] == "asinh", "asinh_", trans[1])
  variable <- paste0(trans[1], "_cont_", trans[2], "_", parameter)
  coeffs_df <- as.data.frame(summary(model, vcov = vcov)$coeftable)
  coeffs_df$name <- rownames(coeffs_df)
  coeffs_df %>%
    #filter only coefficients of the continuous exposure (treatment) variable, including lags, leads, and binned endpoints
    filter(grepl(pattern = "cont", name)) %>%
    #label each coefficient with a number (-12 to 12) in order to graph on a numeric x-axis
    mutate(horizon = 
             case_when(
               name == variable ~ 0,
               #lag terms to 1 to 11 (right side of graph)
               grepl("l\\(", name) ~ as.numeric(sub(".*, (\\d+)\\)", "\\1", name)),
               #lead terms to -11 to -2 (left side of graph)
               grepl("f\\(", name) ~ -as.numeric(sub(".*, (\\d+)\\)", "\\1", name)),
               #binned endpoints to -12 and 12
               name == paste0("lead_", variable, "_end") ~ -12,
               name == paste0("lag_", variable, "_end") ~ 12),
           Model = label) %>%
    rename(SE = 'Std. Error')
  }

#run regression for continuous exposure model, excluding first lead and including endpoint bins
main_cont_model <- feols(cont_model("price_log", c("asinh", "log"), "0001", 1:11, 2:11, TRUE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)
#print cumulative lag and cumulative lead coefficient estimates for the continuous exposure model
cont_coeffs(main_cont_model, c("asinh", "log"), "0001", countyyear_cluster_vcov(main_cont_model), c("l", "f"), lags = 1:11, leads = 2:11)

#check robustness with other specifications:
#run model with no endpoint binning
main_cont_model2 <- feols(cont_model("price_log", c("asinh", "log"), "0001", 1:12, 1:12, FALSE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)
cont_coeffs(main_cont_model2, c("asinh", "log"), "0001", countyyear_cluster_vcov(main_cont_model2), c("l", "f"), lags = 1:12, leads = 1:12)
endpoint_outputs <- bind_rows(cont_output(main_cont_model, c("asinh", "log"), "0001", countyyear_cluster_vcov(main_cont_model), "with Endpoint Binning"),
                              cont_output(main_cont_model2, c("asinh", "log"), "0001", countyyear_cluster_vcov(main_cont_model2), "without Endpoint Binning"))

#run regressions for each of the five parameters and place outputs into single data frame
parameter_models <- list()
parameters <- c("0004", "0002", "0001", "00005", "000025")
for (i in parameters) {
  parameter_models[[paste0("cont_", i, "_model")]] <- feols(cont_model("price_log", c("asinh", "log"), i, 1:11, 2:11, TRUE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)}
parameter_labels <- c("Parameter = 0.0004", "Parameter = 0.0002", "Parameter = 0.0001", "Parameter = 0.00005", "Parameter = 0.000025")
parameter_outputs <- pmap_dfr(
  list(model = parameter_models, parameter = parameters, labels = parameter_labels),
  ~ cont_output(model = ..1, trans = c("asinh", "log"), parameter = ..2, vcov = countyyear_cluster_vcov(..1), label = ..3))

#add endpoints for horizons of 9 and 15 to check robustness of alternative event study structure
panel_all <- panel_all %>%
  group_by(CITY) %>%
  mutate(lag_asinh_cont_log_0001_end_9 = lag_endpoint(asinh_cont_log_0001, 9),
         lead_asinh_cont_log_0001_end_9 = lead_endpoint(asinh_cont_log_0001, 9),
         lag_asinh_cont_log_0001_end_15 = lag_endpoint(asinh_cont_log_0001, 15),
         lead_asinh_cont_log_0001_end_15 = lead_endpoint(asinh_cont_log_0001, 15)) %>%
  ungroup() %>%
  panel(panel.id = ~CITY + month_year)
cont_model_9a <- feols(cont_model("price_log", c("asinh", "log"), "0001", 1:8, 2:8, TRUE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)
cont_coeffs(cont_model_9a, c("asinh", "log"), "0001", countyyear_cluster_vcov(cont_model_9a), "l", 1:8)
cont_model_15a <- feols(cont_model("price_log", c("asinh", "log"), "0001", 1:14, 2:14, TRUE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)
cont_coeffs(cont_model_15a, c("asinh", "log"), "0001", countyyear_cluster_vcov(cont_model_15a), "l", 1:14)
#repeat regressions without endpoint binning
cont_model_9b <- feols(cont_model("price_log", c("asinh", "log"), "0001", 1:9, 1:9, FALSE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)
cont_coeffs(cont_model_9b, c("asinh", "log"), "0001", countyyear_cluster_vcov(cont_model_9b), "l", 1:9)
cont_model_15b <- feols(cont_model("price_log", c("asinh", "log"), "0001", 1:15, 1:15, FALSE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)
cont_coeffs(cont_model_15b, c("asinh", "log"), "0001", countyyear_cluster_vcov(cont_model_15b), "l", 1:15)

#define dynamic function to run bin models
bin_model <- function(dep, bin, distances, lags, leads, endpoints, controls, fe) {
  variables <- paste0(bin, "_", distances)
  lag_terms <- outer(lags, variables, FUN = function(x, y) {paste0("l(", y, ", ", x, ")")}) %>% as.vector()
  lead_terms <- outer(leads, variables, FUN = function(x, y) {paste0("f(", y, ", ", x, ")")}) %>% as.vector()
  main_terms <- if_else(endpoints == TRUE, 
                        paste(c(lag_terms, variables, lead_terms, paste0("lag_binary_end_", max(lags) + 1), paste0("lead_binary_end_", max(leads) + 1)), collapse = " + "),
                        paste(c(lag_terms, variables, lead_terms), collapse = " + "))
  controls_terms <- paste(controls, collapse = " + ")
  fixed_effects <- paste(fe, collapse = " + ")
  equation <- paste0(dep, " ~ ", main_terms, " + ", controls_terms, " | ", fixed_effects)
  return(as.formula(equation))}

#define functions to set up coefficient hypothesis testing:
#creates list of terms being tested
bin_model_terms <- function(horizon, bin, distance, lags = NULL, leads = NULL) {
  lead_terms <- paste0("`f(", bin, "_", distance, ", ", leads, ")`")
  lag_terms <- paste0("`l(", bin, "_", distance, ", ", lags, ")`")
  if (horizon == "l") {
    return(lag_terms)
  } else if (horizon == "f") {
    return(lead_terms)
  } else {
    return(NULL)}}
#places list into an equation string
bin_hyp_string <- function(horizon, bin, distance, lags = NULL, leads = NULL) {
  paste(paste(bin_model_terms(horizon, bin, distance, lags = lags, leads = leads), collapse = " + "), "= 0")}

#runs joint nullity hypothesis tests
#where 'distances' can be a vector (c("0_10", "10_50", "50_100")) or each term individually and horizon can be 'l', 'f', or both
bin_coeffs <- function(model, bin, vcov, horizons, distances, lags = NULL, leads = NULL) {
  for (distance in distances) {
    for (horizon in horizons) {
      print(
        hypotheses(
          model,
          bin_hyp_string(horizon, bin, distance, lags = lags, leads = leads),
          vcov))}}}

#puts cumulative coefficients into data frame
bin_output <- function(model, bin, vcov, horizons, distances, lags = NULL, leads = NULL, label) {
  coeffs <- list()
  #adds each of the three distances to the data frame
  for (horizon in horizons) {
    for (distance in distances) {
      coeffs[[paste0(horizon, "_", distance)]] <- hypotheses(model, bin_hyp_string(horizon, bin, distance, lags = lags, leads = leads), vcov) %>%
        mutate(horizon = if (horizon == "l") "Lag" else "Lead",
               distance_bin = distance)}}
  bind_rows(coeffs) %>%
    #use the middle of the bin for plot
    mutate(
      distance_mean = case_when(
        distance_bin == "0_10" ~ 5,
        distance_bin == "10_50" ~ 30,
        distance_bin == "50_100" ~ 75),
      Bin = case_when(
        distance_bin == "0_10" ~ "Close",
        distance_bin == "10_50" ~ "Medium",
        distance_bin == "50_100" ~ "Far"),
      Model = label)}

#run main binary model
main_binary_model <- feols(bin_model("price_log", "binary", c("0_10", "10_50", "50_100"), 1:3, 1, TRUE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)

#prints the cumulative lag coefficient of the binary model for each distance bin
bin_coeffs(main_binary_model, "binary", countyyear_cluster_vcov(main_binary_model), c("l","f"), c("0_10", "10_50", "50_100"), lags = 1:3, leads = 1)
#create data frame presenting cumulative coefficients
binary_output <- bin_output(main_binary_model, "binary", countyyear_cluster_vcov(main_binary_model), c("l","f"), c("0_10", "10_50", "50_100"), "Binary Model", lags = 1:3, leads = 1)

#compares outputs with two-way cluster standard errors to outputs with conley spatial errors
vcov_binary_outputs <- bind_rows(
  bin_output(main_binary_model, "binary", countyyear_cluster_vcov(main_binary_model), "l", c("0_10", "10_50", "50_100"), lags = 1:3, leads = 1, "Two-Way Cluster: County-Year"),
  bin_output(main_binary_model, "binary", conley_vcov(main_binary_model, 100), "l", c("0_10", "10_50", "50_100"), lags = 1:3, leads = 1, "Spatial: Conley (100km Cutoff)"),
  bin_output(main_binary_model, "binary", conley_vcov(main_binary_model, 200), "l", c("0_10", "10_50", "50_100"), lags = 1:3, leads = 1, "Spatial: Conley (200km Cutoff)"),
  bin_output(main_binary_model, "binary", conley_vcov(main_binary_model, 400), "l", c("0_10", "10_50", "50_100"), lags = 1:3, leads = 1, "Spatial: Conley (400km Cutoff)"))
#adds Driscoll-Kraay standard errors to the above comparison
vcov_binary_outputs_dk <- vcov_binary_outputs %>% bind_rows(
  bin_output(main_binary_model, "binary", vcov_DK(main_binary_model, lag = 12), "l", c("0_10", "10_50", "50_100"), lags = 1:3, leads = 1, "Mix: Driscoll-Kraay"))

#run regression for days bin exposure model
main_days_model <- feols(bin_model("price_log", "days", c("0_10", "10_50", "50_100"), 1:3, 1, TRUE, c("temp", "ppt", "pdsi"), c("CITY", "month_year")), data = panel_all)
#prints the cumulative lag coefficient of the days model for each distance bin
bin_coeffs(main_days_model, "days", countyyear_cluster_vcov(main_days_model), "f", c("0_10", "10_50", "50_100"), leads = 1)

#plots the cumulative lag coefficient estimates for days bin model
days_output <- bin_output(main_days_model, "days", countyyear_cluster_vcov(main_days_model), c("l","f"), c("0_10", "10_50", "50_100"), lags = 1:3, leads = 1, "Days Bins Model")
  
#creates data frame presenting cumulative lag coefficient estimates for continuous and bin models
coeff_summary <- bind_rows(
  hypotheses(main_cont_model, cont_hyp_string("l", c("asinh", "log"), "0001", lags = 1:11, leads = 2:11), countyyear_cluster_vcov(main_cont_model)) %>%
    summarise('Cumulative Lag Estimate' = first(estimate), 'Std. Error' = first(std.error)) %>% mutate(Model = "Continuous Exposure"),
  hypotheses(main_binary_model, bin_hyp_string("l", "binary", "0_10", lags = 1:3, leads = 1), countyyear_cluster_vcov(main_binary_model)) %>%
    summarise('Cumulative Lag Estimate' = first(estimate), 'Std. Error' = first(std.error)) %>% mutate(Model = "Binary Exposure", Bin = "0 to 10km"),
  hypotheses(main_binary_model, bin_hyp_string("l", "binary", "10_50", lags = 1:3, leads = 1), countyyear_cluster_vcov(main_binary_model)) %>%
    summarise('Cumulative Lag Estimate' = first(estimate), 'Std. Error' = first(std.error)) %>% mutate(Model = "Binary Exposure", Bin = "10 to 50km"),
  hypotheses(main_binary_model, bin_hyp_string("l", "binary", "50_100", lags = 1:3, leads = 1), countyyear_cluster_vcov(main_binary_model)) %>%
    summarise('Cumulative Lag Estimate' = first(estimate), 'Std. Error' = first(std.error)) %>% mutate(Model = "Binary Exposure", Bin = "50 to 100km"),
  hypotheses(main_days_model, bin_hyp_string("l", "days", "0_10", lags = 1:3, leads = 1), countyyear_cluster_vcov(main_days_model)) %>%
    summarise('Cumulative Lag Estimate' = first(estimate), 'Std. Error' = first(std.error)) %>% mutate(Model = "Days Exposure", Bin = "0 to 10km"),
  hypotheses(main_days_model, bin_hyp_string("l", "days", "10_50", lags = 1:3, leads = 1), countyyear_cluster_vcov(main_days_model)) %>%
    summarise('Cumulative Lag Estimate' = first(estimate), 'Std. Error' = first(std.error)) %>% mutate(Model = "Days Exposure", Bin = "10 to 50km"),
  hypotheses(main_days_model, bin_hyp_string("l", "days", "50_100", lags = 1:3, leads = 1), countyyear_cluster_vcov(main_days_model)) %>%
    summarise('Cumulative Lag Estimate' = first(estimate), 'Std. Error' = first(std.error)) %>% mutate(Model = "Days Exposure", Bin = "50 to 100km")) %>% select(Model, Bin, 'Cumulative Lag Estimate', 'Std. Error')

saveRDS(main_cont_output, here("regression_outputs", "main_cont_output.rds"))
saveRDS(endpoint_outputs, here("regression_outputs", "endpoint_outputs.rds"))
saveRDS(parameter_outputs, here("regression_outputs", "parameter_outputs.rds"))
saveRDS(binary_output, here("regression_outputs", "main_binary_output.rds"))
saveRDS(vcov_binary_outputs, here("regression_outputs", "vcov_binary_outputs.rds"))
saveRDS(vcov_binary_outputs_dk, here("regression_outputs", "vcov_binary_outputs_dk.rds"))
saveRDS(days_output, here("regression_outputs", "main_days_output.rds"))
saveRDS(coeff_summary, here("regression_outputs", "coeff_summary.rds"))
