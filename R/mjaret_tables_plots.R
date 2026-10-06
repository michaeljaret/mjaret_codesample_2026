#for data cleaning and manipulation
library(tidyverse)
#for graphs and coefficient plots
library(ggplot2)
library(patchwork)
library(scales)
library(ggcorrplot)
#for feols output tables
library(fixest)
#for LaTeX table extraction
library(xtable)
#for sourcing data
library(here)

#load data frames
panel_all <- readRDS(here("intermediate_data", "panel_all.rds"))
main_cont_output <- readRDS(here("regression_outputs", "main_cont_output.rds"))
main_cont_output2 <- readRDS(here("regression_outputs", "main_cont_output2.rds"))
endpoint_outputs <- readRDS(here("regression_outputs", "endpoint_outputs.rds"))
parameter_outputs <- readRDS(here("regression_outputs", "parameter_outputs.rds"))
main_binary_output <- readRDS(here("regression_outputs", "main_binary_output.rds"))
vcov_binary_outputs <- readRDS(here("regression_outputs", "vcov_binary_outputs.rds"))
vcov_binary_outputs_dk <- readRDS(here("regression_outputs", "vcov_binary_outputs_dk.rds"))
main_days_output <- readRDS(here("regression_outputs", "main_days_output.rds"))
coeff_summary <- readRDS(here("regression_outputs", "coeff_summary.rds"))

panel_tex <- unpanel(panel_all)
panel_tex <- panel_tex %>% select(starts_with("price"), asinh_cont_log_0001, starts_with(c("binary", "days")), ppt, temp, pdsi)
#rename column names for presentation
names(panel_tex)
colnames(panel_tex) <- c("ZHVI", "Ln(ZHVI)", "Continuous", "Binary - Close", "Binary - Medium", "Binary - Far", "Days - Close", "Days - Medium", "Days - Far", "Precipitation (in.)", "Temperature (Deg. F.)", "PDSI")

#create data frame with key summary statistics
summary_panel <- bind_rows(
  panel_tex %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Mean"),
  panel_tex %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Median"),
  panel_tex %>%
    summarise(across(everything(), ~ sd(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Std. Dev."),
  panel_tex %>%
    summarise(across(everything(), ~ min(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Min"),
  panel_tex %>%
    summarise(across(everything(), ~ max(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Max"),
  panel_tex %>%
    summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
    mutate(statistic = "Obs.")) %>%
  select(statistic, everything()) %>%
  pivot_longer(cols = -statistic, names_to = "Variable") %>%
  pivot_wider(names_from = statistic, values_from = value)
#output to LaTeX
print(xtable(summary_panel), file = here("tables_plots/tables", "summary_panel.tex"), booktabs = TRUE, include.rownames = FALSE)

#create summary statistics table for continuous exposure variable
cont_exp_tex <- panel_all %>% select(starts_with(c("cont_", "asinh_")))
names(cont_exp_tex)
colnames(cont_exp_tex) <- c("0.0004, None", "0.0002, None", "0.0001, None", "0.00005, None", "0.000025, None",
                            "0.0004, Logged Magnitude Only", "0.0002, Logged Magnitude Only", "0.0001, Logged Magnitude Only", "0.00005, Logged Magnitude Only", "0.000025, Logged Magnitude Only",
                            "0.0004, IHS Only", "0.0002, IHS Only", "0.0001, IHS Only", "0.00005, IHS Only", "0.000025, IHS Only",
                            "0.0004, Both Transformations", "0.0002, Both Transformations", "0.0001, Both Transformations", "0.00005, Both Transformations", "0.000025, Both Transformations")
cont_exp_summary <- bind_rows(
  cont_exp_tex %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Mean"),
  cont_exp_tex %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Median"),
  cont_exp_tex %>%
    summarise(across(everything(), ~ sd(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Std. Dev."),
  cont_exp_tex %>%
    summarise(across(everything(), ~ min(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Min"),
  cont_exp_tex %>%
    summarise(across(everything(), ~ max(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Max")) %>%
  select(statistic, everything()) %>%
  pivot_longer(cols = -statistic, names_to = "Parameter, Transformation") %>%
  pivot_wider(names_from = statistic, values_from = value)
#output to LaTeX
print(xtable(cont_exp_summary), file = here("tables_plots/tables", "summary_cont_exp.tex"), booktabs = TRUE, include.rownames = FALSE)

#plots coefficients of primary continuous exposure model on event horizon
main_cont_plot <- ggplot(data = main_cont_output, aes(x = horizon, y = Estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = 0, color = "red", alpha = 0.1, linewidth = 5) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE)) + 
  scale_x_continuous(breaks = seq(-12, 12, 2), limits = c(-12.1, 12.1)) + 
  labs(x = "Exposure Horizon", y = "Coefficient Estimate") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line())
ggsave(here("tables_plots/plots", "main_cont_plot.pdf"), plot = main_cont_plot, device = "pdf")

#plots coefficients of primary continuous model with and without endpoints
endpoint_plot <- ggplot(data = endpoint_outputs, aes(x = horizon, y = Estimate, color = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = 0, color = "red", alpha = 0.1, linewidth = 5) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE), position = position_dodge(width = 0.4)) + 
  scale_x_continuous(breaks = seq(-12, 12, 2)) + 
  labs(x = "Exposure Horizon", y = "Coefficient Estimate") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(), legend.position = "bottom", plot.caption = element_text(face = "italic", hjust = 0))
ggsave(here("tables_plots/plots", "endpoint_plot.pdf"), plot = endpoint_plot, device = "pdf")

#plot cumulative lag and lead coefficients of binary and days models
binary_lag_plot <- ggplot(data = binary_output %>% filter(horizon == "Lag"), aes(x = distance_mean, y = estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), colour = "red") +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  scale_y_continuous(limits = c(-0.04, 0.04)) +
  labs(x = "Distance from Fire (km)", y = "Cumulative Lagged Coefficient Estimate", title = "Binary Model") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(), axis.title = element_text(size = 10), plot.title = element_text(size = 11))
ggsave(here("tables_plots/plots", "binary_lag_plot.pdf"), plot = binary_lag_plot, device = "pdf")

days_lag_plot <- ggplot(data = days_output %>% filter(horizon == "Lag"), aes(x = distance_mean, y = estimate, color = "blue")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),  colour = "blue") +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  scale_y_continuous(labels = label_number(), limits = c(-0.00175, 0.00175)) +
  labs(x = "Distance from Fire (km)", y = "Cumulative Lagged Coefficient Estimate", title = "Days Model") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(), axis.title = element_text(size = 10), plot.title = element_text(size = 11))
ggsave(here("tables_plots/plots", "days_lag_plot.pdf"), plot = days_lag_plot, device = "pdf")

binary_lead_plot <- ggplot(data = binary_output %>% filter(horizon == "Lead"), aes(x = distance_mean, y = estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), colour = "red") +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  scale_y_continuous(limits = c(-0.04, 0.04)) +
  labs(x = "Distance from Fire (km)", y = "Cumulative Leaded Coefficient Estimate", title = "Binary Model") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(), axis.title = element_text(size = 10), plot.title = element_text(size = 11))
ggsave(here("tables_plots/plots", "binary_lead_plot.pdf"), plot = binary_lead_plot, device = "pdf")

days_lead_plot <- ggplot(data = days_output %>% filter(horizon == "Lead"), aes(x = distance_mean, y = estimate, color = "blue")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),  colour = "blue") +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  scale_y_continuous(labels = label_number(), limits = c(-0.00175, 0.00175)) +
  labs(x = "Distance from Fire (km)", y = "Cumulative Leaded Coefficient Estimate", title = "Days Model") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(), axis.title = element_text(size = 10), plot.title = element_text(size = 11))
ggsave(here("tables_plots/plots", "days_lead_plot.pdf"), plot = days_lead_plot, device = "pdf")

#plot binary and days models together
bins_lead_plot <- binary_lead_plot + days_lead_plot + plot_layout(axes = "collect")
ggsave(here("tables_plots/plots", "bins_lead_plot.pdf"), plot = bins_lead_plot, device = "pdf", width = 5.93, height = 3, units = "in")

bins_lag_plot <- binary_lag_plot + days_lag_plot + plot_layout(axes = "collect")
ggsave(here("tables_plots/plots", "bins_lag_plot.pdf"), plot = bins_lag_plot, device = "pdf", width = 5.93, height = 3, units = "in")

#plots coefficients of binary bins exposure model with varying standard error specifications
vcov_binary_plot <- ggplot(data = vcov_binary_outputs, aes(x = distance_mean, y = estimate, color = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 12)) +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  labs(x = "Distance from Fire (km)", y = "Coefficient Estimate") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(), plot.caption = element_text(face = "italic", hjust = 0), legend.position = "bottom")
ggsave(here("tables_plots/plots", "binary_vcovs.pdf"), plot = vcov_binary_plot, device = "pdf")

#plots coefficients of binary bins exposure model with Driscoll-Kraay standard errors
vcov_binary_plot_dk <- ggplot(data = vcov_binary_outputs_dk %>% filter(!Model == "Spatial: Conley (200km Cutoff)" & !Model == "Spatial: Conley (400km Cutoff)"), aes(x = distance_mean, y = estimate, color = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 10)) +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  labs(x = "Distance from Fire (km)", y = "Coefficient Estimate") +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(), plot.caption = element_text(face = "italic", hjust = 0), legend.position = "bottom")
ggsave(here("tables_plots/plots", "binary_vcovs_dk.pdf"), plot = vcov_binary_plot_dk, device = "pdf")

#plots coefficients of binary bins exposure model

#create results table downloadable in LaTeX
etable(main_cont_model, main_binary_model, main_days_model, 
       vcov = ~COUNTY + year,
       fitstat = c('n', 'ar2', 'war2'),
       file = here("tables_plots/tables", "primary_results.tex"),
       tex = TRUE)

#create data frame presenting cumulative coefficient results for both bin models
bin_output_table <- bind_rows(binary_output, days_output) %>%
  mutate(Coefficient = paste("Cumulative", horizon, "Sum", sep = " ")) %>%
  select(Model, Coefficient, Bin, estimate, std.error) %>%
  rename("Estimate" = estimate, "Std. Error" = std.error) %>%
  pivot_longer(cols = c(Estimate, 'Std. Error'), names_to = "Statistic", values_to = "Value") %>%
  pivot_wider(id_cols = c(Coefficient, Bin, Statistic), names_from = Model, values_from = Value)
#print LaTeX table
print(xtable(bin_output_table, digits = 5), file = here("tables_plots/tables", "bin_table.tex"), booktabs = TRUE, include.rownames = FALSE)

#print LaTeX table for the cumulative lag coefficients of the three primary models
print(xtable(coeff_summary, digits = 5), file = here("tables_plots/tables", "coeff_summary.tex"), booktabs = TRUE, include.rownames = FALSE)

#define function to create correlation matrix by residualizing on fixed effects
horizon_terms <- c(paste0("Lag_", 1:3), "Exposure", paste0("Lead_", 1:3))
horizon_corr <- function(exposure_var, title, lab_size, text_size) {
  panel_horizon <- panel_all %>% unpanel() %>% 
    rename(Exposure = all_of(exposure_var)) %>% #renames exposure variable for presentation
    select(CITY, month_year, Exposure) %>%
    group_by(CITY) %>%
    mutate(
      Lag_1 = lag(Exposure, 1),
      Lag_2 = lag(Exposure, 2),
      Lag_3 = lag(Exposure, 3),
      Lead_1 = lead(Exposure, 1),
      Lead_2 = lead(Exposure, 2),
      Lead_3 = lead(Exposure, 3)) %>%
    filter(if_all(all_of(horizon_terms), ~ !is.na(.))) #creates lag and lead columns
  #creates residualized data frame after running each term against 1 with fixed effects
  resid_horizon <- sapply(horizon_terms, function(i) {
    resid(feols(as.formula(paste0(i, " ~ 1 | CITY + month_year")), 
                data = panel_horizon, panel.id = ~CITY + month_year, notes = FALSE))}) %>% 
    as.data.frame() %>%
    select(Lag_3, Lag_2, Lag_1, everything()) #re-orders matrix for presentation
  #creates correlation matrix and returns plot
  cor_mat <- round(cor(resid_horizon, use = "pairwise.complete.obs"), 2)
  return(ggcorrplot(cor_mat, title = title, lab = TRUE, lab_size = lab_size, lab_col = "black") + 
           theme(axis.text.x = element_text(size = text_size), axis.text.y = element_text(size = text_size)))
}

cont_corrplot <- horizon_corr("asinh_cont_log_0001", "", lab_size = 4, text_size = 10)
ggsave(here("tables_plots/plots", "cont_corrplot.pdf"), plot = cont_corrplot, device = "pdf")

binary_close_corrplot <- horizon_corr("binary_0_10", "Binary - Close", 1.3, 6)
binary_medium_corrplot <- horizon_corr("binary_10_50", "Binary - Medium", 1.3, 6)
binary_far_corrplot <- horizon_corr("binary_50_100", "Binary - Far", 1.3, 6)
#combines all three plots
binary_distances_corrplot <- binary_close_corrplot + binary_medium_corrplot + binary_far_corrplot + plot_layout(guides = "collect")
ggsave(here("tables_plots/plots", "binary_distances_corrplot.pdf"), plot = binary_distances_corrplot, device = "pdf", width = 5.93, height = 2.5, units = "in")

days_close_corrplot <- horizon_corr("days_0_10", "Days - Close", 1.3, 6)
days_medium_corrplot <- horizon_corr("days_10_50", "Days - Medium", 1.3, 6)
days_far_corrplot <- horizon_corr("days_50_100", "Days - Far", 1.3, 6)
days_distances_corrplot <- days_close_corrplot + days_medium_corrplot + days_far_corrplot + plot_layout(guides = "collect")
ggsave(here("tables_plots/plots", "days_distances_corrplot.pdf"), plot = days_distances_corrplot, device = "pdf", width = 5.93, height = 2.5, units = "in")

#creates correlation matrix and plot without event horizon for bin variables
bin_corr <- function(cols, title) {
  panel_bin <- panel_all %>% unpanel() %>%
    select(CITY, month_year, all_of(cols)) %>%
    group_by(CITY)
  
  resid_bin <- sapply(cols, function(i) {
    resid(feols(as.formula(paste0(i, " ~ 1 | CITY + month_year")), 
                data = panel_bin, panel.id = ~CITY + month_year, notes = FALSE))}) %>%
    as.data.frame()
  
  cor_mat <- round(cor(resid_bin, use = "pairwise.complete.obs"), 2)
  #change names for presentation
  rownames(cor_mat) <- c("Close", "Medium", "Far")
  colnames(cor_mat) <- c("Close", "Medium", "Far")
  return(ggcorrplot(cor_mat, title = title, lab = TRUE, lab_size = 3, lab_col = "black") + 
           theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)))
}

binary_all_corrplot <- bin_corr(c("binary_0_10", "binary_10_50", "binary_50_100"), "Binary Bins")
days_all_corrplot <- bin_corr(c("days_0_10", "days_10_50", "days_50_100"), "Days Bins")
#combines both plots
bins_all_corrplot <- binary_all_corrplot + days_all_corrplot + plot_layout(guides = "collect")
ggsave(here("tables_plots/plots", "bins_all_corrplot.pdf"), plot = bins_all_corrplot, device = "pdf", width = 5.93, height = 3, units = "in")
