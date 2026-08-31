#for data cleaning and manipulation
library(tidyverse)
#for coefficient plots
library(ggplot2)
#for feols output tables
library(fixest)
#for LaTeX table extraction
library(xtable)
#for sourcing data
library(here)

#load necessary data
panel_all <- readRDS(here("source_data", "panel_all.rds"))
cont_0001_asinh_output <- readRDS(here("source_data", "cont_0001_asinh_output.rds"))
cont_0001_asinh_output2 <- readRDS(here("source_data", "cont_0001_asinh_output2.rds"))
endpoint_outputs <- readRDS(here("source_data", "endpoint_outputs.rds"))
parameter_outputs <- readRDS(here("source_data", "parameter_outputs.rds"))
binary_output <- readRDS(here("source_data", "binary_output.rds"))
vcov_outputs <- readRDS(here("source_data", "vcov_outputs.rds"))
vcov_outputs_dk <- readRDS(here("source_data", "vcov_outputs_dk.rds"))
days_output <- readRDS(here("source_data", "days_output.rds"))

#create data frame with key summary statistics
summary_panel <- bind_rows(
  panel_all %>% 
    select(starts_with(c("price", "cont", "asinh", "binary", "days")), ppt, temp) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Mean"),
  panel_all %>% 
    select(starts_with(c("price", "cont", "asinh", "binary", "days")), ppt, temp) %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Median"),
  panel_all %>% 
    select(starts_with(c("price", "cont", "asinh", "binary", "days")), ppt, temp) %>%
    summarise(across(everything(), ~ sd(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Std. Dev."),
  panel_all %>% 
    select(starts_with(c("price", "cont", "asinh", "binary", "days")), ppt, temp) %>%
    summarise(across(everything(), ~ min(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Min"),
  panel_all %>% 
    select(starts_with(c("price", "cont", "asinh", "binary", "days")), ppt, temp) %>%
    summarise(across(everything(), ~ max(.x, na.rm = TRUE))) %>%
    mutate(statistic = "Max"),
  panel_all %>% 
    select(starts_with(c("price", "cont", "asinh", "binary", "days")), ppt, temp) %>%
    summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
    mutate(statistic = "Obs.")) %>%
  select(statistic, price, price_log, everything()) %>%
  pivot_longer(cols = -statistic, names_to = "Variable") %>%
  pivot_wider(names_from = statistic, values_from = value)
#output to LaTeX
print(xtable(summary_panel), file = here("tables_plots", "summary_stats.tex"), booktabs = TRUE, include.rownames = FALSE)

#plots coefficients of primary continuous exposure model on event horizon
cont_0001_asinh_plot <- ggplot(data = cont_0001_asinh_output, aes(x = horizon, y = Estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = 0, color = "red", alpha = 0.1, linewidth = 5) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE)) + 
  scale_x_continuous(breaks = seq(-12, 12, 2)) + 
  labs(x = "Exposure Horizon", y = "Coefficient Estimate", title = "Continuous Exposure Model Results") +
  theme_minimal() +
  theme(panel.background = element_blank())
ggsave(here("tables_plots", "cont_0001_asinh_plot.pdf"), plot = cont_0001_asinh_plot, device = "pdf")

#plots coefficients of primary continuous model with and without endpoints
endpoint_plot <- ggplot(data = endpoint_outputs, aes(x = horizon, y = Estimate, color = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = 0, color = "red", alpha = 0.1, linewidth = 5) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE), position = position_dodge(width = 0.4)) + 
  scale_x_continuous(breaks = seq(-12, 12, 2)) + 
  labs(x = "Exposure Horizon", y = "Coefficient Estimate", title = "Model Results With and Without Endpoint Binning",
       caption = str_wrap("Coefficients are presented with 95% confidence intervals. The event horizon is on the x-axis, with the lefthand side showing anticipation trends and the righthand side showing the effects after a wildfire occurs.", 90)) +
  theme_minimal() +
  theme(panel.background = element_blank(), legend.position = "bottom", plot.caption = element_text(face = "italic", hjust = 0))
ggsave(here("tables_plots", "endpoint_plot.pdf"), plot = endpoint_plot, device = "pdf")

#plots coefficients of continuous model with varying exponential decay parameters
parameter_plot <- ggplot(data = parameter_outputs, aes(x = horizon, y = Estimate, color = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = 0, color = "red", alpha = 0.1, linewidth = 5) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE), position = position_dodge(width = 0.3)) + 
  scale_x_continuous(breaks = seq(-12, 12, 2)) + 
  labs(x = "Exposure Horizon", y = "Coefficient Estimate", title = "Model Results With Varying Exponential Decay Parameters",
       caption = str_wrap("Coefficients are presented with 95% confidence intervals. The event horizon is on the x-axis, with the lefthand side showing anticipation trends and the righthand side showing the effects after a wildfire occurs.", 90)) +
  theme_minimal() +
  theme(panel.background = element_blank(), legend.position = "bottom", plot.caption = element_text(face = "italic", hjust = 0))
ggsave(here("tables_plots", "parameter_plot.pdf"), plot = parameter_plot, device = "pdf")

#plots coefficients of binary bins exposure model
binary_plot <- ggplot(data = binary_output %>% filter(horizon == "Lag"), aes(x = distance_mean, y = estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  labs(x = "Distance from Fire (km)", y = "Coefficient Estimate", title = "Binary Bin Exposure Model Results", caption = 
         str_wrap("Cumulative Lag Coefficients (Lags 1 through 11) are presented with 95% confidence intervals for each of the three distance bins (0-10km, 10-50km, and 50-100km)", 90)) +
  theme_minimal() +
  theme(panel.background = element_blank(), plot.caption = element_text(face = "italic", hjust = 0))
ggsave(here("tables_plots", "binary_plot.pdf"), plot = binary_plot, device = "pdf")

#plots coefficients of binary bins exposure model with varying standard error specifications
vcov_plot <- ggplot(data = vcov_outputs, aes(x = distance_mean, y = estimate, color = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 12)) +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  labs(x = "Distance from Fire (km)", y = "Coefficient Estimate", title = "Binary Model Results With Varying Standard Error Specs", caption = 
         str_wrap("Cumulative Lag Coefficients (Lags 1 through 11) are presented with 95% confidence intervals for each of the three distance bins (0-10km, 10-50km, and 50-100km)", 90)) +
  theme_minimal() +
  theme(panel.background = element_blank(), plot.caption = element_text(face = "italic", hjust = 0), legend.position = "bottom")
ggsave(here("tables_plots", "binary_vcovs.pdf"), plot = vcov_plot, device = "pdf")

#plots coefficients of binary bins exposure model with Driscoll-Kraay standard errors
vcov_plot_dk <- ggplot(data = vcov_outputs_dk, aes(x = distance_mean, y = estimate, color = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 10)) +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100)) + 
  labs(x = "Distance from Fire (km)", y = "Coefficient Estimate", title = "Binary Model Results With Varying Standard Error Specs", caption = 
         str_wrap("Cumulative Lag Coefficients (Lags 1 through 11) are presented with 95% confidence intervals for each of the three distance bins (0-10km, 10-50km, and 50-100km)", 90)) +
  theme_minimal() +
  theme(panel.background = element_blank(), plot.caption = element_text(face = "italic", hjust = 0), legend.position = "bottom")
ggsave(here("tables_plots", "binary_vcovs_dk.pdf"), plot = vcov_plot_dk, device = "pdf")

#create results table downloadable in LaTeX
etable(cont_0001_asinh_model, binary_model, days_model, 
       vcov = ~COUNTY + year,
       fitstat = c('n', 'ar2', 'war2'),
       file = here("tables_plots", "primary_results.tex"),
       tex = TRUE)

#create data frame presenting cumulative coefficient results for both bin models
bin_table <- bind_rows(binary_output, days_output) %>%
  mutate(Bin = paste0(str_split_fixed(distance_bin, "_", 2)[,1], "km to ", str_split_fixed(distance_bin, "_", 2)[,2], "km"),
         Coefficient = paste("Cumulative", horizon, "Sum", sep = " ")) %>%
  select(Model, Coefficient, Bin, estimate, std.error) %>%
  rename("Estimate" = estimate, "Std. Error" = std.error) %>%
  pivot_longer(cols = c(Estimate, 'Std. Error'), names_to = "Statistic", values_to = "Value") %>%
  pivot_wider(id_cols = c(Coefficient, Bin, Statistic), names_from = Model, values_from = Value)
#print LaTeX table
print(xtable(bin_table), file = here("tables_plots", "bin_table.tex"), booktabs = TRUE, include.rownames = FALSE)
