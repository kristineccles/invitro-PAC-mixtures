################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: April 17th, 2024
# Note: Must run mixture_modeling.R
# Import files are very large ~ 4GB

#################################################

# Load libraries
library(drc)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggh4x)
library(ggpubr)

FIXED_C = 0 #lower limit

#################################################
#### Import Data ####

file_paths <- list.files(path = "mix_pred_boot", pattern = "*.csv", full.names = TRUE)
data_list <- map(file_paths, read.csv)
data_unlist <- na.omit(as.data.frame(bind_rows(data_list)))

# Create the `id` column
data_unlist$id <- paste(data_unlist$iteration, data_unlist$group, data_unlist$mix, data_unlist$method)

# Apply `is.finite()` only to numeric columns
numeric_cols <- sapply(data_unlist, is.numeric)
data_unlist <- data_unlist[rowSums(sapply(data_unlist[, numeric_cols], is.finite)) == ncol(data_unlist[, numeric_cols]), ]


# Define function to run individual model analyses
group_params <- function(data, fixed_c_value) {
  
  # Initialize vectors for storing results
  iterVector <- vector(mode = "numeric")
  groupVector <- vector(mode = "character")
  slopeVector <- vector(mode = "numeric")
  topVector <- vector(mode = "numeric")
  ec50Vector <- vector(mode = "numeric")
  
  # Initialize an empty data frame to store output data
  output_data <- data.frame(groupVector = character(),
                            slopeVector = numeric(),
                            topVector = numeric(),
                            ec50Vector = numeric(),
                            stringsAsFactors = FALSE)
  
  # Get unique id
  unique_itter <- unique(data$id)
  
  # Loop through each id
  for (i in unique_itter) {
    
    # Subset data for the current id
    sub_data <- subset(data, id == i)
    
    # Apply the selected model
    i_model <- drm(y ~ x, data = sub_data, fct = LL.4(fixed = c(NA, fixed_c_value, NA, NA),
                                                      names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    
    # Extract coefficients
    coefficient <- i_model$coefficients
    
    # Store results for this iteration
    output_data <- rbind(output_data, data.frame(
      groupVector = paste(unique(sub_data$group), unique(sub_data$mix), unique(sub_data$method)),
      slopeVector = coefficient[1],
      topVector = coefficient[2],
      ec50Vector = coefficient[3]
    ))
    print(output_data)
  }
  
  # Calculate summary statistics for the individual model
  measured_summary <-  output_data %>%
    group_by(group = groupVector)%>%
    summarize(slope_mean = quantile(slopeVector * -1, 0.5),  # Adjust for negative slope in hill models
              slope_lower = quantile(slopeVector * -1, 0.025),
              slope_upper = quantile(slopeVector * -1, 0.975),
              top = quantile(topVector, 0.5),
              top_lower = quantile(topVector, 0.025),
              top_upper = quantile(topVector, 0.975),
              EC50 = quantile(ec50Vector, 0.5),
              EC50_lower = quantile(ec50Vector, 0.025),
              EC50_upper = quantile(ec50Vector, 0.975)
    ) %>%
    as.data.frame()
  
  print(measured_summary)
}

# Example call with  dataset and specific FIXED_C constant
# Assuming data is defined and FIXED_C is a constant
output_params <- group_params(data_unlist, FIXED_C)

# Set up for plotting 
output_edit<- output_params %>%
  separate(group, into = c("group", "mixture", "model"), sep = " ", remove = FALSE)%>%
  as.data.frame()
output_edit$mixture<- gsub("ED10", "EXP1", output_edit$mixture)
output_edit$mixture<- gsub("ED50", "EXP2", output_edit$mixture)
output_edit$id <- paste(output_edit$group, output_edit$mixture, output_edit$model)
output_edit$id2 <- paste( output_edit$mixture, output_edit$model)
output_edit$group <- factor(output_edit$group, levels = c("adj", "unadj"), 
                            labels = c("% Contribution equals 100", "% Contribution does not equal 100"))

output_edit$active <- "Active Chemicals"

output_edit_all <- output_edit
output_edit_all$active <- "All Chemicals"

#Format measured data
mix_model_coeff <- tidy(mix_model, conf.int = TRUE)%>%
  dplyr::select(term, curve, estimate, conf.low, conf.high)%>% 
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high))%>%
  mutate(curve = curve,
    slope_mean = estimate_Slope*-1,
    slope_lower = conf.low_Slope*-1,
    slope_upper = conf.high_Slope*-1,
    top = `estimate_Upper Limit`,
    top_lower = `conf.low_Upper Limit`,
    top_upper = `conf.high_Upper Limit`,
    EC50 = estimate_ED50,
    EC50_lower = conf.low_ED50,
    EC50_upper = conf.high_ED50,
    .keep = "none")%>%
  
  as.data.frame()
mix_model_coeff


# Prepare Naming for Plotting
mix_model_coeff_edit<- mix_model_coeff %>%
  separate(curve, into = c("mixture", "adj", "active"), sep = " ", remove = FALSE)%>%
  as.data.frame()
mix_model_coeff_edit <- subset(mix_model_coeff_edit, mixture != "inVitro")
mix_model_coeff_edit$mixture<- gsub("ED10", "EXP1", mix_model_coeff_edit$mixture)
mix_model_coeff_edit$mixture<- gsub("ED50", "EXP2", mix_model_coeff_edit$mixture)
mix_model_coeff_edit$id <- paste(mix_model_coeff_edit$mixture, mix_model_coeff_edit$model)
mix_model_coeff_edit$id2 <- paste( mix_model_coeff_edit$mixture, mix_model_coeff_edit$model)
mix_model_coeff_edit$model <- "Measured"

mix_model_coeff_edit$active <- factor(mix_model_coeff_edit$active, levels = c("Active", "All"),
                        labels = c("Active Chemicals", "All Chemicals"))

#### Plot: Active ####

active_slope<- ggplot()+
  #modeled
  geom_point(data = subset(output_edit, group == "% Contribution equals 100"), aes(x = slope_mean, y = model, group = id, color = model), size = 2)+
  geom_errorbar(data = subset(output_edit, group == "% Contribution equals 100"), aes(x = slope_mean, y = model, xmin=slope_lower, xmax=slope_upper, group = id, 
                                        color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_model_coeff_edit, active == "Active Chemicals"), aes(x = slope_mean, y = model, group = active, color = model,), size = 2)+
  geom_errorbar(data = subset(mix_model_coeff_edit, active == "Active Chemicals"), aes(x = slope_mean, y = model, xmin=slope_lower, xmax=slope_upper, group = active, 
                                                 color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_nested(group ~ mixture+active )+
  scale_color_manual(name = "Group",
                     values =  c(DA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Slope")
active_slope

active_top<- ggplot()+
  #modeled
  geom_point(data = subset(output_edit, group == "% Contribution equals 100"), aes(x = top, y = model, group = id, color = model), size = 2)+
  geom_errorbar(data = subset(output_edit, group == "% Contribution equals 100"), aes(x = top, y = model, xmin=top_lower, xmax=top_upper, group = id, 
                                        color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_model_coeff_edit, active == "Active Chemicals"), aes(x = top, y = model, group = active, color = model,), size = 2)+
  geom_errorbar(data = subset(mix_model_coeff_edit, active == "Active Chemicals"), aes(x = top, y = model, xmin=top_lower, xmax=top_upper, group = active, 
                                                                                       color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_nested(group ~ mixture+active )+
  scale_color_manual(name = "Group",
                     values =  c(DA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Top (%)")
active_top

active_EC50<- ggplot()+
  #modeled
  geom_point(data = subset(output_edit, group == "% Contribution equals 100"), aes(x = log10(EC50), y = model, group = id, color = model), size = 2)+
  geom_errorbar(data = subset(output_edit, group == "% Contribution equals 100"), aes(x = log10(EC50), y = model, xmin=log10(EC50_lower), xmax=log10(EC50_upper), group = id, 
                                        color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_model_coeff_edit, active == "Active Chemicals"), aes(x = log10(EC50), y = model, group = active, color = model,), size = 2)+
  geom_errorbar(data = subset(mix_model_coeff_edit, active == "Active Chemicals"), aes(x = log10(EC50), y = model, xmin=log10(EC50_lower), xmax=log10(EC50_upper), group = active, 
                                                                                       color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_nested(group ~ mixture+active )+
  scale_color_manual(name = "Group",
                     values =  c(DA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Log10 EC50")
active_EC50

combined_CI_active <- ggarrange(active_slope, active_top, active_EC50,
                           ncol = 3,
                           vjust =1,
                           nrow = 1,
                           labels = c("D", "E", "F"),
                           common.legend = TRUE,
                           legend = "bottom")
combined_CI_active

ggsave("combined_CI_plot_active.jpg", combined_CI_active,  height =6, width =14)


#### Plot: All ####
# Edit Mixtures 
mix_model_coeff_edit2 <- mix_model_coeff_edit[!grepl("Total_Dose_uM", mix_model_coeff_edit$curve), ]

all_slope<- ggplot()+
  #modeled
  geom_point(data = output_edit_all, aes(x = slope_mean, y = model, group = id, color = model), size = 2)+
  geom_errorbar(data = output_edit_all, aes(x = slope_mean, y = model, xmin=slope_lower, xmax=slope_upper, group = id, 
                                        color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_model_coeff_edit2, active == "All Chemicals"), aes(x = slope_mean, y = model, group = active, color = model), size = 2)+
  geom_errorbar(data = subset(mix_model_coeff_edit2, active == "All Chemicals"), aes(x = slope_mean, y = model, xmin=slope_lower, xmax=slope_upper, group = active, 
                                                                                       color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_nested(group ~ mixture+active)+
  scale_color_manual(name = "Group",
                     values =  c("DA" = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Slope")
all_slope

all_top<- ggplot()+
  #modeled
  geom_point(data = output_edit_all, aes(x = top, y = model, group = id, color = model), size = 2)+
  geom_errorbar(data = output_edit_all, aes(x = top, y = model, xmin=top_lower, xmax=top_upper, group = id, 
                                        color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_model_coeff_edit2, active == "All Chemicals"), aes(x = top, y = model, group = active, color = model), size = 2)+
  geom_errorbar(data = subset(mix_model_coeff_edit2, active == "All Chemicals"), aes(x = top, y = model, xmin=top_lower, xmax=top_upper, group = active, 
                                                                                       color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_nested(group ~ mixture+active)+
  scale_color_manual(name = "Group",
                     values =  c("DA" = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Top (%)")
all_top

all_EC50<- ggplot()+
  #modeled
  geom_point(data = output_edit_all, aes(x = log10(EC50), y = model, group = id, color = model), size = 2)+
  geom_errorbar(data = output_edit_all, aes(x = log10(EC50), y = model, xmin=log10(EC50_lower), xmax=log10(EC50_upper), group = id, 
                                        color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_model_coeff_edit2, active == "All Chemicals" | curve == "Active_Dose_uM"), aes(x = log10(EC50), y = model, group = active, color = model,), size = 2)+
  geom_errorbar(data = subset(mix_model_coeff_edit2, active == "All Chemicals" | curve == "Active_Dose_uM"), aes(x = log10(EC50), y = model, xmin=log10(EC50_lower), xmax=log10(EC50_upper), group = active, 
                                                                                       color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_nested(group ~ mixture+active)+
  scale_color_manual(name = "Group",
                     values =  c("DA" = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Log10 EC50")
all_EC50


combined_CI_all <- ggarrange(all_slope, all_top, all_EC50,
                                ncol = 3,
                                nrow = 1,
                                vjust =1,
                                labels = c("G", "H", "I"),
                                common.legend = TRUE,
                                legend = "bottom")
combined_CI_all

#### Plot Together ####


boxplot_combined <- ggarrange(all_slope, all_top, all_EC50,
                                active_slope, active_top, active_EC50,
                                ncol = 3,
                                nrow = 2,
                                heights = c(1.5,1),
                                labels = c("D", "E", "F", "G", "H", "I"),
                                common.legend = TRUE,
                                legend = "bottom")

boxplot_combined


combined_plot_final <- ggarrange(combined_plot_all, boxplot_combined,
                             ncol = 1,
                             nrow = 2, 
                             widths = c(1,1.5),
                             common.legend = TRUE,
                             legend = "bottom")
combined_plot_final

ggsave("Figure 4.jpg", combined_plot_final,  height =16, width =14)


# compare just the EC50s

combined_EC50 <- ggarrange(all_EC50, active_EC50,
                                 ncol = 1,
                                 nrow = 2, 
                                heights = c(1.5, 1),
                                labels = c("A", "B"),
                                 common.legend = TRUE,
                                 legend = "bottom")
combined_EC50

ggsave("combined_EC50.jpg", combined_EC50,  height =12, width =5)

