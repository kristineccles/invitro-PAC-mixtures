################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: April 17th, 2024
# Note: Must run 03-mixture_modeling.R and 08-plot_model_params.R
# Import files are very large ~ 4GB
#################################################

# Load libraries
library(drc)
library(viridis)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(ggpmisc)
library(tidyverse)
library(reshape2)
library(tcplfit2)
library(tibble)

FIXED_C = 0

set.seed (4540)
#################################################
#### Import Data ####

file_paths <- list.files(path = "mix_pred_boot", pattern = "*.csv", full.names = TRUE)
data_list <- map(file_paths, read.csv)
data_unlist <- na.omit(as.data.frame(bind_rows(data_list)))
data_unlist[data_unlist == Inf | data_unlist == -Inf] <- NA
data_unlist <- na.omit(data_unlist)

# Create the `id` columns
data_unlist$id <- paste(data_unlist$iteration, data_unlist$group, data_unlist$mix, data_unlist$method)

model_bmd_calc <- list()

for (chem in unique(data_unlist$id)) {
  
  # Subset the data for this chemical
  chem_data <- (data_unlist[data_unlist$id == chem, ])
  #chem_data <- subset(data_unlist, id %in% chem)
  
  # Define BMR setup
  row <- list(
    conc = chem_data$x,
    resp = chem_data$y,
    bmed = 0,        # Baseline centered at 0
    cutoff = 0,   # cut off for hit call
    onesd = 10,       #10% increase = BMR but need to set bmr_scale = 1
    assay = "assay",
    name = chem
  )
  
  # Run Hill model + BMC calculation
  res <- concRespCore(
    row,
    fitmodels = c("hill"),
    conthits = TRUE,    # This enables BMC calculation
    aicc = FALSE,
    bidirectional = FALSE,
    errfun = "dnorm",
    bmr_scale = 1
  )
  
  # Store result
  model_bmd_calc[[chem]] <- res[,c("name", "assay", "bmd", "ac10", "bmdl", "bmdu", "tp")]
}

model_bmd_calc

# Extract BMD summary info from each chemical's result
model_melt <- melt(model_bmd_calc)

# Use dcast to reshape
tidy_model<- dcast(model_melt, L1 ~ variable, value.var = "value")
tidy_model <- tidy_model %>%
  separate(L1, into = c("itter", "group", "mix", "method"), sep = " ")
tidy_model$id2 <- paste(tidy_model$group, tidy_model$mix, tidy_model$method)

model_summary <- na.omit(tidy_model) %>%
  group_by(id2) %>%
  summarise(
    bmdl  = quantile(bmd, 0.025),
    bmd_mean = quantile(bmd, 0.50),
    bmdu = quantile(bmd, 0.975),
    ac10_p5  = quantile(ac10, 0.025),
    ac10_p50 = quantile(ac10, 0.50),
    ac10_p95 = quantile(ac10, 0.975)
  )

model_bmd_coeff_edit<- model_summary %>%
  separate(id2, into = c("group", "mixture","model"), sep = " ", remove = FALSE)%>%
  as.data.frame()


#Change Labels
model_bmd_coeff_edit$mixture<- gsub("ED10", "EXP1", model_bmd_coeff_edit$mixture)
model_bmd_coeff_edit$mixture<- gsub("ED50", "EXP2", model_bmd_coeff_edit$mixture)
model_bmd_coeff_edit$group <- factor(model_bmd_coeff_edit$group, levels = c("adj", "unadj"), 
                            labels = c("% Contribution equals 100", "% Contribution does not equal 100"))

model_bmd_coeff_edit$id <- paste(model_bmd_coeff_edit$mixture, model_bmd_coeff_edit$model)

# Create the plot
ggplot(tidy_model, aes(x = ac10, y = bmd, color = tp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95) +
  labs(title = "AC10 vs BMD10", x = "AC10",y = "BMD") +
  theme_minimal()

###################################################################################
#### Measured Data ####
#mix_bmd_calc <- as.data.frame(bmd(mix_model, bmr = 0.10, def = "relative", backg = 0, display = TRUE))
#mix_bmd_calc$curve <- rownames(mix_bmd_calc)

#mix_ec10_calc <- ED(mix_model, c(10))
#plot(mix_model)
# the two lines give the same result- BMC is not being calculated properly since all the tops 

# Calculate using tcplfit2
mix_bmd_calc2 <- list()

for (chem in unique(df_mix$grouping)) {
  
  # Subset the data for this chemical
  chem_data <- df_mix[df_mix$grouping == chem, ]
  
  # Define BMR setup
  row <- list(
    conc = chem_data$Dose_uM,
    resp = chem_data$MaxResp,
    bmed = 0,        # Baseline centered at 0
    cutoff = 0,   # cut off for hit call
    onesd = 10,       #10% increase = BMR but need to set bmr_scale = 1
    assay = "assay",
    name = chem
  )
  
  # Run Hill model + BMC calculation
  res <- concRespCore(
    row,
    fitmodels = c("hill"),
    conthits = TRUE,    # This enables BMC calculation
    aicc = FALSE,
    bidirectional = FALSE,
    errfun = "dnorm",
    bmr_scale = 1
  )
  
  # Store result
  mix_bmd_calc2[[chem]] <- res[,c("name", "assay", "bmd", "ac10", "bmdl", "bmdu", "tp")]
}

# Extract BMD summary info from each chemical's result
mix_melt <- melt(mix_bmd_calc2)

# Use dcast to reshape
tidy_mix <- dcast(mix_melt, L1 ~ variable, value.var = "value")

# Create the plot
ggplot(tidy_mix, aes(x = ac10, y = bmd, color = tp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95) +
  labs(title = "AC10 vs BMD10", x = "AC10",y = "BMD") +
  theme_minimal()

mix_bmd_coeff <- tidy_mix%>%
  mutate(curve = L1,
         bmd_mean = bmd,
         bmd_lower = bmdl,
         bmd_upper = bmdu,
         .keep = "none")%>%
  as.data.frame()
mix_bmd_coeff

mix_bmd_coeff_edit<- tidy_mix %>%
  separate(L1, into = c("mixture", "active", "adj"), sep = "  ", remove = FALSE)%>%
  as.data.frame()
mix_bmd_coeff_edit$model <- "Measured"

# remove invitro
mix_bmd_edit <- subset(mix_bmd_coeff_edit, mixture != "inVitro")

#Change Labels
mix_bmd_edit$mixture<- gsub("ED10", "EXP1", mix_bmd_edit$mixture)
mix_bmd_edit$mixture<- gsub("ED50", "EXP2", mix_bmd_edit$mixture)

mix_bmd_edit$active <- factor(mix_bmd_edit$active, levels = c("Active", "All"),
                             labels = c("Active Chemicals", "All Chemicals"))
mix_bmd_edit$id <- paste(mix_bmd_edit$mixture, mix_bmd_edit$model)
###########################################################################################################
#### Plot ####
active_bmd<- ggplot()+
  #modeled
  geom_point(data = subset(model_bmd_coeff_edit, group == "% Contribution equals 100"), aes(x = log10(bmd_mean), y = model, group = model, color = model), size = 2)+
  geom_errorbar(data = subset(model_bmd_coeff_edit, group == "% Contribution equals 100"), aes(x = log10(bmd_mean), y = model, xmin=log10(bmdl), xmax=log10(bmdu), group = id, 
                                                                                 color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_bmd_edit, active == "Active Chemicals"), aes(x = log10(bmd), y = model, group = active, color = model,), size = 2)+
  geom_errorbar(data = subset(mix_bmd_edit, active == "Active Chemicals"), aes(x = log10(bmd), y = model, xmin=log10(bmdl), xmax=log10(bmdu), group = active, 
                                                                                       color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_grid(group ~ mixture+active, drop = TRUE)+
  scale_color_manual(name = "Group",
                     values =  c(DA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Log10 BMD10")
active_bmd

mix_bmd_edit2 <- mix_bmd_edit[!grepl("Total_Dose_uM", mix_bmd_edit$adj), ]

all_bmd<- ggplot()+
  #modeled
  geom_point(data = model_bmd_coeff_edit, aes(x = log10(bmd_mean), y = model, group = id, color = model), size = 2)+
  geom_errorbar(data = model_bmd_coeff_edit, aes(x = log10(bmd_mean), y = model, xmin=log10(bmdl), xmax=log10(bmdu), group = id, 
                                                 color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  
  #measured
  geom_point(data = subset(mix_bmd_edit2, active == "All Chemicals" & adj == "Active_Dose_uM"), aes(x = log10(bmd), y = model, group = active, color = model,), size = 2)+
  geom_errorbar(data = subset(mix_bmd_edit2, active == "All Chemicals" & adj == "Active_Dose_uM"), aes(x = log10(bmd), y = model, xmin=log10(bmdl), xmax=log10(bmdu), group = active, 
                                         color = model), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  facet_grid(group ~ mixture+active, drop = TRUE)+
  scale_color_manual(name = "Group",
                     values =  c(DA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Log10 BMD10")
all_bmd


combined_bmd <- ggarrange(all_bmd, active_bmd,
                                ncol = 1,
                                nrow = 2, 
                                heights = c(1.5, 1),
                                common.legend = TRUE,
                                labels = c("A", "B"),
                                legend = "bottom")
combined_bmd

ggsave("combined_bmd10.jpg", combined_bmd,  height =12, width =5)





