################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: April 27th, 2023

# Notes: need to run mixtures_modeling.R AND mixtures_dr.R
# Order to run - 5
#################################################

#----------------------------------------------------------------------
#---------- Compare Predicted vs Measured ---------
#---------------------------------------------------------------------- 

#plot together
MIN_LOG_DOSE <- (-4)  # How to show zero dose values on log scale
MAX_LOG_DOSE <- 3  # Just for scaling

# read in data and prep for plotting
GCA_adj <- read.csv("GCA_df_adj.csv")
GCA_unadj <- read.csv("GCA_df_noadj.csv")
GCA_df <- rbind(GCA_adj,GCA_unadj)
GCA_df$active <- "Active"
GCA_df1 <- GCA_df
GCA_df1$active <- "All"
GCA_df <- rbind(GCA_df, GCA_df1)
GCA_df$active <- factor(GCA_df$active, levels = c("Active", "All"),
                       labels = c("Active Chemicals", "All Chemicals"))

IA_adj <- read.csv("IA_df_adj.csv")
IA_unadj <- read.csv("IA_df_noadj.csv")
IA_df <- rbind(IA_adj,IA_unadj)
IA_df$active <- "Active"
IA_df1 <- IA_df
IA_df1$active <- "All"
IA_df <- rbind(IA_df, IA_df1)
IA_df$active <- factor(IA_df$active, levels = c("Active", "All"),
                        labels = c("Active Chemicals", "All Chemicals"))

DA_adj <- read.csv("DA_df_adj.csv")
DA_unadj <- read.csv("DA_df_noadj.csv")
DA_df <- rbind(DA_adj,DA_unadj)
DA_df$active <- "Active"
DA_df1 <- DA_df
DA_df1$active <- "All"
DA_df <- rbind(DA_df, DA_df1)
DA_df$active <- factor(DA_df$active, levels = c("Active", "All"),
                       labels = c("Active Chemicals", "All Chemicals"))

measured_df <- read.csv("predict_df.csv")

# Modify the labels
measured_df$active <- factor(measured_df$active, levels = c("Active", "All"),
                       labels = c("Active Chemicals", "All Chemicals"))

DA_df$category <- factor(DA_df$category, levels = c("Adjusted", "Unadjusted"), 
                        labels = c("% Contribution Adjusted", "% Contribution Unadjusted"))
IA_df$category <- factor(IA_df$category, levels = c("Adjusted", "Unadjusted"), 
                        labels = c("% Contribution Adjusted", "% Contribution Unadjusted"))
GCA_df$category <- factor(GCA_df$category, levels = c("Adjusted", "Unadjusted"), 
                          labels = c("% Contribution Adjusted", "% Contribution Unadjusted"))

#################################################################################################
#### Plots ####
# EM
#Active_Dose_uM
#Total_Dose_uM
compare_EM <- ggplot()+
  geom_line(data = subset(measured_df, mixture == "EM" & dosecalc =="Active_Dose_uM"), aes(x= log10(x), y = Prediction, color = "Measured in Vitro"), linewidth = 1)+
  geom_ribbon(data=subset(measured_df, mixture == "EM" & dosecalc =="Active_Dose_uM"), aes(x=log10(x), y=Prediction, ymin=Lower, ymax=Upper), fill = "#5A5A5A", alpha=0.2) +
  
  geom_line(data = subset(DA_df, mix_ratio == "EM"), aes(x= log10(mean), y = y, color = "CA"), linewidth = 1)+
  geom_ribbon(data=subset(DA_df, mix_ratio == "EM"), aes(x=log10(mean), y = y, xmin=log10(x_lower), 
                                                                 xmax=log10(x_upper)), fill = "#FDE725FF", alpha=0.2) +
  
  geom_line(data = subset(GCA_df, mix_ratio == "EM"), aes(x= log10(x), y = y, color = "GCA"), linewidth = 1)+
  geom_ribbon(data=subset(GCA_df, mix_ratio == "EM"), aes(x=log10(x), y=y, ymin=y_lower, ymax=y_upper), fill = "#7AD151FF",alpha=0.2) +
  
  geom_line(data = subset(IA_df, mix_ratio == "EM"), aes(x= log10(x), y = mean, color = "IA"), linewidth = 1)+
  geom_ribbon(data=subset(IA_df, mix_ratio == "EM"), aes(x=log10(x), y= mean, ymin=y_lower, ymax=y_upper), fill = "#2A788EFF",alpha=0.2) +
  
  theme_bw()+
  ylim(0, 90)+
  facet_grid(category~active)+
  labs(y="% Max MeBio Response", x= "Log10 Concentration (uM)", color = "Mixture", fill = "Mixture")+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured in Vitro" = "#5A5A5A"),
                     labels = c("CA Predict", "GCA Predict", "IA Predict", "Measured"))+
  theme(legend.position="right",
        axis.title.y=element_blank())
compare_EM
ggsave("compare_pred_meas_EM.jpg",compare_EM,  height = 8, width = 8)

# ED10
compare_ED10 <- ggplot()+
  geom_line(data = subset(measured_df, mixture == "ED10" & dosecalc =="Active_Dose_uM"), aes(x= log10(x), y = Prediction, color = "Measured in Vitro"), linewidth = 1)+
  geom_ribbon(data=subset(measured_df, mixture == "ED10" & dosecalc =="Active_Dose_uM"), aes(x=log10(x), y=Prediction, ymin=Lower, ymax=Upper), fill = "#5A5A5A", alpha=0.2) +
  
  geom_line(data = subset(DA_df, mix_ratio == "ED10"), aes(x= log10(mean), y = y, color = "CA"), linewidth = 1)+
  geom_ribbon(data=subset(DA_df, mix_ratio == "ED10"), aes(x=log10(mean), y = y, xmin=log10(x_lower), 
                                                         xmax=log10(x_upper)), fill = "#FDE725FF", alpha=0.2) +
  
  geom_line(data = subset(GCA_df, mix_ratio == "ED10"), aes(x= log10(x), y = y, color = "GCA"), linewidth = 1)+
  geom_ribbon(data=subset(GCA_df, mix_ratio == "ED10"), aes(x=log10(x), y=y, ymin=y_lower, ymax=y_upper), fill = "#7AD151FF",alpha=0.2) +
  
  geom_line(data = subset(IA_df, mix_ratio == "ED10"), aes(x= log10(x), y = mean, color = "IA"), linewidth = 1)+
  geom_ribbon(data=subset(IA_df, mix_ratio == "ED10"), aes(x=log10(x), y=mean, ymin=y_lower, ymax=y_upper), fill = "#2A788EFF",alpha=0.2) +
  
  theme_bw()+
  ylim(0, 90)+
  facet_grid(category~active)+
  labs(y="% Max MeBio Response", x= "Log10 Concentration (uM)", color = "Mixture", fill = "Mixture")+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured in Vitro" = "#5A5A5A"),
                     labels = c("CA Predict", "GCA Predict", "IA Predict", "Measured"))+
  theme(legend.position="right",
        axis.title.y=element_blank())
compare_ED10
ggsave("compare_pred_meas_ED10.jpg",compare_ED10,  height = 10, width = 10)

# ED50
compare_ED50 <- ggplot()+
  geom_line(data = subset(measured_df, mixture == "ED50" & dosecalc =="Active_Dose_uM"), aes(x= log10(x), y = Prediction, color = "Measured in Vitro"), linewidth = 1)+
  geom_ribbon(data=subset(measured_df, mixture == "ED50" & dosecalc =="Active_Dose_uM"), aes(x=log10(x), y= Prediction, ymin=Lower, ymax=Upper), fill = "#5A5A5A", alpha=0.2) +
  
  geom_line(data = subset(DA_df, mix_ratio == "ED50"), aes(x= log10(mean), y = y, color = "CA"), linewidth = 1)+
  geom_ribbon(data=subset(DA_df, mix_ratio == "ED50"), aes(x=log10(mean), y = y, xmin=log10(x_lower), 
                                                           xmax=log10(x_upper)), fill = "#FDE725FF", alpha=0.2) +
  
  geom_line(data = subset(GCA_df, mix_ratio == "ED50"), aes(x= log10(x), y = y, color = "GCA"), linewidth = 1)+
  geom_ribbon(data=subset(GCA_df, mix_ratio == "ED50"), aes(x=log10(x), y=y, ymin=y_lower, ymax=y_upper), fill = "#7AD151FF",alpha=0.2) +
  
  geom_line(data = subset(IA_df, mix_ratio == "ED50"), aes(x= log10(x), y = mean, color = "IA"), linewidth = 1)+
  geom_ribbon(data=subset(IA_df, mix_ratio == "ED50"), aes(x=log10(x), y=mean, ymin=y_lower, ymax=y_upper), fill = "#2A788EFF",alpha=0.2) +
  
  theme_bw()+
  facet_grid(category~active)+
  ylim(0, 90)+
  labs(y="% Max MeBio Response", x= "Log10 Concentration (uM)", color = "Mixture", fill = "Mixture")+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured in Vitro" = "#5A5A5A"),
                     labels = c("CA Predict", "GCA Predict", "IA Predict", "Measured"))+
  theme(legend.position="right", axis.title.y=element_blank())
compare_ED50
ggsave("compare_pred_meas_ED50.jpg",compare_ED50,  height = 8, width = 8)


combined_plot_all <- ggarrange(compare_EM, compare_ED10, compare_ED50,
                           ncol = 3,
                           vjust =1,
                           labels = "AUTO",
                           #widths = c(1.5, 1),
                           common.legend = TRUE, 
                           legend = "bottom")
combined_plot_all


ggsave("compare_dr.jpg", combined_plot_all,  height =5, width =10)

