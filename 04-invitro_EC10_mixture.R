################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: April 27th, 2023

# Notes: need to run mixtures_modeling.R AND mixtures_dr.R
# Order to run - 4
#################################################

# load libraries
library(drc)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(sjPlot)
library(cowplot)
library(data.table)
library(reshape2)
library(ggpubr)
library(broom)
library(tidyr)
library(viridis)
library(truncnorm)

#### Independent Action ####
# Calculate the array of predicted Y values for fixed X values
# For each x value, predicted y of mixture = 1 - PRODUCT(1 - predictedY(weightFactor*x))
#y = d/(1 + exp(b *(log(x) - e)))
#predY <- 1-prod(1/(1+exp(b*(log(x*fract)- log(e)))))

IA_invitro100_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$Invitro10_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_invitro100_CI <- cbind(x,  melt(IA_invitro100_CI))
colnames(unlist_invitro100_CI) <- c("x", "y_pred", "itteration")

IA_CI_invitro100_final <- unlist_invitro100_CI%>%
  group_by(x)%>%
  summarise(mean = median(y_pred),
            lower = quantile(y_pred,0.025, na.rm =TRUE),
            upper = quantile(y_pred, 0.975, na.rm =TRUE))%>%
  as.data.frame()

IA_CI_invitro100_final$mix_ratio <- "Invitro100"
IA_CI_invitro100_final

# Combined dataframes
IA_df_invitro <- IA_CI_invitro100_final
IA_df_invitro$Method <- "IA"

#rescale
IA_df_invitro$mean <- scales::rescale(IA_df_invitro$mean*100, to=c(0,invitro_top))
IA_df_invitro$lower <- scales::rescale(IA_df_invitro$lower*100, to=c(0,invitro_top))
IA_df_invitro$upper <- scales::rescale(IA_df_invitro$upper*100, to=c(0,invitro_top))

#plot
p1 <- ggplot()+
  geom_line(data = IA_df_invitro, aes(x= log10(x), y = mean, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=IA_df_invitro, aes(x=log10(x), y = mean, ymin=lower, ymax=upper, fill = mix_ratio), alpha=0.2) +
  theme_bw()+
  labs(x="Log10 Dose (uM)", y="% Max MeBio Response", color = "Mixing Ratio \nIndependent Action", 
       fill = "Mixing Ratio \nIndependent Action") +
  scale_fill_viridis(discrete= TRUE)+
  scale_color_viridis(discrete= TRUE)+
  theme_bw() +
  theme(legend.title = element_text(face="bold", size=11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 16, hjust=0.5))
p1

#### 2B: Concentration Addition ####
# x = (e^(e *b)* (d/y - 1))^(1/b)

DA_ED10_invitro_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$Invitro10_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope )*log((1-y)/y))))
    })})

unlist_DA_ED10_invitro_CI <- cbind(y,  melt(DA_ED10_invitro_CI))
colnames(unlist_DA_ED10_invitro_CI) <- c("y", "x_pred", "itteration")

DA_CI_ED10_invitro_final <- unlist_DA_ED10_invitro_CI%>%
  group_by(y)%>%
  summarise(mean = median(x_pred),
            lower = quantile(x_pred,0.025, na.rm =TRUE),
            upper = quantile(x_pred, 0.975, na.rm =TRUE))%>%
  as.data.frame()
DA_CI_ED10_invitro_final$mix_ratio <- "Invitro100"

# Combined dataframes
DA_df_invitro <- (DA_CI_ED10_invitro_final)
DA_df_invitro$Method <- "CA"

#rescale
DA_df_invitro$y <- scales::rescale(DA_df_invitro$ y*100, to=c(0,invitro_top))

#plot
p2 <- ggplot()+
  geom_line(data = DA_df_invitro, aes(x= log10(mean), y = y, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=DA_df_invitro, aes(x=log10(mean), y = y, xmin=log10(lower), xmax=log10(upper), fill = mix_ratio), alpha=0.2) +
  theme_bw()+
  labs(x="Log10 Dose (uM)", y="% Max MeBio Response", color = "Mixing Ratio \nDose Addition", 
       fill = "Mixing Ratio \nDose Addition") +
  theme_bw() +
  scale_fill_viridis(discrete= TRUE)+
  scale_color_viridis(discrete= TRUE)+
  theme(legend.title = element_text(face="bold", size=11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 16, hjust=0.5))
p2

#### 2C: Generalized Concentration Addition ####
#add mixing ratios
B1bootmat_invitro <- left_join(B1bootmat, individual_model_coeff2[,c("curve", "Invitro10_percent")],
                       by= c("chemical" = "curve"), keep=FALSE )
#make new dataframe for CI generation
B1bootmat_list_invitro <- split(B1bootmat_invitro, f=B1bootmat_invitro$itter)

# GCA Predictions
#Emax = (MAX1(C1/ED501) + MAX2(C2/ED502))/(1+(C1/ED501)+(C2/EC2))
#Emix = MAX1(C1/ED501) + MAX2(C2/ED502) + MAX1(C1/ED501)/1+(C1/ED501)+(C2/ED502) +(C3/ED503)

GCA_ED10_invitro_CI <- 
  lapply(1:length(B1bootmat_list_invitro), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      (sum((B1bootmat_list_invitro[[j]]$Top)*((x*B1bootmat_list_invitro[[j]]$Invitro10_percent)/B1bootmat_list_invitro[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list_invitro[[j]]$Invitro10_percent)/B1bootmat_list_invitro[[j]]$ED50)))
    })})

unlist_GCA_ED10_invitro_CI <- cbind(x,  melt(GCA_ED10_invitro_CI))
colnames(unlist_GCA_ED10_invitro_CI) <- c("x", "y_pred", "itteration")

GCA_CI_ED10_invitro_final <- unlist_GCA_ED10_invitro_CI%>%
  group_by(x)%>%
  summarise(mean = median(y_pred),
            lower = quantile(y_pred,0.025, na.rm =TRUE),
            upper = quantile(y_pred, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_ED10_invitro_final$mix_ratio <- "Invitro100"

# Combined dataframes
GCA_df_invitro <- (GCA_CI_ED10_invitro_final)
GCA_df_invitro$Method <- "GCA"

#plot
p3 <- ggplot()+
  geom_line(data = GCA_df_invitro, aes(x= log10(x), y = mean, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=GCA_df_invitro, aes(x=log10(x), y=mean, ymin=lower, ymax=upper, fill = mix_ratio), alpha=0.2) +
  theme_bw()+
  labs(x="Log10 Dose (uM)", y="% Max MeBio Response", color = "Mixing Ratio \nGCA", 
       fill = "Mixing Ratio \nGCA") +
  scale_fill_viridis(discrete= TRUE)+
  scale_color_viridis(discrete= TRUE)+
  theme_bw() +
  theme(legend.title = element_text(face="bold", size=11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 16, hjust=0.5))
p3

#----------------------------------------------------------------------
#---------- PART 3: Curve Fitting for Measured Mixtures  ---------
#---------------------------------------------------------------------- 
invitro_mix <- read.csv("mix_conc_df 5-16-2024.csv")
invitro_mix <- subset(invitro_mix, Chemicalname == "inVitro")
invitro_mix$Active_Dose_uM <- 10^(invitro_mix$DoseActivesOnly_logM)*1e6

mix_model<- drm(MaxResp~Active_Dose_uM, data=invitro_mix,
                type = "continuous",
                fct=LL.4(fixed=c(NA, 0 , NA, NA), 
                         names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                na.action = na.omit)
plot(mix_model)

#get coefficients
mix_model_CI <- tidy(mix_model, conf.int = TRUE)

#reorganize
mix_model_coeff <- mix_model_CI%>%
  dplyr::select(term, curve, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)%>%
  as.data.frame
mix_model_coeff

# add SE
mix_slope_SE <- as.data.frame(subset(mix_model_CI, term == "Slope")$std.error)
colnames(mix_slope_SE) <- "SE_slope"
mix_ED50_SE <- as.data.frame(subset(mix_model_CI, term == "ED50")$std.error)
colnames(mix_ED50_SE) <- "SE_ED50"
mix_UpperLimit_SE <- as.data.frame(subset(mix_model_CI, term == "Upper Limit")$std.error)
colnames(mix_UpperLimit_SE) <- "SE_UpperLimit"
mix_model_coeff_wCI <- cbind(mix_model_coeff, mix_slope_SE, mix_ED50_SE, mix_UpperLimit_SE)


ED10_Invitro<- drm(MaxResp~Active_Dose_uM, data=invitro_mix,
                   type = "continuous",
                   lowerl = c(-Inf,NA, 0),
                   upperl = c(0, 100, Inf),
                   fct=LL.4(fixed=c(NA, FIXED_C , NA, NA), 
                            names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

predict_ED10_Invitro <- as.data.frame(cbind(x, mixture ="ED10_Invitro", (predict(ED10_Invitro, newdata = x,
                                                                                 interval = "confidence", level = 0.95))))

invitro_predict_df <- predict_ED10_Invitro

invitro_predict_df[ invitro_predict_df<0 ] <- 0

summary_predict <- invitro_predict_df%>%
  group_by(mixture)%>%
  summarise(max_pred = max(Prediction),
            max_lower = max(Lower),
            max_upper = max(Upper))
summary_predict

#Plot
measured_mix_plot <- ggplot()+
  geom_line(data = invitro_predict_df, aes(x= log10(x), y = Prediction, color = mixture))+
  geom_ribbon(data = invitro_predict_df, aes(x=log10(x), y = Prediction, ymin=Lower, ymax=Upper, 
                                     fill = mixture), alpha=0.2) +
  theme_bw()+
  #geom_point(data= invitro_mix, aes(x = log10(Dose_uM), y = MaxResp, color = CURVE),stat = "summary", fun = mean) +
  labs(y="% Max MeBio Response", x= "Log10 Concentration (uM)", color = "Mixture", fill = "Mixture")
measured_mix_plot

#----------------------------------------------------------------------
#---------- PART 4: Compare Predicted vs Measured ---------
#---------------------------------------------------------------------- 

#plot together
MIN_LOG_DOSE <- (-5)  # How to show zero dose values on log scale
MAX_LOG_DOSE <- 3  # Just for scaling

# Plots

compare_ED10_100 <- ggplot()+
  
  geom_line(data = subset(DA_df_invitro, mix_ratio == "Invitro100"), aes(x= log10(mean), y = y, color = "CA"), linewidth = 1)+
  geom_ribbon(data=subset(DA_df_invitro, mix_ratio == "Invitro100"), aes(x=log10(mean), y = y, xmin=log10(lower), 
                                                           xmax=log10(upper)), fill = "#FDE725FF", alpha=0.2) +
  
  geom_line(data = subset(GCA_df_invitro, mix_ratio == "Invitro100"), aes(x= log10(x), y = mean, color = "GCA"), linewidth = 1)+
  geom_ribbon(data=subset(GCA_df_invitro, mix_ratio == "Invitro100"), aes(x=log10(x), y=mean, ymin=lower, ymax=upper), fill = "#7AD151FF",alpha=0.2) +
  
  geom_line(data = subset(IA_df_invitro, mix_ratio == "Invitro100"), aes(x= log10(x), y = mean, color = "IA"), linewidth = 1)+
  geom_ribbon(data=subset(IA_df_invitro, mix_ratio == "Invitro100"), aes(x=log10(x), y=mean, ymin=lower, ymax=upper), fill = "#2A788EFF",alpha=0.2) +
  
  geom_line(data = subset(invitro_predict_df), aes(x= log10(x), y = Prediction, color = "Measured in Vitro"), linewidth = 1)+
  geom_ribbon(data=subset(invitro_predict_df), aes(x=log10(x), y=Prediction, ymin=Lower, ymax=Upper), fill = "#5A5A5A", alpha=0.2) +
  xlim(c(-5, 5))+
  theme_bw()+
  labs(y="% Max MeBio Response", x= "Log10 Concentration (uM)", color = "Mixture", fill = "Mixture")+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured in Vitro" = "#5A5A5A"),
                     labels = c("CA Predict", "GCA Predict", "IA Predict", "Measured"))+
  theme(legend.position="bottom")
compare_ED10_100

invitro_pie<- ggplot(individual_model_coeff2, aes(x="", y=Invitro10_percent, fill=curve))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6))+
  labs(fill = "Chemical")
invitro_pie

#### Make point and SE plot ####
# tops are calculated in 02

unlist_invitro100_CI$y <- scales::rescale(unlist_invitro100_CI$y_pred*100, to=c(0,invitro_top))
dr_IA_invitro <- drm(y ~ x, data = unlist_invitro100_CI, 
             type = "continuous",
             lowerl = c(-Inf, 0, 0),
             upperl = c(0, 100, Inf),
             fct = LL.4(fixed = c(NA, FIXED_C, NA, NA), 
                        names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
fit_IA_invitro <- tidy(dr_IA_invitro)%>%
  as.data.frame()
fit_IA_invitro$method <- "IA"

unlist_DA_ED10_invitro_CI$y <- scales::rescale(unlist_DA_ED10_invitro_CI$y*100, to=c(0,invitro_top))
dr_CA_invitro <- drm(y ~ x_pred, data = unlist_DA_ED10_invitro_CI, 
                     type = "continuous",
                     #lowerl = c(-Inf, 0, 0),
                     #upperl = c(0, 100, Inf),
                     fct = LL.4(fixed = c(NA, FIXED_C, NA, NA), 
                                names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
fit_CA_invitro <- tidy(dr_CA_invitro)%>%
  as.data.frame()
fit_CA_invitro$method <- "CA"

dr_GCA_invitro <- drm(y_pred ~ x, data = unlist_GCA_ED10_invitro_CI, 
                     type = "continuous",
                     #lowerl = c(-Inf, 0, 0),
                     #upperl = c(0, 100, Inf),
                     fct = LL.4(fixed = c(NA, FIXED_C, NA, NA), 
                                names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
fit_GCA_invitro <- tidy(dr_GCA_invitro)%>%
  as.data.frame()
fit_GCA_invitro$method <- "GCA"

fit_measured_invitro <- tidy(ED10_Invitro)%>%
  as.data.frame()
fit_measured_invitro$method <- "Measured"

#combine
fit_compiled_invitro <- rbind(fit_GCA_invitro, fit_CA_invitro, fit_IA_invitro, fit_measured_invitro)
# make slope positive

fit_compiled_invitro$estimate[fit_compiled_invitro$term == "Slope"] <- 
  abs(fit_compiled_invitro$estimate[fit_compiled_invitro$term == "Slope"])

fit_compiled_invitro$lower <- fit_compiled_invitro$estimate-fit_compiled_invitro$std.error*1.96
fit_compiled_invitro$upper <- fit_compiled_invitro$estimate+fit_compiled_invitro$std.error*1.96

plot_all_invitro<- ggplot(data = fit_compiled_invitro, aes(x = (estimate), y = method, colour = method))+
  geom_errorbar(data = fit_compiled_invitro, aes(xmin=lower, xmax=upper), 
                width=.2,position=position_dodge(.9))+
  geom_point()+
  facet_grid(.~term, scales = "free")+
  theme_bw()+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"),
                     labels = c("CA Predict", "GCA Predict", "IA Predict", "Measured"))+
  ylab(" ")
plot_all_invitro

combined_plot <- ggarrange(compare_ED10_100, plot_all_invitro, invitro_pie,
                           ncol = 3,
                           vjust =1,
                           labels = "AUTO",
                           widths = c(1, 1.5, 0.75),
                           common.legend = TRUE,
                           legend = "bottom")
combined_plot

ggsave("invitro_combined_plot.jpg", combined_plot,  height =4, width =12)


