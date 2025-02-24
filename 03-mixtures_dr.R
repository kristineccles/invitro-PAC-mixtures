################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: June 14, 2024
# Note: Order to run - 3
# must have created folder in project called "mixture_fit"
#################################################

# load the data
mixture_df <- read.csv("mix_conc_df_edit.csv")
mixture_df <- subset(mixture_df, DoseActivesOnly_logM<0)

#convert to uM
mixture_df$Active_Dose_uM <- 10^(mixture_df$DoseActivesOnly_logM)*1e6
mixture_df$Total_Dose_uM <- 10^(mixture_df$Dose_LogM)*1e6

df_mix<- mixture_df %>%
  dplyr::select(-c(DoseActivesOnly_logM, Dose_LogM))%>%
  pivot_longer(cols = c(Active_Dose_uM, Total_Dose_uM), names_to = "Dose_type", values_to = "Dose_uM")
df_mix <- na.omit(df_mix)
df_mix$grouping <-paste0(df_mix$Chemicalname,"  ", df_mix$Mixture, "  ",df_mix$Dose_type)

# Set Up
#set up the data frame
NUM_PTS <- 1000
MIN_LOGX <- (-5)
MAX_LOGX <- 5

xVec <- 1:NUM_PTS
x <- as.data.frame(10 ^ (MIN_LOGX + ((MAX_LOGX-MIN_LOGX)*xVec/NUM_PTS)))
colnames(x) <- "x"

MAX_Y <- 1
yVec <- vector(mode="numeric", length = 1000 - 1)
y <- as.data.frame(MAX_Y * 1:1000 / 1000)
colnames(y) <-"y"

FIXED_C = 0

#set up parameters for CI generation
MCiter <- 1000

#set seed for reproducibility
set.seed(4565)

#----------------------------------------------------------------------
#---------- Curve Fitting for Measured Mixtures  ---------
#---------------------------------------------------------------------- 
unique_df_mix <- unique(df_mix$grouping)
for(i in unique_df_mix) {
  sub_data <- subset(df_mix, grouping == i)
  plot_chem<- drm(MaxResp~Dose_uM, data=sub_data,
                  type = "continuous",
                  fct=LL.4(fixed=c(NA, 0 , NA, NA), 
                           names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
  select<- as.data.frame(cbind(plot_chem$coefficients, chemical = unlist(i)))
  select$parameter <- rownames(select)
  write.csv(select, paste0("mixture_fit/",i,"_LL4_fit.csv"), row.names = FALSE)
  
  # Predict
  predict_i <- as.data.frame(cbind(x, chemical = i, 
                                   (predict(plot_chem, newdata = x,
                                            se.fit = TRUE, 
                                            interval = "confidence", level = 0.95))))
  write.csv(predict_i, paste0("mixture_fit/",i,"_LL4_predict.csv"), row.names = FALSE)
  
  #plot
  jpeg(paste0("mixture_fit/", i,"_rplot.jpg"))
  plot(plot_chem, type = c("bars"))
  title(i)
  dev.off()
}

mix_model<- drm(MaxResp~Dose_uM, data=df_mix, curveid = grouping,
                type = "continuous",
                fct=LL.4(fixed=c(NA, 0 , NA, NA), 
                         names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                na.action = na.omit)
plot(mix_model)

#get coefficients
mix_model_CI <- tidy(mix_model, conf.int = TRUE)
#get EC10
ED(mix_model, 10)
write.csv(mix_model_CI[,1:6], "mix_model_CI.csv", row.names = FALSE)

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
mix_model_coeff_wCI2 <- mix_model_coeff_wCI%>%
  separate(curve, c("mixture", "active", "dosecalc"), "  ")
#write to csv
write.csv(mix_model_coeff_wCI2, "mixture_fit/mix_model_coeff.csv", row.names = FALSE)

# Mixture Predictions with 95% CI
# split dataframe
df_mix_split <- split(df_mix, df_mix$grouping)
n = names(df_mix_split)
mix_pred <- sapply(n, USE.NAMES = TRUE, simplify = FALSE, FUN = function(i) {
  
  model<- drm(df_mix_split[[i]]$MaxResp~df_mix_split[[i]]$Dose_uM, data=df_mix_split[[i]], curveid = df_mix_split[[i]]$grouping,
              type = "continuous",
              lowerl = c(-Inf,0, 0),
              upperl = c(0, 100, Inf),
              fct=LL.4(fixed=c(NA, FIXED_C , NA, NA), 
                       names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
  response <- predict(model, newdata = x, interval = "confidence", level = 0.95)})

predict_df <- cbind(x, melt(mix_pred, id=c("Prediction", "Lower", "Upper")))
predict_df<- pivot_wider(predict_df, names_from = Var2, values_from = c(value))%>%
  as.data.frame()

# split based on unique character, x
predict_df <- predict_df%>%
  separate(L1, into = c("mixture", "active", "dosecalc"), sep = "\\s+", extra = "merge")
predict_df

predict_df <- predict_df%>%
  mutate(across(c(Prediction, Lower, Upper), ~ ifelse(. < 0, 0, .)))

write.csv(predict_df, "predict_df.csv")

#remove invitro before plotting
predict_df <- subset(predict_df, ! mixture == "inVitro" )

predict_df$active <- factor(predict_df$active, levels = c("Active", "All"),
                            labels = c("Active Chemicals", "All Chemicals"))

predict_df$mixture <- factor(predict_df$mixture, levels = c( "EM","ED10", "ED50"),
                            labels = c("EM", "EXP1", "EXP2"))

# Plots
measured_mix_plot <- ggplot()+
  geom_line(data = predict_df, aes(x= log10(x), y = Prediction, color = active, linetype = dosecalc), linewidth =0.75)+
  geom_ribbon(data = predict_df, aes(x=log10(x), y = Prediction, ymin=Lower, ymax=Upper, fill = active, linetype = dosecalc), alpha=0.2) +
  theme_bw()+
  facet_grid(vars(mixture))+
  scale_fill_manual(values = c("#2A788EFF", "#7AD151FF"))+
  scale_color_manual(values = c("#2A788EFF", "#7AD151FF"))+
  scale_linetype_manual(values = c("solid", "dashed"),labels = c("Active Dose (uM)","Total Dose (uM)" ))+
  labs(y="% Max MeBio Response", x= "Log10 Concentration (uM)", color = "Mixture", fill = "Mixture",
       linetype = "Dose Calculation")+
  ylim(0, 100)
measured_mix_plot
ggsave("measured_mix_plot.jpg",measured_mix_plot,  height = 6, width = 5)

################################################################
#### Parameters Only ####
#remove in vitro
mix_model_coeff_wCI2 <- subset(mix_model_coeff_wCI2, ! mixture  == "inVitro")
mix_model_coeff_wCI2$group <- paste(mix_model_coeff_wCI2$mixture, mix_model_coeff_wCI2$active)

mix_model_coeff_wCI2 <-  subset(mix_model_coeff_wCI2, dosecalc == "Active_Dose_uM")

#Change Labels
mix_model_coeff_wCI2$group<- gsub("ED10", "EXP1", mix_model_coeff_wCI2$group)
mix_model_coeff_wCI2$group<- gsub("ED50", "EXP2", mix_model_coeff_wCI2$group)
mix_model_coeff_wCI2$mixture<- gsub("ED10", "EXP1", mix_model_coeff_wCI2$mixture)
mix_model_coeff_wCI2$mixture<- gsub("ED50", "EXP2", mix_model_coeff_wCI2$mixture)
mix_model_coeff_wCI2$group <- factor(mix_model_coeff_wCI2$group, levels = rev(sort(unique(mix_model_coeff_wCI2$group))))


# make plot of just EC50
measured_mix_ED50<- ggplot(data = mix_model_coeff_wCI2, aes(x = log10(ED50), y = group, color = mixture))+
  geom_point()+
  #facet_grid(vars(dosecalc), scales = "free")+
  scale_color_viridis(discrete= TRUE)+
  geom_errorbar(data = mix_model_coeff_wCI2, aes(xmin=log10(ED50-1.96*SE_ED50), 
                                                 xmax=log10(ED50+1.96*SE_ED50), color = mixture), 
                width=.2, position=position_dodge(.9))+
  theme_bw()+
  labs(color ="Mixture")+
  ylab("Mixture")+
  xlab("Log 10 EC50 (uM)")
measured_mix_ED50

measured_mix_top<- ggplot(data = mix_model_coeff_wCI2, aes(x = (`Upper Limit`), y = group, color = mixture))+
  geom_point()+
  #facet_grid(vars(dosecalc), scales = "free")+
  scale_color_viridis(discrete= TRUE)+
  geom_errorbar(data = mix_model_coeff_wCI2, aes(xmin=(`Upper Limit`-1.96*SE_UpperLimit), 
                 xmax=(`Upper Limit`+1.96*SE_UpperLimit), color = mixture), 
                width=.2, position=position_dodge(.9))+
  theme_bw()+
  labs(color ="Mixture")+
  ylab("Mixture")+
  xlab("Top of the Curve (%)")
measured_mix_top

measured_mix_slope<- ggplot(data = mix_model_coeff_wCI2, aes(x = (Slope*-1), y = group, color = mixture))+
  geom_point()+
  scale_color_viridis(discrete= TRUE)+
  geom_errorbar(data = mix_model_coeff_wCI2, aes(xmin=((Slope-1.96*SE_slope)*-1), 
                                                 xmax=((Slope+1.96*SE_slope)*-1), color = mixture), 
                width=.2, position=position_dodge(.9))+
  theme_bw()+
  labs(color ="Mixture")+
  ylab("Mixture")+
  xlab("Slope")
measured_mix_slope


combined_plot_mix <- ggarrange(measured_mix_ED50, measured_mix_top,measured_mix_slope,
                           ncol = 1,
                           vjust =3,
                           labels = NULL,
                           common.legend = TRUE,
                           legend = "bottom")
combined_plot_mix

combined_plot_mix2 <- ggarrange(measured_mix_plot,combined_plot_mix,
                               ncol = 2,
                               labels = "AUTO",
                               common.legend = TRUE,
                               legend = "bottom")
combined_plot_mix2

ggsave("Fig2.jpg", combined_plot_mix2,  height =8, width =8)

