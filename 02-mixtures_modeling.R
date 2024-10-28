################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: April 17th, 2024
# Note: Order to run -2
# Must run individual_chemical_dr.R first
#################################################

# load libraries
library(drc)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(sjPlot)
library(cowplot)
library(data.table)
library(reshape2)
library(ggpubr)
library(broom)
library(viridis)
library(MCMCglmm)
library(dplyr)

# Set Up
NUM_PTS <- 1000
MIN_LOGX <- (-5)
MAX_LOGX <- 5   

MAX_Y <- 1
yVec <- vector(mode="numeric", length = 1000 - 1)
y <- as.data.frame(MAX_Y * 1:1000 / 1000)
colnames(y) <-"y"

#set up parameters for CI generation
MCiter <- 10000

#set seed for reproducibility
set.seed(8789)


#----------------------------------------------------------------------
#---------- Refit Curve for DA and IA  ---------
#----------------------------------------------------------------------
#need to have the same top and bottom
#only fit for reduced list 
df$include <- df$Chemicalname %in%  individual_coeff_final$curve 
df_reduced <- subset(df, include == TRUE)

#### Fit Parameters with Selected Model Hill Model ####
# fit for mixtures modeling
individual_model2<- drm(MaxResp~Dose_uM, data=df_reduced, curveid = Chemicalname,
                       type = "continuous",
                       lowerl = c(-Inf, 0),
                       upperl = c(0, Inf),
                       fct=LL.4(fixed=c(NA, FIXED_C , FIXED_D, NA), 
                                names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(individual_model2)
# Quick Plot of curve fits
plot(individual_model2)

# Extract coefficients
individual_model_wCI2 <- tidy(individual_model2, conf.int = TRUE)
individual_model_coeff2 <- individual_model_wCI2 %>%
  select(term, curve, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)

# add SE
slope_SE2 <- as.data.frame(subset(individual_model_wCI2, term == "Slope")$std.error)
colnames(slope_SE2) <- "SE_slope"
ED50_SE2 <- as.data.frame(subset(individual_model_wCI2, term == "ED50")$std.error)
colnames(ED50_SE2) <- "SE_ED50"
UpperLimit_SE2 <- as.data.frame(subset(individual_model_wCI2, term == "Upper Limit")$std.error)
colnames(UpperLimit_SE2) <- "SE_UpperLimit"
individual_model_coeff2 <- cbind(individual_model_coeff2,slope_SE2, ED50_SE2)
#write.csv(individual_model_coeff, "individual_model_coeff.csv", row.names = FALSE)

# Bootstrap for CI
coeff_final_list2 <- split(individual_model_coeff2, f = individual_model_coeff2$curve)
slope_boot <- sapply(names(coeff_final_list2), function(x) rtruncnorm(n = MCiter, a = -Inf, b = 0, mean = coeff_final_list2[[x]]$Slope, sd = coeff_final_list2[[x]]$SE_slope/sqrt(2)))
ED50_boot <- sapply(names(coeff_final_list2), function(x) rtruncnorm(n = MCiter, a = 0, b = Inf, mean = coeff_final_list2[[x]]$ED50, sd = coeff_final_list2[[x]]$SE_ED50/sqrt(2)))

# Add Mixing Ratios
MIX_FRACTIONS <- na.omit(read.csv("mixing_fractions.csv"))
individual_model_coeff2 <- left_join(individual_model_coeff2, MIX_FRACTIONS, by = c("curve" = "Chemical"))
individual_model_coeff2 <- na.omit(individual_model_coeff2)

#### Bootstrap for CI ####

coeff_final_list2 <- split(individual_model_coeff2, f=individual_model_coeff2$curve)
n <- names(coeff_final_list2)

slope_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtruncnorm(n = MCiter, a = -Inf, b = 0, mean = coeff_final_list2[[x]]$Slope, sd = coeff_final_list2[[x]]$SE_slope/sqrt(2))}, 
  simplify = FALSE,USE.NAMES = TRUE)
slope_boot_melt <- melt(slope_boot)

ED50_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtruncnorm(n = MCiter, a = 0, b = Inf, mean = coeff_final_list2[[x]]$ED50, sd = coeff_final_list2[[x]]$SE_ED50/sqrt(2))}, 
  simplify = FALSE,USE.NAMES = TRUE)

ED50_boot_melt <- melt(ED50_boot)

bootmat <- as.data.frame(cbind(ED50_boot_melt[,1], slope_boot_melt))
colnames(bootmat) <- cbind("ED50","slope", "chemical")
#add iteration number for each chemical
bootmat$itter <- 1:MCiter

#Mixing Ratios
MIX_FRACTIONS <- na.omit(read.csv("mixing_fractions.csv"))

#check
sum(individual_model_coeff2$EM_percent)
sum(individual_model_coeff2$ED10_percent)
sum(individual_model_coeff2$ED50_percent)

#adjust percentage so it adds to 100%
individual_model_coeff2$EM_percent_100 <- individual_model_coeff2$EM_mM/sum(individual_model_coeff2$EM_mM)
sum(individual_model_coeff2$EM_percent_100)
individual_model_coeff2$ED10_percent_100 <- individual_model_coeff2$EP10_mM/sum(individual_model_coeff2$EP10_mM)
sum(individual_model_coeff2$ED10_percent_100)
individual_model_coeff2$ED50_percent_100 <- individual_model_coeff2$EP50_mM/sum(individual_model_coeff2$EP50_mM)
sum(individual_model_coeff2$ED50_percent_100)

# get weighted top of curves - must use model fit w/o  a top
coeff_fractions <- left_join(individual_coeff_final, MIX_FRACTIONS, by= c("curve" = "Chemical"), keep=FALSE )
coeff_fractions$EM_percent_100 <- coeff_fractions$EM_mM/sum(coeff_fractions$EM_mM)
coeff_fractions$ED10_percent_100 <- coeff_fractions$EP10_mM/sum(coeff_fractions$EP10_mM)
coeff_fractions$ED50_percent_100 <- coeff_fractions$EP50_mM/sum(coeff_fractions$EP50_mM)

invitro_top <- weighted.mean(coeff_fractions$`Upper Limit`, coeff_fractions$Invitro10_percent)
EM_top <- weighted.mean(coeff_fractions$`Upper Limit`, coeff_fractions$EM_percent_100)
ED10_top <- weighted.mean(coeff_fractions$`Upper Limit`, coeff_fractions$ED10_percent_100)
ED50_top <- weighted.mean(coeff_fractions$`Upper Limit`, coeff_fractions$ED50_percent_100)

bootmat <- left_join(bootmat, individual_model_coeff2[,c("curve", 
                                                         "EM_percent_100", "EM_percent",
                                                        "ED10_percent_100", "ED10_percent", 
                                                        "ED50_percent_100", "ED50_percent", "Invitro10_percent")],
                     by= c("chemical" = "curve"), keep=FALSE )
#make new dataframe for CI generation
bootmat_list <- split(bootmat, f=bootmat$itter)

# Refit with slope of -1 for GCA
individual_model_b1<- drm(MaxResp~Dose_uM, data=df_reduced, curveid= Chemicalname, 
                          type = "continuous", 
                          lowerl = c(0, 0),
                          upperl = c(100, Inf),
                          fct=LL.4(fixed=c(FIXED_B, FIXED_C , NA, NA), 
                                   names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

# Quick Plot of curve fits
plot(individual_model_b1)

# get coefficients
B1individual_model_wCI <- tidy(individual_model_b1, conf.int = TRUE)
B1individual_model_wCI

#reorganize
B1individual_model_coeff <- B1individual_model_wCI%>%
  dplyr::select(term, curve, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)%>%
  as.data.frame
B1individual_model_coeff

# add std. error
B1ED50_SE <- as.data.frame(subset(B1individual_model_wCI, term == "ED50")$std.error)
colnames(B1ED50_SE) <- "SE_ED50"

B1UpperLimit_SE <- as.data.frame(subset(B1individual_model_wCI, term == "Upper Limit")$std.error)
colnames(B1UpperLimit_SE) <- "SE_UpperLimit"

B1individual_model_coeff <- cbind(B1individual_model_coeff,B1ED50_SE, B1UpperLimit_SE)

#Bootstrap for CIs
B1coeff_final_list <- split(B1individual_model_coeff, f=B1individual_model_coeff$curve)
n <- names(B1coeff_final_list)

B1Top_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtnorm(n = MCiter, lower = 0, upper = Inf, mean = B1coeff_final_list[[x]]$`Upper Limit`, sd = B1coeff_final_list[[x]]$SE_UpperLimit/sqrt(2))}, 
  simplify = FALSE,USE.NAMES = TRUE)
B1Top_boot_melt <- melt(B1Top_boot)

B1ED50_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtnorm(n = MCiter, lower = 0, upper = Inf, mean =  B1coeff_final_list[[x]]$ED50, sd = B1coeff_final_list[[x]]$SE_ED50/sqrt(2))}, 
  simplify = FALSE,USE.NAMES = TRUE)

B1ED50_boot_melt <- melt(B1ED50_boot)

B1bootmat <- as.data.frame(cbind(ED50_boot_melt[,1], B1Top_boot_melt))
colnames(B1bootmat) <- cbind("ED50","Top", "chemical")
#add itteration number for each chemical
B1bootmat$itter <- 1:MCiter

#add mixing ratios
B1bootmat <- left_join(B1bootmat, individual_model_coeff2[,c("curve", "EM_percent","EM_percent_100", 
                                                             "ED10_percent","ED10_percent_100",
                                                             "ED50_percent", "ED50_percent_100" )],
                       by= c("chemical" = "curve"), keep=FALSE )
#make new dataframe for CI generation
B1bootmat_list <- split(B1bootmat, f=B1bootmat$itter)

#----------------------------------------------------------------------
#---------- Predict Mixture Responses  ---------
#----------------------------------------------------------------------

#### Unadjusted Percent ####
#### Unadjusted Independent Action ####
IA_EM_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$EM_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_IA_EM_CI <- cbind(x,  melt(IA_EM_CI))
colnames(unlist_IA_EM_CI) <- c("x", "y", "itteration")
unlist_IA_EM_CI$group <- "unadj"
unlist_IA_EM_CI$mix <- "EM"
unlist_IA_EM_CI$method <- "IA"
unlist_IA_EM_CI<- unlist_IA_EM_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EM_top))))
write.csv(unlist_IA_EM_CI, "mix_pred_boot/unlist_EM_CI_unadj.csv", row.names = FALSE)

IA_CI_EM_final <- unlist_IA_EM_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()

head(IA_CI_EM_final)
IA_CI_EM_final$mix_ratio <- "EM"

IA_ED10_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$ED10_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_IA_ED10_CI <- cbind(x,  melt(IA_ED10_CI))
colnames(unlist_IA_ED10_CI) <- c("x", "y", "itteration")
unlist_IA_ED10_CI$group <- "unadj"
unlist_IA_ED10_CI$mix <- "ED10"
unlist_IA_ED10_CI$method <- "IA"
unlist_IA_ED10_CI <- unlist_IA_ED10_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED10_top))))
write.csv(unlist_IA_ED10_CI, "mix_pred_boot/unlist_ED10_CI_unadj.csv", row.names = FALSE)

IA_CI_ED10_final <- unlist_IA_ED10_CI%>%
  group_by(x)%>%
  dplyr::summarize(mean = quantile(y,0.5, na.rm =TRUE),
            y_lower = quantile(y, 0.025, rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
head(IA_CI_ED10_final)
IA_CI_ED10_final$mix_ratio <- "ED10"

IA_ED50_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$ED50_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_IA_ED50_CI <- cbind(x,  melt(IA_ED50_CI))
colnames(unlist_IA_ED50_CI) <- c("x", "y", "itteration")
unlist_IA_ED50_CI$group <- "unadj"
unlist_IA_ED50_CI$mix <- "ED50"
unlist_IA_ED50_CI$method <- "IA"
unlist_IA_ED50_CI <- unlist_IA_ED50_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED50_top))))
write.csv(unlist_IA_ED50_CI, "mix_pred_boot/unlist_ED50_CI_unadj.csv", row.names = FALSE)

IA_CI_ED50_final <- unlist_IA_ED50_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y, 0.5, na.rm =TRUE),
            y_lower = quantile(y, 0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
head(IA_CI_ED50_final)
IA_CI_ED50_final$mix_ratio <- "ED50"

# Combined dataframes
IA_df <- rbind(IA_CI_EM_final, IA_CI_ED10_final, IA_CI_ED50_final)
IA_df$Method <- "IA"

IA_df$category <- "Unadjusted"
write.csv(IA_df, "IA_df_noadj.csv", row.names = FALSE)

#plot
p1 <- ggplot()+
  geom_line(data = IA_df, aes(x= log10(x), y = mean, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=IA_df, aes(x=log10(x), y=mean, ymin=y_lower, ymax=y_upper, fill = mix_ratio), alpha=0.2) +
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

#### Unadjusted Concentration Addition ####
# x = (e^(e *b)* (d/y - 1))^(1/b)

DA_EM_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$EM_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope)*log((1-y)/y))))
    })})

unlist_DA_EM_CI <- cbind(y,  melt(DA_EM_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EM_top))))
  
colnames(unlist_DA_EM_CI) <- c("y", "x", "itteration")
unlist_DA_EM_CI$group <- "unadj"
unlist_DA_EM_CI$mix <- "EM"
unlist_DA_EM_CI$method <- "DA"
write.csv(unlist_DA_EM_CI, "mix_pred_boot/unlist_DA_EM_CI_unadj.csv", row.names = FALSE)

DA_CI_EM_final <- unlist_DA_EM_CI%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EM_top))))%>%
  as.data.frame()
DA_CI_EM_final$mix_ratio <- "EM"

DA_ED10_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$ED10_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope)*log((1-y)/y))))
    })})

unlist_DA_ED10_CI <- cbind(y,  melt(DA_ED10_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED10_top))))
colnames(unlist_DA_ED10_CI) <- c("y", "x", "itteration")
unlist_DA_ED10_CI$group <- "unadj"
unlist_DA_ED10_CI$mix <- "ED10"
unlist_DA_ED10_CI$method <- "DA"
write.csv(unlist_DA_ED10_CI, "mix_pred_boot/unlist_DA_ED10_CI_unadj.csv", row.names = FALSE)

DA_CI_ED10_final <- unlist_DA_ED10_CI%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED10_top))))%>%
  as.data.frame()
DA_CI_ED10_final$mix_ratio <- "ED10"

DA_ED50_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$ED50_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope )*log((1-y)/y))))
    })})

unlist_DA_ED50_CI <- cbind(y,  melt(DA_ED50_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED50_top))))
colnames(unlist_DA_ED50_CI) <- c("y", "x", "itteration")
unlist_DA_ED50_CI$group <- "unadj"
unlist_DA_ED50_CI$mix <- "ED50"
unlist_DA_ED50_CI$method <- "DA"
write.csv(unlist_DA_ED50_CI, "mix_pred_boot/unlist_DA_ED50_CI_unadj.csv", row.names = FALSE)

DA_CI_ED50_final <- unlist_DA_ED50_CI%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED50_top))))%>%
  as.data.frame()
DA_CI_ED50_final$mix_ratio <- "ED50"

# Combined dataframes
DA_df <- rbind(DA_CI_EM_final, DA_CI_ED10_final, DA_CI_ED50_final)
DA_df$Method <- "CA"
DA_df$category <- "Unadjusted"
write.csv(DA_df, "DA_df_noadj.csv", row.names = FALSE)

#plot
p2noadj <- ggplot()+
  geom_line(data = DA_df, aes(x= log10(mean), y = y, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=DA_df, aes(x=log10(mean), y = y, xmin=log10(x_lower), xmax=log10(x_upper), fill = mix_ratio), alpha=0.2) +
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
p2noadj

#### Unadjusted Generalized Concentration Addition ####

# GCA Predictions
#Emax = (MAX1(C1/ED501) + MAX2(C2/ED502))/(1+(C1/ED501)+(C2/EC2))
#Emix = MAX1(C1/ED501) + MAX2(C2/ED502) + MAX1(C1/ED501)/1+(C1/ED501)+(C2/ED502) +(C3/ED503)

GCA_EM_CI <- 
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      (sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$EM_percent)/B1bootmat_list[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list[[j]]$EM_percent)/B1bootmat_list[[j]]$ED50)))
    })})

unlist_GCA_EM_CI <- cbind(x,  melt(GCA_EM_CI))
colnames(unlist_GCA_EM_CI) <- c("x", "y", "itteration")
unlist_GCA_EM_CI$group <- "unadj"
unlist_GCA_EM_CI$mix <- "EM"
unlist_GCA_EM_CI$method <- "GCA"
write.csv(unlist_GCA_EM_CI, "mix_pred_boot/unlist_GCA_EM_CI_unadj.csv", row.names = FALSE)

GCA_CI_EM_final <- unlist_GCA_EM_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_EM_final$mix_ratio <- "EM"

GCA_ED10_CI <- 
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      (sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$ED10_percent)/B1bootmat_list[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list[[j]]$ED10_percent)/B1bootmat_list[[j]]$ED50)))
    })})

unlist_GCA_ED10_CI <- cbind(x,  melt(GCA_ED10_CI))
colnames(unlist_GCA_ED10_CI) <- c("x", "y", "itteration")
unlist_GCA_ED10_CI$group <- "unadj"
unlist_GCA_ED10_CI$mix <- "ED10"
unlist_GCA_ED10_CI$method <- "GCA"
write.csv(unlist_GCA_ED10_CI, "mix_pred_boot/unlist_GCA_ED10_CI_unadj.csv", row.names = FALSE)

GCA_CI_ED10_final <- unlist_GCA_ED10_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_ED10_final$mix_ratio <- "ED10"

GCA_ED50_CI <- 
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      (sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$ED50_percent)/B1bootmat_list[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list[[j]]$ED50_percent)/B1bootmat_list[[j]]$ED50)))
    })})

unlist_GCA_ED50_CI <- cbind(x,  melt(GCA_ED50_CI))
colnames(unlist_GCA_ED50_CI) <- c("x", "y", "itteration")
unlist_GCA_ED50_CI$group <- "unadj"
unlist_GCA_ED50_CI$mix <- "ED50"
unlist_GCA_ED50_CI$method <- "GCA"
write.csv(unlist_GCA_ED50_CI, "mix_pred_boot/unlist_GCA_ED50_CI_unadj.csv", row.names = FALSE)

GCA_CI_ED50_final <- unlist_GCA_ED50_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_ED50_final$mix_ratio <- "ED50"

# Combined dataframes
GCA_df <- rbind(GCA_CI_EM_final, GCA_CI_ED10_final, GCA_CI_ED50_final)
GCA_df$Method <- "GCA"
GCA_df$category <- "Unadjusted"
GCA_df$y <- GCA_df$mean
GCA_df$y_lower <- GCA_df$lower
GCA_df$y_upper <- GCA_df$upper
write.csv(GCA_df, "GCA_df_noadj.csv", row.names = FALSE)

#plot
p3noadj <- ggplot()+
  geom_line(data = GCA_df, aes(x= log10(x), y = y, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=GCA_df, aes(x=log10(x), y=y, ymin=y_lower, ymax=y_upper, fill = mix_ratio), alpha=0.2) +
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
p3noadj

#### Adjusted Percent 100 ####
#### Adj Independent Action ####
IA_EM_CI_100 <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$EM_percent_100) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_EM_CI_100 <- cbind(x,  melt(IA_EM_CI_100))
colnames(unlist_EM_CI_100) <- c("x", "y", "itteration")
unlist_EM_CI_100$group <- "adj"
unlist_EM_CI_100$mix <- "EM"
unlist_EM_CI_100$method <- "IA"
write.csv(unlist_EM_CI_100, "mix_pred_boot/unlist_EM_CI_adj.csv", row.names = FALSE)

IA_CI_EM_final_100 <- unlist_EM_CI_100%>%
  group_by(x)%>%
  summarise(mean = quantile(y, 0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  mutate(across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, EM_top))))%>%
  as.data.frame()
IA_CI_EM_final_100$mix_ratio <- "EM"

IA_ED10_CI_100 <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$ED10_percent_100) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_ED10_CI_100 <- cbind(x,  melt(IA_ED10_CI_100))
colnames(unlist_ED10_CI_100) <- c("x", "y", "itteration")
unlist_ED10_CI_100$group <- "adj"
unlist_ED10_CI_100$mix <- "ED10"
unlist_ED10_CI_100$method <- "IA"
write.csv(unlist_ED10_CI_100, "mix_pred_boot/unlist_ED10_CI_adj.csv", row.names = FALSE)

IA_CI_ED10_final_100 <- unlist_ED10_CI_100%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  mutate(across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, ED10_top))))%>%
  as.data.frame()
IA_CI_ED10_final_100$mix_ratio <- "ED10"

IA_ED50_CI_100 <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$ED50_percent_100) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_ED50_CI_100 <- cbind(x,  melt(IA_ED50_CI_100))
colnames(unlist_ED50_CI_100) <- c("x", "y", "itteration")
unlist_ED50_CI_100$group <- "adj"
unlist_ED50_CI_100$mix <- "ED50"
unlist_ED50_CI_100$method <- "IA"
write.csv(unlist_ED50_CI_100, "mix_pred_boot/unlist_ED50_CI_adj.csv", row.names = FALSE)

IA_CI_ED50_final_100 <- unlist_ED50_CI_100%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  mutate(across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, ED50_top))))%>%
  as.data.frame()
IA_CI_ED50_final_100$mix_ratio <- "ED50"

# Combined data frames
IA_df_100 <- rbind(IA_CI_EM_final_100, IA_CI_ED10_final_100, IA_CI_ED50_final_100)
IA_df_100$Method <- "IA"
IA_df_100$category <- "Adjusted"
write.csv(IA_df_100, "IA_df_adj.csv", row.names = FALSE)

#plot
p1adj <- ggplot()+
  geom_line(data = IA_df_100, aes(x= log10(x), y = mean, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=IA_df_100, aes(x=log10(x), y= mean, ymin=y_lower, ymax=y_upper, fill = mix_ratio), alpha=0.2) +
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
p1adj

#### Adj Concentration Addition ####
# x = (e^(e *b)* (d/y - 1))^(1/b)

DA_EM_CI_100 <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$EM_percent_100*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope)*log((1-y)/y))))
    })})

unlist_DA_EM_CI_100 <- cbind(y,  melt(DA_EM_CI_100))
colnames(unlist_DA_EM_CI_100) <- c("y", "x", "itteration")
unlist_DA_EM_CI_100$group <- "adj"
unlist_DA_EM_CI_100$mix <- "EM"
unlist_DA_EM_CI_100$method <- "DA"
write.csv(unlist_DA_EM_CI_100, "mix_pred_boot/unlist_DA_EM_CI_adj.csv", row.names = FALSE)

DA_CI_EM_final_100 <- unlist_DA_EM_CI_100%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EM_top))))%>%
  as.data.frame()
DA_CI_EM_final_100$mix_ratio <- "EM"

DA_ED10_CI_100 <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$ED10_percent_100*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope)*log((1-y)/y))))
    })})

unlist_DA_ED10_CI_100 <- cbind(y,  melt(DA_ED10_CI_100))
colnames(unlist_DA_ED10_CI_100) <- c("y", "x", "itteration")
unlist_DA_ED10_CI_100$group <- "adj"
unlist_DA_ED10_CI_100$mix <- "ED10"
unlist_DA_ED10_CI_100$method <- "DA"
write.csv(unlist_DA_ED10_CI_100, "mix_pred_boot/unlist_DA_ED10_CI_adj.csv", row.names = FALSE)

DA_CI_ED10_final_100 <- unlist_DA_ED10_CI_100%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED10_top))))%>%
  as.data.frame()
DA_CI_ED10_final_100$mix_ratio <- "ED10"

DA_ED50_CI_100 <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$ED50_percent_100*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope )*log((1-y)/y))))
    })})

unlist_DA_ED50_CI_100 <- cbind(y,  melt(DA_ED50_CI_100))
colnames(unlist_DA_ED50_CI_100) <- c("y", "x", "itteration")
unlist_DA_ED50_CI_100$group <- "adj"
unlist_DA_ED50_CI_100$mix <- "ED50"
unlist_DA_ED50_CI_100$method <- "DA"
write.csv(unlist_DA_ED50_CI_100, "mix_pred_boot/unlist_DA_ED50_CI_adj.csv", row.names = FALSE)

DA_CI_ED50_final_100 <- unlist_DA_ED50_CI_100%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, ED50_top))))%>%
  as.data.frame()
DA_CI_ED50_final_100$mix_ratio <- "ED50"

# Combined dataframes
DA_df_100 <- rbind(DA_CI_EM_final_100, DA_CI_ED10_final_100, DA_CI_ED50_final_100)
DA_df_100$Method <- "CA"
DA_df_100$category <- "Adjusted"
write.csv(DA_df_100, "DA_df_adj.csv", row.names = FALSE)

#plot
p2adj <- ggplot()+
  geom_line(data = DA_df_100, aes(x= log10(mean), y = y, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data = DA_df_100, aes(x=log10(mean), y = y, xmin=log10(x_lower), xmax=log10(x_upper), fill = mix_ratio), alpha=0.2) +
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
p2adj

#### Adj Generalized Concentration Addition ####
# GCA Predictions

GCA_EM_CI_100 <- 
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$EM_percent_100)/B1bootmat_list[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list[[j]]$EM_percent_100)/B1bootmat_list[[j]]$ED50))
    })})

unlist_GCA_EM_CI_100 <- cbind(x,  melt(GCA_EM_CI_100))
colnames(unlist_GCA_EM_CI_100) <- c("x", "y", "itteration")
unlist_GCA_EM_CI_100$group <- "adj"
unlist_GCA_EM_CI_100$mix <- "EM"
unlist_GCA_EM_CI_100$method <- "GCA"
write.csv(unlist_GCA_EM_CI_100, "mix_pred_boot/unlist_GCA_EM_CI_adj.csv", row.names = FALSE)


GCA_CI_EM_final_100 <- unlist_GCA_EM_CI_100%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_EM_final_100$mix_ratio <- "EM"

GCA_ED10_CI_100 <- 
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$ED10_percent_100)/B1bootmat_list[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list[[j]]$ED10_percent_100)/B1bootmat_list[[j]]$ED50))
    })})

unlist_GCA_ED10_CI_100 <- cbind(x,  melt(GCA_ED10_CI_100))
colnames(unlist_GCA_ED10_CI_100) <- c("x", "y", "itteration")
unlist_GCA_ED10_CI_100$group <- "adj"
unlist_GCA_ED10_CI_100$mix <- "ED10"
unlist_GCA_ED10_CI_100$method <- "GCA"
write.csv(unlist_GCA_ED10_CI_100, "mix_pred_boot/unlist_GCA_ED10_CI_adj.csv", row.names = FALSE)

GCA_CI_ED10_final_100 <- unlist_GCA_ED10_CI_100%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_ED10_final_100$mix_ratio <- "ED10"

GCA_ED50_CI_100 <- 
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$ED50_percent_100)/B1bootmat_list[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list[[j]]$ED50_percent_100)/B1bootmat_list[[j]]$ED50))
    })})

unlist_GCA_ED50_CI_100 <- cbind(x,  melt(GCA_ED50_CI_100))
colnames(unlist_GCA_ED50_CI_100) <- c("x", "y", "itteration")
unlist_GCA_ED50_CI_100$group <- "adj"
unlist_GCA_ED50_CI_100$mix <- "ED50"
unlist_GCA_ED50_CI_100$method <- "GCA"
write.csv(unlist_GCA_ED50_CI_100, "mix_pred_boot/unlist_GCA_ED50_CI_adj.csv", row.names = FALSE)

GCA_CI_ED50_final_100 <- unlist_GCA_ED50_CI_100%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_ED50_final_100$mix_ratio <- "ED50"

# Combined dataframes
GCA_df_100 <- rbind(GCA_CI_EM_final_100, GCA_CI_ED10_final_100, GCA_CI_ED50_final_100)
GCA_df_100$Method <- "GCA"

GCA_df_100$category <- "Adjusted"
GCA_df_100$y <- GCA_df_100$mean
GCA_df_100$y_lower <- GCA_df_100$lower
GCA_df_100$y_upper <- GCA_df_100$upper

write.csv(GCA_df_100, "GCA_df_adj.csv", row.names = FALSE)

#plot
p3adj <- ggplot()+
  geom_line(data = GCA_df_100, aes(x= log10(x), y = y, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=GCA_df_100, aes(x=log10(x), y=y, ymin=y_lower, ymax=y_upper, fill = mix_ratio), alpha=0.2) +
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
p3adj

