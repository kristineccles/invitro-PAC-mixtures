################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: November 8th, 2023

# Notes: need to run 01-individual_chemical_dr.R and 02-mixtures-modeling.R
# Order to run - 4
#################################################

#### Independent Action ####
# Calculate the array of predicted Y values for fixed X values
# For each x value, predicted y of mixture = 1 - PRODUCT(1 - predictedY(weightFactor*x))

IA_invitro100_CI <- 
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$Invitro10_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_IA_invitro_CI <- cbind(x,  melt(IA_invitro100_CI))
colnames(unlist_IA_invitro_CI) <- c("x", "y", "iteration")
unlist_IA_invitro_CI$group <- "IA"

unlist_IA_invitro_CI <- unlist_IA_invitro_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, invitro_top))))
summary(unlist_IA_invitro_CI)

IA_CI_invitro100_final <- unlist_IA_invitro_CI%>%
  group_by(x)%>%
  summarise(mean = median(y),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()

IA_CI_invitro100_final$mix_ratio <- "Invitro100"
IA_CI_invitro100_final

# Combined dataframes
IA_df_invitro <- IA_CI_invitro100_final
IA_df_invitro$Method <- "IA"

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

unlist_DA_ED10_invitro_CI <- cbind(y,  melt(DA_ED10_invitro_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, invitro_top))))
colnames(unlist_DA_ED10_invitro_CI) <- c("y", "x", "iteration")
unlist_DA_ED10_invitro_CI$group <- "CA"

DA_CI_ED10_invitro_final <- unlist_DA_ED10_invitro_CI%>%
  group_by(y)%>%
  summarise(mean = median(x),
            lower = quantile(x,0.025, na.rm =TRUE),
            upper = quantile(x, 0.975, na.rm =TRUE))%>%
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
GCA_ED10_invitro_CI <- 
  lapply(1:length(B1bootmat_list_invitro), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      (sum((B1bootmat_list_invitro[[j]]$Top)*((x*B1bootmat_list_invitro[[j]]$Invitro10_percent)/B1bootmat_list_invitro[[j]]$ED50))/ 
         (1+sum((x*B1bootmat_list_invitro[[j]]$Invitro10_percent)/B1bootmat_list_invitro[[j]]$ED50)))
    })})

unlist_GCA_ED10_invitro_CI <- cbind(x,  melt(GCA_ED10_invitro_CI))
colnames(unlist_GCA_ED10_invitro_CI) <- c("x", "y", "iteration")
unlist_GCA_ED10_invitro_CI$group <- "GCA"

ggplot(unlist_GCA_ED10_invitro_CI, aes(y=y, x=log10(x)))+
  geom_point()

GCA_CI_ED10_invitro_final <- unlist_GCA_ED10_invitro_CI%>%
  group_by(x)%>%
  summarise(mean = median(y),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
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
invitro_mix <- read.csv("mix_conc_df_edit.csv")
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

invitro_predict_df <- as.data.frame(cbind(x, mixture ="Invitro", (predict(mix_model, newdata = x,
                                                                                 interval = "confidence", level = 0.95))))

invitro_predict_df[ invitro_predict_df<0 ] <- 0

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
#---------- PART 4: Plot to Compare Predicted vs Measured ---------
#---------------------------------------------------------------------- 

#### Compare CR ####
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

#### Compare Pie plot ####
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

#### SE and Point Comparison ####
# format the data for comparison plot

# Extract model coefficients and confidence intervals
coef_summary <- summary(mix_model)$coef  # Extracts coefficient summary
confint_summary <- confint(mix_model)    # Extracts confidence intervals

# Extract values and combine them into a formatted data frame
Measured_invitro <- data.frame(
  group = "Measured",
  slope_mean = coef_summary[1,1]*-1,
  slope_lower = confint_summary[1, 2]*-1,
  slope_upper = confint_summary[1, 1]*-1,
  top = coef_summary[2, 1],
  top_lower = confint_summary[2, 1],
  top_upper = confint_summary[2, 2],
  EC50 = coef_summary[3, 1],
  EC50_lower = confint_summary[3, 1],
  EC50_upper = confint_summary[3, 2]
)

# Define function to run individual model analyses
model_params <- function(data, fixed_c_value) {
  
  # Initialize vectors for storing results
  iterVector <- vector(mode = "numeric")
  groupVector <- vector(mode = "character")
  slopeVector <- vector(mode = "numeric")
  topVector <- vector(mode = "numeric")
  ec50Vector <- vector(mode = "numeric")
  
  # Get unique iterations
  unique_itter <- unique(data$iteration)
  
  # Loop through each iteration
  for (i in unique_itter) {
    
    # Subset data for the current iteration
    sub_data <- subset(data, iteration == i)
    
    # Apply the selected model
    i_model <- drm(y ~ x, data = sub_data, fct = LL.4(fixed = c(NA, fixed_c_value, NA, NA),
                                                      names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    
    # Extract coefficients
    coefficient <- i_model$coefficients
    
    # Store results
    iterVector[i] <- i
    groupVector[i] <- unique(data$group)
    slopeVector[i] <- coefficient[1]
    topVector[i] <- coefficient[2]
    ec50Vector[i] <- coefficient[3]
    
    output_data <- data.frame(iterVector, groupVector, slopeVector, topVector, ec50Vector)
    write.csv(output_data, paste0(unique(data$group), "_param_compare.csv"), row.names = FALSE)
  }
  
  # Calculate summary statistics for the individual model
  measured_summary <-  output_data%>%
    summarize(group = unique(groupVector),
      slope_mean = quantile(slopeVector * -1, 0.5),  # Adjust for negative slope in hill models
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

# Example call with data dataset and specific FIXED_C constant
# Assuming data is defined and FIXED_C is a constant
IA_invitro <- model_params(unlist_IA_invitro_CI, FIXED_C)
DA_invitro <- model_params(unlist_DA_ED10_invitro_CI, FIXED_C)
GCA_invitro <- model_params(unlist_GCA_ED10_invitro_CI, FIXED_C)

# Plot
invitro_combined <- rbind(IA_invitro,DA_invitro, GCA_invitro, Measured_invitro)

invitro_slope<- ggplot()+
  geom_point(data = invitro_combined, aes(x = slope_mean, y = group, color = group), size = 2)+
  geom_errorbar(data = invitro_combined, aes(x = slope_mean, y = group, xmin=slope_lower, xmax=slope_upper, color = group), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Slope")
invitro_slope

invitro_top<- ggplot()+
  geom_point(data = invitro_combined, aes(x = top, y = group, color = group), size = 2)+
  geom_errorbar(data = invitro_combined, aes(x = top, y = group, xmin=top_lower, xmax=top_upper, color = group), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Top of the Curve (%)")
invitro_top

invitro_EC50<- ggplot()+
  geom_point(data = invitro_combined, aes(x = (EC50), y = group, color = group), size = 2)+
  geom_errorbar(data = invitro_combined, aes(x = (EC50), y = group, xmin=(EC50_lower), xmax=(EC50_upper), color = group), 
                width=.1, size = 1, position=position_dodge(.9))+
  theme_bw()+
  scale_color_manual(name = "Group",
                     values =  c(CA = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"))+
  labs( y = "Method", x = "Log 10 EC50 (uM)")
invitro_EC50

invitro_com_plot <- ggarrange(invitro_EC50, invitro_slope, invitro_top,
                           ncol = 3,
                           vjust =1,
                           labels = c("C", "D", "E"),
                           common.legend = TRUE,
                           legend = "bottom")
invitro_com_plot

ggsave("invitro_parameters_plot.jpg", invitro_com_plot,  height =4, width =12)

####################################################################
invitro_comb1 <- ggarrange(compare_ED10_100, invitro_pie,
                              ncol = 2,
                              vjust =1,
                              labels = c("A", "B"),
                              common.legend = FALSE,
                              legend = "bottom")

invitro_final_plot <- ggarrange(invitro_comb1, invitro_com_plot, 
                              ncol = 1,
                              vjust =1,
                              labels = NULL,
                              common.legend = TRUE,
                              legend = "bottom")
invitro_final_plot

ggsave("invitro_parameters_plot.jpg", invitro_final_plot,  height =8, width =10)


