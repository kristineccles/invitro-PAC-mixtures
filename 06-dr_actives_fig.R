################################################
# Figure 1 
# Written By: Kristin Eccles
# Date: June 3, 2022
# Must run individual_chemical_dr.R and mixture_modeling.R first
#################################################

#### Plot just the active chemicals ####
# Refit with slope of -1 for GCA
individual_model_free<- drm(MaxResp~Dose_uM, data=df_reduced, curveid= Chemicalname, 
                            type = "continuous", 
                            lowerl = c(-Inf, 0, 0),
                            upperl = c(0, 100, Inf),
                            fct=LL.4(fixed=c(NA, FIXED_C , NA, NA), 
                                     names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

# Quick Plot of curve fits
plot(individual_model_free)

# get coefficients
freeindividual_model_wCI <- tidy(individual_model_free, conf.int = TRUE)
freeindividual_model_wCI

#reorganize
freeindividual_model_coeff <- freeindividual_model_wCI%>%
  dplyr::select(term, curve, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)%>%
  as.data.frame
freeindividual_model_coeff

# add std. error
freeED50_SE <- as.data.frame(subset(freeindividual_model_wCI, term == "ED50")$std.error)
colnames(freeED50_SE) <- "SE_ED50"

freeUpperLimit_SE <- as.data.frame(subset(freeindividual_model_wCI, term == "Upper Limit")$std.error)
colnames(freeUpperLimit_SE) <- "SE_UpperLimit"

freeUpperSlope_SE <- as.data.frame(subset(freeindividual_model_wCI, term == "Slope")$std.error)
colnames(freeUpperSlope_SE) <- "SE_Slope"

freeindividual_model_coeff <- cbind(freeindividual_model_coeff,freeED50_SE, freeUpperLimit_SE,freeUpperSlope_SE )

#Bootstrap for CIs
freecoeff_final_list <- split(freeindividual_model_coeff, f=freeindividual_model_coeff$curve)
n <- names(freecoeff_final_list)

freeTop_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtruncnorm(n = MCiter, a = 0, b = Inf, mean = freecoeff_final_list[[x]]$`Upper Limit`, sd = freecoeff_final_list[[x]]$SE_UpperLimit/sqrt(2))}, 
  simplify = FALSE,USE.NAMES = TRUE)
freeTop_boot_melt <- melt(freeTop_boot)

freeED50_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtruncnorm(n = MCiter, a = 0, b = Inf, mean = freecoeff_final_list[[x]]$ED50, sd = freecoeff_final_list[[x]]$SE_ED50/sqrt(2))}, 
  simplify = FALSE,USE.NAMES = TRUE)
freeED50_boott_melt <- melt(freeED50_boot)

freeSlope_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtruncnorm(n = MCiter, a = 0, b = Inf, mean = freecoeff_final_list[[x]]$Slope, sd = -freecoeff_final_list[[x]]$SE_Slope/sqrt(2))}, 
  simplify = FALSE,USE.NAMES = TRUE)

freeSlope_boot_melt <- melt(freeSlope_boot)

freebootmat <- as.data.frame(cbind(ED50_boot_melt[,1], freeSlope_boot_melt[,1],freeTop_boot_melt))
colnames(freebootmat) <- cbind("ED50","Slope", "Top", "chemical")


# Calculate dose response for each row and combine into a long data frame
response_data <- do.call(rbind, lapply(1:nrow(freebootmat), function(i) {
  row <- freebootmat[i, ]
  resp <- row$Top/(1+ (exp(row$Slope * (log(x) - (row$ED50)))))
  data.frame(x = x$x, resp = resp, Chemical = row$chemical)
}))

colnames(response_data) <- c("x", "y", "Chemical")
summary_df <- response_data%>%
  group_by(Chemical,  x)%>%
  summarise(mean = mean(y),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()

#edit labsl
chem_labels <- c(
  "Benz(j)aceanthrylene"     = "Benz[j]aceanthrylene",
  "Dibenz(a,h)anthracene"    = "Dibenz[a,h]anthracene",
  "Benzo(a)pyrene"     = "Benzo[a]pyrene",
  "Benzo(b)fluoranthene"     = "Benzo[b]fluoranthene",
  "Indeno(1,2,3-cd)pyrene"  = "Indeno[1,2,3-cd]pyrene",
  "Benzo(k)fluoranthene"     = "Benzo[k]fluoranthene"
)

# Plot the dose response curves
individ_chems <- ggplot() +
  geom_line(data = summary_df, aes(x = log10(x), y = mean, color = Chemical), linewidth = 1) +
  geom_ribbon(data = summary_df, aes(x = log10(x), y = mean, ymin = lower, ymax = upper, fill = Chemical), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE,
                      breaks = names(chem_labels),
                      labels = unname(chem_labels)) +
  scale_fill_viridis(discrete = TRUE,
                     breaks = names(chem_labels),
                     labels = unname(chem_labels)) +
  labs(x = expression("Log"["10"] ~ " Dose (" * mu * "M)"), y = "% Max MeBio Response") +
  theme_minimal() +
  theme(legend.position = "bottom")
individ_chems
ggsave(individ_chems, file="individual_curve_dr.jpg", height = 5, width = 5)
