################################################
# Mixture Contributions
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: April 27th, 2023
#################################################

# Pie Chart 

#Mixing Ratios
mixture_df <- read.csv("mixing_fractions_pie.csv")
mixture_df_reduced <- mixture_df[complete.cases(mixture_df$InvitroEC10), ]

mixture_df_reduced$EM_percent_100 <- mixture_df_reduced$EM_mM/sum(mixture_df_reduced$EM_mM)
mixture_df_reduced$ED10_percent_100 <- mixture_df_reduced$EP10_mM/sum(mixture_df_reduced$EP10_mM)
mixture_df_reduced$ED50_percent_100 <- mixture_df_reduced$EP50_mM/sum(mixture_df_reduced$EP50_mM)

library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)

#all chemicals
em_plt <- ggplot(mixture_df, aes(x="", y=EM_percent, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
em_plt

ec10_plt <- ggplot(mixture_df, aes(x="", y=ED10_percent, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
ec10_plt

ec50_plt <- ggplot(mixture_df, aes(x="", y=ED50_percent, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
ec50_plt

##################################################
#all chemicals but only actives
em_plt6 <- ggplot(mixture_df_reduced, aes(x="", y=EM_percent_100, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
em_plt6

ec10_plt6 <- ggplot(mixture_df_reduced, aes(x="", y=ED10_percent_100, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
ec10_plt6

ec50_plt6 <- ggplot(mixture_df_reduced, aes(x="", y=ED50_percent_100, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
ec50_plt6

#actives only
NAem_plt <- ggplot(mixture_df_reduced, aes(x="", y=NA_EM, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
NAem_plt

NAec10_plt <- ggplot(mixture_df_reduced, aes(x="", y=NA_EP10, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
NAec10_plt

NAec50_plt <- ggplot(mixture_df_reduced, aes(x="", y=NA_EP50, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
NAec50_plt

NA_invitro_plt <- ggplot(mixture_df, aes(x="", y=Invitro10_percent, fill=Chemical))+
  geom_bar(stat = "identity")+
  # Make it circular
  coord_polar("y", start=0)+
  theme_minimal()+
  scale_fill_viridis(discrete = TRUE)+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
NA_invitro_plt

combinedh <- ggarrange(em_plt, ec10_plt, ec50_plt, 
                       em_plt6, ec10_plt6, ec50_plt6,
                       NAem_plt, NAec10_plt, NAec50_plt,
          labels="AUTO", 
          ncol = 3, nrow = 3,
          common.legend = TRUE,
          legend = "bottom")
combinedh

ggsave("pie_plot_combined_h.jpg", 
       combinedh, 
       dpi = 300, 
       width = 8,
       height = 8,
       units = "in")

