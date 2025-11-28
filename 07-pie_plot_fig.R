################################################
# Mixture Contributions
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: April 27th, 2023
# need to run 06-dr_actives_fig for individ_chems
#################################################

chem_labels <- c(
  "Benz(j)aceanthrylene"     = "Benz[j]aceanthrylene",
  "Dibenz(a,h)anthracene"    = "Dibenz[a,h]anthracene",
  "Benzo(a)pyrene"     = "Benzo[a]pyrene",
  "Benzo(b)fluoranthene"     = "Benzo[b]fluoranthene",
  "Indeno(1,2,3-cd)pyrene"  = "Indeno[1,2,3-cd]pyrene",
  "Benzo(k)fluoranthene"     = "Benzo[k]fluoranthene"
)


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

# Helper: a viridis fill scale that relabels legend entries,
# leaving any unmapped chemicals unchanged.
fill_scale_relabel <- scale_fill_viridis(
  discrete = TRUE,
  labels = function(l) dplyr::recode(l, !!!chem_labels, .default = l)
)

# ---- Theme helper to avoid repetition ----
pie_theme <- theme_minimal() +
  theme(
    plot.title   = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks   = element_blank()
  )

# ============================
# All chemicals (three pies)
# ============================
em_plt <- ggplot(mixture_df, aes(x = "", y = EM_percent, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

ec10_plt <- ggplot(mixture_df, aes(x = "", y = ED10_percent, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

ec50_plt <- ggplot(mixture_df, aes(x = "", y = ED50_percent, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

# ==========================================
# All chemicals but only actives (normalized)
# ==========================================
em_plt6 <- ggplot(mixture_df_reduced, aes(x = "", y = EM_percent_100, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

ec10_plt6 <- ggplot(mixture_df_reduced, aes(x = "", y = ED10_percent_100, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

ec50_plt6 <- ggplot(mixture_df_reduced, aes(x = "", y = ED50_percent_100, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

# =================
# Actives only pies
# =================
NAem_plt <- ggplot(mixture_df_reduced, aes(x = "", y = NA_EM, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

NAec10_plt <- ggplot(mixture_df_reduced, aes(x = "", y = NA_EP10, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

NAec50_plt <- ggplot(mixture_df_reduced, aes(x = "", y = NA_EP50, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

NA_invitro_plt <- ggplot(mixture_df_reduced, aes(x = "", y = Invitro10_percent, fill = Chemical)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  pie_theme +
  fill_scale_relabel

# 3x3 grid of pies
combinedh <- ggarrange(
  em_plt, ec10_plt, ec50_plt,
  em_plt6, ec10_plt6, ec50_plt6,
  NAem_plt, NAec10_plt, NAec50_plt,
  labels = c("B", "C", "D", "E", "F", "G", "H", "I", "J"),
  ncol = 3, nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)

# ---- Bring in individual DR panel from 06 (assumes `individ_chems` exists) ----
fig1 <- ggarrange(
  individ_chems, combinedh,
  labels = c("A", " "),
  ncol = 1, nrow = 2,
  heights = c(0.75, 1),
  common.legend = FALSE
)
fig1
# ---- Save ----
ggsave("Fig1.jpg", fig1, dpi = 300, width = 6.5, height = 10, units = "in")
