################################################################################
### Decision for dialysis versus conservative management
### PART 10 - Sensitivity analysis for unmeasured confounding
################################################################################

# set-up
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)
setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/")
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Results/"

# load libraries
library(tidyverse)

################################################################################
### Shared parameters
################################################################################
rr_confounded   <- 0.57
rr_unconfounded <- 1
n               <- 200

################################################################################
### 3D plot
################################################################################
RR_CD       <- seq(1, 10, length.out = n)
P_C1        <- seq(0, 50, length.out = n)
P_C0_matrix <- outer(
  RR_CD, P_C1 / 100,
  function(RR_CD, P_C1) {
    P_C0 <- (P_C1 * (RR_CD - 1) + 1 - rr_confounded) / (rr_confounded * (RR_CD - 1))
    P_C0[P_C0 < 0 | P_C0 > 1] <- NA
    return(P_C0 * 100)
  }
)

# Color mapping: blue -> purple -> green4 along the P_C0 (z) axis
nbcol         <- n
color_palette <- colorRampPalette(c("blue", "purple", "green4"))(nbcol)

z_facet <- (P_C0_matrix[-1, -1] +
              P_C0_matrix[-1, -ncol(P_C0_matrix)] +
              P_C0_matrix[-nrow(P_C0_matrix), -1] +
              P_C0_matrix[-nrow(P_C0_matrix), -ncol(P_C0_matrix)]) / 4
facet_col <- color_palette[cut(z_facet, nbcol)]

png(
  file.path(results_path, paste0("Supplemental/Figure_M_unmeasured_confounding_3D_Plot.png")),
  width  = 1200,
  height = 1000,
  res    = 150
)
par(mar = c(2, 2, 4, 2))
persp(
  RR_CD, P_C1, P_C0_matrix,
  theta     = 310,
  phi       = 20,
  expand    = 0.8,
  col       = facet_col,
  lwd       = 0.2,
  ticktype  = "detailed",
  border    = "black",
  xlab      = "Confounder-outcome strength (RR_CD)",
  ylab      = "Prevalence dialysis (%, P_C1)",
  zlab      = "Prevalence conservative management (%, P_C0)",
  main      = "A. Varying RR_CD, P_C0, and P_C1",
  cex.main  = 1.2,
  cex.axis  = 0.8,
  cex.lab   = 1
)
dev.off()

################################################################################
### 2D plot: fixed P_C0 — confounder strength (RR_CD) vs required P_C1
################################################################################
P_C0            <- c(20, 40, 60, 80)
RR_CD           <- seq(1, 10, length.out = n)

P_C0_labels     <- paste0("P_C0 = ", P_C0)
my_colours      <- setNames(c("#1D9E75", "#378ADD", "#D85A30", "#8B0000"), P_C0_labels)

df_fixed_P_C0 <- map_dfr(P_C0, function(P_C0) {
  P_C1 <- (rr_confounded / rr_unconfounded * (P_C0 / 100 * (RR_CD - 1) + 1) - 1) / (RR_CD - 1)
  tibble(
    RR_CD       = RR_CD,
    P_C1        = P_C1 * 100,
    P_C0_fixed    = factor(paste0("P_C0 = ", P_C0), levels = P_C0_labels)
  )
}) %>%
  filter(P_C1 >= 0, P_C1 <= 100)

ref_df_fixed_P_C0 <- tibble(
  P_C0_fixed   = factor(P_C0_labels, levels = P_C0_labels),
  yintercept = P_C0
)

plot_fixed_P_C0 <- ggplot(df_fixed_P_C0, aes(x = RR_CD, y = P_C1, colour = P_C0_fixed)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_x_continuous(limits = c(1, 10), breaks = seq(1, 10, 1)) +
  scale_colour_manual(values = my_colours) +
  labs(
    x      = "Confounder strength (RR_CD)",
    y      = "Prevalence in dialysis group (%, P_C1)",
    title  = "B. Fixed prevalence in conservative management group",
    subtitle = paste0("RR_confounded = ", rr_confounded, ", RR_unconfounded = ", rr_unconfounded),
    colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", colour = "gray80"),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14)
  )
plot_fixed_P_C0

ggsave(
  file.path(results_path, "Supplemental/Figure_M_unmeasured_confounding_fix_P_C0.png"),
  plot   = plot_fixed_P_C0,
  width  = 8,
  height = 6,
  dpi    = 150
)

################################################################################
### 2D plot: fixed RR_CD — conservative management group prevalence (P_C0) vs required P_C1
################################################################################
RR_CD         <- c(2, 3, 4, 5)
P_C0          <- seq(1, 99, length.out = n)

RR_CD_labels  <- paste0("RR_CD = ", RR_CD)
my_colours2   <- setNames(c("#1D9E75", "#378ADD", "#D85A30", "#8B0000"), RR_CD_labels)

df_fixed_RR_CD <- map_dfr(RR_CD, function(rr_cd) {
  P_C1 <- (rr_confounded / rr_unconfounded * (P_C0 / 100 * (rr_cd - 1) + 1) - 1) / (rr_cd - 1)
  tibble(
    P_C0_fixed  = P_C0,
    P_C1 = P_C1 * 100,
    RR_CD = factor(paste0("RR_CD = ", rr_cd), levels = RR_CD_labels)
  )
}) %>%
  filter(P_C1 >= 0, P_C1 <= 100)

plot_fixed_RR_CD <- ggplot(df_fixed_RR_CD, aes(x = P_C0_fixed, y = P_C1, colour = RR_CD)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_colour_manual(values = my_colours2) +
  labs(
    x        = "Prevalence in conservative management group (%, P_C0)",
    y        = "Prevalence in dialysis group (%, P_C1)",
    title    = "C. Fixed strength between unmeasured confounder and outcome",
    subtitle = paste0("RR_confounded = ", rr_confounded, ", RR_unconfounded = ", rr_unconfounded),
    colour   = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", colour = "gray80"),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14)
  )
plot_fixed_RR_CD

ggsave(
  file.path(results_path, "Supplemental/Figure_M_unmeasured_confounding_fix_RR_CD.png"),
  plot   = plot_fixed_RR_CD,
  width  = 8,
  height = 6,
  dpi    = 150
)

# Example below Table
RR_CD <- 3
P_C0 <- 60
P_C1 <- (rr_confounded / rr_unconfounded * (P_C0 / 100 * (RR_CD - 1) + 1) - 1) / (RR_CD - 1)
cat("Example", "\n",
    "RR_CD:", RR_CD, "\n", # 3
    "P_C0 :", P_C0, "\n",  # 60
    "P_C1 :", P_C1 * 100)  # 12.7

# Examples in text
P_C0 <- rep(c(40, 60, 80), 4)
P_C1 <- c(0, 20, 40, 10, 30, 50, 20, 40, 60, 30, 50, 70)
gap <- P_C0 - P_C1
P_C0_C1 <- P_C1 / P_C0
k <- rr_confounded / rr_unconfounded
RR_CD <- 1 + (1 - k) / (k * P_C0 / 100 - P_C1 / 100)
data.frame(P_C0, P_C0, gap, P_C0_C1, RR_CD)
