################################################################################
### Decision for dialysis versus conservative management
### PART 11 - Illustration relationship between RD, RR, RMST, and HR
################################################################################

# set-up
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/"
)
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Results/"

# load function to plot
source("Code/utils/plots.R")

# Define baseline risk
x_lim <- c(0.4, 0.9)
mortality_risk <- seq(x_lim[1], x_lim[2], by = 0.01)

tau <- 1  # time horizon (set to match the period your baseline risk covers)

# ---- RMST difference function (exponential survival assumption) ----
dRMST <- function(mortality_risk, hr, tau = 1) {
  lambda      <- -log(1 - mortality_risk) / tau
  rmst_ctrl   <- (1 - exp(-lambda * tau)) / lambda
  lambda_t    <- lambda * hr
  rmst_treat  <- (1 - exp(-lambda_t * tau)) / lambda_t
  return(rmst_treat - rmst_ctrl)
}

# ---- CONSTANT HAZARD RATIO (HR = 0.4) ----
hr_const <- 0.4
data_constant_HR <- data.frame(
  mortality_risk = mortality_risk,
  HR            = hr_const,
  RR            = (1 - (1 - mortality_risk)^hr_const) / mortality_risk,
  RD            = (1 - (1 - mortality_risk)^hr_const) - mortality_risk,
  Scenario      = paste0("Constant HR (", hr_const, ")")
)
data_constant_HR$dRMST <- dRMST(data_constant_HR$mortality_risk, data_constant_HR$HR, tau)

# ---- INCREASING HAZARD RATIO (0.2 to 0.8) ----
hr_start <- 0.2
hr_end   <- 0.8
data_inc_HR <- data.frame(mortality_risk = mortality_risk)
data_inc_HR$HR       <- hr_start + (hr_end - hr_start) * (data_inc_HR$mortality_risk - min(mortality_risk)) / (max(mortality_risk) - min(mortality_risk))
data_inc_HR$RR       <- (1 - (1 - data_inc_HR$mortality_risk)^data_inc_HR$HR) / data_inc_HR$mortality_risk
data_inc_HR$RD       <- (1 - (1 - data_inc_HR$mortality_risk)^data_inc_HR$HR) - data_inc_HR$mortality_risk
data_inc_HR$Scenario <- paste("Increasing HR from", hr_start, "to", hr_end)
data_inc_HR$dRMST    <- dRMST(data_inc_HR$mortality_risk, data_inc_HR$HR, tau)

# ---- CONSTANT RISK DIFFERENCE (RD = -0.3) ----
rd_const    <- -0.3
br_filtered <- mortality_risk[mortality_risk > abs(rd_const)]
data_constant_RD <- data.frame(
  mortality_risk = br_filtered,
  RD            = rd_const,
  RR            = (br_filtered + rd_const) / br_filtered,
  HR            = log(1 - (br_filtered + rd_const)) / log(1 - br_filtered),
  Scenario      = paste0("Constant RD (", rd_const * 100, ")")
)
data_constant_RD$dRMST <- dRMST(data_constant_RD$mortality_risk, data_constant_RD$HR, tau)

# ---- Y-axis limits ----
y_lim_rd   <- c(-0.4, 0)
y_lim_rmst <- c(0.3, 0)

# ---- Row 1: Constant HR ----
p_RD_constant_HR    <- plot_metric(data_constant_HR, "RD", "RD", x_lim = x_lim, y_lim = y_lim_rd)
p_RR_constant_HR    <- plot_metric(data_constant_HR, "RR", "RR", x_lim = x_lim)
p_dRMST_constant_HR <- plot_metric(data_constant_HR, "dRMST", "dRMST", x_lim = x_lim, y_lim = y_lim_rmst)
p_HR_constant_HR    <- plot_metric(data_constant_HR, "HR", "HR", x_lim = x_lim)

# ---- Row 2: Increasing HR ----
p_RD_inc_HR    <- plot_metric(data_inc_HR, "RD", "RD", x_lim = x_lim, y_lim = y_lim_rd)
p_RR_inc_HR    <- plot_metric(data_inc_HR, "RR", "RR", x_lim = x_lim)
p_dRMST_inc_HR <- plot_metric(data_inc_HR, "dRMST", "dRMST", x_lim = x_lim, y_lim = y_lim_rmst)
p_HR_inc_HR    <- plot_metric(data_inc_HR, "HR", "HR", x_lim = x_lim)

# ---- Row 3: Constant RD ----
p_RD_constant_RD    <- plot_metric(data_constant_RD, "RD", "RD", "Mortality Risk", x_lim = x_lim, y_lim = y_lim_rd)
p_RR_constant_RD    <- plot_metric(data_constant_RD, "RR", "RR", "Mortality Risk", x_lim = x_lim)
p_dRMST_constant_RD <- plot_metric(data_constant_RD, "dRMST", "dRMST", "Mortality Risk", x_lim = x_lim, y_lim = y_lim_rmst)
p_HR_constant_RD    <- plot_metric(data_constant_RD, "HR", "HR", "Mortality Risk", x_lim = x_lim)

# ---- Combine into 3x4 grid ----
final_plot <- ggpubr::ggarrange(
  p_RD_constant_HR,
  p_dRMST_constant_HR,
  p_RR_constant_HR,
  p_HR_constant_HR,
  p_RD_inc_HR,
  p_dRMST_inc_HR,
  p_RR_inc_HR,
  p_HR_inc_HR,
  p_RD_constant_RD,
  p_dRMST_constant_RD,
  p_RR_constant_RD,
  p_HR_constant_RD,
  nrow = 3,
  ncol = 4
)

# ---- Save ----
ggplot2::ggsave(
  plot     = final_plot,
  filename = file.path(results_path, "Supplemental/Figure_M_HTE_simulation.png"),
  width    = 10,
  height   = 5,
  dpi      = 600
)
