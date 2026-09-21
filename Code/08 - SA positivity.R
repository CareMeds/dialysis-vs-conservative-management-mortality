################################################################################
### Decision for dialysis versus conservative management
### PART 8 - Sensitivity analysis for positivity
################################################################################

# set-up
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/"
)
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Results/"

# load libraries
library(patchwork) # combine figures
library(data.table)
library(foreach)   # parallel computation
library(doRNG)     # handle parallel seeds
set.seed(1)        # set seed for parallel backend

# load functions
source("Code/utils/weighting.R")
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/compute_absolute_relative_risks.R")

# load data
load("Data/cohort_with_prob.Rdata")

# Save untrimmed data
baseline_untrimmed <- baseline
rm(baseline)

################################################################################
### Sensitivity anaysis for non-positivity, investigate non-overlapping region of PS
################################################################################
# Define main survival model formula
survival_model <- survival::Surv(time2event_death_2y, event_death_2y) ~ trt

trim_meths <- c("unweighted", "IPTW", "overlap", "Crump", "Stürmer", "Walker")

# Initialize summary table
summary_table <- data.frame(
  analysis_name = character(),
  N = numeric(),
  Excluded = numeric(),
  Nonoverlap = numeric(),
  Imbalance = numeric(),
  RD = numeric(),
  RD_lower = numeric(),
  RD_upper = numeric(),
  dRMST = numeric(),
  dRMST_lower = numeric(),
  dRMST_upper = numeric(),
  HR = numeric(),
  HR_lower = numeric(),
  HR_upper = numeric(),
  stringsAsFactors = FALSE
)

for (trim_meth in trim_meths) {
  cat("Trimming method:", trim_meth, "\n")
  # Use untrimmed data for common_range, IPTW, and overlap methods
  if (trim_meth %in% c("Crump", "Stürmer", "Walker")) {
    # Apply trimming and reweight IPTW
    baseline_updated <- trim_propensity_scores(
      data = baseline_untrimmed,
      id_name = id_name,
      trt_var = trt_var,
      w_meth = "IPTW",
      model_PS = model_PS,
      catvar = catvar,
      contvar = contvar,
      trim_meth = trim_meth
    )$overlap$data
  } else{
    baseline_updated <- baseline_untrimmed
  }
  
  # Identify non-overlap patients even after trimming
  nonoverlap <- trim_propensity_scores(
    data = baseline_updated,
    id_name = id_name,
    trt_var = trt_var,
    w_meth = "unweighted",
    model_PS = model_PS,
    catvar = catvar,
    contvar = contvar,
    trim_meth = "common_range"
  )$nonoverlap
  
  # define sample and weights
  if (trim_meth == "unweighted") {
    # extract weights
    weights_meth <- rep(1, nrow(baseline_updated))
    
    # use untrimmed sample
    baseline <- copy(baseline_untrimmed)
  } else if (trim_meth == "IPTW" | trim_meth == "overlap") {
    # extract weights
    weights_meth <- baseline_untrimmed[[paste0("sw_", trim_meth)]]
    
    # use untrimmed sample
    baseline <- copy(baseline_untrimmed)
  } else if (trim_meth == "Crump" |
             trim_meth == "Stürmer" |
             trim_meth == "Walker") {
    # extract re-estimated weights
    weights_meth <- baseline_updated$w
    
    # use trimmed and re-estimated PS
    baseline <- copy(baseline_updated)
  }
  
  # create baseline table of non overlap
  table_one <- create_baseline_table(
    data = baseline,
    id_name = "LOPNR",
    weights = weights_meth,
    vars = listvar,
    categoricalVars = catvar,
    continuousVars = contvar,
    IQRVars = non_normal_vars,
    treatmentColumn = trt_var,
    treatmentLabel = treatment_label,
    controlLabel = control_label,
    tableCaption = ""
  )
  
  # Create PS distribution plot
  unweighted_fig_ps <- create_ps_distribution_plot(
    data = baseline,
    PS_varname = "ps",
    trt_varname = trt_var,
    PS_title_hist = trim_meth,
    PS_title_scaled_hist = trim_meth,
    xlab = ifelse(trim_meth == tail(trim_meths, 1), TRUE, FALSE),
    ylab = TRUE,
    palette = manual_colors[1:2]
  )
  assign(paste0("fig_unweighted_ps_", trim_meth), unweighted_fig_ps)
  
  weighted_fig_ps <- create_ps_distribution_plot(
    data = baseline,
    PS_varname = "ps",
    trt_varname = trt_var,
    weights = weights_meth,
    PS_title_hist = trim_meth,
    PS_title_scaled_hist = trim_meth,
    xlab = ifelse(trim_meth == tail(trim_meths, 1), TRUE, FALSE),
    ylab = TRUE,
    palette = manual_colors[1:2]
  )
  assign(paste0("fig_weighted_ps_", trim_meth), weighted_fig_ps)
  
  # compute estimates
  out_est <- compute_estimates_with_CI(
    data = baseline,
    id_name = id_name,
    unit = unit,
    horizon = horizon,
    model_PS = model_PS,
    event_var = outcome_var,
    competing_event_var = competing_events_var,
    time2event_var = time2outcome_var,
    trt_var = trt_var,
    w_meth = ifelse(
      trim_meth == "Crump" |
        trim_meth == "Stürmer" |
        trim_meth == "Walker",
      "IPTW",
      trim_meth
    ),
    trim_meth = ifelse(
      trim_meth == "Crump" |
        trim_meth == "Stürmer" |
        trim_meth == "Walker",
      trim_meth,
      "no_trimming"
    ),
    catvar = catvar,
    contvar = contvar,
    n_bootstraps = n_bootstraps,
    bootstrap_seed = 1
  )
  
  # Append summary statistics for current trimming method
  summary_table <- rbind(
    summary_table,
    data.frame(
      analysis_name = trim_meth,
      N = nrow(baseline),
      Excluded = nrow(baseline_untrimmed) - nrow(baseline),
      Nonoverlap = nrow(nonoverlap),
      Imbalance = sum(table_one$smd_table > 0.1),
      RD = out_est$RD * 100,
      RD_lower = out_est$RD_lower * 100,
      RD_upper = out_est$RD_upper * 100,
      dRMST = out_est$dRMST,
      dRMST_lower = out_est$dRMST_lower,
      dRMST_upper = out_est$dRMST_upper,
      HR = out_est$HR,
      HR_lower = out_est$HR_lower,
      HR_upper = out_est$HR_upper
    )
  )
}

################################################################################
### Combine histograms and scaled histograms of propensity score distributions
### across different trimming methods
################################################################################
# Define the layout:
# Each letter represents a plot area.
# # represents an empty space.
# We define a 4-column grid.
layout_design <- "
  ABCD
  ##EF
  GHIJ
  KLMN
  OPQR
"

# Combine the plots into a list to pass to wrap_plots
plot_list <- list(
  # Row 1
  fig_unweighted_ps_unweighted$hist + ggplot2::labs(title = "1A. Total population"),
  fig_unweighted_ps_unweighted$scaled_hist + ggplot2::labs(title = "1B. Total population"),
  fig_weighted_ps_IPTW$hist + ggplot2::labs(title = "2A. IPTW"),
  fig_weighted_ps_IPTW$scaled_hist + ggplot2::labs(title = "2B. IPTW"),
  
  # Row 2 (Starts at column 3 based on design ##EF)
  fig_weighted_ps_overlap$hist + ggplot2::labs(title = "3A. Overlap weighting"),
  fig_weighted_ps_overlap$scaled_hist + ggplot2::labs(title = "3B. Overlap weighting"),
  
  # Row 3
  fig_unweighted_ps_Crump$hist + ggplot2::labs(title = "4A. Crump trimming"),
  fig_unweighted_ps_Crump$scaled_hist + ggplot2::labs(title = "4B. Crump trimming"),
  fig_weighted_ps_Crump$hist + ggplot2::labs(title = "5A. Crump trimming"),
  fig_weighted_ps_Crump$scaled_hist + ggplot2::labs(title = "5B. Crump trimming"),
  
  # Row 4
  fig_unweighted_ps_Stürmer$hist + ggplot2::labs(title = "6A. Stürmer trimming"),
  fig_unweighted_ps_Stürmer$scaled_hist + ggplot2::labs(title = "6B. Stürmer trimming"),
  fig_weighted_ps_Stürmer$hist + ggplot2::labs(title = "7A. Stürmer trimming"),
  fig_weighted_ps_Stürmer$scaled_hist + ggplot2::labs(title = "7B. Stürmer trimming"),
  
  # Row 5
  fig_unweighted_ps_Walker$hist + ggplot2::labs(title = "8A. Walker trimming"),
  fig_unweighted_ps_Walker$scaled_hist + ggplot2::labs(title = "8B. Walker trimming"),
  fig_weighted_ps_Walker$hist + ggplot2::labs(title = "9A. Walker trimming"),
  fig_weighted_ps_Walker$scaled_hist + ggplot2::labs(title = "9B. Walker trimming")
)

# Generate final plot
final_plot <- patchwork::wrap_plots(plot_list, design = layout_design) +
  patchwork::plot_layout(guides = "collect")

# Save
ggplot2::ggsave(
  plot = final_plot,
  filename = paste0(results_path, "Supplemental/Figure_S_positivity.png"),
  width = 20,
  height = 20,
  dpi = 300
)

# create forest plots of subgroup analysis
summary_table <- summary_table[-1, ]
summary_table$analysis_name <- c(
  "Main analysis",
  "Overlap weighting",
  "Crump trimming",
  "Stürmer trimming",
  "Walker trimming"
)
nonoverlap_forest <- create_forest_plot_all_measures(dt = setDT(summary_table),
                                                     print_metrics = c("N", "Nonoverlap"))
ggplot2::ggsave(
  plot = nonoverlap_forest$combined_plot,
  filename = paste0(results_path, "Supplemental/Figure_S_positivity_forest.png"),
  width = 11,
  height = 8,
  dpi = 300
)
