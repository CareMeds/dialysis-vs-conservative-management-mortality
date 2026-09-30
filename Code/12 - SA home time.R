################################################################################
### Decision for dialysis versus conservative management
### PART 12 - Sensitivity analysis: multi-state home time model
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
library(data.table)
library(mstate)
library(survival)
library(openxlsx)
library(patchwork)

# load functions (create_weights, for the IPTW bootstrap only)
source("Code/utils/competing_risk.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/weighting.R")
source("Code/utils/tables.R")
source("Code/utils/plots.R")

# load data
load("Data/cohort_with_prob.Rdata")
load("Data/long_cohort_hosp_krt.Rdata")

# set parameters
M <- 10
M_bootstrap <- 2

# TRUE = rerun the bootstrap and save the results to the Data folder;
# FALSE = load the previously saved bootstrap results from the Data folder
recompute_bootstrap <- FALSE
bootstrap_file <- "Data/bootstrap_results_home_time.Rdata"

# assumptions for the in-center days (see in_center_days in competing_risk.R; applied per state episode)
in_center_params <- list(
  hd_sessions_per_week = 3,   # in-center HD treatments per week
  pd_training_days     = 5,   # PD training at the start of PD
  pd_visit_interval    = 30   # days between PD outpatient check-ups
)

################################################################################
### Fixed modeling choices (same for the real analysis and every bootstrap replicate)
################################################################################
# every transition with data
transitions <- list(
  c("At home", "Hospitalization"),
  c("At home", "HD"),
  c("At home", "PD"),
  c("At home", "Death"),
  c("Hospitalization", "At home"),
  c("Hospitalization", "HD"),
  c("Hospitalization", "PD"),
  c("Hospitalization", "Death"),
  c("HD", "Hospitalization"),
  c("HD", "PD"),
  c("HD", "Death"),
  c("PD", "HD"),
  c("PD", "Hospitalization"),
  c("PD", "Death")
)

# transitions fit on dialysis-arm patients only (trt == 1); all others are stratified by trt x hospitalizations
dialysis_only <- c(
  "HD -> PD",
  "PD -> HD"
)

# one row per patient per day, day 0 to 729 (2-year horizon)
time_grid <- seq(0, horizon - 1, by = 1)

################################################################################
### Build baseline_long once (resampled by LOPNR in each bootstrap replicate)
################################################################################
bl <- build_baseline_long(long_cohort_hosp,
                          long_cohort_krt,
                          baseline,
                          transitions,
                          dialysis_only,
                          horizon)

write.xlsx(
  bl$transition_n_dt,
  file = paste0(results_path, "Supplemental/Table_M_n_per_transition.xlsx")
)

# states in reporting order (not alphabetical); sets the column order of tables
# (state_cols_in_center) and figures. States absent from the data are dropped
state_order <- c("At home", "Hospitalization", "HD", "PD", "Death")
state_cols <- state_order[state_order %in% unique(bl$baseline_long$state)]

# Columns added per patient: in-center days, and days / in-center days from the
# start of dialysis (first HD or PD row) onward; 0 for patients who never start
in_center_cols <- c("In-center", "Days after dialysis", "In-center after dialysis")

# columns summarised in the estimate tables: the states plus in-center days
estimate_cols <- c(state_cols, in_center_cols)

# set cap for hospitalizations (e.g., 7)
cap_hosp <- bl$cap_hosp

################################################################################
### Run the real analysis
################################################################################
main <- fit_and_simulate(
  baseline              = baseline,
  baseline_long         = bl$baseline_long,
  last_obs_dt           = bl$last_obs_dt,
  observed_part         = bl$observed_part,
  non_censored_part     = bl$non_censored_part,
  cap_hosp              = cap_hosp,
  model_PS              = model_PS,
  id_name               = id_name,
  trt_var               = trt_var,
  horizon               = horizon,
  M                     = M,
  transitions           = transitions,
  dialysis_only         = dialysis_only,
  time_grid             = time_grid,
  state_cols            = state_cols,
  in_center_cols        = in_center_cols,
  in_center_params      = in_center_params
)
days_per_patient         <- main$days_per_patient
days_per_patient_imputed <- main$days_per_patient_imputed
state_prob_dt            <- main$state_prob_dt
iptw_dt                  <- main$iptw_dt
imputations              <- main$imputations

################################################################################
### Bootstrap CI for the confounding+censoring-corrected estimate
################################################################################
# each replicate resamples patients and repeats the whole estimation (weights, Cox models,
# imputations); failed replicates are skipped. Replicates only give the percentile CI,
# the point estimate is the real (non-bootstrapped) one

# do this for all bootstraps - or load the saved results of a previous run
if (recompute_bootstrap) {
  bootstrap_results <- lapply(
    seq_len(n_bootstraps),
    run_one_bootstrap,
    baseline              = baseline,
    baseline_long         = bl$baseline_long,
    last_obs_dt           = bl$last_obs_dt,
    observed_part         = bl$observed_part,
    non_censored_part     = bl$non_censored_part,
    cap_hosp              = cap_hosp,
    model_PS              = model_PS,
    id_name               = id_name,
    trt_var               = trt_var,
    horizon               = horizon,
    M_bootstrap           = M_bootstrap,
    transitions           = transitions,
    dialysis_only         = dialysis_only,
    time_grid             = time_grid,
    state_cols            = state_cols,
    in_center_cols        = in_center_cols,
    in_center_params      = in_center_params
  )
  
  # save the raw results (incl. failed replicates, so the failure count
  # below stays correct after loading) together with the settings used
  bootstrap_settings <- list(
    n_bootstraps = n_bootstraps,
    M_bootstrap  = M_bootstrap,
    created      = Sys.time()
  )
  save(bootstrap_results, bootstrap_settings, file = bootstrap_file)
} else {
  load(bootstrap_file)
  
  # failure count below is based on the loaded replicates
  n_bootstraps <- length(bootstrap_results)
}

# capture number of failed bootstraps
failed <- vapply(bootstrap_results, is.null, logical(1))
cat(sprintf(
  "%d of %d bootstrap replicates failed and were skipped.\n",
  sum(failed),
  n_bootstraps
))
bootstrap_results <- bootstrap_results[!failed]

# percentile CI per (trt, state, adjustment), merged with the real
# (non-bootstrapped) point estimate for each of the 4 adjustment levels
bootstrap_dt <- rbindlist(lapply(bootstrap_results, `[[`, "final"), idcol = "b")
estimate_ci_dt <- bootstrap_dt[, .(lower = quantile(value, 0.025, na.rm = TRUE),
                                   upper = quantile(value, 0.975, na.rm = TRUE)), by = .(trt, state, adjustment)]

point_long <- compute_all_estimates(days_per_patient_imputed,
                                    days_per_patient,
                                    iptw_dt,
                                    estimate_cols,
                                    in_center_cols)
setnames(point_long, "value", "estimate")
estimate_ci_dt <- merge(point_long, estimate_ci_dt, by = c("trt", "state", "adjustment"))
setorder(estimate_ci_dt, trt, state)

################################################################################
### Estimate tables
################################################################################
# excludes Death in reporting
state_cols_in_center <- c(setdiff(state_cols, "Death"), "In-center")

# CI table: Dialysis / Conservative management / Difference per adjustment block
# (Difference is taken within each replicate in compute_all_estimates()).
# build_ci_table_dt() also adds the columns: In-center (% of days alive),
# In-center after dialysis (days), In-center (% of days after starting dialysis).
# The percentages are ratios of means computed in compute_all_estimates(), within
# each bootstrap replicate, so their CIs are percentile CIs of the ratio itself.
# The after-dialysis columns are dialysis arm only (empty for CM and Difference).
diff_label <- "Difference"
Bootstrap_CI_dt <- build_ci_table_dt(estimate_ci_dt, state_cols_in_center)
Bootstrap_CI_with_diff_dt <- copy(Bootstrap_CI_dt)
setnames(Bootstrap_CI_with_diff_dt, "In-center", "In-center (days)")

# the bootstrap CI first, then "Pooled", then one sheet per imputation's
# own unpooled estimates (m1...m10) - 12 tabs total
estimate_sheets <- c(
  list(
    Bootstrap_CI = Bootstrap_CI_with_diff_dt,
    Pooled = build_estimate_table_dt(
      days_per_patient_imputed,
      days_per_patient,
      iptw_dt,
      estimate_cols,
      state_cols_in_center
    )
  ),
  setNames(
    lapply(imputations, function(imp)
      build_estimate_table_dt(
        imp$days_per_patient_imputed,
        days_per_patient,
        iptw_dt,
        estimate_cols,
        state_cols_in_center
      )),
    paste0("m", seq_along(imputations))
  )
)

write.xlsx(
  estimate_sheets,
  file = paste0(
    results_path,
    "Supplemental/Table_S_home_time_estimates.xlsx"
  ),
  rowNames = FALSE
)

################################################################################
### Figure: state occupancy over time (observed + simulated), by trt
################################################################################
# pooled (multiply imputed) mean cumulative days
# time-varying CI for the left (cumulative) panel only: percentile CI per
# (time, trt, state) across the bootstrap replicates' full curves
state_prob_bootstrap_dt <- rbindlist(lapply(bootstrap_results, `[[`, "state_prob_dt"), idcol = "b")
state_prob_ci_dt <- state_prob_bootstrap_dt[, .(
  lower = quantile(mean_cumulative_days, 0.025),
  upper = quantile(mean_cumulative_days, 0.975)
), by = .(time, trt, state)]
state_prob_dt <- merge(state_prob_dt, state_prob_ci_dt, by = c("time", "trt", "state"))

# stacking order: "At home" at bottom, "Death" at top - flip state_levels
# if your ggplot2 version stacks the other way
state_levels <- c("At home", "Hospitalization", "HD", "PD", "Death")
state_prob_dt[, state := factor(state, levels = state_levels)]

# reuses the project's manual_colors if in scope; previously assigned 6
# colors to 5 names - fixed here to exactly 5, Death using "#FF7F00"
state_colors <- setNames(c(manual_colors[c(3, 2, 4, 1)], "#FF7F00"), state_levels)

################################################################################
### Figure: state probability (left) and state occupancy (right), by trt
################################################################################
# only one panel keeps its legend (colors are shared)
p_cum_dialysis    <- make_state_panel(state_prob_dt = state_prob_dt, 
                                      state_colors = state_colors, 
                                      trt_value = 1, 
                                      type = "cumulative", 
                                      show_legend = TRUE) +
  ggplot2::labs(title = "Choose dialysis")
p_stack_dialysis  <- make_state_panel(state_prob_dt = state_prob_dt, 
                                      state_colors = state_colors, 
                                      trt_value = 1, 
                                      type = "stacked", 
                                      show_legend = FALSE)
p_cum_cm          <- make_state_panel(state_prob_dt = state_prob_dt, 
                                      state_colors = state_colors, 
                                      trt_value = 0, 
                                      type = "cumulative",
                                      show_legend = FALSE) +
  ggplot2::labs(title = "Choose conservative management")
p_stack_cm        <- make_state_panel(state_prob_dt = state_prob_dt, 
                                      state_colors = state_colors, 
                                      trt_value = 0, 
                                      type = "stacked",
                                      show_legend = FALSE)

################################################################################
### Annotate the state-probability panels
################################################################################
# annotation includes In-center (a subset of other states, so text only, not a curve)
# rows after the header of the "confounding and censoring" block: Dialysis, then CM
final_block_row <- which(Bootstrap_CI_dt$trt == "Adjusted for confounding and censoring")
p_cum_dialysis <- annotate_state_probability(
  p              = p_cum_dialysis,
  estimate_ci_dt = estimate_ci_dt,
  trt_value      = 1,
  state_cols     = state_cols_in_center,
  state_prob_dt  = state_prob_dt
)
p_cum_cm <- annotate_state_probability(
  p              = p_cum_cm,
  estimate_ci_dt = estimate_ci_dt,
  trt_value      = 0,
  state_cols     = state_cols_in_center,
  state_prob_dt  = state_prob_dt
)

combined_plot <- (p_cum_dialysis | p_stack_dialysis) /
  (p_cum_cm | p_stack_cm) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "center"
  )

ggplot2::ggsave(
  plot = combined_plot,
  filename = paste0(results_path, "Supplemental/Figure_S_home_time.png"),
  width = 10,
  height = 8,
  dpi = 300
)
