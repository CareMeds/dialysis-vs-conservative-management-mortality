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

# Columns added per patient: in-center days, total days at home (days alive - in-center
# days), and days / in-center days from the start of dialysis (first HD or PD row)
# onward; 0 for patients who never start. "Started dialysis" (0/1) and the total home /
# in-center days counted only for starters feed the "among those who started" columns
in_center_cols <- c("In-center", "Time at home", "Days after dialysis", "In-center after dialysis",
                    "Started dialysis", "Home if started", "In-center if started",
                    "Hospital if started")

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
  state_cols            = state_cols,
  in_center_cols        = in_center_cols,
  in_center_params      = in_center_params
)
days_per_patient         <- main$days_per_patient
days_per_patient_imputed <- main$days_per_patient_imputed
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

# total days at home (days alive - in-center days) and its % of days alive, taken
# within each replicate. compute_all_estimates() now returns them; this only derives
# them from the stored replicates of an older saved bootstrap file, so it does not
# need to be rerun
if (!"At home %" %in% bootstrap_dt$state) {
  bootstrap_dt <- bootstrap_dt[state != "Time at home"]
  bootstrap_dt <- add_days_at_home(bootstrap_dt,
                                   alive_cols = setdiff(estimate_cols, c("Death", in_center_cols)),
                                   by_cols = c("b", "trt", "adjustment"))
}

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
# states in the figure annotation (days; excludes Death in reporting)
state_cols_in_center <- c(setdiff(state_cols, "Death"), "In-center", "Time at home")

# The Excel tables only report three measures, in MONTHS: total time at home, total time
# in-center (incl. hospital days) and total time in hospital. Days -> months with
# days_per_month days per month (as in the figure and the RMST). "Time at home" = days
# alive - in-center days (also counts days at home on dialysis).
days_per_month   <- 30.44  # 365.25 / 12 days per month
table_state_cols <- c("Time at home", "In-center", "Hospitalization")
table_labels     <- c(
  "Time at home"    = "Total months at home",
  "In-center"       = "Total months in-center",
  "Hospitalization" = "Total months in hospital"
)

# renames whichever of the labelled columns are present (all at once, by position)
rename_state_cols <- function(dt, labels) {
  dt  <- copy(dt)
  idx <- match(names(labels), names(dt))
  ok  <- !is.na(idx)
  setnames(dt, idx[ok], unname(labels[ok]))
  dt
}

# CI table: Dialysis / Conservative management / Difference per adjustment block
# (Difference is taken within each replicate in compute_all_estimates()); estimates and
# CIs converted from days to months, 1 decimal. extra_cols = FALSE: no percentage columns
# Patients who started dialysis (amongst those who chose dialysis) get their own row with
# the same measures as the arms.
# Row order per adjustment block: choose conservative management, choose dialysis, the
# difference choose dialysis vs. choose conservative management, then start dialysis
table_trt_labels <- c(
  `0`   = "Choose conservative management",
  `1`   = "Choose dialysis",
  diff  = "Choose dialysis vs. choose conservative management",
  start = "Start dialysis (amongst those who chose dialysis)"
)

estimate_ci_months_dt <- copy(estimate_ci_dt)
estimate_ci_months_dt[state %in% table_state_cols,
                      `:=`(estimate = estimate / days_per_month,
                           lower    = lower    / days_per_month,
                           upper    = upper    / days_per_month)]
diff_label <- "Difference"

# check of the table: for every measure and adjustment the point estimates (months) per row
# group, and the difference must equal choose dialysis - choose conservative management.
# Printed so the numbers can be compared with the Excel table
table_check <- dcast(estimate_ci_months_dt[state %in% table_state_cols],
                     adjustment + state ~ trt, value.var = "estimate")
print(table_check[, .(adjustment, state,
                      choose_CM = round(`0`, 2), choose_dialysis = round(`1`, 2),
                      difference = round(diff, 2), start_dialysis = round(start, 2))])
stopifnot(all(abs(table_check$diff - (table_check$`1` - table_check$`0`)) < 1e-8))

Bootstrap_CI_dt <- build_ci_table_dt(estimate_ci_months_dt, table_state_cols,
                                     trt_labels = table_trt_labels,
                                     digits = 1, extra_cols = FALSE)
Bootstrap_CI_with_diff_dt <- rename_state_cols(Bootstrap_CI_dt, table_labels)

# the bootstrap CI first (its point estimates are the pooled estimates), then one sheet per
# imputation's own unpooled estimates (m1...m10) - 11 tabs total
estimate_sheets <- c(
  list(Bootstrap_CI = Bootstrap_CI_with_diff_dt),
  setNames(
    lapply(imputations, function(imp)
      rename_state_cols(build_estimate_table_dt(
        imp$days_per_patient_imputed,
        days_per_patient,
        iptw_dt,
        estimate_cols,
        table_state_cols,
        days_per_month = days_per_month
      ), table_labels)),
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
### Colours for the states (used by the bar figure)
################################################################################
# the bar figure only uses three colours, the first three of manual_colors: green for time at
# home, blue for in-center (the Hospitalization colour) and red for death
state_levels <- c("At home", "Hospitalization", "Death")
state_colors <- setNames(manual_colors[c(3, 2, 1)], state_levels)

################################################################################
### Figure 2: months at home / in-center / death over 2 years, by trt
################################################################################
# Three bars: conservative management, choose dialysis, and an extra bar that splits the time
# at home of the dialysis arm into before / after the start of dialysis, with dashed lines
# linking the blocks so their sizes can be compared. Each of the first two bars: mean total time
# at home + mean total time in-center + death (the remaining time) = 2 years = 24 months.
# Confounding+censoring-adjusted estimates with bootstrap 95% CIs. Days are converted to months
# with days_per_month (30.44 = 365.25 / 12); death = 24 months - months alive. The figure is
# built by create_time_bar_figure() in plots.R
p_time_bars <- create_time_bar_figure(
  estimate_ci_dt   = estimate_ci_dt,
  bootstrap_dt     = bootstrap_dt,
  days_per_patient = days_per_patient,
  state_cols       = state_cols,
  state_colors     = state_colors,
  adjustment_level = "Adjusted for confounding and censoring",
  days_per_month   = days_per_month,
  total_months     = 24
)

p_time_bars
ggplot2::ggsave(plot = p_time_bars,
                filename = paste0(results_path, "Main/Figure_2.pdf"),
                width = 10, height = 5, dpi = 300)
ggplot2::ggsave(plot = p_time_bars,
                filename = paste0(results_path, "Main/Figure_2.png"),
                width = 10, height = 5, dpi = 300)