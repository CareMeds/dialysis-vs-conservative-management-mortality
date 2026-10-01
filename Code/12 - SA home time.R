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
# n_bootstraps <- 10

# TRUE = rerun the bootstrap and save the results to the Data folder;
# FALSE = load the previously saved bootstrap results from the Data folder
recompute_bootstrap <- TRUE
bootstrap_file <- paste0("Data/bootstrap_results_home_time_", n_bootstraps, ".Rdata")

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

# states in reporting order (not alphabetical). States absent from the data are dropped
state_order <- c("At home", "Hospitalization", "HD", "PD", "Death")
state_cols <- state_order[state_order %in% unique(bl$baseline_long$state)]

# Columns added per patient (see add_in_center_per_patient() in competing_risk.R): in-center
# days, total days at home (days alive - in-center days), their components (home without /
# on dialysis, in-center hospitalized / on dialysis), and days / in-center days / hospital
# days from the start of dialysis (first HD or PD row) onward; 0 for patients who never start
# (used for the before / after split in the figure)
in_center_cols <- c("In-center", "Time at home",
                    "Home without dialysis", "Home on dialysis",
                    "In-center hospitalized", "In-center dialysis",
                    "Days after dialysis", "In-center after dialysis", "Hospital after dialysis")

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

# a bootstrap file saved before the home / in-center split columns existed lacks them
stopifnot("saved bootstrap lacks the home / in-center split: set recompute_bootstrap <- TRUE" =
            "Home on dialysis" %in% bootstrap_dt$state)

estimate_ci_dt <- bootstrap_dt[, .(lower = quantile(value, 0.025, na.rm = TRUE),
                                   upper = quantile(value, 0.975, na.rm = TRUE)), by = .(trt, state, adjustment)]

point_long <- compute_all_estimates(days_per_patient_imputed,
                                    days_per_patient,
                                    iptw_dt,
                                    estimate_cols)
setnames(point_long, "value", "estimate")
estimate_ci_dt <- merge(point_long, estimate_ci_dt, by = c("trt", "state", "adjustment"))
setorder(estimate_ci_dt, trt, state)

# Share of the in-center time in the dialysis arm that is due to dialysis (HD sessions / PD
# visits) versus hospitalization, for the text of the manuscript. Ratio of the (weighted) mean
# months, per adjustment level; the 95% CI is the percentile CI over the bootstrap replicates,
# with the share taken within each replicate. The two shares add up to 100%
ic_share <- function(dt, value_col, by_cols) {
  d <- dt[trt == "1" & state %in% c("In-center", "In-center hospitalized", "In-center dialysis")]
  d[, v := get(value_col)]
  d[, {
    ic <- v[state == "In-center"]
    .(hospitalized = 100 * v[state == "In-center hospitalized"] / ic,
      dialysis     = 100 * v[state == "In-center dialysis"] / ic)
  }, by = by_cols]
}
ic_share_point <- ic_share(estimate_ci_dt, "estimate", "adjustment")
ic_share_boot  <- ic_share(bootstrap_dt, "value", c("b", "adjustment"))
stopifnot(all(abs(ic_share_point$hospitalized + ic_share_point$dialysis - 100) < 1e-8))

ic_share_ci <- rbindlist(lapply(c("hospitalized", "dialysis"), function(part) {
  boot <- ic_share_boot[, .(lower = quantile(get(part), 0.025, na.rm = TRUE),
                            upper = quantile(get(part), 0.975, na.rm = TRUE)),
                        by = adjustment]
  merge(ic_share_point[, .(adjustment, estimate = get(part))], boot, by = "adjustment")[
    , part := part]
}))

cat("\nShare of in-center time in the dialysis arm (% of in-center months)\n")
for (adj in unique(ic_share_ci$adjustment)) {
  cat("\n", adj, "\n", sep = "")
  for (pt in c("dialysis", "hospitalized")) {
    r <- ic_share_ci[adjustment == adj & part == pt]
    cat(sprintf("  %-13s %.1f%% (95%% CI %.1f, %.1f)\n",
                paste0(pt, ":"), r$estimate, r$lower, r$upper))
  }
}

################################################################################
### Estimate tables
################################################################################
# The Excel tables only report three measures, in MONTHS: total time at home, total time
# in-center (incl. hospital days) and total time in hospital. Days -> months with
# days_per_month days per month (as in the figure and the RMST). "Time at home" = days
# alive - in-center days (also counts days at home on dialysis).
days_per_month   <- 30.44  # 365.25 / 12 days per month
table_state_cols <- c("Time at home", "In-center", "Hospitalization")

# CI table: Dialysis / Conservative management / Difference per adjustment block
# (Difference is taken within each replicate in compute_all_estimates()); estimates and
# CIs converted from days to months, 1 decimal
estimate_ci_months_dt <- copy(estimate_ci_dt)
split_state_cols <- c("Home without dialysis", "Home on dialysis",
                      "In-center hospitalized", "In-center dialysis")
estimate_ci_months_dt[state %in% c(table_state_cols, split_state_cols),
                      `:=`(estimate = estimate / days_per_month,
                           lower    = lower    / days_per_month,
                           upper    = upper    / days_per_month)]

# check of the table: for every measure and adjustment the point estimates (months) per row
# group, and the difference must equal choose dialysis - choose conservative management.
# Printed so the numbers can be compared with the Excel table
table_check <- dcast(estimate_ci_months_dt[state %in% table_state_cols],
                     adjustment + state ~ trt, value.var = "estimate")
print(table_check[, .(adjustment, state,
                      choose_CM = round(`0`, 2), choose_dialysis = round(`1`, 2),
                      difference = round(diff, 2),
                      before_dialysis = round(before, 2), after_dialysis = round(after, 2))])
stopifnot(all(abs(table_check$diff - (table_check$`1` - table_check$`0`)) < 1e-8))
# time after the start is part of the total for choosing dialysis
stopifnot(all(table_check$after <= table_check$`1` + 1e-8))
# before + after the start of dialysis must add up to the total for choosing dialysis
stopifnot(all(abs(table_check$before + table_check$after - table_check$`1`) < 1e-8))

# the components must add up to their totals, for every group and adjustment
split_check <- dcast(estimate_ci_months_dt[trt %in% c("0", "1", "diff")],
                     adjustment + trt ~ state, value.var = "estimate")
stopifnot(all(abs(split_check[["Home without dialysis"]] + split_check[["Home on dialysis"]] -
                    split_check[["Time at home"]]) < 1e-8))
stopifnot(all(abs(split_check[["In-center hospitalized"]] + split_check[["In-center dialysis"]] -
                    split_check[["In-center"]]) < 1e-8))

# number of patients per arm for the column headers
n_dt <- days_per_patient[, .(N = uniqueN(LOPNR)), by = .(trt = as.character(trt))]

Bootstrap_CI_with_diff_dt <- build_home_time_table_dt(estimate_ci_months_dt, n_dt,
                                                      horizon_months = 24, digits = 1)

# one sheet: the bootstrap CI table (its point estimates are the pooled estimates)
write.xlsx(
  list(Bootstrap_CI = Bootstrap_CI_with_diff_dt),
  file = paste0(
    results_path,
    "Supplemental/Table_S_home_time_estimates_", 
    n_bootstraps,
    ".xlsx"
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

ggplot2::ggsave(plot = p_time_bars,
                filename = paste0(results_path, "Main/Figure_2.pdf"),
                width = 10, height = 5, dpi = 300)
