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
source("Code/utils/data_manipulation.R")
source("Code/utils/weighting.R")

# load data
load("Data/cohort_with_prob.Rdata")
# baseline: id_name, model_PS, w_meths, horizon, etc.
load("Data/long_cohort_hosp_krt.Rdata")
# long_cohort_hosp: LOPNR, hosp_start, hosp_stop,
# decision_date, DODSDAT, followup_end, time2event_death_2y
# long_cohort_krt: LOPNR, krt_modality, krt_startdate, krt_stopdate
# (one row per dialysis modality period - a patient can have more than one
# if they switch between HD and PD)

# Check number of patients
table(baseline$trt)
length(unique(long_cohort_hosp$LOPNR))
length(unique(long_cohort_krt$LOPNR))

# TODO: move functions to utils
# TODO: CI bootstrapping
# TODO: make it multiple imputation? (m=10), then average only the number of days spent in states

# DISCUSS Ilarie: number of infinite draw
# DISCUSS Ilaria: number of ties
# DISCUSS Ilaria: simulated vs observed event times

# DISCUSS Marie: 4 conservative management start PD; 23 conservative management start HD
# DISCUSS Marie: no KRT modality switches (HD -> PD; PD -> HD)
# DISCUSS Marie: those that have one-day hospitalizations on decision dates are not considered at home, those that have one-day hospitalization after the decision date are considered hospitalizations
# DISCUSS Marie: "At home" -> HD/PD starts that don't go through an inpatient admission
# DISCUSS Marie: we split hospitalization when someone starts KRT in that hospitalization

################################################################################
### Prepare the data
################################################################################
# one row per patient: static info that never changes over time
patient_static_dt <- unique(long_cohort_hosp[, .(LOPNR, trt, decision_date, DODSDAT, followup_end)])

# hospitalization episodes, one row per episode
hosp_dt <- long_cohort_hosp[, .(LOPNR, hosp_start, hosp_stop)]

# a same-day hospitalization (hosp_start == hosp_stop) has a real
# 1-day duration - extend hosp_stop by a day so the normal breakpoint chain
# below creates a proper interval for it (whatever comes next automatically
# starts a day later too, via tstop's lead-shift), EXCEPT when it coincides
# with the patient's own decision_date: there, a same-instant admission and
# the start of follow-up would otherwise tie in later day-level analyses,
# so those are left untouched and stay invisible to the breakpoint chain
hosp_dt <- patient_static_dt[, .(LOPNR, decision_date)][hosp_dt, on = "LOPNR"]
hosp_dt[hosp_start == hosp_stop &
          hosp_start != decision_date, hosp_stop := hosp_stop + 1]
hosp_dt[, decision_date := NULL]

# number hospitalization episodes per patient (1st, 2nd, 3rd, ...)
setorder(hosp_dt, LOPNR, hosp_start)
hosp_dt[, hosp_num := seq_len(.N), by = LOPNR]

# dialysis (KRT) periods; multiple rows per patient if modality switches (see file 02)
krt_dt <- long_cohort_krt[, .(LOPNR, krt_modality, krt_startdate, krt_stopdate)]

# number KRT periods (increments on every switch, unlike dialysis_num further down)
setorder(krt_dt, LOPNR, krt_startdate)
krt_dt[, krt_num := seq_len(.N), by = LOPNR]

# stack every candidate breakpoint into one table (every krt_startdate is its own breakpoint)
break_points_dt <- rbindlist(
  list(
    patient_static_dt[, .(LOPNR, break_point = decision_date)],
    patient_static_dt[!is.na(DODSDAT) &
                        DODSDAT <= followup_end, .(LOPNR, break_point = DODSDAT)],
    patient_static_dt[, .(LOPNR, break_point = followup_end)],
    krt_dt[, .(LOPNR, break_point = krt_startdate)],
    hosp_dt[, .(LOPNR, break_point = hosp_start)],
    hosp_dt[, .(LOPNR, break_point = hosp_stop)]
  )
)

# sort chronologically so consecutive rows form [tstart, tstop) intervals
setorder(break_points_dt, LOPNR, break_point)

# drop exact duplicate (LOPNR, break_point) pairs (e.g. same-day coincidence)
break_points_dt <- unique(break_points_dt)

# bring back the static patient-level columns onto every breakpoint row
break_points_dt <- patient_static_dt[break_points_dt, on = "LOPNR"]

# ensure no breakpoints after follow-up end are added
break_points_dt <- break_points_dt[break_point <= followup_end]

# tstop = next breakpoint (lead); last breakpoint per patient gets NA
break_points_dt[, tstop := shift(break_point, type = "lead"), by = LOPNR]

# dates needing their own zero-width row: death only. Same-day
# hospitalizations are handled above (extended to a real 1-day interval,
# or left invisible if they coincide with decision_date) rather than
# needing a separate zero-width row here; same-day KRT periods still do.
zero_width_dates <- unique(rbindlist(list(
  patient_static_dt[!is.na(DODSDAT) &
                      DODSDAT <= followup_end, .(LOPNR, break_point = DODSDAT)],
  krt_dt[krt_startdate == krt_stopdate, .(LOPNR, break_point = krt_startdate)] # two patients start KRT right at the end of follow-up
)))
zero_len_rows <- break_points_dt[zero_width_dates, on = .(LOPNR, break_point), nomatch = NULL]

# collapse to zero width
zero_len_rows[, tstop := break_point]

# add the duplicates back in
break_points_dt <- rbindlist(list(break_points_dt, zero_len_rows))

# drop trailing rows with no zero-width duplicate
break_points_dt <- break_points_dt[!is.na(tstop)]

# zero-width row sorts right before its pair
setorder(break_points_dt, LOPNR, break_point, tstop)

# rename break_point -> tstart now that it's paired with tstop
setnames(break_points_dt, "break_point", "tstart")

# attach the hospitalization episode (if any) each interval falls inside
match_idx <- hosp_dt[break_points_dt, on = .(LOPNR, hosp_start <= tstart, hosp_stop >= tstop), which = TRUE]
break_points_dt[, `:=`(
  hosp_start = hosp_dt$hosp_start[match_idx],
  hosp_stop  = hosp_dt$hosp_stop[match_idx],
  hosp_num   = hosp_dt$hosp_num[match_idx]
)]

# attach the active KRT modality (if any) to each interval - same
# containment join as hosp_num above. TX is intentionally not treated as
# dialysis below (fcase only matches "HD"/"PD").
match_idx <- krt_dt[break_points_dt, on = .(LOPNR, krt_startdate <= tstart, krt_stopdate >= tstop), which = TRUE]
break_points_dt[, `:=`(
  krt_modality  = krt_dt$krt_modality[match_idx],
  krt_start     = krt_dt$krt_startdate[match_idx],
  krt_stop      = krt_dt$krt_stopdate[match_idx],
  krt_num       = krt_dt$krt_num[match_idx]
)]

# find, per patient, which hospitalization (if any) their first-ever
# HD/PD start fell inside - the one admission allowed to switch to HD/PD
first_krt_dt <- krt_dt[krt_modality %in% c("HD", "PD"), .(krt_startdate = min(krt_startdate)), by = LOPNR]
initiating_hosp_dt <- hosp_dt[first_krt_dt, on = .(LOPNR, hosp_start <= krt_startdate, hosp_stop >= krt_startdate), .(LOPNR, initiating_hosp_num = hosp_num)]
initiating_hosp_dt <- unique(initiating_hosp_dt[!is.na(initiating_hosp_num)])
break_points_dt[initiating_hosp_dt, initiating_hosp_num := i.initiating_hosp_num, on = "LOPNR"]

################################################################################
### Assign states
################################################################################
# ---- raw / unnumbered state category, priority order: ----
# 1. Death overrides everything
# 2/3. HD/PD, but only within the initiating hospitalization; every other
#      admission stays "Hospitalization" even if a modality switch occurs
# 4. HD/PD for non-hospitalized intervals
# 5. At home (default)
states_dt <- break_points_dt[, state := fcase(
  !is.na(DODSDAT) &
    tstart == DODSDAT,
  "Death",
  !is.na(hosp_num) &
    !is.na(initiating_hosp_num) &
    hosp_num == initiating_hosp_num &
    !is.na(krt_modality) &
    krt_modality == "HD",
  "HD",
  !is.na(hosp_num) &
    !is.na(initiating_hosp_num) &
    hosp_num == initiating_hosp_num &
    !is.na(krt_modality) &
    krt_modality == "PD",
  "PD",
  !is.na(hosp_start),
  "Hospitalization",
  !is.na(krt_modality) &
    krt_modality == "HD",
  "HD",
  !is.na(krt_modality) &
    krt_modality == "PD",
  "PD",
  default = "At home"
)]

################################################################################
### Merge consecutive same-state rows split by an unrelated breakpoint
### (e.g. krt_startdate splitting a hospitalization)
################################################################################
# number each run of consecutive same-state rows per patient
states_dt[, merge_state_id := rleid(state), by = LOPNR]

# carry forward every column except the grouping keys and tstart/tstop
merge_cols <- setdiff(names(states_dt),
                      c("LOPNR", "tstart", "tstop", "merge_state_id"))

# collapse each run into one row. Most columns just take the last row's
# value. hosp_start/hosp_stop/hosp_num are the exception: if the initiating
# hospitalization gets recoded to HD/PD, its row may sit in the MIDDLE of
# the run rather than at the end, so "last row" would wipe them to NA -
# grab them from wherever they actually are instead (only one row per run
# can ever have them).
preserve_cols <- c("hosp_start", "hosp_stop", "hosp_num")
states_dt <- states_dt[, c(list(tstart = min(tstart), tstop = max(tstop)),
                           Map(function(x, nm)
                             if (nm %in% preserve_cols)
                               x[!is.na(x)][1]
                             else
                               x[.N], .SD, names(.SD))), by = .(LOPNR, merge_state_id), .SDcols = merge_cols]
states_dt[, merge_state_id := NULL]

################################################################################
### Number the states
################################################################################
# numbers consecutive runs of TRUE in flag_col as 1, 2, 3...; FALSE rows get NA
number_runs <- function(dt, flag_col, episode_col, id_col = "LOPNR") {
  run_id_col <- paste0(episode_col, "_run_id")
  
  # temporary run-id for every consecutive stretch of the flag
  dt[, (run_id_col) := rleid(get(flag_col)), by = id_col]
  
  # one row per TRUE run - these are what get numbered
  runs <- unique(dt[get(flag_col) == TRUE, c(id_col, run_id_col), with = FALSE])
  setorderv(runs, c(id_col, run_id_col))
  
  # number those runs 1, 2, 3... in chronological order, per patient
  runs[, (episode_col) := seq_len(.N), by = id_col]
  
  # join episode number back; FALSE rows stay NA (absent from runs)
  dt[runs, (episode_col) := get(paste0("i.", episode_col)), on = c(id_col, run_id_col)]
  
  # drop the temporary run-id column
  dt[, (run_id_col) := NULL]
  invisible(dt)
}

# ---- number dialysis episodes (HD and PD combined into one counter) ----
states_dt[, is_dialysis := state %in% c("HD", "PD")]
number_runs(states_dt, flag_col = "is_dialysis", episode_col = "dialysis_num")

# ---- number "at home" episodes ----
states_dt[, is_home := state == "At home"]
number_runs(states_dt, flag_col = "is_home", episode_col = "home_num")

# find which dialysis episode (if any) was active right before each
# hospitalization started, for the "Hospitalization after Nth dialysis" label
pre_hosp_dialysis <- unique(states_dt[is_dialysis == TRUE, .(LOPNR, hosp_start = tstop, dialysis_num)])
states_dt[pre_hosp_dialysis, pre_hosp_dialysis_num := i.dialysis_num, on = .(LOPNR, hosp_start)]

################################################################################
### Inspect episode-count distributions to choose sensible caps
################################################################################
# max episode number reached per patient, for each episode type
max_hosp_per_patient     <- states_dt[!is.na(hosp_num), .(max_n = max(hosp_num)), by = LOPNR]
max_dialysis_per_patient <- states_dt[!is.na(dialysis_num), .(max_n = max(dialysis_num)), by = LOPNR]
max_home_per_patient     <- states_dt[!is.na(home_num), .(max_n = max(home_num)), by = LOPNR]

# cat("Hospitalization episodes per patient:\n")
# print(table(max_hosp_per_patient$max_n))
# print(quantile(max_hosp_per_patient$max_n, c(0.75, 0.9, 0.95)))
# 
# cat("\nDialysis episodes per patient:\n")
# print(table(max_dialysis_per_patient$max_n))
# print(quantile(max_dialysis_per_patient$max_n, c(0.75, 0.9, 0.95)))
# 
# cat("\nAt-home episodes per patient:\n")
# print(table(max_home_per_patient$max_n))
# print(quantile(max_home_per_patient$max_n, c(0.75, 0.9, 0.95)))

# remove helper columns
states_dt[, `:=`
          (
            is_home = NULL,
            home_num = NULL,
            pre_hosp_dialysis_num = NULL,
            dialysis_num = NULL,
            is_dialysis = NULL
          )]

################################################################################
### Attach baseline covariates to build the final long-format dataset
################################################################################
# patients with zero hospitalizations are missing from states_dt entirely -
# add them back as one "At home" run (+ "Death" if they died)
missing_id <- setdiff(baseline$LOPNR, states_dt$LOPNR)

missing_dt <- baseline[LOPNR %in% missing_id, .(LOPNR, decision_date = visit_date, followup_end, event_death_2y)]

missing_states_dt <- rbindlist(list(missing_dt[, .(
  LOPNR,
  decision_date,
  tstart = decision_date,
  tstop  = followup_end,
  state = "At home",
  home_num = 1L
)], missing_dt[event_death_2y == 1, .(
  LOPNR,
  decision_date,
  tstart = followup_end,
  tstop  = followup_end,
  state = "Death"
)]), fill = TRUE)

states_dt <- rbindlist(list(states_dt, missing_states_dt), fill = TRUE)
setorder(states_dt, LOPNR, tstart)

# running count of hospitalizations before this row: carry hosp_num forward
# from the previous row (locf). Placed after missing_states_dt so patients
# with zero hospitalizations correctly get 0 throughout.
states_dt[, prev_nr_hosp := nafill(shift(hosp_num), type = "locf"), by = LOPNR]

# no admissions yet
states_dt[is.na(prev_nr_hosp), prev_nr_hosp := 0L]

# broadcast baseline covariates onto every row of states_dt; shared column
# names (trt, decision_date) keep baseline's version
states_only_cols <- setdiff(names(states_dt), names(baseline))
baseline_long <- baseline[states_dt[, c("LOPNR", states_only_cols), with = FALSE], on = "LOPNR"]

################################################################################
### Build cause-specific Cox models, with history of hospitalization
################################################################################
# destination state after each episode; baseline_long is already sorted
# by LOPNR, tstart (nothing after the last setorder() reorders it)
baseline_long[, next_state := shift(state, type = "lead"), by = LOPNR]

# time at risk in this episode - clock resets at every state entry
baseline_long[, time_in_state := as.numeric(tstop - tstart)]

# baseline cumulative hazard per stratum, at the model's reference covariate level
get_basehaz <- function(model) {
  # extract baseline hazard from coxph object
  bh <- as.data.table(basehaz(model, centered = FALSE))
  
  # order on strata and time
  setorder(bh, strata, time)
  
  # return baseline hazard dt
  bh
}

# cause-specific Cox model for one origin -> destination transition, plus
# everything simulate_cause_time() needs to simulate from it: the
# baseline hazard per stratum, and the linear predictor for every
# trt x prev_nr_hosp_capped combination the model can see. Both are
# computed once here, right after fitting, instead of being recomputed
# every time a patient's trajectory is simulated.
fit_cause_specific_cox <- function(dt,
                                   origin,
                                   destination,
                                   cap,
                                   no_trt_transitions,
                                   weight_col = "sw_IPTW") {
  # episodes starting in `origin`
  origin_dt <- dt[state == origin]
  
  # other destinations = censored
  origin_dt[, event := fifelse(!is.na(next_state) &
                                 next_state == destination, 1L, 0L)]
  
  # non-parametric stratum by capped number of previous hospitalization
  origin_dt[, prev_nr_hosp_capped := pmin(prev_nr_hosp, cap)]
  
  # include "trt" in the model unless this exact transition is listed as
  # not estimable (no_trt_transitions is a character vector of
  # "origin -> destination" strings, same naming convention as cox_models)
  include_trt <- !paste(origin, "->", destination) %in% no_trt_transitions
  
  # coxph only allows one strata() term, so when trt is included it's
  # combined with prev_nr_hosp_capped into one joint stratum column
  # (prev_nr_hosp_capped itself is kept as-is above, still needed elsewhere)
  if (include_trt) {
    origin_dt[, trt_prev_nr_hosp_capped := paste(trt, prev_nr_hosp_capped, sep = "_")]
    strata_var <- "trt_prev_nr_hosp_capped"
  } else {
    strata_var <- "prev_nr_hosp_capped"
  }
  
  # build formula for cox fit
  formula <- if (include_trt) {
    Surv(time_in_state, event) ~ strata(trt_prev_nr_hosp_capped)
  } else {
    Surv(time_in_state, event) ~ strata(prev_nr_hosp_capped)
  }
  
  # IPTW weights make each stratum's baseline hazard a causal estimate, not just descriptive
  model <- coxph(formula, data = origin_dt, weights = origin_dt[[weight_col]])
  
  # extract baseline hazard of cox model, keyed by whichever stratum column
  # this model actually used
  bh <- get_basehaz(model)
  bh[, (strata_var) := sub(paste0(strata_var, "="), "", strata)]
  bh_by_stratum <- split(bh[, .(time, hazard)], bh[[strata_var]])
  
  # data element named after origin, e.g. "At home" -> "at_home_dt"
  data_name <- paste0(gsub(" ", "_", tolower(origin)), "_dt")
  c(
    list(
      model = model,
      bh_by_stratum = bh_by_stratum,
      include_trt = include_trt
    ),
    setNames(list(origin_dt), data_name)
  )
}

# every transition from the STEPS list above (HD->PD, PD->HD excluded)
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
  c("HD", "At home"),
  c("HD", "Death"),
  c("PD", "Hospitalization"),
  c("PD", "At home"),
  c("PD", "Death")
)

# a few conservative-management (trt=0) patients do end up on HD/PD, but
# almost none of them are ever observed going home from there - trt has
# essentially no events in one arm for these two transitions, so it isn't
# estimable and is left out of the formula just for these
no_trt_transitions <- c("HD -> At home", "PD -> At home")

# fit cox models for all transitions listed above
cap_hosp <- quantile(max_hosp_per_patient$max_n, 0.9)
cox_models <- setNames(
  lapply(transitions, function(state)
    fit_cause_specific_cox(
      dt = baseline_long,
      origin = state[1],
      destination = state[2],
      cap = cap_hosp,
      no_trt_transitions = no_trt_transitions
    )),
  sapply(transitions, function(state)
    paste(state[1], "->", state[2]))
)

# example: "At home" -> "Hospitalization"
summary(cox_models[["At home -> Hospitalization"]]$model)

# simulate one event time for a single cause-specific model (inverts H(t)):
# H(t) = H0(t) * exp(linear predictor). Reads the baseline hazard and
# linear predictor precomputed by fit_cause_specific_cox() (bh_by_stratum,
# lp_by_trt) instead of recomputing basehaz()/predict() here.
simulate_cause_time <- function(fit,
                                trt_chr,
                                hosp_capped,
                                is_first_event = FALSE,
                                current_event_time = NULL,
                                log_env = NULL,
                                transition_name = NULL) {
  # extract baseline hazard for trt and number of previous hospitalizations
  hosp_chr <- as.character(hosp_capped)
  stratum_key <- if (fit$include_trt) {
    paste(trt_chr, hosp_chr, sep = "_")
  } else {
    hosp_chr
  }
  bh_s <- fit$bh_by_stratum[[stratum_key]]
  
  # invert: smallest t where H(t) crosses -log(u); never crossing = no event observed
  # draw from uniform distribution to get draw of cumulative hazard
  U <- runif(1)
  # set threshold on cumulative hazard scale
  if (is_first_event) {
    # cumulative hazard already accrued by current_event_time (step
    # function: value at the last observed time <= current_event_time, or
    # 0 if current_event_time is before any observed event in this stratum)
    bh_entry <- bh_s$hazard[which(bh_s$time <= current_event_time)]
    bh_entry <- if (length(bh_entry) == 0) {
      0
    } else {
      bh_entry[length(bh_entry)]
    }
    # S(t | T > t_entry) = S(t) / S(t_entry) = exp(-(H(t) - H(t_entry)))
    # setting U = S(t | T > t_entry):
    # -log(U) = H(t) - H(t_entry)  =>  H(t) = -log(U) + H(t_entry)
    threshold <- -log(U) + bh_entry
  } else {
    threshold <- -log(U)
  }
  
  # Note: bh_s$hazard is actually the cumulative hazard
  result <- if (bh_s$hazard[length(bh_s$hazard)] < threshold) {
    # set to Inf if cumulative hazard never reaches the threshold - this
    # happens when the drawn threshold exceeds the maximum cumulative
    # hazard ever observed for this transition/stratum, i.e. the baseline
    # hazard curve is being asked to extrapolate beyond its observed range
    Inf
  } else {
    # smallest t where H(t) crosses the threshold
    bh_s$time[which(bh_s$hazard >= threshold)[1]]
  }
  
  # Log diagnostics if a logging environment was supplied: count every
  # candidate draw, and separately count draws that came back Inf (i.e. the
  # hazard threshold was never reached within the observed data, so no event
  # time could be assigned), broken down by transition. This is just for
  # checking afterwards how often each transition had to extrapolate beyond
  # its observed range - it doesn't affect the simulated trajectories.
  if (!is.null(log_env)) {
    log_env$n_candidates <- log_env$n_candidates + 1L
    if (is.infinite(result)) {
      log_env$n_inf <- log_env$n_inf + 1L
      prev <- log_env$n_inf_by_transition[[transition_name]]
      log_env$n_inf_by_transition[[transition_name]] <- if (is.null(prev)) 1L else prev + 1L
    }
  }
  
  result
}

# simulate every competing destination from `origin`, keep the smallest (soonest)
simulate_next_event <- function(cox_models,
                                origin,
                                trt_chr,
                                hosp_capped,
                                is_first_event = FALSE,
                                current_event_time = NULL,
                                log_env = NULL) {
  # select the cox models with possible transitions
  candidates <- names(cox_models)[startsWith(names(cox_models), paste0(origin, " -> "))]
  
  # simulate time to event for each possible transition
  times <- vapply(candidates, function(nm) {
    simulate_cause_time(
      fit                = cox_models[[nm]],
      trt_chr            = trt_chr,
      hosp_capped        = hosp_capped,
      is_first_event     = is_first_event,
      current_event_time = current_event_time,
      log_env            = log_env,
      transition_name    = nm
    )
  }, numeric(1))
  
  # a "tie" = two or more candidate transitions land on the exact same
  # finite simulated time. This is only meaningful when the winning time is
  # finite - if the minimum is Inf, every tied candidate is Inf too, which
  # just means no transition happens before the horizon, not a genuine
  # competing-event tie. which.min() below silently picks whichever tied
  # candidate happens to come first in cox_models, so this tally just
  # counts how often that (currently arbitrary) tie-break is invoked.
  min_time <- min(times)
  n_tied   <- sum(times == min_time)
  if (!is.null(log_env)) {
    log_env$n_steps <- log_env$n_steps + 1L
    if (is.finite(min_time) && n_tied > 1) {
      log_env$n_ties <- log_env$n_ties + 1L
    }
  }
  
  # return the event time and next state for soonest simualted event
  list(time        = min_time,
       destination = sub(paste0(origin, " -> "), "", names(which.min(times))))
}

# example: a real patient's actual current state and covariates
test_patient <- baseline_long[LOPNR == 9116134][.N]
test_patient[, .(LOPNR, tstart, tstop, state, trt, prev_nr_hosp)]
test_stratum <- pmin(test_patient$prev_nr_hosp, cap_hosp)
set.seed(1)
simulate_next_event(
  cox_models  = cox_models,
  origin      = test_patient$state,
  trt_chr     = as.character(test_patient$trt),
  hosp_capped = test_stratum
)
test_current_event_time <- as.numeric(test_patient$tstop - test_patient$tstart)
simulate_next_event(
  cox_models         = cox_models,
  origin             = test_patient$state,
  trt_chr            = as.character(test_patient$trt),
  hosp_capped        = test_stratum,
  is_first_event     = TRUE,
  current_event_time = test_current_event_time
)

################################################################################
### Simulate remaining trajectory for administratively censored patients
################################################################################
# define follow-up days
baseline[, followup_days := as.numeric(followup_end - visit_date)]

# horizon = 730
baseline[, horizon_days_2y := horizon]

# censored = alive at end of raw follow-up, but follow-up stopped before the
# full 2-year horizon for a reason other than death (admin cutoff, loss to
# follow-up, etc.)
baseline[, censor_reason := fcase(
  event_death_2y == 1,
  "death",
  followup_days >= horizon_days_2y,
  "reached_horizon",
  default = "censored"
)]
table(baseline$censor_reason)

# select LOPNR of those patients that were censored
censored_id <- baseline[censor_reason == "censored", LOPNR]

# select data from patients that were censored
last_obs_dt <- baseline_long[order(LOPNR, tstart)][LOPNR %in% censored_id, .SD[.N], by = LOPNR]

# add follow-up days and horizon_days_2y to baseline_long
last_obs_dt <- merge(last_obs_dt, baseline[, .(LOPNR, followup_days, horizon_days_2y)], by = "LOPNR")

# time already spent in the current state as of censoring (days since
# entering current_state, not since decision_date) - needed so the first
# simulated transition is left-truncated at this elapsed time instead of
# incorrectly resetting the clock to zero for a state the patient has
# already been in for a while
last_obs_dt[, time_in_current_state := followup_days - as.numeric(tstart - decision_date)]

# simulate one censored patient's remaining trajectory, from their censoring
# point (start_time, start_state) forward to horizon_days_2y or simulated death.
# time units throughout are days since decision_date (not calendar dates),
# to match simulate_cause_time()'s bh_s$time scale
simulate_patient_trajectory <- function(id,
                                        start_state,
                                        start_time,
                                        time_in_current_state,
                                        prev_nr_hosp,
                                        trt,
                                        horizon_days_2y,
                                        cox_models,
                                        cap,
                                        log_env = NULL,
                                        max_transitions = 200) {
  # initialize state, time, number of hospitalizations and simulated rows.
  # current_time = "now", the point we are simulating forward from
  # (the censoring time on the first iteration). entry_time = the absolute
  # time (days since decision_date) at which current_state was actually
  # entered - equal to current_time for every state except the very first
  # one, where the patient was picked up mid-sojourn at censoring.
  current_state <- start_state
  current_time  <- start_time
  entry_time    <- start_time - time_in_current_state
  current_hosp  <- prev_nr_hosp
  trt_chr <- as.character(trt)  # match trt once, not one data.table per step
  sim_rows <- vector("list", max_transitions)
  n_rows <- 0L
  
  # loop until end of follow-up or simualted death
  for (i in seq_len(max_transitions)) {
    hosp_capped <- pmin(current_hosp, cap)
    
    # simulate next transition (precomputed basehaz/lp - see fit_cause_specific_cox).
    # Only the FIRST simulated transition is left-truncated at the time
    # already spent in current_state at censoring (is_first_event = TRUE):
    # the patient is picked up mid-sojourn, so the draw must condition on
    # having already survived without an event up to that point, rather
    # than resetting the clock to zero as if they had just entered the
    # state. Every subsequent transition (i > 1) genuinely does start a
    # fresh state entry, so the clock correctly resets to zero for those.
    step <- simulate_next_event(
      cox_models = cox_models,
      origin = current_state,
      trt_chr = trt_chr,
      hosp_capped = hosp_capped,
      is_first_event = (i == 1),
      current_event_time = if (i == 1) time_in_current_state else NULL,
      log_env = log_env
    )
    
    # step$time is measured on the state's own clock, i.e. time since
    # entry_time (not since current_time) - so the absolute event time is
    # entry_time + step$time. For i > 1, entry_time == current_time, so
    # this matches the original (correct) behavior exactly.
    new_time <- entry_time + step$time
    
    # stop if the next event has infinite time or falls beyond the horizon
    if (is.infinite(step$time) || new_time >= horizon_days_2y) {
      # no further transition simulated before the horizon: patient stays
      # in current_state for the rest of follow-up, then stop
      n_rows <- n_rows + 1L
      sim_rows[[n_rows]] <- data.table(
        LOPNR = id,
        tstart = current_time,
        tstop = horizon_days_2y,
        state = current_state
      )
      break
    }
    
    # continue if next event is before horizon
    n_rows <- n_rows + 1L
    sim_rows[[n_rows]] <- data.table(
      LOPNR = id,
      tstart = current_time,
      tstop = new_time,
      state = current_state
    )
    
    # add a zero-width Death row (matches missing_states_dt's convention)
    # before stopping, so simulated_dt actually records that death occurred
    if (step$destination == "Death") {
      n_rows <- n_rows + 1L
      sim_rows[[n_rows]] <- data.table(
        LOPNR = id,
        tstart = new_time,
        tstop = new_time,
        state = "Death"
      )
      break
    }
    
    # only a genuine "Hospitalization" destination increments the
    # running count used to stratify the next model - matches prev_nr_hosp
    # in states_dt, which is driven by hosp_num
    if (step$destination == "Hospitalization")
      current_hosp <- current_hosp + 1
    
    current_state <- step$destination
    current_time  <- new_time
    entry_time    <- new_time  # every subsequent state is entered fresh
  }
  
  rbindlist(sim_rows[seq_len(n_rows)])
}

# simulate every censored patient's remaining trajectory. Map() walks
# last_obs_dt's own columns directly instead of looping over censored_id and
# re-filtering last_obs_dt (last_obs_dt[LOPNR == id]) for every single
# patient, which rescans the whole table each time

# Report, across the whole simulation, how many candidate transition draws
# came back Inf (i.e. the cumulative hazard never reached the simulated
# threshold), and how many times two or more candidates tied on the exact
# same finite simulated time - see the check right after simulated_dt is
# built below
sim_time_log <- new.env()
sim_time_log$n_candidates <- 0L
sim_time_log$n_inf <- 0L
sim_time_log$n_inf_by_transition <- list()
sim_time_log$n_ties <- 0L
sim_time_log$n_steps <- 0L

set.seed(1)
simulated_list <- Map(
  simulate_patient_trajectory,
  id                     = last_obs_dt$LOPNR,
  start_state            = last_obs_dt$state,
  start_time             = last_obs_dt$followup_days,
  time_in_current_state  = last_obs_dt$time_in_current_state,
  prev_nr_hosp           = last_obs_dt$prev_nr_hosp,
  trt                    = last_obs_dt$trt,
  horizon_days_2y        = last_obs_dt$horizon_days_2y,
  MoreArgs               = list(cox_models = cox_models, cap = cap_hosp, log_env = sim_time_log)
)
simulated_dt <- rbindlist(simulated_list)
simulated_dt[, days_in_row := as.numeric(tstop - tstart)]
simulated_dt <- merge(simulated_dt, unique(baseline_long[, .(LOPNR, trt)]), by = "LOPNR")

################################################################################
### Check: how many simulated candidate event times were Infinite
################################################################################
# An Inf draw means the simulated cumulative-hazard threshold was never
# reached within the observed range of that transition's baseline hazard -
# i.e. the model is being asked to extrapolate beyond the longest sojourn
# time ever actually observed for that transition. This is expected
# occasionally (that candidate destination simply doesn't happen before a
# competing one does), but a high rate for a specific transition would
# suggest that transition's baseline hazard is sparse/short and simulated
# trajectories relying on it may be unreliable.
cat(sprintf(
  "Simulated %d candidate transition draws in total across all censored patients; %d (%.1f%%) were Inf.\n",
  sim_time_log$n_candidates,
  sim_time_log$n_inf,
  100 * sim_time_log$n_inf / sim_time_log$n_candidates
))

inf_by_transition_dt <- data.table(
  transition  = names(sim_time_log$n_inf_by_transition),
  n_inf       = unlist(sim_time_log$n_inf_by_transition, use.names = FALSE)
)
setorder(inf_by_transition_dt, -n_inf)
inf_by_transition_dt

################################################################################
### Check: how many simulated transitions resulted in a tie
################################################################################
# A tie = two or more candidate transitions landed on the exact same
# finite simulated time within a single step. which.min() in
# simulate_next_event() currently breaks ties by silently picking
# whichever candidate happens to come first in cox_models, rather than at
# random or weighted by hazard - this just quantifies how often that
# (currently arbitrary) tie-break actually gets invoked.
cat(sprintf(
  "Simulated %d transition steps in total across all censored patients; %d (%.1f%%) had a tie.\n",
  sim_time_log$n_steps,
  sim_time_log$n_ties,
  100 * sim_time_log$n_ties / sim_time_log$n_steps
))

# Check if the event times of the new dataset matches the one in the real
# data (per transition): compare simulated time-in-state (for transitions
# that actually occurred during simulation) against the real observed
# time-in-state (event == 1) that each cause-specific model was fit on
setorder(simulated_dt, LOPNR, tstart)
simulated_dt[, next_state := shift(state, type = "lead"), by = LOPNR]
transition_check_dt <- rbindlist(lapply(names(cox_models), function(nm) {
  origin_destination <- strsplit(nm, " -> ")[[1]]
  origin      <- origin_destination[1]
  destination <- origin_destination[2]
  
  sim_times <- simulated_dt[state == origin &
                              next_state == destination, days_in_row]
  
  data_name  <- paste0(gsub(" ", "_", tolower(origin)), "_dt")
  real_times <- cox_models[[nm]][[data_name]][event == 1, time_in_state]
  
  # KS test on the two time-in-state distributions - low p-value flags a
  # transition where simulated and real event times diverge meaningfully
  ks_p <- if (length(sim_times) >= 2 && length(real_times) >= 2) {
    suppressWarnings(ks.test(sim_times, real_times)$p.value)
  } else {
    NA_real_
  }
  
  data.table(
    transition  = nm,
    n_sim       = length(sim_times),
    n_real      = length(real_times),
    mean_sim    = if (length(sim_times) > 0)
      mean(sim_times)
    else
      NA_real_,
    mean_real   = if (length(real_times) > 0)
      mean(real_times)
    else
      NA_real_,
    median_sim  = if (length(sim_times) > 0)
      median(sim_times)
    else
      NA_real_,
    median_real = if (length(real_times) > 0)
      median(real_times)
    else
      NA_real_,
    ks_p_value  = ks_p
  )
}))
transition_check_dt

# combine each censored patient's OBSERVED days (up to censoring) with their
# SIMULATED days (censoring -> horizon/simulated death) into one completed
# trajectory, in the same day-since-decision_date units as simulated_dt
observed_part <- baseline_long[LOPNR %in% censored_id, .(
  LOPNR,
  tstart = as.numeric(tstart - decision_date),
  tstop  = as.numeric(tstop - decision_date),
  state,
  trt,
  source = "observed"
)]
observed_part[, days_in_row := tstop - tstart]

simulated_part <- simulated_dt[, .(LOPNR, tstart, tstop, state, trt, source = "simulated", days_in_row)]

# non-censored patients (died, or reached the horizon) already have a
# complete observed trajectory - no simulation needed for them
non_censored_part <- baseline_long[!(LOPNR %in% censored_id), .(
  LOPNR,
  tstart = as.numeric(tstart - decision_date),
  tstop  = as.numeric(tstop - decision_date),
  state,
  trt,
  source = "observed"
)]
non_censored_part[, days_in_row := tstop - tstart]

# whole cohort: non-censored as observed + censored (observed part + simulated part)
# days_in_row is carried through from each part rather than recomputed here
full_dt <- rbindlist(list(non_censored_part, observed_part, simulated_part))
setorder(full_dt, LOPNR, tstart)

################################################################################
### Unadjusted vs confounding-only vs censoring-only vs both corrected
################################################################################
# simple descriptive tally, not adjusted for confounding or censoring
baseline_long[, days_in_row := as.numeric(tstop - tstart)]

# total days per state, across the whole cohort
baseline_long[, .(total_days = sum(days_in_row)), by = state]

# days per patient per state, raw (censored at followup_end, unweighted)
days_per_patient <- dcast(
  baseline_long,
  LOPNR + trt ~ state,
  value.var = "days_in_row",
  fun.aggregate = sum,
  fill = 0
)

# total days per patient per state, whole cohort, using simulated
# completions in place of administrative censoring
days_per_patient_imputed <- dcast(
  full_dt,
  LOPNR + trt ~ state,
  value.var = "days_in_row",
  fun.aggregate = sum,
  fill = 0
)

# one IPTW weight per patient (same weight column the Cox models use)
iptw_dt <- unique(baseline_long[, .(LOPNR, sw_IPTW)])

state_cols <- setdiff(names(days_per_patient), c("LOPNR", "trt"))
weighted_state_means <- function(dt, weight_dt, state_cols) {
  dt <- merge(dt, weight_dt, by = "LOPNR")
  dt[, lapply(.SD, weighted.mean, w = sw_IPTW), by = trt, .SDcols = state_cols]
}

# full 2x2 factorial: corrected for confounding (yes/no) x censoring (yes/no)
# 1. neither: raw days, unweighted, censored at followup_end
estimate_unadjusted <- days_per_patient[, lapply(.SD, mean), by = trt, .SDcols = state_cols]

# 2. confounding only: IPTW-weighted, still censored at followup_end
estimate_confounding_only <- weighted_state_means(days_per_patient, iptw_dt, state_cols)

# 3. censoring only: unweighted, simulated completions
estimate_censoring_only <- days_per_patient_imputed[, lapply(.SD, mean), by = trt, .SDcols = state_cols]

# 4. both: IPTW-weighted, simulated completions
estimate_confounding_and_censoring <- weighted_state_means(days_per_patient_imputed, iptw_dt, state_cols)

# write all four estimates to one workbook, one sheet each. Rounded to
# whole days for this table only - the underlying estimate_* objects keep
# full precision (and keep the Death column), since
# estimate_confounding_and_censoring is also used elsewhere (at 1 decimal,
# in months, and without Death - see state_levels below) to annotate the
# figure. table_state_cols excludes Death here so the table only reports
# days spent alive in each of the other four states.
table_state_cols <- setdiff(state_cols, "Death")

round_state_cols <- function(dt, cols) {
  dt <- copy(dt)[, c("trt", cols), with = FALSE]
  dt[, (cols) := lapply(.SD, round, digits = 0), .SDcols = cols]
  dt
}

write.xlsx(
  rbind(
    data.table(trt = "Unadjusted"),
    round_state_cols(estimate_unadjusted, table_state_cols),
    data.table(trt = "Adjusted for confounding only"),
    round_state_cols(estimate_confounding_only, table_state_cols),
    data.table(trt = "Adjusted for censoring only"),
    round_state_cols(estimate_censoring_only, table_state_cols),
    data.table(trt = "Adjusted for confounding and censoring"),
    round_state_cols(estimate_confounding_and_censoring, table_state_cols),
    fill = TRUE
  ),
  file = paste0(results_path, "Supplemental/Table_S_home_time_estimates.xlsx")
)

################################################################################
### Figure: state occupancy over time (observed + simulated), by trt
################################################################################
# analog of file 9's stacked-area panel, generalized from 3 events to all 5
# states, and covering the full 2-year horizon for everyone via full_dt's
# simulated completions rather than stopping at each patient's censoring point

# one row per patient per day, day 0 to day 729 (730 days = the 2-year
# horizon; excluding day 730 avoids double-counting the last day)
time_grid <- seq(0, horizon - 1, by = 1)

setkey(full_dt, LOPNR, tstart)
grid_dt <- CJ(LOPNR = unique(full_dt$LOPNR), time = time_grid)

# state each patient was in on each day (roll = TRUE carries the last
# known state forward until the next transition, incl. through death)
daily_state_dt <- full_dt[grid_dt, on = .(LOPNR, tstart = time), roll = TRUE]
setnames(daily_state_dt, "tstart", "time")

# attach IPTW weights so we count weighted, not raw, patients per state/day
daily_state_dt <- merge(daily_state_dt, iptw_dt, by = "LOPNR")

# total IPTW weight per arm - denominator to turn weighted sums into means, i.e.,
# total_weight_trt = sum(w_i)
total_weight_by_trt <- unique(daily_state_dt[, .(LOPNR, trt, sw_IPTW)])[, .(total_weight = sum(sw_IPTW)), by = trt]

# weighted number of patients in each state, per day (snapshot, not cumulative), i.e.,
# weighted_n(t) = Σᵢ wᵢ × 1[patient i is in state s on day t]
state_prob_dt <- daily_state_dt[, .(weighted_n = sum(sw_IPTW)), by = .(time, trt, state)]
setorder(state_prob_dt, trt, state, time)

# running total of weighted patient-days spent in each state so far, i.e., 
# cumulative_weighted_days(T) = Σₜ₌₀ᵀ weighted_n(t) = Σᵢ wᵢ × (Σₜ₌₀ᵀ 1[patient i in state s on day t ])= Σᵢ wᵢ × xᵢ
state_prob_dt[, cumulative_weighted_days := cumsum(weighted_n), by = .(trt, state)]

# left panel: normalize by arm weight -> weighted mean days per patient (comparable
# across arms, same scale as estimate_confounding_and_censoring)
state_prob_dt <- merge(state_prob_dt, total_weight_by_trt, by = "trt")
state_prob_dt[, mean_cumulative_days := cumulative_weighted_days / total_weight]

# sum across all 5 states, per day, as the denominator for proportions
state_prob_dt[, total_cumulative_days := sum(cumulative_weighted_days), by = .(trt, time)]

# right panel: each state's share of cumulative time so far - sums to 1.0 at every day
state_prob_dt[, prop_cumulative_days := cumulative_weighted_days / total_cumulative_days]

# stacking order: "At home" at the bottom (best outcome), "Death" at the
# top - flip state_levels if your ggplot2 version stacks the other way.
# Death IS shown in the figure's curves (unlike the text annotation table -
# see annotation_state_levels below, which excludes it there only).
state_levels <- c("At home", "Hospitalization", "HD", "PD", "Death")
state_prob_dt[, state := factor(state, levels = state_levels)]

# reuses the project's manual_colors if it's in scope (it isn't
# sourced by this file currently); otherwise falls back to a plain 5-color
# palette - swap in manual_colors[1:5] directly for project-wide consistency.
# NOTE: this previously assigned 6 colors ("#FF7F00" plus a duplicate
# manual_colors[1]) to 5 names - a pre-existing length mismatch. Fixed here
# to exactly 5 colors for the 5 states (Death now uses "#FF7F00" alone).
state_colors <- setNames(c(manual_colors[c(3, 2, 4, 1)], "#FF7F00"), state_levels)

################################################################################
### Figure: state probability (left) and state occupancy (right), by trt
################################################################################
# one panel per (trt, plot type), built manually (rather than facet_wrap)
# so the layout can be controlled directly: dialysis row on top,
# conservative management below; probability left, occupancy right
make_state_panel <- function(state_prob_dt,
                             state_colors,
                             trt_value,
                             type = c("cumulative", "stacked"),
                             show_legend = TRUE) {
  type <- match.arg(type)
  dt <- state_prob_dt[trt == trt_value]
  
  p <- if (type == "cumulative") {
    ggplot2::ggplot(dt,
                    ggplot2::aes(
                      x = time / 365,
                      y = mean_cumulative_days,
                      color = state
                    )) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_color_manual(values = state_colors, name = "State")
  } else {
    ggplot2::ggplot(dt,
                    ggplot2::aes(
                      x = time / 365,
                      y = prop_cumulative_days,
                      fill = state
                    )) +
      ggplot2::geom_area(
        color = "black",
        linewidth = 0.2,
        alpha = 0.85,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = state_colors, name = "State") +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        breaks = seq(0, 1, by = 0.1)
      )
  }
  
  p <- p +
    ggplot2::scale_x_continuous(breaks = seq(0, 2, by = 0.5)) +
    ggplot2::labs(
      x = "Time (years)",
      y = if (type  == "cumulative")
        "Mean cumulative days per patient"
      else
        "Proportion of cumulative person-time"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# only one panel keeps its legend - probability (color) and occupancy
# (fill) use the same states/colors, so showing both would just duplicate
# the same "State" key twice
p_cum_dialysis    <- make_state_panel(state_prob_dt, state_colors, 1, "cumulative", show_legend = TRUE) +
  ggplot2::labs(title = "Dialysis")
p_stack_dialysis  <- make_state_panel(state_prob_dt, state_colors, 1, "stacked", show_legend = FALSE)
p_cum_cm          <- make_state_panel(state_prob_dt, state_colors, 0, "cumulative", show_legend = FALSE) +
  ggplot2::labs(title = "Conservative management")
p_stack_cm        <- make_state_panel(state_prob_dt, state_colors, 0, "stacked", show_legend = FALSE)

################################################################################
### Annotate the state-probability panels
################################################################################
# annotate one probability panel with the 4-column state/days/months
# table (At home, Hospitalization, HD, PD), using the
# confounding+censoring-corrected estimate throughout. Death is excluded
# from this text annotation only (via annotation_state_levels below) - it
# still appears as a curve/area in the figure itself (state_levels).
annotation_state_levels <- setdiff(state_levels, "Death")

annotate_state_probability <- function(p,
                                       trt_value,
                                       annotation_state_levels,
                                       estimate_confounding_and_censoring) {
  days_vals   <- as.numeric(estimate_confounding_and_censoring[trt == trt_value, ..annotation_state_levels])
  months_vals <- days_vals / 30.5
  
  p +
    ggplot2::annotate(
      "text",
      x = 0,
      y = max(state_prob_dt$mean_cumulative_days),
      hjust = 0,
      vjust = 1,
      size = 3,
      label = paste0("At 2 years\n", paste(annotation_state_levels, collapse = "\n"))
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = max(state_prob_dt$mean_cumulative_days),
      hjust = 0,
      vjust = 1,
      size = 3,
      label = paste0("days\n", paste(fmt(days_vals, 0), collapse = "\n"))
    ) +
    ggplot2::annotate(
      "text",
      x = 0.8,
      y = max(state_prob_dt$mean_cumulative_days),
      hjust = 0,
      vjust = 1,
      size = 3,
      label = paste0("months\n", paste(fmt(months_vals, 1), collapse = "\n"))
    )
}

p_cum_dialysis <- annotate_state_probability(
  p = p_cum_dialysis,
  trt_value = 1,
  annotation_state_levels,
  estimate_confounding_and_censoring
)
p_cum_cm <- annotate_state_probability(
  p = p_cum_cm,
  trt_value = 0,
  annotation_state_levels,
  estimate_confounding_and_censoring
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
  filename = paste0(results_path, "Main/Figure_3.png"),
  width = 10,
  height = 8,
  dpi = 300
)
