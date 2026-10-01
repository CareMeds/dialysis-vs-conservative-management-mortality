################################################################################
### In-center days
################################################################################

# Days spent in-center per state episode: Hospitalization = all days;
# HD = days * sessions/7; PD = training days + one visit per interval.
# in_center_params: named list with hd_sessions_per_week, pd_training_days and
# pd_visit_interval (defined in the analysis script).
in_center_days <- function(state, days, in_center_params) {
  # in-center days by state
  fcase(
    state == "Hospitalization",
    days,
    state == "HD",
    days * in_center_params$hd_sessions_per_week / 7,
    state == "PD",
    pmin(days, in_center_params$pd_training_days) +
      pmax(days - in_center_params$pd_training_days, 0) / in_center_params$pd_visit_interval,
    default = 0
  )
}

# Sums the row-level in_center column per patient onto the days-per-state table.
# Also adds "Time at home" = days alive (all rows except Death) - in-center days, so it
# includes days at home on dialysis (HD/PD outside the in-center sessions/visits), and its
# components / those of the in-center days:
#   "Home without dialysis"  = days in the "At home" state
#   "Home on dialysis"       = HD + PD days - in-center dialysis days
#   "In-center hospitalized" = in-center days in the Hospitalization state (all hospital days)
#   "In-center dialysis"     = in-center days in the HD / PD states (sessions / visits)
# so Home without + Home on dialysis = Time at home, hospitalized + dialysis = In-center.
add_in_center_per_patient <- function(days_dt, rows_dt) {
  # indicate which part is after dialysis
  # (NA_real_ for patients who never start dialysis: min() of nothing is Inf, which
  # data.table coerces to NA when tstart is an integer column, turning those
  # patients' after-dialysis sums, and so the arm means, into NA)
  rows_dt[, dialysis_start := {
    s <- tstart[state %in% c("HD", "PD")]
    if (length(s) > 0) as.numeric(min(s)) else NA_real_
  }, by = LOPNR]
  rows_dt[, after_dialysis := !is.na(dialysis_start) & as.numeric(tstart) >= dialysis_start]
  
  # sum in-center days, days after dialysis and in-center days after dialysis per patient
  ic <- rows_dt[, .(
    `In-center` = sum(in_center),
    `Time at home` = sum(days_in_row[state != "Death"]) - sum(in_center),
    `Home without dialysis` = sum(days_in_row[state == "At home"]),
    `Home on dialysis` = sum(days_in_row[state %in% c("HD", "PD")]) -
      sum(in_center[state %in% c("HD", "PD")]),
    `In-center hospitalized` = sum(in_center[state == "Hospitalization"]),
    `In-center dialysis` = sum(in_center[state %in% c("HD", "PD")]),
    `Days after dialysis` = sum(days_in_row[after_dialysis & state != "Death"]),
    `In-center after dialysis` = sum(in_center[after_dialysis]),
    `Hospital after dialysis` = sum(days_in_row[after_dialysis & state == "Hospitalization"])
  ), by = .(LOPNR, trt)]
  
  # attach the per-patient sums to the days-per-state table
  merge(days_dt, ic, by = c("LOPNR", "trt"), all.x = TRUE)
}

################################################################################
### Data preparation
################################################################################

# Numbers consecutive runs of TRUE in flag_col within each patient.
number_runs <- function(dt, flag_col, episode_col, id_col = "LOPNR") {
  # name of the temporary run-id column
  run_id_col <- paste0(episode_col, "_run_id")
  
  # temporary run-id for every consecutive stretch of the flag
  dt[, (run_id_col) := rleid(get(flag_col)), by = id_col]
  
  # one row per TRUE run - these are what get numbered
  runs <- unique(dt[get(flag_col) == TRUE, c(id_col, run_id_col), with = FALSE])
  
  # sort runs by patient and run id
  setorderv(runs, c(id_col, run_id_col))
  
  # number those runs 1, 2, 3... in chronological order, per patient
  runs[, (episode_col) := seq_len(.N), by = id_col]
  
  # join episode number back; FALSE rows stay NA (absent from runs)
  dt[runs, (episode_col) := get(paste0("i.", episode_col)), on = c(id_col, run_id_col)]
  
  # drop the temporary run-id column
  dt[, (run_id_col) := NULL]
  
  # return the table invisibly (modified by reference)
  invisible(dt)
}

# Episodes and patients per transition and arm; include_trt = FALSE for dialysis_only transitions.
count_transitions_by_trt <- function(baseline_long,
                                     transitions,
                                     dialysis_only) {
  # count episodes and patients for every transition
  transition_n_dt <- rbindlist(lapply(transitions, function(tr) {
    # episodes and unique patients per arm
    origin <- tr[1]
    destination <- tr[2]
    dt <- baseline_long[state == origin &
                          next_state == destination, .(n_episodes = .N, n_patients = uniqueN(LOPNR)), by = trt]
    
    # label the rows with the transition name
    dt[, transition := paste(origin, "->", destination)]
    
    # return counts for this transition
    dt
  }))
  
  # one row per transition, counts per arm as columns
  transition_n_dt <- dcast(
    transition_n_dt,
    transition ~ trt,
    value.var = c("n_episodes", "n_patients"),
    fill = 0
  )
  
  # names of all transitions
  transition_names <- sapply(transitions, function(tr)
    paste(tr[1], "->", tr[2]))
  
  # flag transitions that include trt in the strata
  include_trt_dt <- data.table(
    transition  = transition_names,
    include_trt = !transition_names %in% dialysis_only
  )
  
  # combine counts and include_trt flag
  merge(transition_n_dt, include_trt_dt, by = "transition")
}

# Builds the long-format state-episode data (states, next_state, time in state,
# in-center days) and the pieces used for simulation. Done once per patient
# sample; resampled by LOPNR for each bootstrap replicate.
build_baseline_long <- function(long_cohort_hosp,
                                long_cohort_krt,
                                baseline,
                                transitions,
                                dialysis_only,
                                horizon) {
  # Prepare the data ----------------------------------------------------------#
  # one row per patient: static info that never changes over time
  patient_static_dt <- unique(long_cohort_hosp[, .(LOPNR, trt, decision_date, DODSDAT, followup_end)])
  
  # Hospitalizations ----------------------------------------------------------#
  # hospitalization episodes, one row per episode
  hosp_dt <- long_cohort_hosp[, .(LOPNR, hosp_start, hosp_stop)]
  
  # same-day hospitalizations last 1 day: extend hosp_stop (unless at decision_date)
  hosp_dt <- patient_static_dt[, .(LOPNR, decision_date)][hosp_dt, on = "LOPNR"]
  
  # extend same-day hospitalizations by one day
  hosp_dt[hosp_start == hosp_stop &
            hosp_start != decision_date, hosp_stop := hosp_stop + 1]
  
  # drop the helper column
  hosp_dt[, decision_date := NULL]
  
  # number hospitalization episodes per patient (1st, 2nd, 3rd, ...)
  setorder(hosp_dt, LOPNR, hosp_start)
  
  # number hospitalizations per patient
  hosp_dt[, hosp_num := seq_len(.N), by = LOPNR]
  
  # Start of KRT --------------------------------------------------------------#
  # dialysis (KRT) periods
  krt_dt <- long_cohort_krt[, .(LOPNR, krt_modality, krt_start, krt_stop)]
  
  # number KRT periods (increments on every switch)
  setorder(krt_dt, LOPNR, krt_start)
  
  # number KRT periods per patient
  krt_dt[, krt_num := seq_len(.N), by = LOPNR]
  
  # Break points --------------------------------------------------------------#
  # stack every candidate breakpoint into one table
  break_points_dt <- rbindlist(
    list(
      patient_static_dt[, .(LOPNR, break_point = decision_date)],
      patient_static_dt[!is.na(DODSDAT) &
                          DODSDAT <= followup_end, .(LOPNR, break_point = DODSDAT)],
      patient_static_dt[, .(LOPNR, break_point = followup_end)],
      krt_dt[, .(LOPNR, break_point = krt_start)],
      krt_dt[, .(LOPNR, break_point = krt_stop)],
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
  
  # dates needing a zero-width row: death, KRT starting on the last follow-up day
  zero_width_dates <- unique(rbindlist(list(patient_static_dt[!is.na(DODSDAT) &
                                                                DODSDAT <= followup_end, .(LOPNR, break_point = DODSDAT)], krt_dt[krt_start == krt_stop, .(LOPNR, break_point = krt_start)])))
  # zero-width rows for those dates
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
  
  # Add hospitalization dates to long dt ----------------------------------------#
  # attach the hospitalization episode (if any) each interval falls inside
  match_idx <- hosp_dt[break_points_dt, on = .(LOPNR, hosp_start <= tstart, hosp_stop >= tstop), which = TRUE]
  
  # attach hospitalization info to each interval
  break_points_dt[, `:=`(
    hosp_start = hosp_dt$hosp_start[match_idx],
    hosp_stop  = hosp_dt$hosp_stop[match_idx],
    hosp_num   = hosp_dt$hosp_num[match_idx]
  )]
  
  # Add KRT dates to long dt ----------------------------------------------------#
  # attach the active KRT modality (if any)
  match_idx <- krt_dt[break_points_dt, on = .(LOPNR, krt_start <= tstart, krt_stop >= tstop), which = TRUE]
  
  # attach KRT info to each interval
  break_points_dt[, `:=`(
    krt_modality  = krt_dt$krt_modality[match_idx],
    krt_start     = krt_dt$krt_start[match_idx],
    krt_stop      = krt_dt$krt_stop[match_idx],
    krt_num       = krt_dt$krt_num[match_idx]
  )]
  
  # which hospitalization (if any) a patient's first HD/PD start fell inside
  first_krt_dt <- krt_dt[krt_modality %in% c("HD", "PD"), .(krt_start = min(krt_start)), by = LOPNR]
  
  # hospitalization in which the first HD/PD started
  initiating_hosp_dt <- hosp_dt[first_krt_dt, on = .(LOPNR, hosp_start <= krt_start, hosp_stop >= krt_start), .(LOPNR, initiating_hosp_num = hosp_num)]
  
  # keep one row per patient
  initiating_hosp_dt <- unique(initiating_hosp_dt[!is.na(initiating_hosp_num)])
  
  # add the initiating hospitalization number to the intervals
  break_points_dt[initiating_hosp_dt, initiating_hosp_num := i.initiating_hosp_num, on = "LOPNR"]
  
  # Assign states -------------------------------------------------------------#
  # priority: Death > HD/PD in initiating hospitalization > Hospitalization > HD/PD > At home
  states_dt <- break_points_dt[, state := fcase(
    !is.na(DODSDAT) &
      tstart == DODSDAT,
    "Death",!is.na(hosp_num) &
      !is.na(initiating_hosp_num) &
      hosp_num == initiating_hosp_num &
      !is.na(krt_modality) &
      krt_modality == "HD",
    "HD",!is.na(hosp_num) &
      !is.na(initiating_hosp_num) &
      hosp_num == initiating_hosp_num &
      !is.na(krt_modality) &
      krt_modality == "PD",
    "PD",!is.na(hosp_start),
    "Hospitalization",!is.na(krt_modality) &
      krt_modality == "HD",
    "HD",!is.na(krt_modality) &
      krt_modality == "PD",
    "PD",
    default = "At home"
  )]
  
  # merge consecutive same-state rows split by an unrelated breakpoint
  states_dt[, merge_state_id := rleid(state), by = LOPNR]
  
  # carry forward every column except the grouping keys and tstart/tstop
  merge_cols <- setdiff(names(states_dt),
                        c("LOPNR", "tstart", "tstop", "merge_state_id"))
  
  # collapse each run into one row (keep hosp_* from anywhere in the run)
  preserve_cols <- c("hosp_start", "hosp_stop", "hosp_num")
  # collapse each run of same-state rows into one row
  states_dt <- states_dt[, c(list(tstart = min(tstart), tstop = max(tstop)),
                             Map(function(x, nm)
                               if (nm %in% preserve_cols)
                                 x[!is.na(x)][1]
                               else
                                 x[.N], .SD, names(.SD))), by = .(LOPNR, merge_state_id), .SDcols = merge_cols]
  # drop the run id
  states_dt[, merge_state_id := NULL]
  
  # Number dialysis episodes (HD and PD combined into one counter) ------------#
  states_dt[, is_dialysis := state %in% c("HD", "PD")]
  number_runs(states_dt, flag_col = "is_dialysis", episode_col = "dialysis_num")
  
  # Number "at home" episodes -------------------------------------------------#
  states_dt[, is_home := state == "At home"]
  number_runs(states_dt, flag_col = "is_home", episode_col = "home_num")
  
  # Episode-count distribution ------------------------------------------------#
  # max hospitalization episode number per patient
  max_hosp_per_patient <- states_dt[!is.na(hosp_num), .(max_n = max(hosp_num)), by = LOPNR]
  
  # choose a sensible cap for prev_nr_hosp_capped
  cap_hosp <- as.numeric(quantile(max_hosp_per_patient$max_n, 0.9))
  
  # remove helper columns
  states_dt[, `:=`
            (
              is_home = NULL,
              home_num = NULL,
              dialysis_num = NULL,
              is_dialysis = NULL
            )]
  
  # Attach baseline covariates to build the final long-format dataset ---------#
  # patients without hospitalizations: add back as At home (+ Death if died)
  missing_id <- setdiff(baseline$LOPNR, states_dt$LOPNR)
  
  # patients without any hospitalization or KRT rows
  missing_dt <- baseline[LOPNR %in% missing_id, .(LOPNR, decision_date = visit_date, followup_end, event_death_2y)]
  
  # give them one At home row (and a Death row if they died)
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
  
  # add them to the states table
  states_dt <- rbindlist(list(states_dt, missing_states_dt), fill = TRUE)
  
  # sort chronologically per patient
  setorder(states_dt, LOPNR, tstart)
  
  # running count of hospitalizations before this row: carry hosp_num forward (locf)
  states_dt[, prev_nr_hosp := nafill(shift(hosp_num), type = "locf"), by = LOPNR]
  
  # no admissions yet for these rows
  states_dt[is.na(prev_nr_hosp), prev_nr_hosp := 0L]
  
  # add baseline covariates onto every row of states_dt
  states_only_cols <- setdiff(names(states_dt), names(baseline))
  
  # add baseline covariates to every row
  baseline_long <- baseline[states_dt[, c("LOPNR", states_only_cols), with = FALSE], on = "LOPNR"]
  
  # destination state after each episode
  baseline_long[, next_state := shift(state, type = "lead"), by = LOPNR]
  
  # time at risk in this episode - clock resets at every state entry
  baseline_long[, time_in_state := as.numeric(tstop - tstart)]
  
  # transition counts need no model; computed once here
  transition_n_dt <- count_transitions_by_trt(baseline_long, transitions, dialysis_only)
  
  # Data for simulating censored patients
  # follow-up days from the visit date
  baseline[, followup_days := as.numeric(followup_end - visit_date)]
  
  # horizon = 730
  baseline[, horizon_days_2y := horizon]
  
  # censored = alive at end of follow-up but before the 2-year horizon
  baseline[, censor_reason := fcase(
    event_death_2y == 1,
    "death",
    followup_days >= horizon_days_2y,
    "reached_horizon",
    default = "censored"
  )]
  
  # select LOPNR of those patients that were censored
  censored_id <- baseline[censor_reason == "censored", LOPNR]
  
  # select data from patients that were censored
  last_obs_dt <- baseline_long[order(LOPNR, tstart)][LOPNR %in% censored_id, .SD[.N], by = LOPNR]
  
  # add follow-up days and horizon_days_2y to baseline_long
  last_obs_dt <- merge(last_obs_dt, baseline[, .(LOPNR, followup_days, horizon_days_2y)], by = "LOPNR")
  
  # time already spent in current state at censoring (for left truncation)
  last_obs_dt[, time_in_current_state := followup_days - as.numeric(tstart - decision_date)]
  
  # observed parts don't depend on the simulation: built once
  observed_part <- baseline_long[LOPNR %in% censored_id, .(
    LOPNR,
    tstart = as.numeric(tstart - decision_date),
    tstop  = as.numeric(tstop - decision_date),
    state,
    trt,
    source = "observed"
  )]
  
  # days in each observed row
  observed_part[, days_in_row := tstop - tstart]
  
  # non-censored patients are fully observed: no simulation needed
  non_censored_part <- baseline_long[!(LOPNR %in% censored_id), .(
    LOPNR,
    tstart = as.numeric(tstart - decision_date),
    tstop  = as.numeric(tstop - decision_date),
    state,
    trt,
    source = "observed"
  )]
  
  # days in each non-censored row
  non_censored_part[, days_in_row := tstop - tstart]
  
  # return all pieces needed for simulation
  list(
    baseline_long     = baseline_long,
    cap_hosp          = cap_hosp,
    last_obs_dt       = last_obs_dt,
    observed_part     = observed_part,
    non_censored_part = non_censored_part,
    transition_n_dt   = transition_n_dt
  )
}

################################################################################
### Cause-specific Cox models
################################################################################

# Baseline hazard by stratum and time.
get_basehaz <- function(model) {
  # extract baseline hazard from coxph object
  bh <- as.data.table(basehaz(model, centered = FALSE))
  
  # order on strata and time
  setorder(bh, strata, time)
  
  # return baseline hazard dt
  bh
}

# Cause-specific Cox model for one origin -> destination transition, stratified by
# trt x capped prior hospitalizations. Transitions in dialysis_only are fit on
# trt == 1 only, stratified by hospitalizations alone. coxph warnings are kept
# in $warnings.
fit_cause_specific_cox <- function(dt,
                                   origin,
                                   destination,
                                   cap,
                                   dialysis_only,
                                   weight_col = "sw_IPTW") {
  # prepare data for model fit ------------------------------------------------#
  # episodes starting in `origin`
  origin_dt <- dt[state == origin]
  
  # dialysis-only transitions: keep dialysis-arm episodes (trt constant, dropped from strata)
  is_dialysis_only <- paste(origin, "->", destination) %in% dialysis_only
  
  # only for dialysis-only transitions
  if (is_dialysis_only) {
    # keep dialysis-arm episodes only
    origin_dt <- origin_dt[as.character(trt) == "1"]
  }
  
  # other destinations = censored
  origin_dt[, event := fifelse(!is.na(next_state) &
                                 next_state == destination, 1L, 0L)]
  
  # non-parametric stratum by capped number of previous hospitalization
  origin_dt[, prev_nr_hosp_capped := pmin(prev_nr_hosp, cap)]
  
  # trt is in the strata unless the transition is dialysis-only
  include_trt <- !is_dialysis_only
  
  # one joint stratum column (trt_hosp, or hosp only when trt is excluded)
  if (include_trt) {
    # joint stratum: trt_hosp
    origin_dt[, trt_prev_nr_hosp_capped := paste(trt, prev_nr_hosp_capped, sep = "_")]
  } else {
    # stratum: hosp only
    origin_dt[, trt_prev_nr_hosp_capped := as.character(prev_nr_hosp_capped)]
  }
  
  # fit the Nelson-Aalen estimator --------------------------------------------#
  # build formula for cox fit
  formula <- Surv(time_in_state, event) ~ strata(trt_prev_nr_hosp_capped)
  
  # fit weighted Cox; capture coxph warnings into fit$warnings
  fit_warnings <- character()
  
  # unweighted: simulation uses observed within-stratum rates
  # IPTW is applied only when averaging
  model <- withCallingHandlers(
    coxph(formula, data = origin_dt, ties = "breslow"),
    warning = function(w) {
      # store the warning message
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      # suppress the warning printout
      invokeRestart("muffleWarning")
    }
  )
  
  # extract baseline hazard of cox model, keyed by the joint stratum column ---#
  bh <- get_basehaz(model)
  
  # clean the stratum labels
  bh[, trt_prev_nr_hosp_capped := sub("trt_prev_nr_hosp_capped=", "", strata)]
  
  # baseline hazard per stratum
  bh_by_stratum <- split(bh[, .(time, hazard)], bh$trt_prev_nr_hosp_capped)
  
  # data element named after origin, e.g. "At home" -> "at_home_dt"
  data_name <- paste0(gsub(" ", "_", tolower(origin)), "_dt")
  
  # return model, hazards, flags and episode data -----------------------------#
  c(
    list(
      model = model,
      bh_by_stratum = bh_by_stratum,
      include_trt = include_trt,
      warnings = fit_warnings,
      n_events = sum(origin_dt$event),
      n_zero_event_strata = origin_dt[, sum(event), by = trt_prev_nr_hosp_capped][V1 == 0, .N]
    ),
    setNames(list(origin_dt), data_name)
  )
}

# One row per fitted transition: events, empty strata, and coxph warnings.
summarise_cox_fits <- function(cox_models) {
  # one row per transition
  rbindlist(lapply(names(cox_models), function(nm) {
    # the fit for this transition
    f <- cox_models[[nm]]
    
    # summary of this fit
    data.table(
      transition          = nm,
      include_trt         = f$include_trt,
      n_events            = f$n_events,
      n_zero_event_strata = f$n_zero_event_strata,
      n_warnings          = length(f$warnings),
      warnings            = paste(unique(f$warnings), collapse = " | ")
    )
  }))
}

################################################################################
### Trajectory simulation
################################################################################

# Simulated event time per patient for one cause-specific model, by inverting the
# stratum's cumulative hazard at -log(U) (left-truncated at current_event_time
# for the first event). Inf if never reached or no baseline hazard for the stratum.
simulate_cause_time <- function(fit,
                                trt_chr,
                                hosp_capped,
                                is_first_event = FALSE,
                                current_event_time = NULL) {
  # initialize to simulate event ----------------------------------------------#
  # number of patients in this batch
  n <- length(trt_chr)
  
  # hospitalization count as text
  hosp_chr <- as.character(hosp_capped)
  
  # stratum key per patient
  stratum_key <- if (fit$include_trt) {
    # trt_hosp key
    paste(trt_chr, hosp_chr, sep = "_")
  } else {
    # hosp-only key
    hosp_chr
  }
  
  # one U draw per patient ----------------------------------------------------#
  U <- runif(n)
  
  # simulated times
  result <- numeric(n)
  
  # invert against each stratum's own baseline hazard curve in turn -----------#
  for (s in unique(stratum_key)) {
    # patients in this stratum
    idx <- which(stratum_key == s)
    
    # baseline hazard of this stratum
    bh_s <- fit$bh_by_stratum[[s]]
    
    # no baseline hazard for this stratum: transition never occurs (Inf)
    if (is.null(bh_s) || nrow(bh_s) == 0) {
      # never transitions
      result[idx] <- Inf
      
      # next stratum
      next
    }
    
    # highest cumulative hazard in this stratum
    max_hazard <- bh_s$hazard[length(bh_s$hazard)]
    
    # first event: condition on survival until entry
    if (is_first_event) {
      # hazard already accrued at entry (step function, 0 before first event)
      cet <- current_event_time[idx]
      
      # position of entry time on the hazard curve
      entry_idx <- findInterval(cet, bh_s$time)
      
      # cumulative hazard at entry
      bh_entry <- ifelse(entry_idx == 0, 0, bh_s$hazard[pmax(entry_idx, 1)])
      
      # conditional on surviving to entry: H(t) = -log(U) + H(t_entry)
      threshold <- -log(U[idx]) + bh_entry
    } else {
      # threshold without truncation
      threshold <- -log(U[idx])
    }
    
    # first time H(t) crosses the threshold; Inf if never
    cross_idx <- findInterval(threshold, bh_s$hazard, left.open = TRUE) + 1L
    
    # event time at the threshold, or Inf
    result[idx] <- ifelse(threshold > max_hazard, Inf, bh_s$time[pmin(cross_idx, length(bh_s$hazard))])
  }
  
  # return event times
  result
}

# Simulates all competing transitions out of `origin` for a batch of patients and
# returns each patient's soonest event time and destination.
simulate_next_event <- function(cox_models,
                                origin,
                                trt_chr,
                                hosp_capped,
                                is_first_event = FALSE,
                                current_event_time = NULL) {
  # select the cox models with possible transitions ---------------------------#
  candidates <- names(cox_models)[startsWith(names(cox_models), paste0(origin, " -> "))]
  
  # time to event per candidate transition (one column each) ------------------#
  times_mat <- vapply(candidates, function(nm) {
    # draw times for this candidate transition
    simulate_cause_time(
      fit                = cox_models[[nm]],
      trt_chr            = trt_chr,
      hosp_capped        = hosp_capped,
      is_first_event     = is_first_event,
      current_event_time = current_event_time
    )
  }, numeric(length(trt_chr)))
  
  # restore matrix shape (vapply drops it for a single candidate)
  dim(times_mat) <- c(length(trt_chr), length(candidates))
  
  # per-patient minimum across candidates -------------------------------------#
  row_min <- as.vector(do.call(pmin, asplit(times_mat, 2)))
  
  hit        <- times_mat == row_min
  winner_idx <- max.col(hit * 1, ties.method = "random")
  
  # return the event time and next state for soonest simulated event, per patient
  list(time        = row_min,
       destination = sub(paste0(origin, " -> "), "", candidates[winner_idx]))
}

# Completes each censored patient's trajectory up to horizon_days_2y or death.
# Vectorized: every round advances all active patients by one transition,
# grouped by current state.
simulate_patient_trajectories <- function(id,
                                          start_state,
                                          start_time,
                                          time_in_current_state,
                                          prev_nr_hosp,
                                          trt,
                                          horizon_days_2y,
                                          cox_models,
                                          cap,
                                          max_transitions = 200) {
  # initialize ----------------------------------------------------------------#
  n <- length(id)
  
  # current state per patient
  current_state <- start_state
  
  # current time per patient
  current_time  <- start_time
  
  # time each patient entered the current state
  entry_time    <- start_time - time_in_current_state
  
  # current hospitalization count
  current_hosp  <- prev_nr_hosp
  
  # match trt once, not one data.table per step
  trt_chr <- as.character(trt)
  
  # patients still being simulated
  active <- rep(TRUE, n)
  
  # storage for output rows
  out_rows <- vector("list", max_transitions * 3L)
  
  # number of output blocks
  n_out <- 0L
  
  # loop until every patient has reached the end of follow-up or died ---------#
  for (i in seq_len(max_transitions)) {
    # patients still active
    act_idx <- which(active)
    
    # stop when no one is active
    if (length(act_idx) == 0L)
      break
    
    # left truncation only in the first round
    is_first_event <- (i == 1)
    
    # snapshot: patients that move state this round must not be re-processed
    round_state <- current_state[act_idx]
    
    # process each origin state in turn
    for (st in unique(round_state)) {
      # patients in this state
      grp <- act_idx[round_state == st]
      
      # capped hospitalization count
      hosp_capped_grp <- pmin(current_hosp[grp], cap)
      
      # next transition for all patients in this state; only round 1 is left-truncated
      step <- simulate_next_event(
        cox_models          = cox_models,
        origin              = st,
        trt_chr             = trt_chr[grp],
        hosp_capped         = hosp_capped_grp,
        is_first_event      = is_first_event,
        current_event_time  = if (is_first_event)
          time_in_current_state[grp]
        else
          NULL
      )
      
      # step$time is on the state's own clock (since entry_time)
      new_time <- entry_time[grp] + step$time
      
      # no transition before the horizon: stay in state until the horizon
      finishes <- is.infinite(step$time) |
        (new_time >= horizon_days_2y[grp])
      
      # patients without a further transition before the horizon
      if (any(finishes)) {
        # their indices
        fin_idx <- grp[finishes]
        
        # next output block
        n_out <- n_out + 1L
        
        # final episode until the horizon
        out_rows[[n_out]] <- data.table(
          LOPNR = id[fin_idx],
          tstart = current_time[fin_idx],
          tstop = horizon_days_2y[fin_idx],
          state = st,
          is_first_transition = is_first_event
        )
        
        # no longer active
        active[fin_idx] <- FALSE
      }
      
      # next event before the horizon: record episode, advance (or stop if death)
      if (any(!finishes)) {
        # indices of patients that transition
        cont_idx      <- grp[!finishes]
        
        # time of their transition
        cont_new_time <- new_time[!finishes]
        
        # their destination state
        cont_dest     <- step$destination[!finishes]
        
        # next output block
        n_out <- n_out + 1L
        
        # completed episode
        out_rows[[n_out]] <- data.table(
          LOPNR = id[cont_idx],
          tstart = current_time[cont_idx],
          # round 1: tstop - tstart is the remaining time, not the full episode
          tstop = cont_new_time,
          state = st,
          is_first_transition = is_first_event
        )
        
        # patients that die
        died <- cont_dest == "Death"
        # record death
        if (any(died)) {
          # zero-width Death row (matches missing_states_dt)
          died_idx <- cont_idx[died]
          
          # next output block
          n_out <- n_out + 1L
          
          # zero-width Death row
          out_rows[[n_out]] <- data.table(
            LOPNR = id[died_idx],
            tstart = cont_new_time[died],
            tstop = cont_new_time[died],
            state = "Death",
            is_first_transition = FALSE
          )
          
          # no longer active
          active[died_idx] <- FALSE
        }
        
        # patients still alive
        alive_idx <- cont_idx[!died]
        
        # advance them
        if (length(alive_idx) > 0L) {
          # their transition time
          alive_new_time <- cont_new_time[!died]
          
          # their new state
          alive_dest     <- cont_dest[!died]
          
          # only a Hospitalization destination increments the hospitalization count
          hosp_inc <- alive_dest == "Hospitalization"
          
          # count hospitalizations
          if (any(hosp_inc)) {
            # add one hospitalization
            current_hosp[alive_idx[hosp_inc]] <- current_hosp[alive_idx[hosp_inc]] + 1
          }
          
          # new state
          current_state[alive_idx] <- alive_dest
          
          # new current time
          current_time[alive_idx]  <- alive_new_time
          entry_time[alive_idx]    <- alive_new_time  # every subsequent state is entered fresh
        }
      }
    }
  }
  
  # combine all output blocks
  rbindlist(out_rows[seq_len(n_out)])
}

################################################################################
### Imputation
################################################################################

# Whole-cohort trajectories: observed rows plus one imputation's simulated rows.
# The first simulated row of a censored patient continues the episode that was
# ongoing at censoring (same state), so it is merged into the last observed row.
# days_in_row and in_center are computed AFTER completion, so each state episode
# is counted once (matters for the non-linear PD rule).
build_full_dt_from_simulated <- function(simulated_dt,
                                         observed_part,
                                         non_censored_part,
                                         in_center_params) {
  base_cols <- c("LOPNR", "tstart", "tstop", "state", "trt")
  
  # keep only the raw episode columns; days and in-center are recomputed below
  ncen <- non_censored_part[, ..base_cols][, source := "observed"]
  obs  <- copy(observed_part[, ..base_cols])
  sim  <- simulated_dt[, c(base_cols, "is_first_transition"), with = FALSE]
  
  # first simulated row per patient = remainder of the ongoing episode
  sim_first <- sim[is_first_transition == TRUE]
  sim_rest  <- sim[is_first_transition == FALSE, ..base_cols][, source := "simulated"]
  
  # flag each censored patient's last observed row
  setorder(obs, LOPNR, tstart)
  obs[, is_last := seq_len(.N) == .N, by = LOPNR]
  obs[, source := "observed"]
  
  # attach the first simulated row's info to that patient's rows
  obs[sim_first, on = "LOPNR", `:=`(sim_state = i.state,
                                    sim_tstart = i.tstart,
                                    sim_tstop = i.tstop)]
  
  # sanity checks: same state, and the simulated part starts where the observed part ends
  bad <- obs[is_last & !is.na(sim_tstop) &
               (state != sim_state |
                  abs(tstop - sim_tstart) > 1e-8)]
  if (nrow(bad) > 0)
    warning(
      nrow(bad),
      " censored patients: first simulated row does not continue the last observed row"
    )
  
  # extend the last observed episode to the end of the first simulated row
  obs[is_last &
        !is.na(sim_tstop), `:=`(tstop = sim_tstop, source = "observed+simulated")]
  obs[, c("is_last", "sim_state", "sim_tstart", "sim_tstop") := NULL]
  
  # whole cohort
  full_dt <- rbindlist(list(ncen, obs, sim_rest), use.names = TRUE)
  setorder(full_dt, LOPNR, tstart)
  
  # days and in-center days on the completed episodes
  full_dt[, days_in_row := as.numeric(tstop - tstart)]
  full_dt[, in_center   := in_center_days(state, days_in_row, in_center_params)]
  
  full_dt
}

# One imputation: simulate censored patients, build full trajectories and days per
# patient (incl. in-center).
run_one_imputation <- function(m,
                               last_obs_dt,
                               cox_models,
                               cap_hosp,
                               observed_part,
                               non_censored_part,
                               iptw_dt,
                               in_center_params) {
  # simulate all censored patients' remaining trajectories (vectorized)
  simulated_dt <- simulate_patient_trajectories(
    id                     = last_obs_dt$LOPNR,
    start_state            = last_obs_dt$state,
    start_time             = last_obs_dt$followup_days,
    time_in_current_state  = last_obs_dt$time_in_current_state,
    prev_nr_hosp           = last_obs_dt$prev_nr_hosp,
    trt                    = last_obs_dt$trt,
    horizon_days_2y        = last_obs_dt$horizon_days_2y,
    cox_models             = cox_models,
    cap                    = cap_hosp
  )
  
  # trt is already on last_obs_dt
  simulated_dt <- merge(simulated_dt, last_obs_dt[, .(LOPNR, trt)], by = "LOPNR")
  
  # whole-cohort trajectories
  full_dt_m <- build_full_dt_from_simulated(simulated_dt,
                                            observed_part,
                                            non_censored_part,
                                            in_center_params)
  
  # total days per patient per state, using this replicate's simulated completions
  days_per_patient_imputed_m <- dcast(
    full_dt_m,
    LOPNR + trt ~ state,
    value.var = "days_in_row",
    fun.aggregate = sum,
    fill = 0
  )
  # add in-center days per patient
  days_per_patient_imputed_m <- add_in_center_per_patient(days_per_patient_imputed_m, full_dt_m)
  
  # return days
  list(days_per_patient_imputed = days_per_patient_imputed_m)
}

################################################################################
### Estimation and bootstrap
################################################################################

# Refits IPTW weights and Cox models, runs M imputations and pools them.
fit_and_simulate <- function(baseline,
                             baseline_long,
                             cap_hosp,
                             last_obs_dt,
                             observed_part,
                             non_censored_part,
                             model_PS,
                             id_name,
                             trt_var,
                             horizon,
                             M,
                             transitions,
                             dialysis_only,
                             state_cols,
                             in_center_cols,
                             in_center_params) {
  # local copies (columns are added by reference)
  baseline <- copy(baseline)
  
  # local copy (modified by reference)
  baseline_long <- copy(baseline_long)
  
  # fit cox models
  cox_models <- setNames(
    lapply(transitions, function(state)
      fit_cause_specific_cox(
        dt = baseline_long,
        origin = state[1],
        destination = state[2],
        cap = cap_hosp,
        dialysis_only = dialysis_only
      )),
    sapply(transitions, function(state)
      paste(state[1], "->", state[2]))
  )
  
  # re-derive IPTW weights on this sample (captures weighting uncertainty)
  out_weights <- create_weights(
    data = baseline,
    id_name = id_name,
    trt_var = trt_var,
    model_PS = model_PS,
    w_meth = "IPTW",
    verbose = FALSE
  )
  
  # IPTW weights per patient
  iptw_dt <- data.table(LOPNR = out_weights$data[[id_name]], sw_IPTW = out_weights$data$w)
  
  # drop an old weight column if it already exists
  if ("sw_IPTW" %in% names(baseline_long)) {
    # remove it
    baseline_long[, sw_IPTW := NULL]
  }
  
  # add the new weights
  baseline_long <- merge(baseline_long, iptw_dt, by = "LOPNR")
  
  # Run M imputations and pool
  imputations <- lapply(
    seq_len(M),
    run_one_imputation,
    last_obs_dt       = last_obs_dt,
    cox_models        = cox_models,
    cap_hosp          = cap_hosp,
    observed_part     = observed_part,
    non_censored_part = non_censored_part,
    iptw_dt           = iptw_dt,
    in_center_params  = in_center_params
  )
  
  # pool esitmates ------------------------------------------------------------#
  # pool by averaging the M replicates per patient
  days_per_patient_imputed <- rbindlist(lapply(imputations, `[[`, "days_per_patient_imputed"),
                                        idcol = "m")[, lapply(.SD, mean), by = .(LOPNR, trt), .SDcols = c(state_cols, in_center_cols)]
  
  # raw (unweighted) days per patient per state for the unadjusted rows -------#
  baseline_long[, days_in_row := as.numeric(tstop - tstart)]
  
  # in-center days in each row
  baseline_long[, in_center := in_center_days(state, days_in_row, in_center_params)]
  
  # days per patient per state
  days_per_patient <- dcast(
    baseline_long,
    LOPNR + trt ~ state,
    value.var = "days_in_row",
    fun.aggregate = sum,
    fill = 0
  )
  
  # add in-center days per patient
  days_per_patient <- add_in_center_per_patient(days_per_patient, baseline_long)
  
  # return everything
  list(
    days_per_patient         = days_per_patient,
    days_per_patient_imputed = days_per_patient_imputed,
    cox_models               = cox_models,
    iptw_dt                  = iptw_dt,
    imputations              = imputations
  )
}

# relabel a table with the new IDs
relabel <- function(dt, id_map, id_name) {
  # join the ID map
  out <- merge(
    dt,
    id_map,
    by.x = id_name,
    by.y = "orig_id",
    allow.cartesian = TRUE
  )
  
  # use the new ID
  out[, (id_name) := new_id]
  
  # drop the helper column
  out[, new_id := NULL]
  
  # return the relabeled table
  out
}

# Bootstrap resample of patients with replacement; duplicates get new IDs.
resample_patients <- function(baseline,
                              baseline_long,
                              last_obs_dt,
                              observed_part,
                              non_censored_part,
                              id_name = "LOPNR") {
  # all patient IDs
  orig_ids <- baseline[[id_name]]
  
  # sample IDs with replacement
  sampled_orig_ids <- sample(orig_ids, size = length(orig_ids), replace = TRUE)
  
  # map original to new IDs for duplicates
  id_map <- data.table(orig_id = sampled_orig_ids,
                       new_id = paste0(
                         sampled_orig_ids,
                         "_",
                         ave(sampled_orig_ids, sampled_orig_ids, FUN = seq_along)
                       ))
  
  # return all resampled tables
  list(
    baseline              = relabel(
      dt = copy(baseline),
      id_map = id_map,
      id_name = id_name
    ),
    baseline_long         = relabel(
      dt = copy(baseline_long),
      id_map = id_map,
      id_name = id_name
    ),
    last_obs_dt           = relabel(
      dt = copy(last_obs_dt),
      id_map = id_map,
      id_name = id_name
    ),
    observed_part         = relabel(
      dt = copy(observed_part),
      id_map = id_map,
      id_name = id_name
    ),
    non_censored_part     = relabel(
      dt = copy(non_censored_part),
      id_map = id_map,
      id_name = id_name
    )
  )
}

# One bootstrap replicate; returns NULL if it fails.
run_one_bootstrap <- function(b,
                              baseline,
                              baseline_long,
                              last_obs_dt,
                              observed_part,
                              non_censored_part,
                              cap_hosp,
                              model_PS,
                              id_name,
                              trt_var,
                              horizon,
                              M_bootstrap,
                              transitions,
                              dialysis_only,
                              state_cols,
                              in_center_cols,
                              in_center_params) {
  # resample patients
  resampled <- resample_patients(baseline,
                                 baseline_long,
                                 last_obs_dt,
                                 observed_part,
                                 non_censored_part,
                                 id_name)
  
  # run the analysis, NULL on error
  result <- tryCatch(
    fit_and_simulate(
      baseline              = resampled$baseline,
      baseline_long         = resampled$baseline_long,
      last_obs_dt           = resampled$last_obs_dt,
      observed_part         = resampled$observed_part,
      non_censored_part     = resampled$non_censored_part,
      cap_hosp              = cap_hosp,
      model_PS              = model_PS,
      id_name               = id_name,
      trt_var               = trt_var,
      horizon               = horizon,
      M                     = M_bootstrap,
      transitions           = transitions,
      dialysis_only         = dialysis_only,
      state_cols            = state_cols,
      in_center_cols        = in_center_cols,
      in_center_params      = in_center_params
    ),
    error = function(e) {
      # report the error
      message("Bootstrap replicate failed: ", conditionMessage(e))
      # failed replicate
      NULL
    }
  )
  
  # skip failed replicates
  if (is.null(result))
    return(NULL)
  
  # return estimates and curves
  list(
    final         = compute_all_estimates(
      days_per_patient_imputed = result$days_per_patient_imputed,
      days_per_patient = result$days_per_patient,
      iptw_dt = result$iptw_dt,
      state_cols = c(state_cols, in_center_cols)
    ),
    cox_fit_summary = summarise_cox_fits(result$cox_models)
  )
}

################################################################################
### Summary helpers
################################################################################

# IPTW-weighted mean of state_cols by arm.
weighted_state_means <- function(dt, weight_dt, state_cols) {
  # add weights
  dt <- merge(dt, weight_dt, by = "LOPNR")
  
  # weighted means per arm
  dt[, lapply(.SD, weighted.mean, w = sw_IPTW), by = trt, .SDcols = state_cols]
}

################################################################################
### Multi-state model used to compute time-until-dialysis
################################################################################

# Extract state occupation probabilities
Pprobtrans <- function(bh12, bh13, bh23) {
  # number of time points
  N <- length(bh12)
  out_P11 <- rep(NA_real_, N)
  out_P12 <- rep(NA_real_, N)
  out_P13 <- rep(NA_real_, N)
  P12_1 <- 0
  P13_1 <- 0
  P11_1 <- 1
  # step through time
  for (i in 1:N) {
    P13_2 <- P11_1 * bh13[i] + P12_1 * bh23[i] + P13_1
    P11_2 <- P11_1 * (1 - bh12[i] - bh13[i])
    P12_2 <- P11_1 * bh12[i] + P12_1 * (1 - bh23[i])
    P11_1 <- out_P11[i] <- P11_2
    P12_1 <- out_P12[i] <- P12_2
    P13_1 <- out_P13[i] <- P13_2
  }
  
  # collect results
  df <- data.frame(P11 = out_P11, P12 = out_P12, P13 = out_P13)
  
  # return the probabilities
  return(df)
}

# State probabilities and RMST from an illness-death model (used by other scripts).
state_probabilities <- function(dt, horizon) {
  # transition matrix
  tmat <- mstate::trans.illdeath(c("Decision", "KRT", "Death"))
  
  # set status and time at horizon for KRT and death
  dt[, KRT_event := fifelse(time2event_KRT_inf <= horizon, event_KRT_inf, 0)]
  dt[, KRT_time := fifelse(time2event_KRT_inf <= horizon, time2event_KRT_inf, horizon)]
  dt[, death_event := fifelse(time2event_death_inf <= horizon, event_death_inf, 0)]
  dt[, death_time := fifelse(time2event_death_inf <= horizon,
                             time2event_death_inf,
                             horizon)]
  
  # prepare data for competing risk model
  msdia <- mstate::msprep(
    time = c(NA, "KRT_time", "death_time"),
    status = c(NA, "KRT_event", "death_event"),
    data = dt,
    trans = tmat,
    keep = "sw_IPTW"
  )
  
  # fit competing risk model
  cox_dia <- survival::coxph(
    Surv(Tstart, Tstop, status) ~ strata(trans),
    data = msdia,
    weights = msdia$sw_IPTW,
    method = "breslow"
  )
  
  # Manually calculate state occupation probabilities
  bh <- survival::basehaz(cox_dia, centered = F)
  bh_11 <- bh[bh$strata == "trans=1", ]
  bh_12 <- bh[bh$strata == "trans=2", ]
  bh_13 <- bh[bh$strata == "trans=3", ]
  
  # ensure that baseline hazard is extracted at uniform times
  alltimes <- 0:(horizon + 1)
  
  # fill the NAs of the cumulative hazard with the previous value
  bh_11_allt <- merge(bh_11, data.frame(time = alltimes), all = T) |>
    tidyr::fill(hazard, .direction = "down")
  
  # transition 2 hazard on all times, filled forward
  bh_12_allt <- merge(bh_12, data.frame(time = alltimes), all = T) |>
    tidyr::fill(hazard, .direction = "down")
  
  # transition 3 hazard on all times, filled forward
  bh_13_allt <- merge(bh_13, data.frame(time = alltimes), all = T) |>
    tidyr::fill(hazard, .direction = "down")
  
  # if there is a remaining NA at the iniatial times, that should be a 0
  bh_11_allt$hazard[is.na(bh_11_allt$hazard)] <- 0
  bh_12_allt$hazard[is.na(bh_12_allt$hazard)] <- 0
  bh_13_allt$hazard[is.na(bh_13_allt$hazard)] <- 0
  
  # extract hazard from cumulative hazard
  bh_11_allt$haz <- diff(c(0, bh_11_allt$hazard))
  bh_12_allt$haz <- diff(c(0, bh_12_allt$hazard))
  bh_13_allt$haz <- diff(c(0, bh_13_allt$hazard))
  
  # obtain state probabilities
  state_prob <- Pprobtrans(bh12 = bh_11_allt$haz,
                           bh13 = bh_12_allt$haz,
                           bh23 = bh_13_allt$haz)
  # P11 = event-free; P12 = decision -> KRT; P13 = decision -> death
  
  # Average time spent in the decison state in the next two years
  RMST_11 <- sum(state_prob$P11 * diff(c(alltimes, horizon))) / 30.44
  
  # Average time spent in dialysis in the next two years
  RMST_12 <- sum(state_prob$P12 * diff(c(alltimes, horizon))) / 30.44
  
  # Average time spent alive in the next two years
  RMST_13 <- sum((1 - state_prob$P13) * diff(c(alltimes, horizon))) / 30.44
  
  # return times, probabilities and RMST
  return(list(
    times = alltimes,
    state_prob = state_prob,
    RMST = list(
      RMST_11 = RMST_11,
      RMST_12 = RMST_12,
      RMST_13 = RMST_13
    )
  ))
}