# -----------------------------
# Core function to compute a single outcome for a single window
# -----------------------------
calc_event <- function(cohort_dt,
                       end_follow_up = NULL,
                       censor_at_death = FALSE, # Censor at death or not
                       outcome_dataset = NULL,
                       id_name = NULL,     # Column containing IDs (e.g., LOPNR) 
                       date_col = NULL,    # Column containing the event date (e.g., DODSDAT or INDATUM)
                       code_col = NULL,    # Column containing ICD or event code (e.g., ULORSAK, hdia)
                       codes = NULL,       # Regex pattern to filter events (e.g., "^I21|^I22|^I23")
                       window = NULL,      # Follow-up window in years
                       extra_cols = NULL,  # Extra columns
                       suffix = "") {      # Suffix for output column names (e.g., "_death_5y")
  
  # Make a copy of the cohort to avoid modifying the original data
  dt <- copy(cohort_dt)
  
  # rename ID name temporarily
  setnames(dt, id_name, "ID")
  
  # If an external outcome dataset is supplied (e.g., CV death, MI, Stroke)
  if (!is.null(outcome_dataset)) {
    # Filter rows by ICD/event codes if a code column and regex are supplied
    if (!is.null(code_col) & !is.null(codes)) {
      outcome_dataset <- outcome_dataset[grepl(codes, get(code_col))]
    }
    
    # Keep ID + event date + any extra columns requested
    keep_cols <- c(id_name, date_col, extra_cols)
    outcome_dataset <- outcome_dataset[, ..keep_cols]
    setnames(outcome_dataset, c(id_name, date_col), c("ID", "event_date"))
    
    # Join visit_date from dt (bringing only visit_date)
    outcome_dataset <- dt[, ][outcome_dataset,
                              on = "ID",
                              nomatch = 0]
    
    # Filter events after visit_date and get the earliest event per ID
    setorder(outcome_dataset, ID, event_date)
    earliest_events <- outcome_dataset[event_date > visit_date,
                                       .SD[1, c("event_date", extra_cols), with = FALSE],
                                       by = "ID"]

    # Merge earliest_events back into dt
    dt <- earliest_events[dt, on = "ID"]
  } else {
    # If no external dataset, assume the event date is in the cohort itself (e.g., DODSDAT for all-cause death)
    dt[, event_date := get(date_col)]
  }
  
  # -----------------------------
  # Compute censoring date for the follow-up window
  # -----------------------------
  # Take the minimum of:
  # 1. End of follow-up (end_follow_up)
  # 2. Visit date + window, taking into account leap year dates
  # 3. If censor_at_death is TRUE, censor at date of death (DODSDAT)
  # This ensures that follow-up does not exceed the specified window or study end
  if (censor_at_death) {
    dt[, censor := pmin(DODSDAT,
                        end_follow_up,
                        lubridate::add_with_rollback(visit_date, lubridate::years(window)),
                        na.rm = TRUE)]
  } else {
    # For mortality outcome: do NOT censor at DODSDAT
    dt[, censor := pmin(end_follow_up,
                        lubridate::add_with_rollback(visit_date, lubridate::years(window)),
                        na.rm = TRUE)]
  }
  
  # -----------------------------
  # Define the event indicator
  # -----------------------------
  # Event = 1 if event_date exists, is after visit_date, and occurs before or on the censor date
  # Event = 0 otherwise
  dt[, event := fifelse(!is.na(event_date) &
                          event_date > visit_date &
                          event_date <= censor,
                        1L,
                        0L)]
  
  # -----------------------------
  # Keep only the first event per patient
  # -----------------------------
  # TODO: check if needed
  # If the event didn't occur inside this window, the extra info doesn't apply
  if (!is.null(extra_cols)) {
    dt[event == 0L, (extra_cols) := NA]
  }
  
  # Order by ID, descending event (so 1 comes before 0), and ascending event_date
  setorder(dt, ID, -event, event_date)
  
  # Keep the first row per patient
  dt <- dt[, .SD[1], by = ID]
  
  # -----------------------------
  # Compute event date and time to event
  # -----------------------------
  # If the event occurred, event_dt = event_date; otherwise, event_dt = censor
  dt[, event_dt := fifelse(event == 1L, event_date, censor)]
  
  # Time to event in days, automatically taking into account leap years
  dt[, time2event := as.numeric(event_dt - visit_date)]
  
  # -----------------------------
  # Rename columns to include outcome and window suffix
  # -----------------------------
  keep <- c("ID", "event", "event_dt", "time2event", extra_cols)
  dt <- dt[, ..keep]
  setnames(dt, "ID", id_name) # revert naming
  
  var <- c("event", "event_dt", "time2event", extra_cols)
  setnames(dt, old = var, new = paste0(var, suffix)) # add suffix for clarity
  
  return(dt)
}

# -----------------------------
# Wrapper to compute multiple outcomes for multiple windows
# -----------------------------
add_multiple_outcomes <- function(cohort_dt,
                                  windows,
                                  id_name, 
                                  end_follow_up,
                                  outcomes_list) {
  # Copy the cohort to avoid modifying the original
  dt <- copy(cohort_dt)
  
  # Loop through each follow-up window
  for (win_name in names(windows)) {
    window <- windows[[win_name]]
    
    # Loop through each outcome (death, cvdeath, MI, stroke, etc.)
    for (outcome_name in names(outcomes_list)) {
      outcome <- outcomes_list[[outcome_name]]
      
      # Compute the outcome for this window using the core function
      res <- calc_event(
        cohort_dt = dt,
        end_follow_up = end_follow_up,
        censor_at_death = ifelse(outcome_name == "death", FALSE, TRUE), # only censor all-cause death if it is not the outcome
        outcome_dataset = outcome$dataset, # dataset containing the events
        id_name = id_name,                 # column with the IDs
        date_col = outcome$date_col,       # column with the event date
        code_col = outcome$code_col,       # column with ICD/event code (optional)
        codes = outcome$codes,             # regex pattern to filter events (optional)
        window = window,
        suffix = paste0("_", outcome_name, "_", win_name), # e.g., "_death_5y"
        extra_cols = outcome$extra_cols
      )
      
      # Merge the computed outcome back into the main cohort table
      dt <- merge(dt, res, by = id_name, all.x = TRUE)
    }
  }
  
  return(dt[])
}
