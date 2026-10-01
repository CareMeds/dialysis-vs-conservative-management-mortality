# Define the stand_diff function
# Returns the absolute standardized difference (unsigned), consistent with
# tableone's own SMD convention (|mean diff| / pooled SD for continuous
# variables, a non-negative Mahalanobis-type distance for categorical ones).
# Without abs(), this would return a signed value whenever the control
# group's proportion exceeds the treatment group's for a given category
# level - inconsistent with every other SMD in the table.
stand_diff <- function(pT, pC) {
  d <- abs(pT - pC) / (sqrt((pT * (1 - pT) + pC * (1 - pC)) / 2))
  return(d)
}

# Helper: Extract % value inside parentheses from strings like "35 (12)"
extract_proportion <- function(entry) {
  match <- regmatches(entry, regexpr("\\(([^)]+)\\)", entry))
  if (length(match) > 0) {
    return(as.numeric(sub("\\(([^)]+)\\)", "\\1", match)))
  } else {
    return(NA)
  }
}

# Function to compute SMD for each category of a categorical variable
create_baseline_table <- function(data,
                                  id_name,
                                  weights = NULL,
                                  vars,
                                  categoricalVars,
                                  continuousVars,
                                  IQRVars = NULL,
                                  treatmentColumn = NULL,
                                  treatmentLabel = NULL,
                                  controlLabel,
                                  treatmentValue = NULL,
                                  controlValue = NULL,
                                  tableCaption,
                                  tableRowLabels = NA) {
  # treatmentValue/controlValue are the ACTUAL values found in treatmentColumn
  # (e.g. 0/1, or "HD"/"PD"). They default to "1"/"0" to preserve behavior
  # for existing calls that stratify on 0/1-coded columns (S, trt_var, etc.),
  # where treatmentLabel/controlLabel are only cosmetic renames applied after
  # the fact. For columns already coded with the display values themselves
  # (like dialysis_type == "HD"/"PD"), pass treatmentValue/controlValue
  # explicitly.
  if (is.null(treatmentValue))
    treatmentValue <- "1"
  if (is.null(controlValue))
    controlValue <- "0"
  # extract IDs
  data <- copy(data)
  data[, ID := get(id_name)]
  
  #-----------------------------
  # Create overall baseline table (weighted or unweighted)
  #-----------------------------
  if (is.null(weights)) {
    table_overall <- tableone::CreateTableOne(
      vars = vars,
      data = data,
      factorVars = categoricalVars,
      includeNA = TRUE
    )
  } else {
    weighted_baseline <- survey::svydesign(ids = ~ ID,
                                           weights = ~ weights,
                                           data = data)
    table_overall <- tableone::svyCreateTableOne(
      vars = vars,
      data = weighted_baseline,
      factorVars = categoricalVars,
      includeNA = TRUE
    )
  }
  
  # Correctly specify decimal arguments (catDigits etc.)
  table_overall <- print(
    table_overall,
    printToggle = FALSE,
    nonnormal = IQRVars,
    noSpaces = TRUE,
    catDigits = 1,
    contDigits = 1,
    pDigits = 3,
    format = "fp"
  )
  
  # Only apply to rows belonging to categorical variables
  cont_rows <- grepl(paste0("^(", paste(continuousVars, collapse = "|"), ")"), trimws(rownames(table_overall)))
  rn <- rownames(table_overall)
  table_overall[!cont_rows, ] <- apply(table_overall[!cont_rows, , drop = FALSE], 2, fix_counts)
  rownames(table_overall) <- rn
  
  # Optional row labels
  if (!is.null(tableRowLabels) && length(tableRowLabels) > 1) {
    row.names(table_overall) <- tableRowLabels
  }
  colnames(table_overall) <- c("Overall")
  
  #-----------------------------
  # Stratified baseline table by treatment group (if provided)
  #-----------------------------
  if (!is.null(treatmentColumn)) {
    if (is.null(weights)) {
      table_stratified <- tableone::CreateTableOne(
        vars = vars,
        data = data,
        factorVars = categoricalVars,
        strata = treatmentColumn
      )
    } else {
      weighted_baseline <- survey::svydesign(ids = ~ ID,
                                             weights = ~ weights,
                                             data = data)
      table_stratified <- tableone::svyCreateTableOne(
        vars = vars,
        data = weighted_baseline,
        factorVars = categoricalVars,
        strata = treatmentColumn
      )
    }
    
    # Correctly specify decimal arguments (catDigits etc.)
    table_stratified <- print(
      table_stratified,
      printToggle = FALSE,
      nonnormal = IQRVars,
      noSpaces = TRUE,
      catDigits = 1,
      contDigits = 1,
      pDigits = 3,
      smd = TRUE
    )
    
    # Only apply to rows belonging to categorical variables
    cont_rows <- grepl(paste0("^(", paste(continuousVars, collapse = "|"), ")"), trimws(rownames(table_stratified)))
    rn <- rownames(table_stratified)
    table_stratified_rounded <- table_stratified
    table_stratified_rounded[!cont_rows, ] <- apply(table_stratified[!cont_rows, , drop = FALSE], 2, fix_counts)
    rownames(table_stratified_rounded) <- rn
    
    # Keep treatment, control, and SMD columns
    missing_vals <- setdiff(c(treatmentValue, controlValue),
                            colnames(table_stratified_rounded))
    if (length(missing_vals) > 0) {
      stop(
        "create_baseline_table: could not find column(s) ",
        paste(sprintf('"%s"', missing_vals), collapse = ", "),
        " in the stratified table produced from treatmentColumn = \"",
        treatmentColumn,
        "\".\n",
        "Available columns are: ",
        paste(sprintf(
          '"%s"', colnames(table_stratified_rounded)
        ), collapse = ", "),
        ".\n",
        "Pass treatmentValue/controlValue matching the ACTUAL values in `",
        treatmentColumn,
        "` ",
        "(treatmentLabel/controlLabel are only used to relabel the output columns)."
      )
    }
    table_stratified <- as.matrix(cbind(table_stratified_rounded[, c(treatmentValue, controlValue)], table_stratified[, "SMD"]))
    
    # Optional row labels
    if (!is.null(tableRowLabels) && length(tableRowLabels) > 1) {
      row.names(table_stratified) <- tableRowLabels
    }
    
    colnames(table_stratified) <- c(treatmentLabel, controlLabel, "SMD")
    
    #-----------------------------
    # Handle SMD values
    #-----------------------------
    # 1. Extract numeric SMDs from the table
    SMDs <- table_stratified[which(table_stratified[, "SMD"] != ""), "SMD"]
    smd_table <- suppressWarnings(as.numeric(SMDs))
    smd_table[which(SMDs == "<0.001")] <- 0
    names(smd_table) <- names(SMDs)
    
    # 2. Compute categorical SMDs manually using extracted proportions
    smd_values <- sapply(1:nrow(table_stratified), function(i) {
      treatment_entry <- table_stratified[i, treatmentLabel]
      control_entry <- table_stratified[i, controlLabel]
      
      pT <- extract_proportion(treatment_entry) / 100
      pC <- extract_proportion(control_entry) / 100
      
      if (!is.na(pT) && !is.na(pC)) {
        return(stand_diff(pT, pC))  # user-defined function
      } else {
        return(NA)
      }
    })
    smd_values <- round(smd_values, 3)
    
    # 3. Merge manual SMDs and existing ones
    existing_smd <- suppressWarnings(as.numeric(table_stratified[, "SMD"]))
    new_smd <- ifelse(is.na(existing_smd), smd_values, existing_smd)
    table_stratified[, "SMD"] <- fmt(new_smd, 3)
    
    # replace NA SMD by empty string
    table_stratified[which(table_stratified == "NA")] <- ""
    
    #-----------------------------
    # Combine overall + stratified tables
    #-----------------------------
    table_1 <- cbind(table_overall, table_stratified)
    colnames(table_1) <- c("Overall", treatmentLabel, controlLabel, "SMD")
    
  } else {
    # If no treatment column, only overall table
    table_1 <- table_overall
    colnames(table_1) <- c("Overall")
    smd_table <- NA
  }
  
  #-----------------------------
  # Pretty print using knitr::kable
  #-----------------------------
  formatted_table <- knitr::kable(table_1,
                                  align = "c",
                                  caption = as.character(tableCaption))
  
  #-----------------------------
  # Return both formatted and raw data
  #-----------------------------
  return(list(
    formatted_table = formatted_table,
    raw_table = table_1,
    smd_table = smd_table
  ))
}

# Insert a labeled header row before each named section of a baseline table
# (blank data cells, with the section name as the row's label), instead of
# manually rbind()-ing hardcoded row-index slices with plain blank
# separator rows.
#
# `section_sizes` is a named integer vector where each name is the section
# label shown in the output (e.g. "Demographics") and each value is the
# number of rows that section occupies in `tbl`. The function verifies
# nrow(tbl) == sum(section_sizes) and fails loudly if they don't match,
# instead of silently inserting rows in the wrong place - which is exactly
# what happens if a variable is added/removed/reordered and the manual
# index ranges (e.g. 1:27, 28:53, ...) aren't updated too.
insert_section_breaks <- function(tbl, section_sizes) {
  stopifnot(
    "insert_section_breaks: nrow(tbl) does not match sum(section_sizes) - a variable was likely added, removed or reordered without updating section_sizes" =
      nrow(tbl) == sum(section_sizes)
  )
  n_sections <- length(section_sizes)
  section_names <- names(section_sizes)
  if (is.null(section_names) || any(section_names == "")) {
    section_names <- paste0("Section ", seq_len(n_sections))
  }
  ends <- cumsum(section_sizes)
  starts <- c(1, head(ends, -1) + 1)
  
  pieces <- vector("list", n_sections)
  for (i in seq_len(n_sections)) {
    pieces[[i]] <- tbl[starts[i]:ends[i], , drop = FALSE]
  }
  
  # a one-row, all-blank slice with the section name as its row label - this
  # is what makes the section name visible in the exported table (rownames
  # become the leftmost label column via write.xlsx(..., rowNames = TRUE)),
  # while still visually separating sections the same way a blank row did
  header_row <- function(name) {
    matrix(rep("", ncol(tbl)),
           nrow = 1,
           dimnames = list(name, colnames(tbl)))
  }
  
  out <- rbind(header_row(section_names[1]), pieces[[1]])
  if (n_sections > 1) {
    for (i in 2:n_sections) {
      out <- rbind(out, header_row(section_names[i]), pieces[[i]])
    }
  }
  out
}

risk_model_table <- function(model_cox,
                             predictor_labels,
                             horizon,
                             digits = 3) {
  # ── Validate labels match model terms ───────────────────────────────────────
  model_terms <- broom::tidy(model_cox) |> dplyr::pull(term)
  
  # ── Coefficient and HR table ────────────────────────────────────────────────
  predictor_rows <- dplyr::left_join(
    broom::tidy(model_cox, exponentiate = FALSE, conf.int = TRUE),
    broom::tidy(model_cox, exponentiate = TRUE, conf.int = TRUE),
    by     = "term",
    suffix = c("_log", "_hr")
  ) |>
    dplyr::mutate(
      Predictor = predictor_labels,
      coef_CI   = fmt_ci(estimate_log, conf.low_log, conf.high_log, digits = digits),
      HR_CI     = fmt_ci(estimate_hr, conf.low_hr, conf.high_hr, digits = digits),
      Wald      = fmt(statistic_log^2)  # z² = Wald chi-square (1 df)
    ) |>
    dplyr::select(Predictor, coef_CI, HR_CI, Wald)
  
  # ── Baseline hazard at horizon ──────────────────────────────────────────────
  bh <- suppressWarnings(survival::basehaz(model_cox))
  h0 <- bh$hazard[bh$time == horizon]
  
  baseline_row <- data.frame(
    Predictor = paste0("Baseline hazard at ", horizon, " years"),
    coef_CI   = sprintf("%.*f", digits, h0),
    HR_CI     = "",
    Wald      = ""
  )
  
  risk_model_table <- rbind(baseline_row, predictor_rows)
  
  return(
    list(
      h0 = h0,
      coef = coefficients(model_cox),
      centers = model_cox$means,
      risk_model_table = risk_model_table
    )
  )
}

# Build a table for results
build_results_table <- function(data,
                                out_est,
                                label,
                                outcome_var,
                                trt_var,
                                control_label,
                                treatment_label,
                                w_meth,
                                unit) {
  results_df <- data.frame(Control = label, Treatment = "")
  rownames(results_df) <- "Outcome"
  colnames(results_df) <- c(control_label, treatment_label)
  
  # Header row for the weighting method
  results_df[ifelse(w_meth == "unweighted",
                    "Unweighted",
                    paste("Weighting", w_meth)), ] <- rep("", 2)
  
  # Sample size
  results_df["Sample size", ] <- c(sum(data[[trt_var]] == 0), sum(data[[trt_var]] == 1))
  
  # Number of events
  results_df["Number of events", ] <- c(sum(data[[outcome_var]] == 1 &
                                              data[[trt_var]] == 0), sum(data[[outcome_var]] == 1 &
                                                                           data[[trt_var]] == 1))
  
  # Absolute risks
  results_df[paste("Risk, % (95% CI)", w_meth), ] <- c(
    fmt_ci(out_est$R0 * 100, out_est$R0_lower * 100, out_est$R0_upper * 100),
    fmt_ci(out_est$R1 * 100, out_est$R1_lower * 100, out_est$R1_upper * 100)
  )
  
  # Risk difference
  results_df[paste("Risk difference, % (95% CI)", w_meth), ] <- c("Reference",
                                                                  fmt_ci(out_est$RD * 100, out_est$RD_lower * 100, out_est$RD_upper * 100))
  
  # Risk ratio
  results_df[paste("Risk ratio (95% CI)", w_meth), ] <- c("Reference",
                                                          fmt_ci(out_est$RR, out_est$RR_lower, out_est$RR_upper, 2))
  
  # RMST
  results_df[paste0("RMST, ", unit, " (95% CI) ", w_meth), ] <- c(
    fmt_ci(out_est$RMST0, out_est$RMST0_lower, out_est$RMST0_upper),
    fmt_ci(out_est$RMST1, out_est$RMST1_lower, out_est$RMST1_upper)
  )
  
  # RMST difference
  results_df[paste0("\u0394RMST, ", unit, " (95% CI) ", w_meth), ] <- c("Reference",
                                                                        fmt_ci(out_est$dRMST, out_est$dRMST_lower, out_est$dRMST_upper))
  
  # Hazard ratio
  results_df[paste("HR (95% CI)", w_meth), ] <- c("Reference",
                                                  fmt_ci(out_est$HR, out_est$HR_lower, out_est$HR_upper, 2))
  
  return(results_df)
}

# Computes all four adjustment-level estimates (unadjusted / confounding
# only / censoring only / both) as one long table: trt, state,
# adjustment, value. Used for the real, non-bootstrapped point estimate and by the
# bootstrap loop (once per replicate, so each adjustment level gets its own CI).
compute_all_estimates <- function(days_per_patient_imputed,
                                  days_per_patient,
                                  iptw_dt,
                                  state_cols) {
  ests <- list(
    "Unadjusted"                             = days_per_patient[, lapply(.SD, mean), by = trt, .SDcols = state_cols],
    "Adjusted for confounding only"          = weighted_state_means(days_per_patient, iptw_dt, state_cols),
    "Adjusted for censoring only"            = days_per_patient_imputed[, lapply(.SD, mean), by = trt, .SDcols = state_cols],
    "Adjusted for confounding and censoring" = weighted_state_means(days_per_patient_imputed, iptw_dt, state_cols)
  )
  long_dt <- rbindlist(lapply(names(ests), function(adj) {
    long <- melt(
      ests[[adj]],
      id.vars = "trt",
      variable.name = "state",
      value.name = "value"
    )
    long[, adjustment := adj]
    long
  }))
  long_dt[, `:=`(trt = as.character(trt), state = as.character(state))]
  
  # Dialysis - Conservative management difference, added as trt "diff".
  # Computed here so every caller (the real analysis and each bootstrap
  # replicate) gets the difference taken within the same replicate
  diff_dt <- dcast(long_dt, state + adjustment ~ trt, value.var = "value")
  diff_dt <- diff_dt[, .(trt = "diff", state, adjustment, value = `1` - `0`)]
  long_dt <- rbind(long_dt, diff_dt, use.names = TRUE)
  
  # Time after the start of dialysis, averaged over ALL patients who chose dialysis
  # (patients who never start contribute 0, so nothing is conditioned on a future event).
  # Extra "trt" group "after", same measures as the arms. Home after the start = days alive
  # after the start - in-center days after the start (includes days at home on HD/PD)
  after_measures_dt <- long_dt[trt == "1", {
    v <- function(s)
      value[state == s]
    .(
      trt   = "after",
      state = c("Time at home", "In-center", "Hospitalization"),
      value = c(
        v("Days after dialysis") - v("In-center after dialysis"),
        v("In-center after dialysis"),
        v("Hospital after dialysis")
      )
    )
  }, by = adjustment]
  
  # Time BEFORE the start of dialysis, same measures, same population (all patients who
  # chose dialysis): total (dialysis arm) - after. Taken within each replicate / adjustment,
  # so before + after add up exactly to the total of "Choose dialysis". Extra "trt" group "before"
  total_dt <- long_dt[trt == "1" &
                        state %in% c("Time at home", "In-center", "Hospitalization")]
  after_dt2 <- after_measures_dt[, .(adjustment, state, after = value)]
  before_measures_dt <- merge(total_dt, after_dt2, by = c("adjustment", "state"))
  before_measures_dt <- before_measures_dt[, .(trt = "before", state, adjustment, value = value - after)]
  rbind(long_dt, after_measures_dt, before_measures_dt, use.names = TRUE)
}

# Table of mean months over the horizon: one block per adjustment level, columns Conservative
# management / Dialysis / Dialysis vs conservative management, rows At home and In-center,
# each followed by its split (without dialysis / on dialysis, or hospitalized / on dialysis).
# The split rows are only filled for the dialysis column: only ~3% of the conservative arm ever
# starts dialysis, so a split for that arm (and so for the difference) would rest on very few
# events (the PD transitions on 1-3 events). Say so in the table footnote.
# estimate_ci_months_dt: long (trt, state, adjustment, estimate, lower, upper) in MONTHS.
# n_dt: data.table(trt, N) with the number of patients per arm (shown in the column headers).
build_home_time_table_dt <- function(estimate_ci_months_dt,
                                     n_dt,
                                     horizon_months = 24,
                                     digits = 1) {
  dt <- copy(estimate_ci_months_dt)
  dt[, formatted := fmt_ci(estimate, lower, upper, digits = digits)]
  
  # row layout: label, state, and whether the row is only filled for the dialysis column
  rows <- data.table(
    label = c(
      "At home",
      "   Without dialysis",
      "   On dialysis",
      "In-center",
      "   Hospitalized",
      "   On dialysis"
    ),
    state = c(
      "Time at home",
      "Home without dialysis",
      "Home on dialysis",
      "In-center",
      "In-center hospitalized",
      "In-center dialysis"
    ),
    dialysis_only = c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
  )
  
  # column headers with the group sizes
  n_fmt <- function(code)
    formatC(n_dt[trt == code, N], format = "d", big.mark = ",")
  col_names <- c(
    first = paste0("Mean months over ", horizon_months, " months"),
    `0`   = paste0("Conservative management (N=", n_fmt("0"), ")"),
    `1`   = paste0("Dialysis (N=", n_fmt("1"), ")"),
    diff  = "Dialysis vs conservative management"
  )
  
  adjustments <- c(
    "Unadjusted",
    "Adjusted for confounding only",
    "Adjusted for censoring only",
    "Adjusted for confounding and censoring"
  )
  
  rbindlist(lapply(adjustments, function(adj) {
    body <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
      cells <- vapply(c("0", "1", "diff"), function(g) {
        if (rows$dialysis_only[i] && g != "1")
          return(NA_character_)
        out <- dt[adjustment == adj &
                    trt == g & state == rows$state[i], formatted]
        if (length(out) == 0)
          NA_character_
        else
          out
      }, character(1))
      data.table(
        first = rows$label[i],
        `0` = cells[["0"]],
        `1` = cells[["1"]],
        diff = cells[["diff"]]
      )
    }))
    # adjustment level as a header row above its block
    rbind(
      data.table(
        first = adj,
        `0` = NA_character_,
        `1` = NA_character_,
        diff = NA_character_
      ),
      body
    )
  }))[, setNames(.SD, col_names[names(.SD)])]
}