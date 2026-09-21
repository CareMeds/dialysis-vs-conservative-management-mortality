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
  if (is.null(treatmentValue)) treatmentValue <- "1"
  if (is.null(controlValue)) controlValue <- "0"
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
  cont_rows <- grepl(
    paste0("^(", paste(continuousVars, collapse = "|"), ")"),
    trimws(rownames(table_overall))
  )
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
    cont_rows <- grepl(
      paste0("^(", paste(continuousVars, collapse = "|"), ")"),
      trimws(rownames(table_stratified))
    )
    rn <- rownames(table_stratified)
    table_stratified_rounded <- table_stratified
    table_stratified_rounded[!cont_rows, ] <- apply(table_stratified[!cont_rows, , drop = FALSE], 2, fix_counts)
    rownames(table_stratified_rounded) <- rn
    
    # Keep treatment, control, and SMD columns
    missing_vals <- setdiff(c(treatmentValue, controlValue), colnames(table_stratified_rounded))
    if (length(missing_vals) > 0) {
      stop(
        "create_baseline_table: could not find column(s) ",
        paste(sprintf('"%s"', missing_vals), collapse = ", "),
        " in the stratified table produced from treatmentColumn = \"", treatmentColumn, "\".\n",
        "Available columns are: ", paste(sprintf('"%s"', colnames(table_stratified_rounded)), collapse = ", "), ".\n",
        "Pass treatmentValue/controlValue matching the ACTUAL values in `", treatmentColumn, "` ",
        "(treatmentLabel/controlLabel are only used to relabel the output columns)."
      )
    }
    table_stratified <- as.matrix(cbind(table_stratified_rounded[, c(treatmentValue, controlValue)],
                                        table_stratified[, "SMD"]))
    
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
    matrix(
      rep("", ncol(tbl)),
      nrow = 1,
      dimnames = list(name, colnames(tbl))
    )
  }
  
  out <- rbind(header_row(section_names[1]), pieces[[1]])
  if (n_sections > 1) {
    for (i in 2:n_sections) {
      out <- rbind(out, header_row(section_names[i]), pieces[[i]])
    }
  }
  out
}

create_table_with_ci <- function(data_absolute_risks,
                                 data_column_headers = c("Treatment (%)",
                                                         "Control (%)",
                                                         "Risk difference (%)",
                                                         "Risk ratio"),
                                 data_decimals = 2,
                                 row_labels_header = "Time (months)",
                                 row_labels = "row_labels",
                                 row_labels_decimals = 0,
                                 table_caption = "Absolute risks of death",
                                 .extension_mean = "_estimate",
                                 .extension_low_CI = "_conf.low",
                                 .extension_high_CI = "_conf.high") {
  # Get column names without mean/CI extensions
  mean_columns <- grep(.extension_mean, colnames(data_absolute_risks), value = TRUE)
  column_names <- gsub(.extension_mean, "", mean_columns)
  
  # For each variable, format mean and confidence intervals
  table_df <- data_absolute_risks
  for (column in column_names) {
    table_df$mean <- table_df[[paste0(column, .extension_mean)]]
    table_df$low_CI <- table_df[[paste0(column, .extension_low_CI)]]
    table_df$high_CI <- table_df[[paste0(column, .extension_high_CI)]]
    
    table_df[, formatted := fmt_ci(mean, low_CI, high_CI, digits = data_decimals)]
    table_df[[column]] <- table_df$formatted
  }
  table_df <- table_df |> select(column_names)
  
  # Add labels
  colnames(table_df) <- data_column_headers
  table_df[[row_labels_header]] <- round(data_absolute_risks[[row_labels]], digits = row_labels_decimals)
  table_df <- table_df[, c(ncol(table_df), 1:(ncol(table_df) - 1))]
  
  # Create table
  formatted_table <- knitr::kable(table_df,
                                  align = "c",
                                  caption = as.character(table_caption))
  
  return(list(raw_table = table_df, formatted_table = formatted_table))
}

# Weighted Variance Function
# Formula: sum(w * (x - mu)^2) / (sum(w) - sum(w^2)/sum(w))
# Note: If all w=1, this reduces to sum(x-mu)^2 / (N-1), which is standard var()
calc_wvar <- function(x, w, wm) {
  sum(w * (x - wm)^2) / (sum(w) - sum(w^2)/sum(w))
}

# Helper: Weighted Covariance Matrix & Means
get_wstats <- function(v, w, levs) {
  if(length(v) == 0) return(NULL)
  
  # Dummy Matrix (rows=obs, cols=levels)
  mat <- t(sapply(v, function(x) as.numeric(levs == x)))
  if(nrow(mat) == 1) mat <- t(mat)
  colnames(mat) <- levs
  
  # Remove last column (k-1 degrees of freedom)
  if (ncol(mat) > 1) mat <- mat[, -ncol(mat), drop = FALSE]
  
  # Weighted Means (Proportions)
  w_props <- colSums(mat * w) / sum(w)
  
  # Center matrix for covariance calc
  mat_centered <- sweep(mat, 2, w_props, "-")
  
  # Weighted Covariance: (X' W X) / (sum(w) - correction)
  mat_weighted <- mat_centered * sqrt(w) 
  CovMat <- crossprod(mat_weighted)
  denom <- sum(w) - (sum(w^2) / sum(w))
  
  return(list(p = w_props, cov = CovMat / denom))
}

calculate_smd <- function(vec1, vec2, w1 = NULL, w2 = NULL) {
  
  # --- 1. Robust Input Handling ---
  # Convert inputs to simple vectors to handle data.table columns or lists
  vec1 <- unlist(as.vector(vec1))
  vec2 <- unlist(as.vector(vec2))
  
  # --- 2. Handle Weights ---
  # If weights are missing (NULL), assign 1 to everyone (Unweighted mode)
  if (is.null(w1)) w1 <- rep(1, length(vec1))
  else w1 <- unlist(as.vector(w1))
  
  if (is.null(w2)) w2 <- rep(1, length(vec2))
  else w2 <- unlist(as.vector(w2))
  
  # --- 3. Clean NAs (Synchronized) ---
  # Remove observations where either Data OR Weight is NA
  valid1 <- !is.na(vec1) & !is.na(w1)
  vec1 <- vec1[valid1]; w1 <- w1[valid1]
  
  valid2 <- !is.na(vec2) & !is.na(w2)
  vec2 <- vec2[valid2]; w2 <- w2[valid2]
  
  # Stop if empty
  if (length(vec1) < 2 || length(vec2) < 2) return(NA)
  
  # --- 4. Logic Switch: Numeric vs Categorical ---
  is_numeric <- is.numeric(vec1) && is.numeric(vec2)
  
  if (is_numeric) {
    # === A. CONTINUOUS VARIABLES ===
    
    # Weighted Mean
    wm1 <- sum(vec1 * w1) / sum(w1)
    wm2 <- sum(vec2 * w2) / sum(w2)
    
    # Weighted Variance
    wv1 <- calc_wvar(vec1, w1, wm1)
    wv2 <- calc_wvar(vec2, w2, wm2)
    
    # Pooled SD
    pooled_sd <- sqrt((wv1 + wv2) / 2)
    
    # SMD
    smd <- abs(wm1 - wm2) / pooled_sd
    
  } else {
    # === B. CATEGORICAL VARIABLES (Mahalanobis) ===
    
    v1_char <- as.character(vec1)
    v2_char <- as.character(vec2)
    all_levs <- sort(unique(c(v1_char, v2_char)))
    
    res1 <- get_wstats(v1_char, w1, all_levs)
    res2 <- get_wstats(v2_char, w2, all_levs)
    
    if (is.null(res1) || is.null(res2)) return(NA)
    
    # Pooled Covariance
    S_pooled <- (res1$cov + res2$cov) / 2
    
    # Mahalanobis Distance
    diff_p <- res1$p - res2$p
    
    # Calculate D^2
    dist_sq <- tryCatch({
      t(diff_p) %*% solve(S_pooled) %*% diff_p
    }, error = function(e) return(NA))
    
    smd <- sqrt(abs(dist_sq))
    smd <- as.numeric(smd)
  }
  
  return(smd)
}

risk_model_table <- function(model_cox,
                             predictor_labels,
                             horizon,
                             digits = 2) {
  # ── Validate labels match model terms ───────────────────────────────────────
  model_terms <- broom::tidy(model_cox) |> dplyr::pull(term)
  
  # ── Coefficient and HR table ────────────────────────────────────────────────
  predictor_rows <- dplyr::left_join(
    broom::tidy(model_cox, exponentiate = FALSE, conf.int = TRUE),
    broom::tidy(model_cox, exponentiate = TRUE,  conf.int = TRUE),
    by     = "term",
    suffix = c("_log", "_hr")
  ) |>
    dplyr::mutate(
      Predictor = predictor_labels,
      coef_CI   = fmt_ci(estimate_log, conf.low_log, conf.high_log, digits = digits),
      HR_CI     = fmt_ci(estimate_hr,  conf.low_hr,  conf.high_hr,  digits = digits),
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
  
  return(list(h0 = h0, coef = coefficients(model_cox), centers = model_cox$means, risk_model_table = risk_model_table))
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
  results_df["Sample size", ] <- c(sum(data[[trt_var]] == 0),
                                   sum(data[[trt_var]] == 1))
  
  # Number of events
  results_df["Number of events", ] <- c(
    sum(data[[outcome_var]] == 1 & data[[trt_var]] == 0),
    sum(data[[outcome_var]] == 1 & data[[trt_var]] == 1)
  )
  
  # Absolute risks
  results_df[paste("Risk, % (95% CI)", w_meth), ] <- c(
    fmt_ci(out_est$R0 * 100, out_est$R0_lower * 100, out_est$R0_upper * 100),
    fmt_ci(out_est$R1 * 100, out_est$R1_lower * 100, out_est$R1_upper * 100)
  )
  
  # Risk difference
  results_df[paste("Risk difference, % (95% CI)", w_meth), ] <- c(
    "Reference",
    fmt_ci(out_est$RD * 100, out_est$RD_lower * 100, out_est$RD_upper * 100)
  )
  
  # Risk ratio
  results_df[paste("Risk ratio (95% CI)", w_meth), ] <- c(
    "Reference",
    fmt_ci(out_est$RR, out_est$RR_lower, out_est$RR_upper, 2)
  )
  
  # RMST
  results_df[paste0("RMST, ", unit, " (95% CI) ", w_meth), ] <- c(
    fmt_ci(out_est$RMST0, out_est$RMST0_lower, out_est$RMST0_upper),
    fmt_ci(out_est$RMST1, out_est$RMST1_lower, out_est$RMST1_upper)
  )
  
  # RMST difference
  results_df[paste0("\u0394RMST, ", unit, " (95% CI) ", w_meth), ] <- c(
    "Reference",
    fmt_ci(out_est$dRMST, out_est$dRMST_lower, out_est$dRMST_upper)
  )
  
  # Hazard ratio
  results_df[paste("HR (95% CI)", w_meth), ] <- c(
    "Reference",
    fmt_ci(out_est$HR, out_est$HR_lower, out_est$HR_upper, 2)
  )
  
  return(results_df)
}