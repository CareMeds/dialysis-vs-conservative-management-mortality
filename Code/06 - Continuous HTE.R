################################################################################
### Decision for dialysis versus conservative management
### PART 6 - Continuous heterogeneous treatment effect estimation
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
library(survival)
library(rms)
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
load("Data/cohort_with_models.Rdata")

################################################################################
### Internal 2-year time-to-event risk model
################################################################################
# perform internal validation 
validate <- FALSE

# make outcome for elig cohort
elig_Surv <- survival::Surv(elig_cohort$time2event_death_2y, elig_cohort$event_death_2y)

# use external predictors identified by Chava
elig_cohort[, log_crp := log(crp + 1)]
baseline[, log_crp := log(crp + 1)]
risk_model_cox <- survival::coxph(
  elig_Surv ~ age + egfr2021 + cancer + dm + ihd +
    vhd + pvd + female + albumin + log_crp,
  data = elig_cohort,
  method = "breslow"
)

# save table
table_risk <- risk_model_table(
  model_cox  = risk_model_cox,
  predictor_labels = c(
    "   Age (per 1 year)",
    "   eGFR (per 1 ml/min/1.73 m²)",
    "   Malignancies",
    "   Diabetes mellitus",
    "   Ischemic heart disease",
    "   Valvular heart disease",
    "   Peripheral vascular disease",
    "   Female vs. male",
    "   Albumin (per 1 g/L)",
    "   CRP (log-transformed, per unit increase)"
  ),
  horizon      = horizon
)

# create plot on hazard scale
dd <- rms::datadist(elig_cohort)
options(datadist = "dd")
risk_model <- rms::cph(
  elig_Surv ~ age + egfr2021 + cancer + dm + ihd +
    vhd + pvd + female + albumin + log(crp + 1),
  data = elig_cohort,
  method = c("breslow"),
  y = TRUE,
  x = TRUE
)
png(
  filename = "Results/Supplemental/Figure_M_HTE_Risk_model.png",
  width = 2000,
  height = 2000,
  res = 300
)
print(plot(Predict(risk_model, fun = exp), ylab = "Hazard Ratio"))
dev.off()

# ################################################################################
# ### Internal validation
# ################################################################################
# --- Initialise storage --------------------------------------------------
if (!validate){
  n_bootstraps_elig <- 1
} else {
  n_bootstraps_elig <- 1000
}
n_iter <- n_bootstraps_elig + 1L          # iteration 1 = original sample

metrics <- list(
  orig = data.frame(
    elig_Intercept = numeric(n_iter),
    elig_Slope     = numeric(n_iter),
    elig_AUC       = numeric(n_iter),
    trt_Intercept  = numeric(n_iter),
    trt_Slope      = numeric(n_iter),
    trt_AUC        = numeric(n_iter)
  ),
  boot = data.frame(
    elig_Intercept = numeric(n_iter),
    elig_Slope     = numeric(n_iter),
    elig_AUC       = numeric(n_iter),
    trt_Intercept  = numeric(n_iter),
    trt_Slope      = numeric(n_iter),
    trt_AUC        = numeric(n_iter)
  )
)

cal_plots <- list()

# --- Main loop -----------------------------------------------------------
for (B in seq_len(n_iter)) {
  
  is_original <- (B == 1L)
  
  # 1. Draw samples -------------------------------------------------------
  if (is_original) {
    boot_elig    <- elig_cohort
    boot_trt_dec <- baseline
  } else {
    boot_elig    <- elig_cohort[sample(nrow(elig_cohort), replace = TRUE), ]
    boot_trt_dec <- baseline[sample(nrow(baseline),       replace = TRUE), ]
  }
  
  # 2. Fit risk model on bootstrap sample ---------------------------------
  boot_risk_model <- rms::cph(
    survival::Surv(time2event_death_2y, event_death_2y) ~
      age + egfr2021 + cancer + dm + ihd + vhd + pvd + female + albumin,
    data   = boot_elig,
    method = "breslow",
    y      = TRUE,
    x      = TRUE
  )
  
  # 3. Evaluate on original sample (apparent for B=1, optimism for B>1) --
  elig_orig <- compute_measures(boot_risk_model,
                                data = elig_cohort,
                                plot = is_original)
  trt_orig  <- compute_measures(boot_risk_model,
                                data = baseline,
                                plot = is_original)
  
  if (is_original) {
    cal_plots$elig <- calibration_plot(elig_orig)
    cal_plots$trt  <- calibration_plot(trt_orig)
  }
  
  # 4. Evaluate on bootstrap sample (needed only for optimism correction) -
  if (!is_original) {
    elig_boot <- compute_measures(boot_risk_model, data = boot_elig)
    trt_boot  <- compute_measures(boot_risk_model, data = boot_trt_dec)
  }
  
  # 5. Store results ------------------------------------------------------
  metrics$orig[B, ] <- c(extract_metrics(elig_orig),
                         extract_metrics(trt_orig))
  
  if (!is_original) {
    metrics$boot[B, ] <- c(extract_metrics(elig_boot),
                           extract_metrics(trt_boot))
  }
}

# --- Compute optimism (rows 2:n = bootstrap iterations) -----------------
# metrics$boot[1, ] is never filled (all zeros) — optimism correctly uses [-1, ]
optimism   <- colMeans(metrics$orig[-1, ] - metrics$boot[-1, ])
apparent   <- metrics$orig[1, ]
orig_boots <- metrics$orig[-1, ]

# --- Build and save plots -----------------------------------------------
for (cohort in list(list(prefix = "elig", label = "elig"),
                    list(prefix = "trt",  label = "trt"))) {
  
  annotated <- annotate_cal_plot(
    cal_plot_obj  = cal_plots[[cohort$prefix]],
    apparent_vals = cohort_apparent(cohort$prefix),
    boot_vals     = cohort_boots(cohort$prefix),
    optimism_vals = cohort_optimism(cohort$prefix)
  )
  
  save_cal_plot(
    cal_plot_obj   = cal_plots[[cohort$prefix]],
    annotated_plot = annotated,
    filename       = paste0("Figure_M_HTE_calibration_", cohort$label, ".png")
  )
}

################################################################################
### Explore linearity of age
################################################################################
fit_age_linear <- survival::coxph(
  survival::Surv(time2event_death_2y, event_death_2y) ~
    trt * age,
  data = baseline,
  method = c("breslow"),
  weights = baseline$sw_IPTW,
  x = TRUE,
  y = TRUE
)

# save table
ITE_model_age <- risk_model_table(
  model_cox  = fit_age_linear,
  predictor_labels = c(
    "   Dialysis",
    "   Age",
    "   Dialysis * age"
  ),
  horizon      = horizon
)

################################################################################
### Explore linearity of risk
################################################################################
baseline$lp_risk <- predict(risk_model, newdata = baseline)
fit_lp_linear <- survival::coxph(
  survival::Surv(time2event_death_2y, event_death_2y) ~
    trt * lp_risk,
  data = baseline,
  method = c("breslow"),
  weights = baseline$sw_IPTW,
  x = TRUE,
  y = TRUE
)

# save table
ITE_model_lp <- risk_model_table(
  model_cox  = fit_lp_linear,
  predictor_labels = c(
    "   Dialysis",
    "   Linear predictor",
    "   Dialysis * Linear predictor"
  ),
  horizon      = horizon
)

# save to xlsx
summmarized_table <- rbind(
  c("Risk model", "", "", ""),
  table_risk$risk_model_table,
  c("ITE model for age", "", "", ""),
  ITE_model_age$risk_model_table,
  c("ITE model for linear predictor", "", "", ""),
  ITE_model_lp$risk_model_table
)
openxlsx::write.xlsx(summmarized_table, 
                     file = paste0(results_path, "Supplemental/Table_S_HTE_risk_HR.xlsx"), 
                     rowNames = FALSE)

################################################################################
### Continuous HTE
################################################################################
# probabilities between 40% and 90% mortality risk
p_range <- c(0.2, 0.4, 0.9)

# set names of metrics
names_metrics <- c("RD", "dRMST", "RR", "HR") 

# set subgroups
data_sets <- c("baseline", "baseline[Davies_score >= 2]")

for (nr_analysis in 1:2) {
  cat("Perform analysis on",
      ifelse(nr_analysis == 1, "full data\n", "Davies >= 2 data\n"))
  # determine LP and survival probability using risk model
  analysis_data <- eval(parse(text = data_sets[nr_analysis]))
  analysis_data$lp_risk <- predict(risk_model, newdata = analysis_data)
  analysis_data$pred_risk <- PredictionTools::fun.event(h0 = table_risk$h0, 
                                                        lp = analysis_data$lp_risk)
  
  # compute estimates across 100 risk points
  for (B in 1:(n_bootstraps + 1)) {
    if (B == 1) {
      # Original sample
      bootstrap <- analysis_data
    } else{
      # create bootstrap sample
      bootstrap <- analysis_data[sample(1:nrow(analysis_data), replace = TRUE), ]
    }
    
    # re-estimate weights
    bootstrap_reestimated <- create_weights(
      data = bootstrap,
      id_name = id_name,
      model_PS = model_PS,
      w_meth = "IPTW",
      catvar = catvar,
      contvar = contvar,
      verbose = FALSE
    )
    bootstrap$sw_IPTW <- bootstrap_reestimated$data$w
    
    # check SMDs
    if (B == 1) {
      table_one_weighted <- create_baseline_table(
        data = bootstrap_reestimated$data,
        id_name = "LOPNR",
        weights = bootstrap_reestimated$data$w,
        vars = listvar,
        categoricalVars = catvar,
        continuousVars = contvar,
        IQRVars = non_normal_vars,
        treatmentColumn = trt_var,
        treatmentLabel = treatment_label,
        controlLabel = control_label,
        tableCaption = paste("Subgroup", 2)
      )
      cat("Number of SMDs > 0.1",
          sum(table_one_weighted$smd_table > 0.1),
          "\n")
    }
    
    # compute HTE across age
    show_test <- ifelse(B == 1, TRUE, FALSE)
    age_df <- compute_HTE(
      data = bootstrap,
      unit = unit,
      horizon = horizon,
      event_var = outcome_var,
      time2event_var = time2outcome_var,
      effect_modifier = "age",
      effect_modifier_range = seq(65, 95, 1),
      add_interaction = TRUE,
      test_relative_HTE = show_test
    )
    
    # compute HTE across sex
    show_test <- ifelse(B == 1, TRUE, FALSE)
    sex_df <- compute_HTE(
      data = bootstrap,
      unit = unit,
      horizon = horizon,
      event_var = outcome_var,
      time2event_var = time2outcome_var,
      effect_modifier = "female",
      effect_modifier_range = c(1, 2), # as.numeric() creates 1 and 2
      add_interaction = TRUE,
      test_relative_HTE = show_test
    )
    
    # compute HTE across predicted risk
    # get corresponding linear predictor
    lp_range <- log(-log(1 - p_range) / table_risk$h0)
    pred_risk_df <- compute_HTE(
      data = bootstrap,
      unit = unit,
      horizon = horizon,
      event_var = outcome_var,
      time2event_var = time2outcome_var,
      effect_modifier = "lp_risk",
      effect_modifier_range = sort(c(
        seq(lp_range[1], lp_range[3], length.out = 99), lp_range[2]
      )),
      add_interaction = TRUE,
      test_relative_HTE = show_test
    )
    
    # save results
    if (B == 1) {
      # report effect modifier using probabilities and not linear predictor
      pred_risk_df$effect_modifier_range <- PredictionTools::fun.event(h0 = table_risk$h0,
                                                                       lp = pred_risk_df$effect_modifier_range)
      summary(pred_risk_df$effect_modifier_range)
      
      # create seperate dt for each estimate
      dfs <- list(age = age_df, pred_risk = pred_risk_df, sex = sex_df)
      for (df_name in names(dfs)) {
        # set colnames of dt to effect_modifier_range
        for (metric in names_metrics) {
          var_name <- paste0(df_name, "_", metric)
          assign(var_name, dfs[[df_name]][, c("effect_modifier_range", metric)])
        }
        
        # set p-value
        p_raw <- unique(dfs[[df_name]]$p_for_HTE)
        p_fmt <- ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))
        assign(paste0("p_value_", df_name), p_fmt)
      }
    } else{
      # bootstrapped sample
      age_RD[, paste0("boot_", B - 1)] <- age_df$RD
      age_RR[, paste0("boot_", B - 1)] <- age_df$RR
      age_dRMST[, paste0("boot_", B - 1)] <- age_df$dRMST
      age_HR[, paste0("boot_", B - 1)] <- age_df$HR
      
      sex_RD[, paste0("boot_", B - 1)] <- sex_df$RD
      sex_RR[, paste0("boot_", B - 1)] <- sex_df$RR
      sex_dRMST[, paste0("boot_", B - 1)] <- sex_df$dRMST
      sex_HR[, paste0("boot_", B - 1)] <- sex_df$HR
      
      pred_risk_RD[, paste0("boot_", B - 1)] <- pred_risk_df$RD
      pred_risk_RR[, paste0("boot_", B - 1)] <- pred_risk_df$RR
      pred_risk_dRMST[, paste0("boot_", B - 1)] <- pred_risk_df$dRMST
      pred_risk_HR[, paste0("boot_", B - 1)] <- pred_risk_df$HR
    }
  }
  
  # add nr_analysis to dt
  for (df_name in names(dfs)) {
    for (metric in names_metrics) {
      var_name <- paste0(df_name, "_", metric)
      assign(paste0(var_name, "_", nr_analysis), get(var_name))
    }
    assign(paste0("p_value_", df_name, "_", nr_analysis), get(paste0("p_value_", df_name)))
  }
}

# compute figures and tables
HTE_table_age <- data.frame()
HTE_table_sex <- data.frame()
HTE_table_pred_risk <- data.frame()
for (nr_analysis in 1:2) {
  for (effect_modifier in c("age", "sex", "pred_risk")) {
    # create histogram of variable stratified by treatment
    analysis_data <- eval(parse(text = data_sets[nr_analysis]))
    analysis_data$lp_risk <- predict(risk_model, newdata = analysis_data)
    analysis_data$pred_risk <- PredictionTools::fun.event(h0 = table_risk$h0, 
                                                          lp = analysis_data$lp_risk)
    
    hist_stratified <- create_histogram_stratified(
      dt = analysis_data,
      var_name = effect_modifier,
      trt_name = trt_var,
      manual_colors = manual_colors
    )
    
    for (measure in names_metrics) {
      # legend of histogram
      if (effect_modifier == "pred_risk" & measure == "dRMST") {
        # add legend title to dRMST plot
        hist_stratified <- hist_stratified +
          ggplot2::theme(legend.text = ggplot2::element_text(size = 18))
      } else {
        # remove legends from other plots
        hist_stratified <- hist_stratified +
          ggplot2::theme(legend.position = "none")
      }
      
      # Retrieve the object from the environment
      df_name <- paste0(effect_modifier, "_", measure, "_", nr_analysis)
      estimates_df <- get(df_name)
      
      # Calculate quantiles
      # Assuming the first column(s) are the bootstrap replicates to be excluded
      estimates_df[, paste0(measure, "_lower")] <- apply(estimates_df[, -c(1:2)],
                                                         1,
                                                         quantile,
                                                         probs = 0.025,
                                                         na.rm = TRUE)
      estimates_df[, paste0(measure, "_upper")] <- apply(estimates_df[, -c(1:2)],
                                                         1,
                                                         quantile,
                                                         probs = 0.975,
                                                         na.rm = TRUE)
      
      # save table
      if (effect_modifier == "age") {
        # every round-number age already evaluated (min to max: seq(65,95,1))
        sel_rows <- order(estimates_df$effect_modifier_range)
        range <- estimates_df[sel_rows, "effect_modifier_range"]
      } else if (effect_modifier == "sex") {
        sel_rows <- c(
          which(estimates_df$effect_modifier_range == 1),
          which(estimates_df$effect_modifier_range == 2)
        )
        range <- estimates_df[sel_rows, "effect_modifier_range"]
      } else if (effect_modifier == "pred_risk") {
        # approximate pred_risk_step_pct steps by rounding the already-computed
        # grid (100 points spread evenly in the linear-predictor scale, not
        # evenly in probability) to the nearest multiple of that step, keeping
        # - for each distinct step value from min to max - the one
        # already-evaluated row closest to that exact value. Steps with no
        # nearby evaluated point are simply absent from the table (a gap),
        # since no new grid points are computed here.
        pred_risk_step_pct <- 5
        pct <- round(estimates_df$effect_modifier_range * 100 / pred_risk_step_pct) * pred_risk_step_pct
        sel_rows <- sapply(sort(unique(pct)), function(target) {
          candidates <- which(pct == target)
          candidates[which.min(abs(estimates_df$effect_modifier_range[candidates] - target / 100))]
        })
        range <- paste0(pct[sel_rows], "%")
      }
      est <- sprintf(ifelse(measure == "RR" |
                              measure == "HR", "%.2f", "%.1f"),
                     estimates_df[sel_rows, measure])
      CI <- paste0(
        "(",
        sprintf(
          ifelse(measure == "RR" | measure == "HR", "%.2f", "%.1f"),
          estimates_df[sel_rows, paste0(measure, "_lower")]
        ),
        ", ",
        sprintf(
          ifelse(measure == "RR" | measure == "HR", "%.2f", "%.1f"),
          estimates_df[sel_rows, paste0(measure, "_upper")]
        ),
        ")"
      )
      
      if (measure == "RD") {
        HTE_tab <- cbind(range, est, CI)
      } else{
        HTE_tab <- cbind(HTE_tab, est, CI)
      }
      
      # make effect plot
      if (effect_modifier != "sex") {
        effect_plot_out <- effect_plot(
          estimates_df = estimates_df,
          effect_modifier = effect_modifier,
          y_middle = ifelse(measure == "HR" | measure == "RR", 1, 0),
          measure = measure,
          y_min_RD = -70,
          y_max_RD = ifelse(nr_analysis == 1, 21, 56),
          y_min_RR = 0,
          y_max_RR = ifelse(nr_analysis == 1, 1.3, 1.8),
          y_min_dRMST = ifelse(nr_analysis == 1, -3, -8),
          y_max_dRMST = 10,
          y_min_HR = 0,
          y_max_HR = ifelse(nr_analysis == 1, 1.3, 1.8),
          show_favor_annotation = (measure == "RD")
        )
        
        # control x-axis of effect plot
        effect_plot_out <- effect_plot_out + 
          ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank()
          )
        if (effect_modifier == "age") {
          effect_plot_out <- effect_plot_out +
            ggplot2::scale_x_continuous(limits = c(60, 100),
                                        breaks = seq(60, 100, 5))
        } else{
          effect_plot_out <- effect_plot_out +
            ggplot2::scale_x_continuous(
              limits = c(0, 1),
              breaks = seq(0, 1, 0.1),
              labels = seq(0, 100, 10) # scales::percent_format(accuracy = 1)
            )
        }
        
        # combine the effect and histogram (x-axis title on the histogram is
        # blanked here since a single shared x-axis label is added under each
        # full row of 4 plots further below)
        combined_effect_hist <- (
          effect_plot_out /
            (hist_stratified + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
        ) +
          plot_layout(heights = c(1, 0.2))
        
        assign(paste0(effect_modifier, "_", measure, "_plot"),
               combined_effect_hist)
      }
    }
    # prefix this block with a subgroup label row (blank data cells) carrying
    # the overall interaction p-value in the last column - one p-value covers
    # the whole block now that age/pred_risk have many rows each, rather than
    # a fixed 2-row block as before
    n_cols <- ncol(HTE_tab)
    subgroup_label <- if (nr_analysis == 1) "Full cohort" else "Davies comorbidity score >=2 subgroup"
    p_val <- get(paste0("p_value_", effect_modifier, "_", nr_analysis))
    header_row <- matrix(c(subgroup_label, rep("", n_cols - 1), p_val), nrow = 1)
    data_rows <- cbind(HTE_tab, rep("", nrow(HTE_tab)))
    HTE_block <- rbind(header_row, data_rows)
    
    if (effect_modifier == "age") {
      HTE_table_age <- rbind(HTE_table_age, HTE_block)
    } else if (effect_modifier == "sex") {
      HTE_table_sex <- rbind(HTE_table_sex, HTE_block)
    } else if (effect_modifier == "pred_risk") {
      HTE_table_pred_risk <- rbind(HTE_table_pred_risk, HTE_block)
    }
  }
  
  # small text-only panels used as row titles
  age_row_title <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Effect modification by age") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 24))
  
  pred_risk_row_title <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Effect modification by 2-year mortality risk (%)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 24))
  
  # small text-only panels used as one shared x-axis label per row
  # (the individual columns' own x-axis titles are blanked above)
  age_x_label <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Age in years") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 18))
  
  pred_risk_x_label <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Predicted 2-year mortality risk (%)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 18))
  
  # save figure
  ggplot2::ggsave(
    plot = (
      age_row_title /
        (age_RD_plot | age_RR_plot | age_dRMST_plot | age_HR_plot) /
        age_x_label /
        pred_risk_row_title /
        (
          pred_risk_RD_plot | 
            pred_risk_RR_plot |
            pred_risk_dRMST_plot | 
            pred_risk_HR_plot
        ) /
        pred_risk_x_label
    ) +
      patchwork::plot_layout(heights = c(0.05, 1, 0.04, 0.05, 1, 0.04)),
    filename = paste0(
      results_path,
      ifelse(
        nr_analysis == 1,
        "Main/Figure_4.pdf",
        "Supplemental/Figure_S_HTE_DCS.png"
      )
    ),
    width = 20,
    height = 15,
    dpi = 300
  )
}

# save tables (age and pred_risk now span their full range; p-values are
# already embedded per-block via the subgroup label row built above, so no
# separate p-value vector needs assembling here)
colnames(HTE_table_age) <- c("Age",
                             "RD",
                             "95% CI",
                             "dRMST",
                             "95% CI",
                             "RR",
                             "95% CI",
                             "HR",
                             "95% CI",
                             "p-value")
colnames(HTE_table_sex) <- c("Sex",
                             "RD",
                             "95% CI",
                             "dRMST",
                             "95% CI",
                             "RR",
                             "95% CI",
                             "HR",
                             "95% CI",
                             "p-value")
colnames(HTE_table_pred_risk) <- c("Predicted risk",
                                   "RD",
                                   "95% CI",
                                   "dRMST",
                                   "95% CI",
                                   "RR",
                                   "95% CI",
                                   "HR",
                                   "95% CI",
                                   "p-value")
openxlsx::write.xlsx(
  HTE_table_age,
  rowNames = FALSE,
  file = paste0(results_path, "Supplemental/Table_S_HTE_age.xlsx")
)
openxlsx::write.xlsx(
  HTE_table_sex,
  rowNames = FALSE,
  file = paste0(results_path, "Other/Table_S_HTE_sex.xlsx")
)
openxlsx::write.xlsx(
  HTE_table_pred_risk,
  rowNames = FALSE,
  file = paste0(results_path, "Supplemental/Table_S_HTE_pred_risk.xlsx")
)

# save variables
save(
  id_name,
  listvar,
  listvar_main,
  catvar,
  contvar,
  non_normal_vars,
  outcome_var,
  time2outcome_var,
  competing_events_var,
  treatment_label,
  control_label,
  baseline,
  model_PS,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  trt_var,
  horizon,
  unit,
  n_bootstraps,
  manual_colors,
  estimates_df,
  metrics,
  table_risk,
  ITE_model_lp,
  ITE_model_age,
  age_RD,
  age_RR,
  age_dRMST,
  age_HR,
  sex_RD,
  sex_RR,
  sex_dRMST,
  sex_HR,
  pred_risk_RD,
  pred_risk_RR,
  pred_risk_dRMST,
  pred_risk_HR,
  file = file.path("Data/cohort_with_prob.Rdata")
)