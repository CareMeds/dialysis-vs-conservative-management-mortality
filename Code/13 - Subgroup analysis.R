################################################################################
### Decision for dialysis versus conservative management
### PART 13 - Subgroup analysis
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

# -----------------------------
# Build a results table for one subgroup, in the same row-by-row
# format as the main results_df table in 05_-_ATE.R (Sample size,
# Number of events, Risk, Risk difference, Risk ratio, RMST, dRMST, HR).
# -----------------------------
for (w_meth in w_meths[1:2]) {
  # Set weights: use 1 for unweighted, otherwise use specified weights
  if (w_meth == "unweighted") {
    weights_meth <- rep(1, nrow(baseline))
  } else {
    weights_meth <- baseline[[paste0("sw_", w_meth)]]
  }
  
  for (subgroup in c("age", "dialysis_type")) {
    if (subgroup == "age") {
      sub_baseline <- list(young = baseline[baseline$age < 80], old = baseline[baseline$age >= 80])
    } else {
      sub_baseline <- list(HD = baseline[trt == 0 | dialysis_type == "HD"], 
                           PD = baseline[trt == 0 | dialysis_type == "PD"])
    }
    
    for (sub_name in names(sub_baseline)) {
      sub_data <- sub_baseline[[sub_name]]
      
      # compute estimates
      out_est <- compute_estimates_with_CI(
        data = sub_data,
        id_name = id_name,
        unit = unit,
        horizon = horizon,
        elig_cohort = elig_cohort,
        model_PS = model_PS,
        event_var = outcome_var,
        competing_event_var = competing_events_var,
        time2event_var = time2outcome_var,
        trt_var = trt_var,
        w_meth = w_meth,
        trim_meth = "no_trimming",
        catvar = catvar,
        contvar = contvar,
        n_bootstraps = n_bootstraps,
        bootstrap_seed = 1
      )
      assign(paste0("out_est_", w_meth, "_", subgroup, "_", sub_name),
             out_est)
      
      # Build a results table for this subgroup, matching file 05's Table_2 format
      sub_table <- build_results_table(
        data = sub_data,
        out_est = out_est,
        label = paste0(subgroup, ": ", sub_name),
        outcome_var = outcome_var,
        trt_var = trt_var,
        control_label = control_label,
        treatment_label = treatment_label,
        w_meth = w_meth,
        unit = unit
      )
      assign(paste0("table_", w_meth, "_", subgroup, "_", sub_name),
             sub_table)
    }
  }
}

openxlsx::write.xlsx(
  rbind(
    cbind(table_unweighted_age_old[-2, ], table_IPTW_age_old[-2, ]),
    cbind(table_unweighted_age_young[-2, ], table_IPTW_age_young[-2, ])
  ),
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/Table_S_HTE_subgroup_age.xlsx")
)

openxlsx::write.xlsx(
  rbind(
    cbind(
      table_unweighted_dialysis_type_HD[-2, ],
      table_IPTW_dialysis_type_HD[-2, ]
    ),
    cbind(
      table_unweighted_dialysis_type_PD[-2, ],
      table_IPTW_dialysis_type_PD[-2, ]
    )
  ),
  rowNames = TRUE,
  file = paste0(
    results_path,
    "Supplemental/Table_S_HTE_subgroup_dialysis_type.xlsx"
  )
)
