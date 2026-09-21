################################################################################
### Decision for dialysis versus conservative management
### PART 3 - Analysis
################################################################################

# set-up
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/")
set.seed(1)
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/"
)
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Results/"

# load libraries
library(data.table)
library(patchwork) # combine figures

# load functions
source("Code/utils/weighting.R")
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/compute_absolute_relative_risks.R")

# load data
load("Data/analysis_data_elig_cohort_new.Rdata")
elig_cohort <- cohort_final
load("Data/analysis_data_cohort_new.Rdata")
baseline <- cohort_final

# Set variables to include in the baseline table
listvar_main <- c(
  "age",
  "female",
  "Davies_score",
  "egfr2021",
  "egfr_cat",
  "sbp",
  "dbp",
  "calcium_total",
  "phosphate",
  "albumin",
  "hb",
  "cancer",
  "ihd",
  "pvd",
  "hf",
  "dm",
  "psycho",
  "hyperten",
  "aki",
  "bblock",
  "calblock",
  "diuretic",
  "rasi",
  "lipid",
  "phosbinder",
  "esa",
  "antiplatelet",
  "anticoag",
  "iron_cat",
  "n_hospital",
  "edu"
)

listvar <- c(
  "age",
  "age_cat",
  "female",
  "Davies_score",
  "Davies_score_cat",
  "region",
  "clinic_level",
  "calendar_year_cat",
  "egfr2021",
  "egfr_cat",
  "sbp",
  "sbp_cat",
  "dbp",
  "dbp_cat",
  "calcium_total",
  "phosphate",
  "albumin",
  "hb",
  "crp",
  "prd_cat",
  "cancer",
  "ihd",
  "pvd",
  "hf",
  "dm",
  "scvd",
  "copd",
  "cirr",
  "psycho",
  "acs",
  "hyperten",
  "vhd",
  "cevd",
  "af",
  "arrh",
  "lung",
  "thrombo",
  "liver",
  "fracture",
  "aki",
  "bblock",
  "calblock",
  "diuretic",
  "rasi",
  "lipid",
  "phosbinder",
  "esa",
  "vitamind",
  "digoxin",
  "vasodilator",
  "antiplatelet",
  "anticoag",
  "iron_cat",
  "n_hospital",
  "n_cvd_hospital",
  "edu"
)

# define which variables are continuous and which are categorical
contvar <- c(
  "age",
  "egfr2021",
  "Davies_score",
  "sbp",
  "dbp",
  "calcium_total",
  "phosphate",
  "albumin",
  "hb",
  "crp",
  "n_hospital",
  "n_cvd_hospital"
)
catvar <- listvar[!listvar %in% contvar]

# define labels
treatment_label <- "Dialysis"
control_label <- "Conservative management"
non_normal_vars <- c("age",
                     "Davies_score",
                     "egfr2021",
                     "n_hospital",
                     "n_cvd_hospital")

################################################################################
### Propensity score model #####################################################
################################################################################
model_PS <- trt ~ rms::pol(age, 2) + age_cat +
  female +
  rms::pol(egfr2021, 2) + egfr_cat +
  Davies_score_cat +
  calendar_year_cat +
  sbp + sbp_cat +
  dbp + dbp_cat +
  calcium_total +
  phosphate +
  albumin +
  hb +
  log(crp + 1) + 
  prd_cat +
  cancer +
  ihd +
  pvd +
  hf +
  dm +
  scvd +
  copd +
  cirr +
  psycho +
  acs +
  hyperten +
  vhd +
  cevd +
  af +
  arrh +
  lung +
  thrombo +
  liver +
  fracture +
  aki +
  bblock +
  calblock +
  diuretic +
  rasi +
  lipid +
  phosbinder +
  esa +
  vitamind +
  digoxin +
  vasodilator +
  antiplatelet +
  anticoag +
  n_hospital +
  n_cvd_hospital +
  edu +
  iron_cat +
  clinic_level +
  region +
  vasodilator:rms::pol(age, 2) +
  region:rms::pol(age, 2) +
  aki:rms::pol(egfr2021, 2)
  
# Create IPTW weights on full cohort
id_name <- "LOPNR"
trt_var <- "trt"
out_weights <- create_weights(
  data = baseline,
  id_name = id_name,
  trt_var = trt_var,
  model_PS = model_PS,
  w_meth = "IPTW",
  verbose = FALSE
)

# add weights to data
baseline <- out_weights$data

# check if SMDs after weighting below 0.1
table_one <- create_baseline_table(
  data = baseline,
  id_name = id_name,
  weights = baseline$w,
  vars = listvar,
  categoricalVars = catvar,
  continuousVars = contvar,
  IQRVars = non_normal_vars,
  treatmentColumn = trt_var,
  treatmentLabel = treatment_label,
  controlLabel = control_label,
  tableCaption = ""
)
cat("Number of SMDs > 0.1: \n")
print(table_one$smd_table[which(table_one$smd_table > 0.1)])

# save coefficients of PS model
coef_PS_overall <- out_weights$coef_ps
PS_df <- data.frame(fmt_ci(
  as.numeric(out_weights$coef_ps),
  as.numeric(out_weights$lower_ps),
  as.numeric(out_weights$upper_ps),
  digits = 2
),
fmt_ci(exp(as.numeric(
  out_weights$coef_ps
)), exp(as.numeric(
  out_weights$lower_ps
)), exp(as.numeric(
  out_weights$upper_ps
))))
rownames(PS_df) <- c(
  "Intercept",
  "Age",
  "Age^2",
  "Age category 70-74 versus 65-69",
  "Age category 75-79 versus 65-69",
  "Age category >=80 versus 65-69",
  "Female",
  "eGFR",
  "eGFR^2",
  "eGFR category 10-14 versus < 10",
  "eGFR category 15-20 versus < 10",
  "Davies score category 2-4 versus < 2",
  "Davies score category >4 versus < 2",
  "Calendar year 2013-2017 versus 2007-2012",
  "Calendar year 2018-2021 versus 2007-2012",
  "SBP",
  "SBP category 120-139 versus <120",
  "SBP category 140-159 versus <120",
  "SBP category >160 versus <120",
  "DBP",
  "DBP category 80-89 versus <80",
  "DBP category 90-99 versus <80",
  "DBP category >100 versus <80",
  "Total calcium, mmol/L",
  "Phosphorus, mmol/L",
  "Albumin, g/L",
  "Haemoglobin , mmol/L",
  "C-reactive protein",
  "Primary kidney disease hypertension versus diabetic nephropathy",
  "Primary kidney disease other versus diabetic nephropathy",
  "Malignancy",
  "Ischemic heart disease",
  "Peripheral vascular disease",
  "Heart failure",
  "Diabetes mellitus",
  "Systemic collagen vascular disease",
  "COPD",
  "Cirrhosis",
  "Psychiatric illness",
  "Acute coronary syndrome",
  "Hypertension",
  "Valvular heart disease",
  "Other cerebrovascular disease",
  "Atrial fibrillation",
  "Other arrhythmias",
  "Other lung disease",
  "Venous thromboembolism",
  "Liver disease",
  "Fracture",
  "Acute kidney injury",
  "Beta blockers",
  "Calcium channel blockers",
  "Diuretics",
  "Renin-angiotensin system inhibitors",
  "Lipid lowering agents",
  "Phosphate binder",
  "Erythropoietin stimulating agents",
  "Vitamin D",
  "Digoxin",
  "Vasodilator",
  "Antiplatelet agents",
  "Anticoagulants",
  "Number of hospitalizations in previous year",
  "Number of cardiovascular hospitalizations in previous year",
  "Education",
  "Iron intravenous versus no iron",
  "Iron per os versus no iron",
  "Clinic level regional versus local",
  "Clinic level academic versus local",
  "Other versus Örebro/Uppsala",
  "Södra versus Örebro/Uppsala",
  "Stockholm versus Örebro/Uppsala",
  "Västra versus Örebro/Uppsala",
  "Age * vasodilator",
  "Age^2 * vasodoliator",
  "Age * Other regions",
  "Age^2 * Other regions",
  "Age * Sodra",
  "Age^2 * Sodra",
  "Age * Stockholm",
  "Age^2 * Stockholm",
  "Age * Vastra",
  "Age^2 * Vastra",
  "eGFR * Acute kidney injury",
  "eGFR^2 * Acute kidney injury"
)
colnames(PS_df) <- c("Coefficients (95% CI)", "HR (95% CI)")
openxlsx::write.xlsx(
  PS_df,
  rowNames = TRUE,
  file = paste0(results_path, "Other/Table_Propensity_Score_Model.xlsx")
)

# variance inflation factor
car::vif(glm(model_PS, family = "binomial", data = baseline), type = "terms")

# C-statistic of propensity score model is 0.86
pROC::auc(pROC::roc(
  baseline$trt,
  predict(glm(
    model_PS, family = "binomial", data = baseline
  )),
  levels = c(0, 1),
  direction = "<"
))

################################################################################
### Compute weights using several techniques ###################################
################################################################################
# compute for all weighting methods
# 1) IPTW (ATE)
# 2) SMR (ATT)
# 3) Fine stratification weights (ATE)
# 4) Fine stratification weights (ATT)
# 5) Overlap weighting
# 6) IPTW and probability of belonging to eligibility cohort (generalizability)
w_meths <- c("unweighted",
             "IPTW",
             "overlap",
             "IPSW",
             "IPSW_IPTW",
             "SMR_ATT",
             "SMR_ATU")
model_S <- update(model_PS, S ~ . - rms::pol(age, 2):vasodilator)
for (w_meth in w_meths[-1]) {
  out_weights <- create_weights(
    data = baseline,
    id_name = id_name,
    trt_var = trt_var,
    elig_cohort = elig_cohort,
    model_PS = model_PS,
    model_S = model_S,
    w_meth = w_meth,
    catvar = catvar,
    contvar = contvar,
    verbose = FALSE
  )
  baseline[[paste0("sw_", w_meth)]] <- as.numeric(out_weights$data$w)
}

# Summarize statistics on weights across methods
summ_weights <- data.table::rbindlist(lapply(w_meths[-1], function(w_meth) {
  # build a small data.table with trt and the corresponding weight
  weights_dt <- baseline[, .(trt, weights = get(paste0("sw_", w_meth)))]
  
  # summarize weights (data.table with Statistic, Overall, Control, Treated)
  summ_dt <- summarize.weights(weights_dt)
  
  # add Method column for this iteration
  summ_dt[, Method := w_meth]
  
  # reorder columns: Method first
  data.table::setcolorder(summ_dt,
                          c("Method", "Statistic", "Overall", "Control", "Treated"))
}), fill = TRUE)

# Write to Excel (one worksheet)
openxlsx::write.xlsx(
  x = summ_weights,
  file = file.path(results_path, "Supplemental/Table_S_weights.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

# set sw_ISPW to 1 for those not in baseline, and to IPSW for those in baseline
idx <- match(elig_cohort$LOPNR, baseline$LOPNR) # gives NA for S = 0 rows
elig_cohort[, sw_IPSW := data.table::fifelse(is.na(idx), 1, baseline$sw_IPSW[idx])]

# save cohort
save(
  id_name,
  listvar,
  listvar_main,
  catvar,
  contvar,
  non_normal_vars,
  treatment_label,
  control_label,
  baseline,
  model_PS,
  model_S,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  trt_var,
  file = file.path("Data/cohort_with_weights.Rdata")
)
