################################################################################
### Decision for dialysis versus conservative management
### PART 7 - Example ITE calculation for random patient
################################################################################

# set-up
rm(list = ls(all.names = TRUE))
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/"
)
set.seed(1)

# load library
library(data.table)

# load data
load("Data/cohort_with_prob.Rdata")

################################################################################
### select patient
################################################################################
nr_patient <- 1
patients_risk <- baseline[, c("age",
                              "egfr2021",
                              "cancer",
                              "dm",
                              "ihd",
                              "vhd",
                              "pvd",
                              "female",
                              "albumin",
                              "log_crp")] 
patient_risk <- as.numeric(as.matrix(patients_risk[nr_patient, ]))

patients_ITE <- baseline[, c("trt", "lp_risk")]
patient_ITE <- as.numeric(as.matrix(patients_ITE[nr_patient, ]))
patient_ITE <- c(patient_ITE, patient_ITE[1] * patient_ITE[2])

################################################################################
### fit risk model
################################################################################
elig_Surv <- survival::Surv(elig_cohort$time2event_death_2y, elig_cohort$event_death_2y)

# extract coefficients
coef_risk <- table_risk$coef

# extract means
centers_risk <- table_risk$centers

# manual calculation linear predictor
manual_lp_risk <- sum(coef_risk * (patient_risk - centers_risk))
manual_lp_risk

# check using predict
risk_model <- survival::coxph(
  elig_Surv ~ age + egfr2021 + cancer + dm + ihd +
    vhd + pvd + female + albumin + log_crp,
  data = elig_cohort,
  method = c("breslow"),
  y = TRUE,
  x = TRUE
)
predict(risk_model, newdata = baseline[nr_patient, ], type = "lp")

# baseline hazard
bh_risk <- survival::basehaz(risk_model)
h0_risk <- bh_risk$hazard[bh_risk$time == horizon]
h0_risk
table_risk$h0

# predicted probability
manual_prob_risk <- 1 - exp(-h0_risk * exp(manual_lp_risk))
manual_prob_risk

# check with function
PredictionTools::fun.event(h0 = h0_risk, lp = manual_lp_risk)

################################################################################
### ITE predicted risk model
################################################################################
# ITE model
coef_ITE <- ITE_model_lp$coef

# centers
centers_ITE <- ITE_model_lp$centers

# manual calculation linear predictor
manual_lp_ITE <- sum(coef_ITE * (patient_ITE - centers_ITE))
manual_lp_ITE

# check using predict
ITE_model <- survival::coxph(
  survival::Surv(time2event_death_2y, event_death_2y) ~
    trt * lp_risk,
  data = baseline,
  method = c("breslow"),
  weights = baseline$sw_IPTW,
  x = TRUE,
  y = TRUE
)
predict(ITE_model, newdata = baseline[nr_patient, ], type = "lp")

# baseline hazard
bh_ITE <- survival::basehaz(ITE_model)
h0_ITE <- bh_ITE$hazard[bh_ITE$time == horizon]
h0_ITE
ITE_model_lp$h0

# predicted probability
manual_prob_ITE <- 1 - exp(-h0_ITE * exp(manual_lp_ITE))
manual_prob_ITE

# check with function
PredictionTools::fun.event(h0 = h0_ITE, lp = manual_lp_ITE)

# check with original code
data_1 <- data.frame(trt = as.factor(1), lp_risk = patient_ITE[2])
sf_1 <- survival::survfit(ITE_model, newdata = data_1)
risk_1 <- 1 - summary(sf_1, times = horizon)$surv

# counterfactual patient
patient_ITE_counterfactual <- data.frame(trt = 0,
                                         lp_risk = patient_ITE[2],
                                         `trt * lp_risk` = 0)

# counterfactual lp
manual_lp_ITE_counterfactual <- sum(coef_ITE * (patient_ITE_counterfactual - centers_ITE))
manual_lp_ITE_counterfactual

# check lp
patient_ITE_counterfactual$trt <- as.factor(patient_ITE_counterfactual$trt)
predict(ITE_model, newdata = patient_ITE_counterfactual, type = "lp")
patient_ITE_counterfactual$trt <- 0

# counterfactual predicted probability
manual_prob_ITE_counterfactual <- 1 - exp(-h0_ITE * exp(manual_lp_ITE_counterfactual))
manual_prob_ITE_counterfactual

# check with function
PredictionTools::fun.event(h0 = h0_ITE, lp = manual_lp_ITE_counterfactual)

# check with original code
data_0 <- data.frame(trt = as.factor(0), lp_risk = patient_ITE[2])
sf_0 <- survival::survfit(ITE_model, newdata = data_0)
risk_0 <- 1 - summary(sf_0, times = horizon)$surv

# risk difference
RD <- manual_prob_ITE - manual_prob_ITE_counterfactual
RD

# check
as.numeric(risk_1 - risk_0)

# dRMST
rmst_1 <- as.numeric(summary(sf_1, rmean = horizon)$table["rmean"]) / 30.5
rmst_0 <- as.numeric(summary(sf_0, rmean = horizon)$table["rmean"]) / 30.5

# hazard ratio
HR <- exp(manual_lp_ITE - manual_lp_ITE_counterfactual)
HR

# check
exp(ITE_model_lp$coef["trt1"])

cat(
  "Patient is",
  baseline[nr_patient, ][["age"]],
  "years old \n",
  "eGFR                                 :",
  baseline[nr_patient, ][["egfr2021"]],
  "ml/min/m3 \n",
  "Albumin                              :",
  baseline[nr_patient, ][["albumin"]],
  "g/L \n",
  "CRP mg/L                             :",
  baseline[nr_patient, ][["crp"]],
  "mg/L \n",
  "Log(CRP+1)                           :",
  baseline[nr_patient, ][["log_crp"]],
  "\n",
  "Linear predictor of risk             :",
  paste0(round(coef_risk, 3), " * (", round(patient_risk, 3), " - ", round(centers_risk, 3), ")", collapse = " + "),
  "\n",
  "Linear predictor value               :",
  round(manual_lp_risk, 3),
  "\n",
  "Risk probability                     :",
  paste0("1-exp(-", round(h0_risk, 3), " * exp(", round(manual_lp_risk, 3), "))"),
  "\n",
  "Risk probability                     :",
  round(manual_prob_risk, 3),
  "\n",
  "Linear predictor under dialysis      :",
  paste0(round(coef_ITE, 3), " * (", round(patient_ITE, 3), " - ", round(centers_ITE, 3), ")", collapse = " + "),
  "\n",
  "Linear predictor under dialysis      :",
  round(manual_lp_ITE, 3),
  "\n",
  "Risk probability - dialysis          :",
  paste0("1-exp(-", round(h0_ITE, 3), " * exp(", round(manual_lp_ITE, 3), "))"),
  "\n",
  "Risk probability - dialysis          :",
  round(manual_prob_ITE, 3),
  "\n",
  "Linear predictor under CC            :",
  paste0(round(coef_ITE, 3), " * (", round(patient_ITE_counterfactual, 3), " - ", round(centers_ITE, 3), ")", collapse = " + "),
  "\n",
  "Linear predictor under CC            :",
  round(manual_lp_ITE_counterfactual, 3),
  "\n",
  "Risk probability - CC                :",
  paste0("1-exp(-", round(h0_ITE, 3), " * exp(", round(manual_lp_ITE_counterfactual, 3), "))"),
  "\n",
  "Risk probability - CC                :",
  round(manual_prob_ITE_counterfactual, 3),
  "\n",
  "Risk difference                      :",
  round(RD, 3),
  "\n",
  "RMST under dialysis                  :",
  round(rmst_1, 1),
  "\n",
  "RMST under CC                        :",
  round(rmst_0, 1),
  "\n",
  "Difference in RMST                   :",
  round(rmst_1 - rmst_0, 1),
  "\n",
  "Hazard ratio                         :",
  paste0("exp(", round(manual_lp_ITE, 3), " - ", round(manual_lp_ITE_counterfactual, 3), ")"),
  "\n",
  "Hazard ratio                         :",
  round(HR, 2),
  "\n"
)
