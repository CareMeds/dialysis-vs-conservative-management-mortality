# Effect of Choosing Dialysis versus Conservative Care on Survival in Older Adults with Advanced CKD: Nationwide Target Trial Emulation Study

#### Carolien C.H.M. Maas, PhD [1,2], Xuerui Zhang, MD [1], Juan-Jesus Carrero, PhD [2,3], Antoine Créon, MD, MSc [2], Roosa Lankinen, MD, PhD [2], Ilaria Prosepe, PhD [4], Friedo W. Dekker, PhD [1], Willem Jan W. Bos, MD, PhD [5,6], Marie Evans, MD, PhD [7,8], Edouard L. Fu, PhD [1,2]

[1] Department of Clinical Epidemiology, Leiden University Medical Center, Leiden, The Netherlands<br>
[2] Department of Medical Epidemiology and Biostatistics, Karolinska Institutet, Stockholm, Sweden<br>
[3] Department of Clinical Sciences, Danderyd University Hospital, Stockholm, Sweden<br>
[4] Department of Biomedical Data Sciences, Leiden University Medical Center, Leiden, The Netherlands<br>
[5] Department of Internal Medicine, St. Antonius Ziekenhuis, Nieuwegein, The Netherlands<br>
[6] Department of Internal Medicine, Leiden University Medical Center, Leiden, The Netherlands<br>
[7] Department of Clinical Science, Intervention and Technology, Division of Renal Medicine, Karolinska Institutet, Stockholm, Sweden<br>
[8] Karolinska University Hospital, Stockholm, Sweden

---

## About this repository

This repository contains the R code for a target trial emulation comparing the decision to start **dialysis** versus **conservative management (CM)** in older adults with advanced chronic kidney disease (CKD), using nationwide Swedish registry data. The primary outcome is all-cause mortality over a 2-year horizon. Effects are estimated with inverse probability of treatment weighting (IPTW), with additional weighting for generalizability to the full eligible population, and are reported as risks, risk differences, risk ratios, restricted mean survival time (RMST) and hazard ratios.

### Data availability

The individual-level registry data cannot be shared publicly due to Swedish data protection regulations. The code is provided for transparency and reproducibility; it will not run without access to the underlying data.

## Repository structure

```
Code/
├── Run all.R                                   # runs the full pipeline (00 → 13)
├── 00 - Data preparation.R
├── 01 - Apply eligibility criteria.R
├── 02 - Covariate and outcome derivation.R
├── 03 - Compute weights.R
├── 04 - Descriptives.R
├── 05 - ATE.R
├── 06 - Continuous HTE.R
├── 07 - Example Supplemental Methods.R
├── 08 - SA positivity.R
├── 09 - Competing risk time-to-dialysis.R
├── 10 - SA unmeasured confounding.R
├── 11 - Supplemental Figure HTE.R
├── 12 - SA home time.R
├── 13 - Subgroup analysis.R
└── utils/                                      # helper functions sourced by the scripts
    ├── competing_risk.R
    ├── compute_absolute_relative_risks.R
    ├── data_manipulation.R
    ├── outcome_derivation.R
    ├── plots.R
    ├── tables.R
    └── weighting.R
```

## How to run

1. Set the paths at the top of each script (`setwd()`, `results_path`, and the `source()` paths in `Run all.R`). These currently point to the project folder on the secure research server.
2. The working directory is expected to contain a `Data/` folder (raw and intermediate `.Rdata` files) and a `Results/` folder with `Main/` and `Supplemental/` subfolders.
3. Run `Code/Run all.R` to execute the full pipeline, or run the numbered scripts individually in order. Each script clears the workspace, loads the intermediate data it needs from `Data/`, and writes its tables and figures to `Results/`.

Scripts 00–03 produce the intermediate datasets (e.g., `cohort_with_weights.Rdata`, `cohort_with_prob.Rdata`) that all later scripts depend on.

### R packages

Main packages used: `data.table`, `tidyverse` (`dplyr`, `tidyr`, `ggplot2`, `lubridate`), `survival`, `rms`, `mstate`, `prodlim`, `riskRegression`, `WeightIt`, `sbw`, `survey`, `geepack`, `tableone`, `pROC`, `broom`, `PredictionTools`, `foreach`, `doRNG`, `doParallel`, `openxlsx`, `patchwork`, `cowplot`, `ggpubr`, `ggthemes`, `ggtext`, `gridExtra`, `scales`, `RColorBrewer`, `pacman`.

Bootstrap confidence intervals use 1,000 replicates (`n_bootstraps` in `05 - ATE.R`); seeds are fixed with `set.seed(1)` and `doRNG` for parallel computation.

---

## Analysis scripts

### 00 - Data preparation
Combines the data files for CKD patients into one, ensuring the right encoding for each variable.

### 01 - Apply eligibility criteria
Applies the eligibility criteria:
1. Estimated glomerular filtration rate (eGFR) <20 mL/min/1.73m²
2. Aged ≥65 with a Davies Comorbidity Score* ≥2, or aged ≥80
3. All lab measurements with a maximum one-year look-back
4. No history of kidney transplantation or dialysis
5. No history of HIV or dementia

*Davies comorbidities are malignancies, ischemic heart disease, peripheral vascular disease, heart failure, diabetes mellitus, systemic collagen vascular disease, COPD, cirrhosis, psychiatric illness, and HIV.

First, the time intervals in which patients are eligible are defined. Second, the data table is expanded for each day within the time interval. Third, if the treatment decision falls on an eligible date, the patient is included in the cohort.

### 02 - Covariate and outcome derivation
Obtains information on:
1. Other comorbidities from inpatient and outpatient files
2. Number of hospitalizations from the inpatient file
3. Medications based on ATC codes
4. ESA and iron from the CKD file
5. Primary kidney disease from the CKD file
6. Education on treatment choice from the CKD file
7. Geographical clinic level from the CKD file
8. Nursing home information from the outpatient file
9. Outcomes (e.g., all-cause mortality) from the death file

### 03 - Compute weights
Propensity score (PS) model
1. Fit the PS model
2. Check that the coefficients of the PS model are not too extreme
3. Check that the AUC of the PS model is not too high

Generalizability model
1. Define S = 1 for eligible patients with a recorded treatment decision and S = 0 for eligible patients without a recorded treatment decision
2. Fit a model with S as the outcome and the confounders X as predictors

Compute weights
1. Inverse probability of treatment weights (IPTW)
2. Overlap weights
3. Inverse probability of selection weights (IPSW)
4. Combined selection and treatment weights (IPSW × IPTW)
5. SMR weights for the ATT
6. SMR weights for the ATU

Describe weights
1. Check the minimum and maximum of the weights
2. IPTW: check that the SMD is below 0.1 for all confounders
3. Generalizability: check that the TASMD is below 0.1, comparing those with a recorded treatment decision versus all eligible patients

### 04 - Descriptives
1. Histogram of time until dialysis
2. Table of patient characteristics for the entire eligible cohort before and after selection weighting
3. Table of patient characteristics before and after treatment weighting
4. Love plots, checking that SMD and TASMD are below 0.1 after weighting
5. Summary of the number of decisions

### 05 - ATE
Computes average treatment effects:
1. Risks (Kaplan-Meier estimates)
2. Risk difference
3. Risk ratio
4. RMST (area under the Kaplan-Meier curves)
5. Difference in RMST
6. Hazard ratio

Also reports the ATT and ATU (SMR weights).

### 06 - Continuous HTE
1. Fit a risk model using age, eGFR, malignancies, diabetes mellitus, ischemic heart disease, valvular heart disease, peripheral vascular disease, sex, and albumin.
2. Predict the risk of mortality for each individual.
3. Fit a model with all-cause mortality as the outcome and treatment, predicted mortality risk, and their interaction as predictors.
4. For every level of predicted mortality risk, estimate the risk difference, difference in RMST, and hazard ratio.

### 07 - Example Supplemental Methods
Calculation of the individualized treatment effect for the example patient in Supplemental Materials 1.

### 08 - SA positivity
Sensitivity analysis comparing average treatment effect estimates when applying overlap weighting, Crump trimming, Stürmer trimming, or Walker trimming.

### 09 - Competing risk time-to-dialysis
Some patients choosing dialysis do not immediately initiate dialysis. To estimate how much time patients spend in the dialysis state in the two years following their decision, a multi-state illness-death model is fit for those choosing dialysis, with the states "Alive without dialysis", "Dialysis", and "Death".

### 10 - SA unmeasured confounding
Sensitivity analysis for unmeasured confounding. Shows which combinations of the confounder–outcome risk ratio (RR<sub>CD</sub>) and the prevalence of the unmeasured confounder in each treatment group (P<sub>C1</sub>, P<sub>C0</sub>) would be needed to move the observed risk ratio (0.57) to the null, presented as a 3D surface and two 2D cross-sections.

### 11 - Supplemental Figure HTE
Creates a figure illustrating the mathematical relationship between absolute and relative heterogeneity of treatment effect.

### 12 - SA home time 
Sensitivity analysis estimating how the 2-year period after the treatment decision is spent, using a multi-state model with the states **At home**, **Hospitalization**, **Hemodialysis (HD)**, **Peritoneal dialysis (PD)** and **Death**.

1. Build a long (patient-day) dataset of hospitalizations and kidney replacement therapy, and count the number of patients per transition.
2. Fit cause-specific Cox models for each transition, stratified by treatment and number of prior hospitalizations (capped); HD ↔ PD transitions are fit on the dialysis arm only.
3. Simulate each patient's trajectory after censoring up to 2 years, using multiple imputation (M = 10).
4. Estimate the mean number of days in each state, adjusted for confounding (IPTW) and censoring, and derive **in-center days** (all hospital days, HD sessions at 3 per week, PD training of 5 days plus one outpatient visit every 30 days), overall and after the start of dialysis.
5. Obtain percentile confidence intervals by bootstrapping the whole procedure (weights, Cox models, imputations). Set `recompute_bootstrap <- TRUE` to rerun the bootstrap; otherwise previously saved results are loaded from `Data/bootstrap_results_home_time.Rdata`.
6. Plot cumulative days and stacked state occupancy over time per treatment arm.

Outputs: `Table_S_n_per_transition.xlsx`, `Table_S_home_time_estimates.xlsx` (bootstrap CI, pooled estimate, and one sheet per imputation), and `Figure_S_home_time.png`.

### 13 - Subgroup analysis 
Estimates treatment effects (unweighted and IPTW) within subgroups, reported in the same format as the main results table (sample size, number of events, risk, risk difference, risk ratio, RMST, difference in RMST, hazard ratio):
1. **Age**: <80 versus ≥80 years
2. **Dialysis type**: HD versus PD (each compared against all patients choosing conservative management)

Outputs: `Table_S_HTE_subgroup_age.xlsx` and `Table_S_HTE_subgroup_dialysis_type.xlsx`.

---

## Helper functions (`Code/utils/`)

| File | Contents |
|---|---|
| `data_manipulation.R` | Factor encoding and dummies, age calculation, look-back of past information, CKD-EPI 2021 eGFR, eligibility intervals, ICD code lists and diagnosis dictionary, number formatting |
| `outcome_derivation.R` | Derivation of (multiple) event outcomes and time-to-event within a given window |
| `weighting.R` | Creation of IPTW, overlap, selection and SMR weights; propensity score trimming (Crump, Stürmer, Walker); weight summaries |
| `compute_absolute_relative_risks.R` | Weighted risks, risk differences, risk ratios, RMST and hazard ratios with bootstrap CIs; risk-based HTE |
| `competing_risk.R` | Multi-state modelling: building the long dataset, cause-specific Cox models, trajectory simulation, multiple imputation, bootstrap, state occupancy and in-center days |
| `tables.R` | Baseline tables, (TA)SMDs, results tables, risk model tables, home-time estimate tables |
| `plots.R` | PS distributions, love plots, Kaplan-Meier curves, forest plots, HTE effect plots, calibration plots, state occupancy panels |

## Outputs

All results are written to `Results/Main/` (Figures 1–3, Table 2, Table S1, statistics reported in the text) and `Results/Supplemental/` (supplemental tables `Table_S_*.xlsx`, supplemental figures `Figure_S_*.png`, and figures for the supplemental methods `Figure_M_*.png`).
