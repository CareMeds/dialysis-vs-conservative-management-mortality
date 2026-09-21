################################################################################
### Decision for dialysis versus conservative management
### PART 1 - Apply eligibility criteria
################################################################################

# set-up
rm(list = ls(all.names = TRUE))
set.seed(1)
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/"
)

# load libraries
library(data.table)

# load functions
source("Code/utils/data_manipulation.R")

# load data
load("Data/merged_ckd.Rdata")
load("Data/cleaned/snr_inpatient.Rdata")
load("Data/cleaned/snr_outpatient.Rdata")

# create new diagnoses dictionary
new_diag <- TRUE
if (!new_diag) {
  load("Data/Davies_65_80.Rdata")
  load("Data/hiv_dementia.Rdata")
}

################################################################################
### Add row for xth birthday for those who had visits before and after that birthdate
################################################################################
merged_ckd_bday65 <- add_bday_rows(
  dt = merged_ckd,
  bday_year = 65,
  birthdate_name = "birthdate",
  id_name = "LOPNR",
  date_name = "visit_date",
  death_name = "DODSDAT",
  study_end = end_date
)
merged_ckd_bday80 <- add_bday_rows(
  dt = merged_ckd_bday65,
  bday_year = 80,
  birthdate_name = "birthdate",
  id_name = "LOPNR",
  date_name = "visit_date",
  death_name = "DODSDAT",
  study_end = end_date
)

# fill bday65 and sex for new rows that are added for bday80
for (var in c("bday65", "female", "decision_date1", "decision_modality1")) {
  merged_ckd_ext <- retrieve_past_info(
    dt = merged_ckd_bday80,
    dictionary = merged_ckd_bday65,
    id_name = "LOPNR",
    date_name = "visit_date",
    var_name = var,
    lookback_months = Inf
  )
}

################################################################################
### Eligibility criterion 1: include those visits with eGFR < 20
################################################################################
# obtain eGFR information for missing eGFR up to a year ago
filled_eGFR <- retrieve_past_info(
  dt = merged_ckd_ext,
  dictionary = merged_ckd,
  id_name = "LOPNR",
  date_name = "visit_date",
  var_name = "egfr2021",
  lookback_months = 12
)

# remove rows with missing eGFR, retain relevant columns
egfr_dt <- filled_eGFR[!is.na(egfr2021), c("LOPNR", "visit_date", "birthdate", "DODSDAT", "egfr2021")]

# compute next visit date per patient
# shift() can find lag or lead (forward) values fast within each group
egfr_dt[, next_visit := shift(visit_date, type = "lead"), by = LOPNR]

# filter rows with low eGFR (< 20)
low_egfr <- egfr_dt[egfr2021 < 20, c("LOPNR", "visit_date", "next_visit", "birthdate", "DODSDAT")]

# determine end date of eligibility for each low eGFR episode
low_egfr[, elig_from_1 := visit_date]
low_egfr[, elig_until_1 := eligibility_end(
  visit_date = visit_date,
  next_visit = next_visit,
  date_of_death = DODSDAT,
  study_end = end_date
)]

################################################################################
### Eligibility criterion 2: age >= 65 and Davies >=2 OR age >= 80
################################################################################
# get patient information from merged_ckd
age_dt <- merged_ckd_bday80[low_egfr[, -c("DODSDAT")], on = c("LOPNR", "visit_date")]

################################################################################
### Update eligibility for patients age >= 80 ##################################
################################################################################
age_above_80 <- age_dt[age >= 80]

# patients are eligible from first visit_date at or after 80th birthday, until they die
age_above_80[, elig_until_2 := eligibility_end(
  visit_date = visit_date,
  next_visit = next_visit,
  date_of_death = DODSDAT,
  study_end = end_date
)]
elig_above_80 <- age_above_80[, c(
  "LOPNR",
  "visit_date",
  "birthdate",
  "DODSDAT",
  "female",
  "elig_until_2",
  "egfr2021",
  "age",
  "decision_date1",
  "decision_modality1"
)]

################################################################################
### Update eligibility for patients age < 80 ###################################
################################################################################
# 1. Define cohort for comorbidity lookup
# 1.1 Identify visits eligible for comorbidity lookup:
# - Include only individuals who survived to at least age 65
#   (DODSDAT is missing or occurs after bday65).
# - Keep only visits where patient age is < 80.
age_65_80 <- age_dt[(is.na(DODSDAT) | DODSDAT > bday65) & age < 80]
setkey(age_65_80, LOPNR, visit_date)
IDs_65_80 <- unique(age_65_80[, LOPNR])

# 1.2 Define maximum lookup date per individual:
# - Restrict to patients in age_65_80
# - Use 80th birthday (bday80) as upper bound for to obtain diagnoses
# - Set max_date to 80th birthday or end of follow-up
max_date_dict <- unique(merged_ckd_bday80[(LOPNR %in% IDs_65_80) &
                                            (is.na(DODSDAT) | DODSDAT > bday80), .(LOPNR, max_date = pmin(end_date, bday80))])

# 1.3 Create dictionary of relevant comorbidity diagnoses (Davies score):
# Use diagnoses.dictionary() to extract inpatient and outpatient records
# for individuals in IDs_65_80, limited by max_date_dict.
Davies_comorbidities <- c("cancer",
                          "ihd",
                          "pvd",
                          "hf",
                          "dm",
                          "scvd",
                          "copd",
                          "cirr",
                          "psycho")
if (new_diag) {
  Davies_65_80 <- diagnoses.dictionary(
    inpatient_dt = UT_R_PAR_SV_123160_2023,
    outpatient_dt = UT_R_PAR_OV_123160_2023,
    lopnr_obtain_diag = IDs_65_80,
    comorbidities = Davies_comorbidities,
    max_date_dict = max_date_dict
  )$diagnoses_dt
  save(Davies_65_80, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Data/Davies_65_80.Rdata")
}

# 1.4 Add comorbidities to dt
age_65_80_comorb <- merge(age_65_80,
                          Davies_65_80,
                          by = c("LOPNR", "visit_date"),
                          all.x = TRUE)

# 2. Add eligible comorbidity visit dates within the time window to dataset
# 2.1 Filter diagnoses to study period (2007-01-01 through 2023-12-31)
start_date <- as.IDate("2007-01-01")
diag_filtered <- Davies_65_80[visit_date >= start_date &
                                visit_date <= end_date]

# 2.2 Construct set of new eligible visit dates after initial eligibility
# Start from all diagnoses in diag_filtered for patients in IDs_65_80.
# Exclude diagnoses that occur exactly on an existing visit_date in age_65_80
# to avoid re-using the same visits.
new_visit_dates <- diag_filtered[LOPNR %in% IDs_65_80]

# 2.3 Keep only visits within eligible range
eligible_visits <- new_visit_dates[age_65_80, on = .(LOPNR, visit_date >= elig_from_1, visit_date <= elig_until_1), nomatch = 0, .(
  LOPNR,
  # from new_visit_dates
  visit_date = x.visit_date,
  decision_date1,
  decision_modality1,
  birthdate,
  bday80,
  DODSDAT,
  female,
  cancer,
  ihd,
  pvd,
  hf,
  dm,
  scvd,
  copd,
  cirr,
  psycho
)]

# ensure to not add rows that are already in age_65_80
new_visit_dates_range <- eligible_visits[!age_65_80, on = c("LOPNR", "visit_date")]

# 2.4 Ensure eligibility of new visit dates:
# - Patient age must be between 65 and 79 (inclusive of 65, exclusive of 80).
# - Visit_date must fall before the recorded death date (DODSDAT).
new_visit_dates_range <- calculate_age(dt = new_visit_dates_range,
                                       birthdate_col = "birthdate",
                                       visitdate_col = "visit_date")
new_rows <- new_visit_dates_range[age >= 65 & age < 80 &
                                    (visit_date < DODSDAT |
                                       is.na(DODSDAT))]
setkey(new_rows, LOPNR, visit_date)

# 2.5 only add new rows that have eGFR value < 20 in past year
new_rows[, egfr2021 := NA_real_]
new_rows_eGFR <- retrieve_past_info(
  dt = new_rows,
  dictionary = merged_ckd,
  id_name = "LOPNR",
  date_name = "visit_date",
  var_name = "egfr2021",
  lookback_months = 12
)
new_rows_eGFR <- new_rows_eGFR[!is.na(egfr2021) & egfr2021 < 20]

# 2.6 merge rows from visit dates and new rows where comorbidities might change
age_65_80_full <- rbindlist(list(age_65_80_comorb[, c(
  "LOPNR",
  "visit_date",
  "decision_date1",
  "decision_modality1",
  "birthdate",
  "bday80",
  "DODSDAT",
  "female",
  ..Davies_comorbidities,
  "age",
  "egfr2021"
)], new_rows_eGFR))
setkey(age_65_80_full, LOPNR, visit_date)

# 2.7 retrieve past info for missing comorbidities
for (comorbidity in Davies_comorbidities) {
  age_65_80_comorbidities <- retrieve_past_info(
    dt = age_65_80_full,
    dictionary = age_65_80_full,
    id_name = "LOPNR",
    date_name = "visit_date",
    var_name = comorbidity,
    lookback_months = ifelse(comorbidity == "cancer", 3 * 12, Inf),
    fill_with_zero = TRUE # set remaining NA to zero
  )
}

# Calculate Davies score
age_65_80_comorbidities[, Davies_score := rowSums(.SD, na.rm = TRUE), .SDcols = Davies_comorbidities]

# Include those visits for which patients are 65 <= age <80 and Davies score >= 2
elig_65_80_Davies <- age_65_80_comorbidities[age >= 65 &
                                               Davies_score >= 2]
setkey(elig_65_80_Davies, LOPNR, visit_date)

# 5. Define eligibility time frames
# update next_visit since rows from snr_inpatient and snr_outpatient were added
elig_65_80_Davies[, next_visit := shift(visit_date, type = "lead"), by = "LOPNR"]
elig_65_80_Davies[, elig_until_2 := eligibility_end(
  visit_date = visit_date,
  next_visit = next_visit,
  date_of_death = DODSDAT,
  study_end = end_date
)]

# ensure that elig_until_2 does not run after first age>80 visit
first_80_visits <- elig_above_80[, .(first_80_visit = min(visit_date)), by = "LOPNR"]
elig_65_80_Davies <- first_80_visits[elig_65_80_Davies, on = "LOPNR"]
elig_65_80_Davies[, elig_until_2 := pmin(
  elig_until_2,
  # original eligibility
  first_80_visit - 1,
  # cannot exceed first 80+ visit
  na.rm = TRUE
)]

# 6. combine data where age >= 65 and Davies >=2 OR age >= 80
age_crit_dt <- rbindlist(list(elig_65_80_Davies[, c(
  "LOPNR",
  "visit_date",
  "birthdate",
  "DODSDAT",
  "female",
  "elig_until_2",
  "egfr2021",
  "age",
  "decision_date1",
  "decision_modality1",
  "Davies_score",
  ..Davies_comorbidities
)], elig_above_80), fill = TRUE) # elig_above_80 does not have comorbidities yet
setkey(age_crit_dt, LOPNR, visit_date)

################################################################################
### Eligibility criterion 3: all lab measurements (max 1 year look-back)
################################################################################
# create laboratory measurements dictionary
lab_vars <- c("albumin", "calcium_total", "dbp", "hb", "phosphate", "sbp")
lab_dictionary <- unique(merged_ckd[, .SD, .SDcols = c("LOPNR", "visit_date", lab_vars)])

# get lab measurements that are available at each visit_date
lab_dt <- lab_dictionary[age_crit_dt, on = c("LOPNR", "visit_date")]

# retrieve labs up to one year ago for those that are NA
for (lab_var in lab_vars) {
  lab_dt <- retrieve_past_info(
    dt = lab_dt,
    dictionary = lab_dictionary,
    id_name = "LOPNR",
    date_name = "visit_date",
    var_name = lab_var,
    lookback_months = 12
  )
}

# filter rows where none of the lab_vars are NA
lab_complete_dt <- lab_dt[complete.cases(lab_dt[, lab_vars, with = FALSE])]

################################################################################
### Eligibility criterion 4: no history of kidney transplantation or dialysis
################################################################################
# Extract the date at which KRT is started or the first decision for transplantation or dialysis is made
first_dia_trans <- merged_ckd[!is.na(krt_startdate) |
                                !is.na(decision_tx_date1) |
                                !is.na(decision_date1),
                              .(date_dia_trans = min(krt_startdate, decision_tx_date1,
                                                     decision_date1 + 1, na.rm = TRUE)),
                              by = "LOPNR"]

# Append first date of dialysis or transplantation
trans_dia_dt <- merge(lab_complete_dt,
                      first_dia_trans,
                      by = "LOPNR",
                      all.x = TRUE)

# Remove visit_dates before which transplantation or dialysis was registered
no_trans_dia_dt <- trans_dia_dt[(is.na(date_dia_trans) |
                                   visit_date < date_dia_trans)]

# Update eligibility interval for start or decision of dialysis or transplantation
# Case 1: no history of start or kidney replacing therapy, if so eligible until day before
# Case 2: no history of decision of transplantation, if so eligible until day before
# Case 3: no history of decision of dialysis, if so eligible until decision date
no_trans_dia_dt[, elig_until_4 :=
                  fifelse(
                    !is.na(date_dia_trans) &
                      date_dia_trans > visit_date &
                      date_dia_trans <= elig_until_2,
                    date_dia_trans - 1,
                    elig_until_2
                  )]

################################################################################
### Eligibility criterion 5: no history of hiv or dementia
################################################################################
# 1. obtain dictionary of diagnoses
if (new_diag) {
  hiv_dementia <- diagnoses.dictionary(
    inpatient_dt = UT_R_PAR_SV_123160_2023,
    outpatient_dt = UT_R_PAR_OV_123160_2023,
    lopnr_obtain_diag = unique(no_trans_dia_dt$LOPNR),
    comorbidities = c("hiv", "dementia"),
    max_date_dict = data.table(
      LOPNR = unique(no_trans_dia_dt$LOPNR),
      max_date = end_date
    )
  )$diagnoses_dt
  save(hiv_dementia, 
       file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Data/hiv_dementia.Rdata")
}

# 2. create dt that finds first hiv or dementia event
hiv_dem <- hiv_dementia[LOPNR %in% no_trans_dia_dt$LOPNR &
                          (hiv == 1 | dementia == 1)]
first_hiv_dem <- hiv_dem[, .(date_hiv_dementia = min(visit_date)), by = "LOPNR"]

# 3. append first date to data
hiv_dementia_dt <- merge(no_trans_dia_dt,
                         first_hiv_dem,
                         by = "LOPNR",
                         all.x = TRUE)

# 3. remove visit_dates before which hiv or dementia was registered
no_hiv_dementia_dt <- hiv_dementia_dt[is.na(date_hiv_dementia) |
                                        visit_date < date_hiv_dementia]

# 3.1 Truncate eligibility interval at HIV or dementia diagnosis
# For each row in no_hiv_dementia_dt:
#   - If date_hiv_dementia exists and falls within the current interval
#       (visit_date ≤ date_hiv_dementia ≤ elig_until_4),
#       set elig_until_5 = date_hiv_dementia - 1
#   - Otherwise, keep elig_until_5 = elig_until_4
no_hiv_dementia_dt[, elig_until_5 :=
                     fifelse(
                       !is.na(date_hiv_dementia) &
                         date_hiv_dementia > visit_date &
                         date_hiv_dementia <= elig_until_4,
                       date_hiv_dementia - 1,
                       elig_until_4
                     )]

################################################################################
### Fill comorbidities for those above 80
################################################################################
# Retrieve Davies comorbidities for these patients from inpatient and outpatient records
Davies_80 <- diagnoses.dictionary(
  inpatient_dt = UT_R_PAR_SV_123160_2023,
  outpatient_dt = UT_R_PAR_OV_123160_2023,
  lopnr_obtain_diag = unique(no_hiv_dementia_dt[is.na(Davies_score), LOPNR]),
  comorbidities = Davies_comorbidities,
  max_date_dict = data.table(
    LOPNR = unique(no_hiv_dementia_dt$LOPNR),
    max_date = end_date
  )
)$diagnoses_dt

# fill comorbidities for those aged above 80
for (comorbidity in Davies_comorbidities) {
  age_80_comorbidities <- retrieve_past_info(
    dt = no_hiv_dementia_dt,
    dictionary = Davies_80,
    id_name = "LOPNR",
    date_name = "visit_date",
    var_name = comorbidity,
    lookback_months = ifelse(comorbidity == "cancer", 3 * 12, Inf),
    fill_with_zero = TRUE # set remaining NA to zero
  )
}

# Calculate Davies score
age_80_comorbidities[, Davies_score := rowSums(.SD, na.rm = TRUE), .SDcols = Davies_comorbidities]

################################################################################
### Remove visits after palliative care or nursing home admission
################################################################################
# extract information on palliative care
outpatient <- UT_R_PAR_OV_123160_2023
setDT(outpatient)
first_palliative <- outpatient[LOPNR %in% age_80_comorbidities$LOPNR &
                                 (MVO == "061" | MVO == "020" | 
                                    MVO == "243" | MVO == "246"),
                               c("LOPNR", "INDATUMA"),
                               with = FALSE][, `:=`
                                             (INDATUMA = as.IDate(as.character(INDATUMA), format = "%Y%m%d"))][order(INDATUMA), .SD[1], by = "LOPNR"]
setnames(first_palliative, "INDATUMA", "date_palliative_nursing")

# append first date to data
palliative_dt <- merge(age_80_comorbidities,
                       first_palliative,
                       by = "LOPNR",
                       all.x = TRUE)

# remove visit_dates before which palliative care was registered
no_palliative_dt <- palliative_dt[is.na(date_palliative_nursing) |
                                    visit_date < date_palliative_nursing]

# update eligibility window
no_palliative_dt[, elig_until_6 :=
                   fifelse(
                     !is.na(date_palliative_nursing) &
                       date_palliative_nursing > visit_date &
                       date_palliative_nursing <= elig_until_5,
                     date_palliative_nursing - 1,
                     elig_until_5
                   )]

################################################################################
### Eligibility cohort: expand cohort
################################################################################
# expand data frame by adding a row for each day between visit_date and elig_until
elig_dt_expanded <- no_palliative_dt[, .(date = seq(visit_date, elig_until_6, by = "day")), by = c(
  "LOPNR",
  "visit_date",
  "elig_until_6",
  "birthdate",
  "DODSDAT",
  "female",
  "egfr2021",
  lab_vars,
  "Davies_score",
  Davies_comorbidities,
  "decision_modality1",
  "decision_date1"
)]
elig_dt_expanded[, c("visit_date", "elig_until_6") := NULL]
setnames(elig_dt_expanded, "date", "visit_date")

################################################################################
### Analysis cohort: those who have registered treatment decision
################################################################################
# only keep eligible date at which the initial decision is made
cohort <- elig_dt_expanded[visit_date == decision_date1 &
                             !is.na(decision_modality1)]
nrow(cohort)

# remove transplantations
cohort[, trt := factor(fifelse(decision_modality1 == "Konservativ behandling", 0L, 1L))]
table(cohort$trt, useNA = "ifany")

# define age
cohort <- calculate_age(dt = cohort,
                        birthdate_col = "birthdate",
                        visitdate_col = "visit_date")

################################################################################
### Eligibility cohort
################################################################################
# randomly select a date among eligible dates for patients without a decision
elig_cohort_no_dec <- elig_dt_expanded[!(LOPNR %in% cohort$LOPNR), .SD[sample(.N, 1)], by = "LOPNR"]

# define age
elig_cohort_no_dec <- calculate_age(dt = elig_cohort_no_dec,
                                    birthdate_col = "birthdate",
                                    visitdate_col = "visit_date")

# define trt
elig_cohort_no_dec[, trt := NA]

# combine cohorts
elig_cohort <- rbind(cohort, elig_cohort_no_dec)

# Define S = 1 if in analysis data set and 0 otherwise
elig_cohort <- elig_cohort[, S := ifelse(is.na(trt), 0, 1)]

################################################################################
### Flow chart
################################################################################
flow_chart <- data.frame(
  names = c("Initial number of patients",
            "Excluded patients with eGFR < 20",
            "Excluded patients with age>=65 & Davies>=2 OR age >= 80",
            "Excluded patients with all lab measurements",
            "Excluded patients without transplantation or dialysis",
            "Excluded patients with history of hiv or dementia",
            "Excluded patients without palliative care",
            "Eligible patients",
            "Excluded patients with treatment decision dialysis or CC",
            "Final cohort",
            "Patients who chose conservative management",
            "Patients who chose dialysis"),
  counts = c(  merged_ckd[, uniqueN(LOPNR)],
               merged_ckd[, uniqueN(LOPNR)] - low_egfr[, uniqueN(LOPNR)],
               low_egfr[, uniqueN(LOPNR)] - age_crit_dt[, uniqueN(LOPNR)],
               age_crit_dt[, uniqueN(LOPNR)] - lab_complete_dt[, uniqueN(LOPNR)],
               lab_complete_dt[, uniqueN(LOPNR)] - no_trans_dia_dt[, uniqueN(LOPNR)],
               no_trans_dia_dt[, uniqueN(LOPNR)] - no_hiv_dementia_dt[, uniqueN(LOPNR)],
               no_hiv_dementia_dt[, uniqueN(LOPNR)] - no_palliative_dt[, uniqueN(LOPNR)],
               no_palliative_dt[, uniqueN(LOPNR)],
               no_palliative_dt[, uniqueN(LOPNR)] - cohort[, uniqueN(LOPNR)],
               cohort[, .N],
               cohort[trt == 0, .N],
               cohort[trt == 1, .N]))
openxlsx::write.xlsx(
  x = flow_chart,
  file = "Data/Figures/Flow_chart.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)
  
# save file
save(
  elig_cohort,
  cohort,
  merged_ckd,
  end_date,
  egfr_dt,
  low_egfr,
  age_dt,
  age_65_80_comorb,
  age_crit_dt,
  lab_complete_dt,
  lab_dt,
  trans_dia_dt,
  no_trans_dia_dt,
  palliative_dt,
  no_palliative_dt,
  hiv_dementia_dt,
  no_hiv_dementia_dt,
  elig_dt_expanded,
  file = "Data/new_cohort.Rdata"
)

# check
sort(colSums(is.na(elig_cohort)))
sort(colSums(is.na(cohort)))
