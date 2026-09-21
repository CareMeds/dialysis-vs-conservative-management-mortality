################################################################################
### Decision for dialysis versus conservative management
### PART 2 - Coavariate and outcome derivation
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
pacman::p_load("dplyr", "tidyr", "readr", "lubridate", "stringr")
source("Code/utils/data_manipulation.R")
source("Code/utils/outcome_derivation.R")

# load data
load("Data/new_cohort.Rdata")
load("Data/merged_ckd.Rdata")
load("Data/cleaned/snr_inpatient.Rdata")
load("Data/cleaned/snr_outpatient.Rdata")
load("Data/cleaned/snr_rrt_long.Rdata")
load("Data/cleaned/snr_lmed.Rdata")
load("Data/cleaned/snr_death.Rdata")

# convert inpatient
inpatient <- UT_R_PAR_SV_123160_2023
setDT(inpatient)
inpatient[, INDATUMA := as.IDate(as.character(INDATUMA), format = "%Y%m%d")]
inpatient[, UTDATUMA := as.IDate(as.character(UTDATUMA), format = "%Y%m%d")]

# convert inpatient
outpatient <- UT_R_PAR_OV_123160_2023
setDT(outpatient)

# convert lmed
lmed <- UT_R_LMED_123160_2023
setDT(lmed)
lmed[, EDATUM := as.IDate(EDATUM, format = "%Y%m%d")]

# ID name
id_name <- "LOPNR"

################################################################################
### Add covariates to eligible cohort without treatment decision and analysis cohort
################################################################################
for (cohort_name in c("cohort", "elig_cohort")) {
  cat("Working on", cohort_name, "\n")
  # evaluate for correct cohort
  working_cohort <- eval(parse(text = cohort_name))
  
  ################################################################################
  ### Type of dialysis (HD or PD)
  ################################################################################
  working_cohort[, dialysis_type := fifelse(
    decision_modality1 == "Konservativ behandling", "Conservative management",
    fifelse(
      decision_modality1 %in% c("HD", "Sj\xe4lv-HD", "Hem-HD"), "HD",
      fifelse(
        decision_modality1 %in% c("PD", "Assisterad PD"), "PD",
        NA_character_
      )
    )
  )]
  working_cohort[, dialysis_type := factor(
    dialysis_type,
    levels = c("PD", "HD", "Conservative management")
  )]
  
  ################################################################################
  ### Add comorbidities acs, hyperten, vhd, cevd, af, arrh, lung, thrombo, liver, fracture, aki
  ################################################################################
  other_comorb <- c(
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
    "aki"
  )
  
  # create dictionary of comorbidities from snr_inpatient and snr_outpatient
  in_out_dict <- diagnoses.dictionary(
    inpatient_dt = UT_R_PAR_SV_123160_2023,
    outpatient_dt = UT_R_PAR_OV_123160_2023,
    lopnr_obtain_diag = unique(working_cohort$LOPNR),
    comorbidities = other_comorb,
    max_date_dict = data.table(LOPNR = unique(working_cohort$LOPNR), max_date = end_date)
  )
  if (cohort_name == "cohort") {
    save(in_out_dict, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Data/other_comorb.Rdata")
  }
  
  # merge with cohort
  cohort_other_comorb <- merge(
    working_cohort,
    in_out_dict$diagnoses_dt,
    by = c(id_name, "visit_date"),
    all.x = TRUE
  )
  
  # retrieve past info infinitely after merging with cohort
  for (comorbidity in other_comorb) {
    cohort_other_comorb <- retrieve_past_info(
      dt = cohort_other_comorb,
      dictionary = in_out_dict$diagnoses_dt,
      id_name = id_name,
      date_name = "visit_date",
      var_name = comorbidity,
      lookback_months = Inf,
      fill_with_zero = TRUE # set remaining NA to zero
    )
  }
  
  ################################################################################
  ### Hospitalizations in past year
  ################################################################################
  # obtain inpatient info, only relevant patients
  dia_dt <- inpatient[LOPNR %in% working_cohort$LOPNR, .SD[1], by = c(id_name, "INDATUMA")][, c(id_name, "INDATUMA", "HDIA"), with = FALSE]
  
  # calculate the number of hospitalizations in past year (any + cardiovascular)
  # define start_date from which to look at hospitalizations by taking into account leap years properly
  hospital <- dia_dt[working_cohort[, .(
    id = get(id_name),
    visit_date,
    start_date = lubridate::add_with_rollback(
      visit_date,
      lubridate::period(-12, units = "months"),
      roll_to_first = FALSE
    )
  )], on = .(LOPNR = id, INDATUMA >= start_date, INDATUMA <= visit_date), nomatch = NULL][, .(n_hospital = .N,
                                                                                              n_cvd_hospital = sum(startsWith(HDIA, "I"), na.rm = TRUE)), by = .(LOPNR)]
  
  # merge with cohort
  cohort_hosp <- merge(cohort_other_comorb,
                       hospital,
                       by = id_name,
                       all.x = TRUE)
  
  # replace NA by zero, i.e., if there was no hospitalization in the past year
  cohort_hosp[, `:=` (
    n_hospital = fifelse(is.na(n_hospital), 0L, n_hospital),
    n_cvd_hospital = fifelse(is.na(n_cvd_hospital), 0L, n_cvd_hospital)
  )]
  
  ################################################################################
  ### Medications and iron based on ATC code
  ################################################################################
  lmed_dict <- lmed[LOPNR %in% working_cohort$LOPNR, c(id_name, "EDATUM", "ATC"), with = FALSE][, unique(.SD)]
  
  # create medications data frames
  med_patterns <- list(
    bblock       = "^C07",
    calblock     = "^C08",
    diuretic     = "^C03",
    rasi         = "^C09A|^C09B|^C09C|^C09D",
    lipid        = "^C10",
    esa          = "^B03XA",
    phosbinder   = "^V03AE02|^V03AE03|^V03AE04|^V03AE05|^V03AE06|^V03AE07|^V03AE08",
    vitamind     = "^A11CC",
    digoxin      = "^C01AA05",
    vasodilator  = "^C01D",
    antiplatelet = "^B01AC",
    anticoag     = "^B01AA|^B01AE07|^B01AF|^B01AX05",
    # potasbinder  = "^V01AE01|^V03AE09|^V03AE10",
    iron_po      = "^B03AA|^B03AB|^B03AD|^B03AE",
    iron_iv      = "^B03AC"
  )
  
  # extract medications
  lmed_dt <- create_dummies(
    dt = lmed_dict,
    var_names = names(med_patterns),
    patterns = med_patterns,
    id_name = id_name,
    col_name = "ATC",
    date_name = "EDATUM"
  )
  if (cohort_name == "cohort") {
    save(lmed_dt, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Data/lmed_dt.Rdata")
  }
  
  # append to cohort
  cohort_med <- merge(cohort_hosp,
                      lmed_dt,
                      by = c(id_name, "visit_date"),
                      all.x = TRUE)
  
  # retrieve past info infinitely after merging with cohort
  for (medication in names(med_patterns)) {
    cohort_med <- retrieve_past_info(
      dt = cohort_med,
      dictionary = lmed_dt,
      id_name = id_name,
      date_name = "visit_date",
      var_name = medication,
      lookback_months = 6,
      fill_with_zero = TRUE # set remaining NA to zero
    )
  }
  
  ################################################################################
  ### Combine esa and iron from CKD and medications data
  ################################################################################
  # extract iron
  esa_iron_ckd_dt <- merged_ckd[LOPNR %in% working_cohort$LOPNR &
                                  (!is.na(iron_med) |
                                     !is.na(esa)), c(id_name, "visit_date", "esa", "iron_med", "iron_type", "crp"), with = FALSE][, unique(.SD)]
  
  # iron from medications dt
  esa_iron_med_dt <- cohort_med[!is.na(esa) |
                                  !is.na(iron_iv) |
                                  !is.na(iron_po), c(id_name, "visit_date", "esa", "iron_iv", "iron_po"), with = FALSE][, unique(.SD)]
  
  # combine iron_dt from CKD and medications data
  esa_iron_dt <- merge(
    esa_iron_ckd_dt,
    esa_iron_med_dt,
    by = c(id_name, "visit_date"),
    all.x = TRUE
  )
  esa_iron_dt[, `:=`
              (
                esa = fifelse((!is.na(esa.x) & esa.x == 1) |
                                (!is.na(esa.y) & esa.y == 1), 1, 0),
                iron_iv = fifelse((!is.na(iron_type) &
                                     iron_type == "i.v.") |
                                    (!is.na(iron_iv) &
                                       iron_iv == 1), 1, 0),
                iron_po = fifelse((!is.na(iron_type) &
                                     iron_type == "p.o.") |
                                    (!is.na(iron_po) &
                                       iron_po == 1), 1, 0)
              )]
  esa_iron_dt <- esa_iron_dt[, c(id_name, "visit_date", "esa", "iron_iv", "iron_po"), with = FALSE]
  if (cohort_name == "cohort") {
    save(esa_iron_dt, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Data/esa_iron_dt.Rdata")
  }
  
  # add esa and iron to cohort by one year look back for iron
  cohort_esa_iron <- copy(cohort_med)
  for (var_name in c("esa", "iron_iv", "iron_po", "crp")) {
    cohort_esa_iron[, (var_name) := NA_real_]
    if (var_name == "crp") {
      dict <- esa_iron_ckd_dt
    } else{
      dict <- esa_iron_dt
    }
    cohort_esa_iron <- retrieve_past_info(
      dt = cohort_esa_iron,
      dictionary = dict,
      id_name = id_name,
      date_name = "visit_date",
      var_name = var_name,
      lookback_months = 12,
      fill_with_zero = TRUE # set remaining NA to zero
    )
  }
  
  ################################################################################
  ### Primary kidney disease from snr_ckd data
  ################################################################################
  # obtain information from snr_ckd
  # 0 = no primary kidney disease
  # 1 = Diabetesnefropati
  # 2 = Hyperoni
  # 3 = Other, i.e., Adult polycystisk njursjukdom, Glomerulonefrit, Pyelonefrit, Renovaskular, Uremi UNS
  prd_dt <- merged_ckd[LOPNR %in% working_cohort$LOPNR &
                         !is.na(prd_cat), c(id_name, "visit_date", "prd_cat"), with =
                         FALSE][, unique(.SD)]
  prd_dt[, prd_cat := fifelse(prd_cat == "Diabetesnefropati",
                              1,
                              fifelse(prd_cat == "Hypertoni", 2, 3))]
  
  # add prd to cohort_iron, looking infinitely back
  cohort_prd <- prd_dt[cohort_esa_iron, on = c(id_name, "visit_date"), roll = TRUE]
  
  # define as factor
  cohort_prd[, prd_cat := factor(
    prd_cat,
    levels = c(1L, 2L, 3L),
    labels = c("Diabetesnefropati", "Hypertoni", "Other")
  )]
  
  ################################################################################
  ### Education category
  ################################################################################
  # extract education, only keep earliest education date for each patient
  edu_dt <- merged_ckd[LOPNR %in% working_cohort$LOPNR &
                         !is.na(info_date1) &
                         !is.na(info_type1), c(id_name, "info_date1", "info_type1"), with = FALSE][, `:=`
                                                                                                   (visit_date = as.IDate(info_date1, format = "%m/%d/%Y"),
                                                                                                     edu = 1)][order(visit_date), .SD[1], by = id_name]
  
  # merge with cohort
  cohort_prd[, edu := NA_real_]
  cohort_edu <- retrieve_past_info(
    dt = cohort_prd,
    dictionary = edu_dt,
    id_name = id_name,
    date_name = "visit_date",
    var_name = "edu",
    lookback_months = Inf,
    fill_with_zero = TRUE # set remaining NA to zero
  )
  
  ################################################################################
  ### Clinic level
  ################################################################################
  geo_dt <- merged_ckd[LOPNR %in% working_cohort$LOPNR &
                         (!is.na(clinic) |
                            !is.na(county)), # at least one is non-missing
                       c(id_name, "visit_date", "clinic", "county"), with = FALSE][, unique(.SD)]
  
  # Define mapping as a named vector
  county_to_region <- c(
    "Stockholm"       = "Stockholm",
    "Skane"           = "Sodra",
    "Blekinge"        = "Sodra",
    "Kronoberg"       = "Sodra",
    "Halland"         = "Sodra",
    "Orebro"          = "Orebro.Uppsala",
    "Uppsala"         = "Orebro.Uppsala",
    "Sodermanland"    = "Orebro.Uppsala",
    "Vastmanland"     = "Orebro.Uppsala",
    "Vastra Gotaland" = "Vastra",
    "Varmland"        = "Vastra",
    "Jonkoping"       = "Vastra",
    "Dalarna"         = "Vastra",
    "Gavleborg"       = "Vastra",
    "Kalmar"          = "Other regions",
    # "Sydostra",
    "Ostergotland"    = "Other regions",
    # "Sydostra",
    "Gotland"         = "Other regions",
    # "Sydostra",
    "Norrbotten"      = "Other regions",
    # "Norra",
    "Vasterbotten"    = "Other regions",
    # "Norra",
    "Vasternorrland"  = "Other regions",
    # "Norra",
    "Jamtland"        = "Other regions",
    # "Norra",
    "Ok\xe4nd"        = "Other regions",
    # meaning Unknown, these are all referred to other disciplines, merge with reference
    "Utrikes"         = "Other regions"  # meaning Emigrated, merge with reference
  )
  
  # Add region column based on mapping
  geo_dt[, region := county_to_region[county]]
  
  # create clinic levels
  # local clinics
  clinic_lev1 <- c(
    "Avesta",
    "Bollnas",
    "Eksjo",
    "Gallivare",
    "Hassleholm",
    "Karlshamn",
    "Karlskoga",
    "Koping",
    "Ljungby",
    # "Lyckesele",
    "Lycksele",
    # new
    "Mora",
    "Motala",
    "Nykoping",
    "Pitea",
    "Skelleftea",
    "Solleftea",
    "Varberg/Kungsbacka",
    "Varnamo",
    "Vastervik",
    "Ystad",
    "Angelholm",
    "Ornskoldsvik",
    "Visby",
    "Falkoping",
    "Gbg, Lundby",
    "Trelleborg",
    "Utrikes"     # emigration merged with local
  )
  
  # regional clinics
  clinic_lev2 <- c(
    "Boras",
    "Danderyd",
    "Eskilstuna",
    "Falun",
    # "Gävle",
    "Gavle",
    # new
    "Halmstad",
    "Helsingborg",
    "Jonkoping",
    "Kalmar",
    "Karlskrona",
    "Karlstad",
    "Kristianstad",
    "Norrkoping",
    "Skovde",
    "Sunderby",
    "Sundsvall",
    "Trollhattan, NAL",
    "Vasteras",
    "Vaxjo",
    "Ostersund",
    "Trollhattan",
    "Ej Njurmedicin"  # referred back to other discipline, merged with regional
  )
  
  # academic clinics
  clinic_lev3 <- c(
    # "Gbg SU/Ostra dialysmott",
    "Gbg SU/Ostra",
    # new
    "Gbg, SU/Njurmed",
    "Gbg, SU/Trpl",
    # new
    "Karolinska Njur med",
    "Linkoping",
    "Lund Njurmed",
    # "Malmo, njurmed",
    "Malmo, Heleneholms",
    # new
    "Molndal",
    "Uppsala, med",
    "Uppsala, Trpl",
    "Umea",
    "Orebro",
    "Gbg/Ostra",
    "Huddinge-K Njur med",
    # "Huddinge-K Njur med (Gammal)",
    "Huddinge-K, Trpl",
    # new
    "Malmo",
    "Solna-K Njur med",
    # "Solna-K Njur med (Gammal)",
    "Solna, diaverum",
    # new
    "Nacka",
    # dialysis unit, academic
    "Sodertalje",
    # dialysis unit, academic
    "J\xe4rf\xe4lla"       # dialysis unit, academic
  )
  
  # merge and assign clinic levels, infinite look-back
  cohort_geo <- geo_dt[cohort_edu, on = .(LOPNR, visit_date), roll = TRUE]
  
  # create factor for clinic level
  cohort_geo[, clinic_level := factor(
    fifelse(
      clinic %in% clinic_lev1,
      1,
      fifelse(
        clinic %in% clinic_lev2,
        2,
        fifelse(clinic %in% clinic_lev3, 3, NA)
      )
    ),
    levels = 1:3,
    labels = c("Local", "Regional", "Academic")
  )]
  
  ################################################################################
  ### Nursing home
  ################################################################################
  # keep only nursing home
  nursing_dt <- outpatient[MVO == "020" |
                             MVO == "243" |
                             MVO == "246", c(id_name, "INDATUMA"), with = FALSE][, `:=`
                                                                                 (visit_date = as.IDate(as.character(INDATUMA), format = "%Y%m%d"),
                                                                                   nursing_home = 1)][order(visit_date), .SD[1], by = id_name]
  
  # merge with cohort
  cohort_geo[, nursing_home := NA_real_]
  cohort_nursing <- retrieve_past_info(
    dt = cohort_geo,
    dictionary = nursing_dt,
    id_name = id_name,
    date_name = "visit_date",
    var_name = "nursing_home",
    lookback_months = Inf,
    fill_with_zero = TRUE # set remaining NA to zero
  )
  
  ################################################################################
  ### Calendar year
  ################################################################################
  cohort_year <- cohort_nursing[, calendar_year := year(as.IDate(visit_date))]
  
  ################################################################################
  ### Primary outcome all-cause mortality
  ################################################################################
  # Prepare datasets
  # All-cause death
  death_dt <- cohort_year[, c(id_name, "DODSDAT"), with = FALSE][, unique(.SD)]
  
  # CV death
  cv_death_dt <- merged_ckd[, c(id_name, "DODSDAT", "ULORSAK"), with = FALSE][, unique(.SD)]
  
  # MI and Stroke
  mi_stroke_dt <- merge(inpatient, death_dt, by = id_name, all.x = TRUE)
  
  # KRT date
  krt_dt <- merged_ckd[, c(id_name, "krt_startdate", "krt_modality"), with = FALSE][, unique(.SD)]
  
  # Define outcomes
  outcomes_list <- list(
    # All-cause death: uses DODSDAT from cohort
    death = list(
      # dataset with all-cause death dates
      dataset = death_dt,
      # column with date of death
      date_col = "DODSDAT",
      # no code column needed
      code_col = NULL,
      # no ICD filtering
      codes = NULL
    ),
    
    # Cardiovascular death: uses DODSDAT and ULORSAK from merged CKD dataset
    cvdeath = list(
      # dataset with cause-specific death
      dataset = cv_death_dt,
      # date of death
      date_col = "DODSDAT",
      # ICD code column for cause of death
      code_col = "ULORSAK",
      # pattern to select cardiovascular deaths
      codes = "^I"
    ),
    
    # Myocardial infarction (MI): ICD codes I21, I22, I23 from inpatient dataset
    mi = list(
      # inpatient dataset filtered for cohort
      dataset = mi_stroke_dt,
      # date of diagnoses registration
      date_col = "INDATUMA",
      # diagnosis code column
      code_col = "HDIA",
      # ICD codes for MI
      codes = "^I21|^I22|^I23"
    ),
    
    # Stroke: ICD codes I60-I64 from inpatient dataset
    stroke = list(
      # inpatient dataset filtered for cohort
      dataset = mi_stroke_dt,
      # date of diagnoses registration
      date_col = "INDATUMA",
      # diagnosis code column
      code_col = "HDIA",
      # ICD codes for stroke
      codes = "^I60|^I61|^I62|^I63|^I64"
    ),
    
    # KRT: uses krt_startdate
    KRT = list(
      # dataset with all-cause death dates
      dataset = krt_dt,
      # column with date of death
      date_col = "krt_startdate",
      # no code column needed
      code_col = NULL,
      # no ICD filtering
      codes = NULL,
      # add other columns
      extra_cols = "krt_modality"
    )
  )
  
  # Define windows
  windows <- list("1y" = 1,
                  "2y" = 2,
                  "inf" = 100)
  
  # add outcomes to cohort_dt
  cohort_outcomes <- add_multiple_outcomes(
    cohort_dt = cohort_year,
    windows = windows,
    id_name = id_name,
    end_follow_up = as.IDate("2024-12-31"),
    outcomes_list = outcomes_list
  )
  
  # add non-CV death and MACE outcomes
  for (window in names(windows)) {
    # Non-CV death: death event but not CV death
    cohort_outcomes[, paste0("noncvdeath_", window) :=
                      fifelse(get(paste0("event_death_", window)) == 1 &
                                get(paste0("event_cvdeath_", window)) == 0, 1L, 0L)]
    
    # MACE: composite of CV death, MI, Stroke
    # Event indicator columns
    mace_events <- c(
      paste0("event_cvdeath_", window),
      paste0("event_mi_", window),
      paste0("event_stroke_", window)
    )
    
    # Event date columns
    mace_dates <- c(
      paste0("event_dt_cvdeath_", window),
      paste0("event_dt_mi_", window),
      paste0("event_dt_stroke_", window)
    )
    
    # Compute MACE indicator: 1 if any of the components occurred
    cohort_outcomes[, paste0("event_mace_", window) := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = mace_events]
    
    # Compute MACE date: earliest date among components
    cohort_outcomes[, paste0("event_dt_mace_", window) := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = mace_dates]
    
    # Compute time to MACE from visit date
    cohort_outcomes[, paste0("time2event_mace_", window) :=
                      as.numeric(get(paste0("event_dt_mace_", window)) - visit_date)]
  }
  
  # ---- follow-up end date for the 2-year window ----
  # defined here (rather than locally inside the hospitalization/KRT
  # derivation below) so it's also carried into cohort_final/baseline,
  # available for every patient - including those with no hospitalizations
  # at all, who never appear in long_cohort_hosp/long_cohort_krt
  cohort_outcomes[, followup_end := as.IDate(visit_date + time2event_death_2y)]
  
  ################################################################################
  ### Secondary outcome hospitalization
  ################################################################################
  if (cohort_name == "cohort") {
    # ---- Build hospitalization table for cohort patients only ----
    # Keep one row per unique (patient, start, stop) combination in case of
    # exact duplicate records in the raw inpatient register.
    raw_hosp_dt <- inpatient[LOPNR %in% working_cohort$LOPNR, .SD[1], 
                             by = c(id_name,
                                    "INDATUMA", 
                                    "UTDATUMA")][, c(id_name,
                                                     "INDATUMA",
                                                     "UTDATUMA"), with = FALSE]
    
    # rename to more descriptive column names
    setnames(raw_hosp_dt,
             c("INDATUMA", "UTDATUMA"),
             c("hosp_start", "hosp_stop"))
    
    # ---- Attach visit_date (treatment decision date) to each hospitalization ----
    # Inner join: only hospitalizations for patients that are in cohort are kept.
    hosp_sub <- raw_hosp_dt[cohort_outcomes[, .(LOPNR,
                                                visit_date,
                                                trt,
                                                time2event_death_2y,
                                                DODSDAT,
                                                followup_end)], on = "LOPNR", nomatch = 0]
    
    # ---- Keep only hospitalizations relevant to the post-visit period ----
    # Two cases we want:
    #   a) hospitalization starts strictly after visit_date
    #   b) hospitalization is ongoing at visit_date (started before/at, ends at/after)
    hosp_sub <- hosp_sub[hosp_start > visit_date |
                           (hosp_start <= visit_date &
                              hosp_stop >= visit_date)]
    
    # For hospitalizations ongoing at visit_date, truncate the start date to
    # visit_date itself -- we only care about the part of the stay that happens
    # from the treatment decision onward.
    hosp_sub[hosp_start <= visit_date &
               hosp_stop >= visit_date, 
             hosp_start := visit_date]
    
    # ---- Restrict to 2-year follow-up window from visit_date ----
    hosp_sub <- hosp_sub[hosp_start <= followup_end]              # drop stays starting after window
    
    # ---- Collapse overlapping/touching hospitalizations into single episodes ----
    # Patients can be registered at multiple departments during overlapping (or
    # even touching) periods. We merge these into one continuous "episode" of
    # hospitalization, using the classic sort + cumulative-max-end trick:
    #   - sort by patient and start date
    #   - track the running max hosp_stop seen so far (per patient)
    #   - a new episode begins whenever a hospitalization's start date comes
    #     strictly after that running max stop date (i.e. there's a real gap)
    #   - touching stays (start == previous max stop) are NOT strictly after,
    #     so they get merged into the same episode -- this is intentional
    setorder(hosp_sub, LOPNR, hosp_start, hosp_stop)
    
    hosp_sub[, prev_max_stop := shift(cummax(as.numeric(hosp_stop)),
                                      fill = -Inf), 
             by = LOPNR]
    
    hosp_sub[, episode := cumsum(as.numeric(hosp_start) > prev_max_stop),
             by = LOPNR]
    
    # collapse each episode into a single row spanning its full duration
    setorder(hosp_sub, LOPNR, episode, hosp_stop)
    merge_cols <- setdiff(names(hosp_sub),
                          c("LOPNR", "episode", "hosp_start", "hosp_stop", "prev_max_stop"))
    
    long_cohort_hosp <- hosp_sub[, c(list(hosp_start = min(hosp_start),
                                          hosp_stop  = max(hosp_stop)),
                                     lapply(.SD, function(x) x[.N])),
                                 by = .(LOPNR, episode), .SDcols = merge_cols]
    
    long_cohort_hosp[, episode := NULL]
    setnames(long_cohort_hosp, "visit_date", "decision_date")
    setcolorder(
      long_cohort_hosp,
      c(
        "LOPNR",
        "decision_date",
        "trt",
        "hosp_start",
        "hosp_stop",
        "DODSDAT",
        "followup_end",
        "time2event_death_2y"
      )
    )
    
    # ---- Build KRT (dialysis) modality periods per patient ----
    # snr_rrt can have multiple KRT rows per patient (e.g. a PD -> HD switch).
    # Keep all of them instead of just the first, and turn them into periods:
    # each modality holds from its own krt_startdate until the NEXT switch.
    # Patients with no dialysis at all have no rows here, so they'll simply
    # get NA (never dialysis) when joined onto hosp_sub below.
    raw_krt_dt <- unique(snr_rrt_long[LOPNR %in% working_cohort$LOPNR &
                                        !is.na(krt_startdate) &
                                        !is.na(krt_modality),
                                      c(id_name, "krt_modality", "krt_startdate"), with = FALSE])
    
    # --- UNTRACED -------------------------------------------------------------
    # Assume untraced patients stay on their initial KRT modality
    long_cohort_krt <- raw_krt_dt[krt_modality != "UNTRACED"]
    
    # --- RECOVERED ------------------------------------------------------------
    # RECOVERED always coincides with another modality on the same date - drop it
    long_cohort_krt <- long_cohort_krt[krt_modality != "RECOVERED"]
    
    # --- Same-date conflicts: HD > PD > TX -------------------------------------
    # When multiple modalities are logged on the same date, keep HD if present,
    # else PD, else TX (this also drops any lingering GRAFT FAILURE row, since
    # it always coincides with one of these)
    
    # HD wins if present
    HD_dates <- unique(long_cohort_krt[krt_modality == "HD", .(LOPNR, krt_startdate)])
    long_cohort_krt[HD_dates, has_HD := TRUE, on = .(LOPNR, krt_startdate)]
    long_cohort_krt <- long_cohort_krt[is.na(has_HD) | krt_modality == "HD"]
    long_cohort_krt[, has_HD := NULL]
    
    # else PD wins if present
    PD_dates <- unique(long_cohort_krt[krt_modality == "PD", .(LOPNR, krt_startdate)])
    long_cohort_krt[PD_dates, has_PD := TRUE, on = .(LOPNR, krt_startdate)]
    long_cohort_krt <- long_cohort_krt[is.na(has_PD) | krt_modality == "PD"]
    long_cohort_krt[, has_PD := NULL]
    
    # else TX wins if present
    tx_dates <- unique(long_cohort_krt[krt_modality == "TX", .(LOPNR, krt_startdate)])
    long_cohort_krt[tx_dates, has_tx := TRUE, on = .(LOPNR, krt_startdate)]
    long_cohort_krt <- long_cohort_krt[is.na(has_tx) | krt_modality == "TX"]
    long_cohort_krt[, has_tx := NULL]
    
    # bring in each patient's followup_end (already computed on cohort_outcomes
    # above - a patient can be on dialysis with no hospitalization at all, so
    # this can't be sourced from hosp_sub, which only has hospitalized patients)
    long_cohort_krt <- unique(cohort_outcomes[, .(LOPNR, followup_end)])[long_cohort_krt, on = "LOPNR"]
    
    # drop KRT events starting after this patient's own follow-up window -
    # they can't affect anything we're modeling, and dropping them here means
    # the fallback below is the ONLY thing bounding the true last in-window
    # period, rather than a real switch date sometimes doing it instead
    long_cohort_krt <- long_cohort_krt[krt_startdate <= followup_end]
    
    setorder(long_cohort_krt, LOPNR, krt_startdate)
    
    # the last period per patient has no "next" switch - bound it at
    # followup_end instead of leaving it open (nothing past followup_end is
    # ever used anyway, so this keeps krt_stopdate meaningful/non-NA
    # whenever a patient has any dialysis at all)
    long_cohort_krt[, krt_stopdate := shift(krt_startdate, type = "lead"), by = LOPNR]
    long_cohort_krt[is.na(krt_stopdate), krt_stopdate := followup_end]
    
    # save the long-format hospitalization + KRT tables for later use
    # (multi-state modeling of home/hospital/death transitions, PART 12)
    save(long_cohort_hosp, long_cohort_krt, file = "Data/long_cohort_hosp_krt.Rdata")
  }
  
  ################################################################################
  ### Categorize some covariates
  ################################################################################
  # Categorize and set factors for all variables in one call
  cohort_final <- cohort_outcomes[, `:=`(
    age_cat = cut(
      age,
      breaks = c(0, 70, 75, 80, 10000),
      labels = c("65-69", "70-74", "75-79", ">=80"),
      right = FALSE
    ),
    
    Davies_score_cat = factor(cut(
      Davies_score,
      breaks = c(0, 2, 4, 10),
      labels = c("<2", "2-4", ">4"),
      right = FALSE
    )),
    
    egfr_cat = factor(cut(
      egfr2021,
      breaks = c(0, 10, 15, 20),
      labels = c("<10", "10-14", "15-20"),
      right = FALSE
    )),
    
    iron_cat = factor(
      fifelse(iron_iv == 1, 1, fifelse(iron_po == 1, 2, 0)),
      levels = c(0, 1, 2),
      labels = c("No iron", "IV iron", "PO iron")
    ),
    
    sbp_cat = factor(cut(
      sbp,
      breaks = c(0, 120, 140, 160, 1000),
      labels = c("<120", "120-139", "140-159", ">160"),
      right = FALSE
    )),
    
    dbp_cat = factor(cut(
      dbp,
      breaks = c(0, 80, 90, 100, 1000),
      labels = c("<80", "80-89", "90-99", ">100"),
      right = FALSE
    )),
    
    calendar_year_cat = factor(cut(
      calendar_year,
      breaks = c(0, 2013, 2018, 10000),
      labels = c("2007-2012", "2013-2017", "2018-2021"),
      right = FALSE
    )),
    
    female = factor(female),
    cancer = factor(cancer),
    ihd = factor(ihd),
    pvd = factor(pvd),
    hf = factor(hf),
    dm = factor(dm),
    scvd = factor(scvd),
    copd = factor(copd),
    cirr = factor(cirr),
    psycho = factor(psycho),
    acs = factor(acs),
    hyperten = factor(hyperten),
    vhd = factor(vhd),
    cevd = factor(cevd),
    af = factor(af),
    arrh = factor(arrh),
    lung = factor(lung),
    thrombo = factor(thrombo),
    liver = factor(liver),
    fracture = factor(fracture),
    aki = factor(aki),
    bblock = factor(bblock),
    calblock = factor(calblock),
    diuretic = factor(diuretic),
    rasi = factor(rasi),
    lipid = factor(lipid),
    phosbinder = factor(phosbinder),
    esa = factor(esa),
    vitamind = factor(vitamind),
    digoxin = factor(digoxin),
    vasodilator = factor(vasodilator),
    antiplatelet = factor(antiplatelet),
    anticoag = factor(anticoag),
    edu = factor(edu),
    clinic_level = factor(clinic_level),
    region = factor(region)
  )]
  
  # save cohort
  save(cohort_final, file = file.path(paste0(
    "Data/analysis_data_", cohort_name, "_new.Rdata"
  )))
}
