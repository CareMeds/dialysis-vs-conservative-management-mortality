################################################################################
### Decision for dialysis versus conservative management
### PART 0 - Data preparation
################################################################################

# remove history
rm(list = ls(all.names = TRUE))

# load data
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/"
)
load("Data/cleaned/snr_ckd.Rdata")        # CKD patients
load("Data/cleaned/snr_rrt.Rdata")        # Renal Replacement Therapy (transplantation)
load("Data/cleaned/snr_hdpd.Rdata")       # lab data from hemodialysis or peritoneal dialysis patients
load("Data/cleaned/snr_death.Rdata")      # vital status

# load functions
source("Code/utils/data_manipulation.R")

# set system to English (US), to format dates correctly
Sys.setlocale("LC_TIME", "en_US.UTF-8")

# load libraries
library(data.table)
library(patchwork)

################################################################################
### Select relevant variables from CKD file ####################################
################################################################################
# Convert to data.table
setDT(SOS_CKDFINAL2024)

# Select relevant columns and filter non-missing birthdate
ckd <- SOS_CKDFINAL2024[, c(
  "LOPNR",
  "decision_date1",
  "decision_date2",
  "decision_date3",
  "decision_modality1",
  "decision_modality2",
  "decision_modality3",
  "decision_tx1",
  "decision_tx_date1",
  "info_date1",
  "info_type1",
  "birthdate",
  "deathdate",
  "female",
  "visit_date",
  "visittype",
  "screa",
  "sbp",
  "dbp",
  "calcium_total",
  "phosphate",
  "albumin",
  "hb",
  "esa",
  "iron_med",
  "iron_type",
  "prd_cat",
  "clinic",
  "county",
  "crp"
)]

# convert to date
ckd[, `:=`(
  birthdate = as.IDate(birthdate, format("%d%b%Y")),
  DODSDAT = as.IDate(deathdate, format("%m/%d/%Y")),
  decision_date1 = as.IDate(decision_date1, format("%m/%d/%Y")),
  decision_date2 = as.IDate(decision_date2, format("%m/%d/%Y")),
  decision_date3 = as.IDate(decision_date3, format("%m/%d/%Y")),
  decision_tx_date1 = as.IDate(decision_tx_date1, format("%m/%d/%Y")),
  visit_date = as.IDate(visit_date, format("%m/%d/%Y"))
)]

# only keep unique rows
ckd <- ckd[, unique(.SD)]
nrow(ckd) # N = 403070

# Sort by LOPNR and visit_date
setDT(ckd)
setkey(ckd, LOPNR, visit_date)

################################################################################
### Attach information from other cohorts
################################################################################
# death information
snr_death <- UT_R_DORS_123160_2023[, c("LOPNR", "DODSDAT", "ULORSAK")]
setDT(snr_death)
snr_death[, `:=`(DODSDAT = as.IDate(as.character(DODSDAT), format("%Y%m%d")))]

# renal replacement therapy information
snr_rrt <- unique(SOS_KRTDATA2024[, c(
  "LOPNR",
  "birthdate",
  "clinic",
  "county",
  "deathdate",
  "female",
  "krt_start",
  "krt_startdate",
  "modality_cat",
  "prd_cat"
)])
setDT(snr_rrt)
setnames(snr_rrt, "modality_cat", "krt_modality")
snr_rrt[, `:=`(
  birthdate = as.IDate(birthdate, format("%d%b%Y")),
  DODSDAT = as.IDate(deathdate, format("%m/%d/%Y")),
  krt_startdate = as.IDate(krt_startdate, format("%d%b%Y"))
)]

# some IDs have multiple rows
snr_rrt[, uniqueN(LOPNR)]          # 43276
snr_rrt[, .N, by = "LOPNR"][N > 1] # 21292
snr_rrt[LOPNR == 18316963, ]       # example

# long snr RRT data
snr_rrt_long <- copy(snr_rrt)

# dup_check <- snr_rrt_long[, .N, by = LOPNR][N > 1, LOPNR]
# 
# snr_rrt_long[LOPNR %in% dup_check, .(
#   n_rows        = .N,
#   n_dates       = uniqueN(krt_startdate),
#   n_modalities  = uniqueN(krt_modality)
# ), by = LOPNR][, .(
#   # patients where every duplicate row shares the same date
#   only_same_date       = sum(n_dates == 1),
#   # patients where duplicates span more than one date (real switches)
#   spans_multiple_dates = sum(n_dates > 1)
# )]

save(snr_rrt_long, file = "Data/cleaned/snr_rrt_long.Rdata")

# only keep krt_start == 1 information, this does not delete IDs
snr_rrt <- snr_rrt[krt_start == 1, ]
snr_rrt[, uniqueN(LOPNR)]      # 43276
snr_rrt[, .N, by = "LOPNR"][N > 1] # 1
# example has two krt start dates
# but is removed later on because they are not in primary CKD data set
snr_rrt[LOPNR == 35333149, ]

# dialysis information, note to not use future information
snr_hdpd <- unique(SOS_DIALYSISDATA2024[, c("LOPNR",
                                            "birthdate",
                                            "deathdate",
                                            "female",
                                            "vascular_access",
                                            "visit_date")])
setDT(snr_hdpd)

# Define dates and visittype
snr_hdpd[, `:=`(
  birthdate = as.IDate(birthdate, format("%d%b%Y")),
  DODSDAT = as.IDate(deathdate, format("%m/%d/%Y")),
  visit_date = as.IDate(visit_date, format("%m/%d/%Y"))
)]

# add data on vital status
merged_ckd <- merge(
  ckd,
  snr_death,
  by = "LOPNR",
  all.x = TRUE,
  suffixes = c("", ".death")
)

# add data on renal replacement therapy
merged_ckd <- merge(
  merged_ckd,
  snr_rrt,
  by = "LOPNR",
  all.x = TRUE,
  suffixes = c("", ".rrt")
)

# add data on birthdate and sex
merged_ckd <- merge(
  merged_ckd,
  snr_hdpd,
  by = c("LOPNR", "visit_date"),
  all = TRUE,
  suffixes = c("", ".hdpd")
)
sum(is.na(merged_ckd$birthdate))

# take the minimum of these merged variables
colnames.merge <- c("clinic", "county", "birthdate", "DODSDAT", "female", "prd_cat")
for (colname in colnames.merge) {
  # print(colSums(is.na(merged_ckd[, .SD, .SDcols = patterns(paste0("^", colname))])))
  merged_ckd[, (colname) := do.call(fcoalesce, .SD), .SDcols = patterns(paste0("^", colname))]
  # print(colSums(is.na(merged_ckd[, .SD, .SDcols = patterns(paste0("^", colname))])))
}

# Remove intermediate columns
merged_ckd[, c(grep("\\.hdpd$", colnames(merged_ckd), value = TRUE),
               grep("\\.rrt$", colnames(merged_ckd), value = TRUE)) := NULL]

# Set key for efficient joins and lookups
setkey(merged_ckd, LOPNR, visit_date)

# Check missings
colSums(is.na(merged_ckd))
nrow(merged_ckd) # N = 576292

################################################################################
### Registration errors
################################################################################
# ensure birthdate is known
merged_ckd <- merged_ckd[!is.na(birthdate)]
nrow(merged_ckd) # N = 576292

# ensure visit_date is before date of death or has not died yet
merged_ckd <- merged_ckd[visit_date < DODSDAT | is.na(DODSDAT)]
nrow(merged_ckd) # N = 576247

################################################################################
### Calculate age and eGFR at each visit
################################################################################
merged_ckd <- calculate_age(dt = merged_ckd,
                            birthdate_col = "birthdate",
                            visitdate_col = "visit_date")

# fill serum creatinine
merged_ckd <- retrieve_past_info(
  dt = merged_ckd,
  dictionary = merged_ckd,
  id_name = "LOPNR",
  date_name = "visit_date",
  var_name = "screa",
  lookback_months = 12
)

# calculate eGFR using the 2021 creatinine equation
merged_ckd[, egfr2021 := ckd_epi_2021_cr(screa, age, female)]

# transform hemoglobin g/l into mmol/l
merged_ckd[, hb := hb * 0.6206 / 10]

# Ensure to only include visits up to 2023
end_date <- as.IDate("2024-12-31")
merged_ckd <- merged_ckd[visit_date <= end_date]

# save file
save(merged_ckd, end_date, file = "Data/merged_ckd.Rdata")
