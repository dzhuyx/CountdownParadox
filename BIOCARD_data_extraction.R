###############################################################################
## BIOCARD Data Extraction for Countdown Paradox Analysis
##
## Purpose: Extract and save clean analysis-ready datasets from BIOCARD:
##   1. Survival outcome data (onset.age, censor.age, event indicator)
##   2. Demographics and covariates (sex, education, APOE4)
##   3. Per-subject CSF measurement summary (mean CSF age, n measures, span)
##
## AABCs (age at biomarker-clock event) are computed by the SILA pipeline
## (BIOCARD_csf_sila.R, BIOCARD_plasma_sila.R), not in this script.
##
## Output files (saved to data/ subdirectory):
##   - survival_data.rda          (survival outcome + demographics + covariates)
##   - analysis_data_merged.rda   (per-subject CSF summary + survival + demographics)
##
## Data source directory:
##   BIOCARD/BIOCARD_codes/  (all .xlsx and .csv data files)
##
## Author: Yuxin Zhu
## Date: February 2026
###############################################################################

rm(list = ls())

library(readxl)
library(dplyr)
library(survival)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set CP_PROJECT_ROOT to the analysis root (the folder that holds data/ and results/).")
biocard_dir  <- Sys.getenv("CP_BIOCARD_DIR")
if (biocard_dir == "") stop("CP_BIOCARD_DIR is not set. Set it to the directory containing the raw BIOCARD data files.")
output_dir   <- file.path(project_root, "data")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== BIOCARD Data Extraction for Countdown Paradox ===\n")
cat("Data source:", biocard_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

###############################################################################
## PART 1: Load Raw Data Files
###############################################################################

cat("--- Loading raw data files ---\n")

# Demographics
dat_demo <- read_excel(path = file.path(biocard_dir,
                                         "BIOCARD_Demographics_2024.08.07.xlsx"))
cat("Demographics:", nrow(dat_demo), "rows\n")

# Diagnosis
dat_dx <- read_excel(path = file.path(biocard_dir,
                                       "BIOCARD_DiagnosisData_2024.09.08.xlsx"))
dat_dx$DIAGDATE <- do.call("c", lapply(dat_dx$DIAGDATE, function(x) {
      result <- as.Date(strsplit(as.character(x), " UTC")[[1]],
                        "%Y-%m-%d")
      return(result)
}))
dat_dx$DOB <- do.call("c", lapply(dat_dx$DOB, function(x) {
      result <- as.Date(strsplit(as.character(x), " UTC")[[1]],
                        "%Y-%m-%d")
      return(result)
}))
dat_dx$AD_primary <- 0
dat_dx$AD_primary[dat_dx$PROBADIF == 1 | dat_dx$POSSADIF == 1] <- 1
dat_dx$AD_contrib <- 0
dat_dx$AD_contrib[dat_dx$PROBAD == 1 | dat_dx$POSSAD == 1] <- 1
cat("Diagnosis:", nrow(dat_dx), "rows,",
    length(unique(dat_dx$SUBJECT_ID)), "subjects\n")

# CSF Lumipulse biomarkers
dat_csf <- read_excel(path = file.path(biocard_dir,
                                        "BIOCARD_CSF_Lumipulse_NFL_GFAP_Data_2025.02.24.xlsx"))
dat_csf$CSF_DATE <- do.call("c", lapply(dat_csf$CSF_DATE, function(x) {
      result <- as.Date(strsplit(as.character(x), " UTC")[[1]],
                        "%Y-%m-%d")
      return(result)
}))
# Rename columns
colnames(dat_csf)[colnames(dat_csf) == "NF_Light_pg_mL"] <- "NfL"
colnames(dat_csf)[colnames(dat_csf) == "GFAP_pg_mL"] <- "GFAP"
cat("CSF Lumipulse:", nrow(dat_csf), "rows,",
    length(unique(dat_csf$SUBJECT_ID)), "subjects\n")

# Genetics (APOE)
dat_gen <- read_excel(path = file.path(biocard_dir,
                                        "BIOCARD_Genetics_Data_2023.03.28.xlsx"))
dat_gen_new <- read_excel(path = file.path(biocard_dir,
                                            "New participants_BIOCARD ApoE Genotypes 2023-2024.xlsx"))
dat_gen_new <- merge(dat_gen_new, dat_demo[, c("JHUANONID", "LETTERCODE",
                                               "NIHID", "SUBJECT_ID")],
                     by = "LETTERCODE",
                     all.x = T)
dat_gen <- merge(dat_gen, dat_gen_new,
                 by = c("JHUANONID", "LETTERCODE",
                        "NIHID", "SUBJECT_ID", "APOECODE"),
                 all = T)
# Replace APOECODE values larger than 10 by the number divided by 10
dat_gen$APOECODE[dat_gen$APOECODE > 10] <- dat_gen$APOECODE[dat_gen$APOECODE > 10] / 10
# Create dichotomous APOE4 indicator
dat_gen$APOE4 <- NA
dat_gen$APOE4[dat_gen$APOECODE %in% c(2.2, 2.3, 3.3)] <- 0
dat_gen$APOE4[dat_gen$APOECODE %in% c(2.4, 3.4, 4.4)] <- 1
cat("Genetics:", nrow(dat_gen), "rows\n")

# Create dat_apoe (used later for merging)
dat_apoe <- dat_gen[, c("SUBJECT_ID", "LETTERCODE", "APOECODE", "APOE4")]
colnames(dat_apoe)[colnames(dat_apoe) == "APOE4"] <- "apoe4"

# Cognitive data
dat_cog <- read_excel(path = file.path(biocard_dir,
                                        "BIOCARD_CognitiveData_With_Composite_Scores_2024.09.07.xlsx"))
dat_cog$VISITDATE <- do.call("c", lapply(dat_cog$VISITDATE, function(x) {
      result <- as.Date(strsplit(as.character(x), " UTC")[[1]],
                        "%Y-%m-%d")
      return(result)
}))
dat_cog$MMSE <- dat_cog$C1201D
cat("Cognitive:", nrow(dat_cog), "rows\n\n")


###############################################################################
## PART 2: Create Survival Data from Diagnosis
##         (Proven algorithm from BIOCARD — do not modify)
###############################################################################

cat("--- Creating survival data ---\n")

dat.by.ID <- split(dat_dx, dat_dx$SUBJECT_ID)

# Source shared helper functions (FindJump, GetSurv) — single canonical copy
source(file.path(project_root, "shared_utils.R"))

dat_surv <- do.call(rbind, lapply(dat.by.ID, GetSurv))
dat_surv <- merge(dat_surv, dat_demo,
                  by = c("SUBJECT_ID", "LETTERCODE"))
cat("Survival data created:", nrow(dat_surv), "subjects,",
    sum(dat_surv$d == 1), "events\n\n")


###############################################################################
## PART 3: Apply Exclusion Criteria
###############################################################################

cat("--- Applying exclusion criteria ---\n")

# List A: Withdrawn subjects
listA <- read.csv(file = file.path(biocard_dir, "list_A_122021.csv"))
cat("List A (withdrawn):", nrow(listA), "subjects\n")

# List B: Impaired at baseline
listB <- read_excel(path = file.path(biocard_dir,
                                      "LIST_B_IMPAIRED_AT_BASELINE.09.22.2015.xlsx"))
cat("List B (impaired at baseline):", nrow(listB), "subjects\n")

# Subjects with a single CSF measure (need >= 2 for individual regression)
n_CSF <- dat_csf %>%
      filter(!is.na(AB42AB40)) %>%
      group_by(SUBJECT_ID) %>%
      summarise(n_CSF = n())
list_1CSF <- setdiff(unique(dat_csf$SUBJECT_ID),
                          n_CSF$SUBJECT_ID[n_CSF$n_CSF >= 2])
cat("Single CSF measure:", length(list_1CSF), "subjects\n")

# New BIOCARD enrollees (SUBJECT_ID >= 400)
list_new <- dat_demo$SUBJECT_ID[which(dat_demo$SUBJECT_ID >= 400)]
cat("New enrollees (ID >= 400):", length(list_new), "subjects\n")

# Subjects with onset before CSF baseline — excluded from survival analysis
# This exclusion is correct for the survival analysis
# (subjects who already had MCI before their first CSF draw cannot contribute
# at-risk person-time). It was REMOVED from BIOCARD_csf_sila.R because SILA
# trajectory fitting benefits from all subjects' biomarker data regardless of
# MCI onset timing.
dat_csfo <- dat_csf  # preserve original CSF data
dat_csf_temp <- merge(dat_csf, dat_surv,
                      by = c("SUBJECT_ID", "LETTERCODE"))
dat_csf_temp$CSF_age <- as.numeric(dat_csf_temp$CSF_DATE - dat_csf_temp$DOB) / 365.25
dat_csf_1 <- dat_csf_temp[which(dat_csf_temp$onset.age > dat_csf_temp$CSF_age), ]
ID_onsetbeforebaseline <- setdiff(unique(dat_csfo$SUBJECT_ID),
                                  unique(dat_csf_1$SUBJECT_ID))
cat("Onset before CSF baseline:", length(ID_onsetbeforebaseline), "subjects\n")

# Combine all exclusion IDs
ID_exclude <- unique(c(
      listA$ID,
      listB$STUDY_ID,
      ID_onsetbeforebaseline,
      list_1CSF,
      list_new
))
cat("Total unique excluded:", length(ID_exclude), "subjects\n")

# Apply exclusions to CSF data
# NOTE: Use original CSF data (dat_csfo) for the included subjects
#       This preserves all CSF measures including those before onset
dat_csf <- dat_csfo[dat_csfo$SUBJECT_ID %in%
                          setdiff(unique(dat_csfo$SUBJECT_ID), ID_exclude), ]
cat("Included subjects with CSF:", length(unique(dat_csf$SUBJECT_ID)), "\n\n")


###############################################################################
## PART 4: Per-subject CSF measurement summary
###############################################################################

cat("--- Computing per-subject CSF measurement summary ---\n")

# Compute CSF_age
dat_dob <- dat_dx[-which(duplicated(dat_dx$SUBJECT_ID)),
                  c("SUBJECT_ID", "DOB")]
dat_csf <- merge(dat_csf, dat_dob, by = "SUBJECT_ID")
dat_csf$CSF_age <- as.numeric(dat_csf$CSF_DATE - dat_csf$DOB) / 365.25

# One row per subject with mean CSF age (used as the analysis-cohort skeleton)
csf_subject_summary <- dat_csf %>%
      group_by(SUBJECT_ID, LETTERCODE) %>%
      summarise(mean_CSF_age = mean(CSF_age), .groups = "drop")

cat("CSF subject summary:", nrow(csf_subject_summary), "subjects\n\n")

###############################################################################
## PART 5: Prepare Clean Output Datasets
###############################################################################

cat("--- Preparing output datasets ---\n")

# 5.1: Survival data (restricted to CSF subjects)
dat_surv_clean <- dat_surv[dat_surv$SUBJECT_ID %in% csf_subject_summary$SUBJECT_ID, ]
dat_surv_clean <- dat_surv_clean[, c("SUBJECT_ID", "LETTERCODE",
                                      "baseline.age", "onset.age", "d",
                                      "diag", "DOB", "censor.age",
                                      "SEX", "EDUC")]
dat_surv_clean$Sex_F <- as.numeric(dat_surv_clean$SEX == 2)

# 5.2: APOE data (restricted to CSF subjects)
dat_apoe_clean <- dat_apoe[dat_apoe$SUBJECT_ID %in% csf_subject_summary$SUBJECT_ID, ]
dat_apoe_clean <- dat_apoe_clean[!duplicated(dat_apoe_clean$SUBJECT_ID), ]
dat_apoe_clean <- dat_apoe_clean[, c("SUBJECT_ID", "apoe4")]

# 5.3: Merged analysis dataset
analysis_data <- merge(csf_subject_summary, dat_surv_clean,
                       by = c("SUBJECT_ID", "LETTERCODE"), all.x = TRUE)
analysis_data <- merge(analysis_data, dat_apoe_clean,
                       by = "SUBJECT_ID", all.x = TRUE)

# Standardized education
analysis_data$EDUC_z <- as.numeric(scale(analysis_data$EDUC))

# Number of CSF measures per subject
n_CSF_per_subject <- dat_csf %>%
      group_by(SUBJECT_ID) %>%
      summarise(n_CSF_measures = n(),
                CSF_span_years = as.numeric(max(CSF_DATE) - min(CSF_DATE)) / 365.25)
analysis_data <- merge(analysis_data, n_CSF_per_subject,
                       by = "SUBJECT_ID", all.x = TRUE)

cat("Analysis dataset:", nrow(analysis_data), "subjects\n")
cat("  Events:", sum(analysis_data$d == 1), "\n")
cat("  Censored:", sum(analysis_data$d == 0), "\n\n")


###############################################################################
## PART 6: Save Output Files
###############################################################################

cat("--- Saving output files ---\n")

# Survival data
survival_data <- dat_surv_clean
save(survival_data, file = file.path(output_dir, "survival_data.rda"))
cat("Saved:", file.path(output_dir, "survival_data.rda"), "\n")

# Complete merged analysis dataset
save(analysis_data, file = file.path(output_dir, "analysis_data_merged.rda"))
cat("Saved:", file.path(output_dir, "analysis_data_merged.rda"), "\n")

# Exclusion accounting (for CONSORT flow diagram)
exclusion_counts <- data.frame(
      step = c("Total with any CSF",
               "Excluded: List A (withdrawn)",
               "Excluded: List B (impaired at baseline)",
               "Excluded: Onset before CSF baseline",
               "Excluded: Single CSF measure",
               "Excluded: New enrollees (ID >= 400)",
               "Total excluded (unique)",
               "Included in analysis"),
      n = c(length(unique(dat_csfo$SUBJECT_ID)),
            length(intersect(unique(dat_csfo$SUBJECT_ID), listA$ID)),
            length(intersect(unique(dat_csfo$SUBJECT_ID), listB$STUDY_ID)),
            length(ID_onsetbeforebaseline),
            length(list_1CSF),
            length(intersect(unique(dat_csfo$SUBJECT_ID), list_new)),
            length(ID_exclude),
            length(unique(csf_subject_summary$SUBJECT_ID)))
)
save(exclusion_counts, file = file.path(output_dir, "exclusion_counts.rda"))
write.csv(exclusion_counts, file = file.path(output_dir, "exclusion_counts.csv"),
          row.names = FALSE)
cat("Saved:", file.path(output_dir, "exclusion_counts.csv"), "\n")

# Provenance metadata
provenance <- data.frame(
      script = "BIOCARD_data_extraction.R",
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version = paste0(R.version$major, ".", R.version$minor),
      stringsAsFactors = FALSE
)
write.csv(provenance, file = file.path(output_dir, "data_extraction_provenance.csv"),
          row.names = FALSE)
cat("Saved: data_extraction_provenance.csv\n")

cat("\n=== Extraction complete ===\n")
cat("\nExclusion summary:\n")
print(exclusion_counts, row.names = FALSE)

cat("\nKey variables in analysis_data:\n")
cat("  mean_CSF_age, n_CSF_measures, CSF_span_years  (CSF measurement info)\n")
cat("  onset.age, censor.age, d  (survival outcome)\n")
cat("  baseline.age, Sex_F, EDUC, EDUC_z, apoe4  (covariates)\n")

