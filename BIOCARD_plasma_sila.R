###############################################################################
## BIOCARD Plasma SILA Trajectory Fitting — p-tau181 (NTK/Simoa)
##
## Purpose: Run silaR on BIOCARD plasma p-tau181 longitudinal data (JHU phase
##   only) to estimate trajectory-based age at onset of abnormality (EAOA).
##
## Analyte: PTAU181 (pg/mL), already increasing with pathology, val0 = 0.8
## Phase: JHU only (Blood_DATE >= 2009-01-01)
## Exclusions: List A (withdrawn), List B (impaired at baseline),
##             onset before JHU baseline
##
## Input:  BIOCARD_Blood_NTK_Data_2024.09.18.xlsx
## Output: data/BIOCARD_plasma_SILA_intermediate.rda
##
## Usage:  source() from CountdownParadox_Analysis/ directory
##
## Author: Yuxin Zhu
## Date: March 2026
###############################################################################

rm(list = ls())
library(readxl)
library(dplyr)
library(silaR)
library(tibble)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set it to <project>/CountdownParadox_Analysis.")
out_dir      <- project_root
data_dir     <- file.path(project_root, "data")
biocard_dir  <- Sys.getenv("CP_BIOCARD_DIR")
if (biocard_dir == "") stop("CP_BIOCARD_DIR is not set. Set it to the directory containing the raw BIOCARD data files.")


###############################################################################
## SECTION 1: Load Data
###############################################################################

cat("================================================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("================================================================\n\n")

# --- 1a: NTK blood data --- #
ntk_file <- file.path(biocard_dir, "BIOCARD_Blood_NTK_Data_2024.09.18.xlsx")
dat_blood <- read_excel(ntk_file, sheet = 1)
cat(sprintf("NTK blood data: %d rows, %d subjects\n",
            nrow(dat_blood), length(unique(dat_blood$SUBJECT_ID))))

# --- 1b: DOB from diagnosis data (for exact age calculation) --- #
dat_dx <- read_excel(path = file.path(biocard_dir,
                                       "BIOCARD_DiagnosisData_2024.09.08.xlsx"))
dat_dx$DIAGDATE <- do.call("c", lapply(dat_dx$DIAGDATE, function(x) {
      as.Date(strsplit(as.character(x), " UTC")[[1]], "%Y-%m-%d")
}))
dat_dx$DOB <- do.call("c", lapply(dat_dx$DOB, function(x) {
      as.Date(strsplit(as.character(x), " UTC")[[1]], "%Y-%m-%d")
}))
dat_dob <- dat_dx[!duplicated(dat_dx$SUBJECT_ID), c("SUBJECT_ID", "DOB")]
cat(sprintf("DOB available for %d subjects\n", nrow(dat_dob)))

# --- 1c: Compute exact age at blood draw --- #
dat_blood$Blood_DATE <- as.Date(dat_blood$Blood_DATE)
dat_blood <- merge(dat_blood, dat_dob, by = "SUBJECT_ID", all.x = TRUE)
dat_blood$blood_age <- as.numeric(dat_blood$Blood_DATE - dat_blood$DOB) / 365.25

# Check DOB coverage
n_no_dob <- sum(is.na(dat_blood$DOB))
if (n_no_dob > 0) {
      cat(sprintf("WARNING: %d rows missing DOB — using BIRTHYEAR approximation\n", n_no_dob))
      # Fallback: mid-year approximation
      no_dob <- is.na(dat_blood$blood_age)
      dat_blood$blood_age[no_dob] <-
            as.numeric(format(dat_blood$Blood_DATE[no_dob], "%Y")) -
            dat_blood$BIRTHYEAR[no_dob] + 0.5
}
cat(sprintf("Age range: %.1f to %.1f\n", min(dat_blood$blood_age, na.rm = TRUE),
            max(dat_blood$blood_age, na.rm = TRUE)))

# --- 1d: Demographics --- #
dat_demo <- read_excel(path = file.path(biocard_dir,
                                         "BIOCARD_Demographics_2024.08.07.xlsx"))

# --- 1e: Genetics (APOE) --- #
dat_gen <- read_excel(path = file.path(biocard_dir,
                                        "BIOCARD_Genetics_Data_2023.03.28.xlsx"))
dat_gen_new <- read_excel(path = file.path(biocard_dir,
                                            "New participants_BIOCARD ApoE Genotypes 2023-2024.xlsx"))
dat_gen_new <- merge(dat_gen_new, dat_demo[, c("JHUANONID", "LETTERCODE",
                                               "NIHID", "SUBJECT_ID")],
                     by = "LETTERCODE", all.x = TRUE)
dat_gen <- merge(dat_gen, dat_gen_new,
                 by = c("JHUANONID", "LETTERCODE", "NIHID", "SUBJECT_ID", "APOECODE"),
                 all = TRUE)
dat_gen$APOECODE[dat_gen$APOECODE > 10] <- dat_gen$APOECODE[dat_gen$APOECODE > 10] / 10
dat_gen$apoe4 <- NA
dat_gen$apoe4[dat_gen$APOECODE %in% c(2.2, 2.3, 3.3)] <- 0
dat_gen$apoe4[dat_gen$APOECODE %in% c(2.4, 3.4, 4.4)] <- 1
dat_apoe <- dat_gen[!duplicated(dat_gen$SUBJECT_ID), c("SUBJECT_ID", "apoe4")]


###############################################################################
## SECTION 2: Filter to JHU Phase
###############################################################################

cat("\n================================================================\n")
cat("SECTION 2: FILTER TO JHU PHASE\n")
cat("================================================================\n\n")

jhu_cutoff <- as.Date("2009-01-01")
dat_jhu <- dat_blood[dat_blood$Blood_DATE >= jhu_cutoff, ]
cat(sprintf("JHU phase (>= %s): %d visits, %d subjects\n",
            jhu_cutoff, nrow(dat_jhu), length(unique(dat_jhu$SUBJECT_ID))))
cat(sprintf("  Date range: %s to %s\n", min(dat_jhu$Blood_DATE), max(dat_jhu$Blood_DATE)))


###############################################################################
## SECTION 3: Apply Exclusion Criteria
###############################################################################

cat("\n================================================================\n")
cat("SECTION 3: APPLY EXCLUSION CRITERIA\n")
cat("================================================================\n\n")

# List A: Withdrawn subjects
listA <- read.csv(file = file.path(biocard_dir, "list_A_122021.csv"))
cat(sprintf("List A (withdrawn): %d subjects\n", nrow(listA)))

# List B: Impaired at baseline
listB <- read_excel(path = file.path(biocard_dir,
                                      "LIST_B_IMPAIRED_AT_BASELINE.09.22.2015.xlsx"))
cat(sprintf("List B (impaired at baseline): %d subjects\n", nrow(listB)))

# Subjects with onset before JHU baseline
# Reconstruct survival data
dat_dx$AD_primary <- 0
dat_dx$AD_primary[dat_dx$PROBADIF == 1 | dat_dx$POSSADIF == 1] <- 1
dat_dx$AD_contrib <- 0
dat_dx$AD_contrib[dat_dx$PROBAD == 1 | dat_dx$POSSAD == 1] <- 1

# Source shared helper functions (FindJump, GetSurv) — single canonical copy
source(file.path(project_root, "shared_utils.R"))

dat_surv <- do.call(rbind, lapply(split(dat_dx, dat_dx$SUBJECT_ID), GetSurv))
dat_surv <- merge(dat_surv, dat_demo[, c("SUBJECT_ID", "LETTERCODE", "SEX", "EDUC")],
                  by = c("SUBJECT_ID", "LETTERCODE"), all.x = TRUE)
dat_surv$Sex_F <- as.numeric(dat_surv$SEX == 2)
dat_surv$EDUC_z <- as.numeric(scale(dat_surv$EDUC))
cat(sprintf("Survival data: %d subjects, %d events\n",
            nrow(dat_surv), sum(dat_surv$d == 1)))

# For plasma: JHU baseline = first JHU blood draw
# Onset before JHU baseline means symptom onset before first JHU plasma visit
jhu_first <- dat_jhu %>%
      group_by(SUBJECT_ID) %>%
      summarise(jhu_baseline_age = min(blood_age), .groups = "drop")

surv_jhu <- merge(dat_surv, jhu_first, by = "SUBJECT_ID")
# Exclude subjects whose onset was before their first JHU blood draw
ID_onset_before_jhu <- surv_jhu$SUBJECT_ID[surv_jhu$d == 1 &
                                                 surv_jhu$onset.age < surv_jhu$jhu_baseline_age]
cat(sprintf("Onset before JHU baseline: %d subjects\n", length(ID_onset_before_jhu)))

# New BIOCARD enrollees (SUBJECT_ID >= 400) — exclude for consistency with
# CSF SILA and data_extraction.R
list_new <- dat_demo$SUBJECT_ID[which(dat_demo$SUBJECT_ID >= 400)]
cat(sprintf("New enrollees (ID >= 400): %d subjects\n", length(list_new)))

# Combine exclusions
ID_exclude <- unique(c(
      listA$ID,
      listB$STUDY_ID,
      ID_onset_before_jhu,
      list_new
))

jhu_ids <- unique(dat_jhu$SUBJECT_ID)
n_excl_A <- length(intersect(jhu_ids, listA$ID))
n_excl_B <- length(intersect(jhu_ids, listB$STUDY_ID))
n_excl_onset <- length(ID_onset_before_jhu)
cat(sprintf("\nExclusions among JHU subjects:\n"))
cat(sprintf("  List A: %d\n", n_excl_A))
cat(sprintf("  List B: %d\n", n_excl_B))
cat(sprintf("  Onset before JHU baseline: %d\n", n_excl_onset))
n_excl_new <- length(intersect(jhu_ids, list_new))
cat(sprintf("  New enrollees (ID >= 400): %d\n", n_excl_new))

dat_jhu <- dat_jhu[!(dat_jhu$SUBJECT_ID %in% ID_exclude), ]
cat(sprintf("\nAfter exclusions: %d visits, %d subjects\n",
            nrow(dat_jhu), length(unique(dat_jhu$SUBJECT_ID))))


###############################################################################
## SECTION 4: Prepare SILA Input — PTAU181
###############################################################################

cat("\n================================================================\n")
cat("SECTION 4: PREPARE silaR INPUT — PTAU181\n")
cat("================================================================\n\n")

# Filter to valid PTAU181 measurements
ptau_valid <- dat_jhu[!is.na(dat_jhu$PTAU181) & !is.na(dat_jhu$blood_age), ]

# Subjects with >=2 JHU visits
n_per_subj <- tapply(ptau_valid$SUBJECT_ID, ptau_valid$SUBJECT_ID, length)
ids_ge2 <- as.integer(names(n_per_subj[n_per_subj >= 2]))
ptau_multi <- ptau_valid[ptau_valid$SUBJECT_ID %in% ids_ge2, ]

cat(sprintf("Subjects with valid PTAU181: %d total, %d with >=2 JHU visits\n",
            length(unique(ptau_valid$SUBJECT_ID)), length(ids_ge2)))
cat(sprintf("Total measurements (>=2 visits): %d\n", nrow(ptau_multi)))

# Visit distribution
time_spans <- do.call(rbind, lapply(split(ptau_multi, ptau_multi$SUBJECT_ID), function(x) {
      data.frame(SUBJECT_ID = x$SUBJECT_ID[1],
                 n_visits = nrow(x),
                 span_yr = max(x$blood_age) - min(x$blood_age))
}))
cat(sprintf("Visit counts: %d with 2, %d with 3, %d with 4, %d with 5+\n",
            sum(time_spans$n_visits == 2),
            sum(time_spans$n_visits == 3),
            sum(time_spans$n_visits == 4),
            sum(time_spans$n_visits >= 5)))
cat(sprintf("Time span: mean = %.1f yr, median = %.1f yr, range = [%.1f, %.1f]\n",
            mean(time_spans$span_yr), median(time_spans$span_yr),
            min(time_spans$span_yr), max(time_spans$span_yr)))

# Create silaR input
sila_df <- tibble(
      subid = as.numeric(ptau_multi$SUBJECT_ID),
      age   = ptau_multi$blood_age,
      val   = ptau_multi$PTAU181
)
# Remove duplicate ages within subject
sila_df <- sila_df[!duplicated(sila_df[, c("subid", "age")]), ]

cat(sprintf("\nsilaR input: %d observations from %d subjects\n",
            nrow(sila_df), length(unique(sila_df$subid))))
cat(sprintf("  PTAU181: mean = %.3f, SD = %.3f, range = [%.3f, %.3f]\n",
            mean(sila_df$val), sd(sila_df$val),
            min(sila_df$val), max(sila_df$val)))
cat(sprintf("  Age: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(sila_df$age), sd(sila_df$age),
            min(sila_df$age), max(sila_df$age)))


###############################################################################
## SECTION 5: Train SILA Model
###############################################################################

cat("\n================================================================\n")
cat("SECTION 5: TRAIN SILA MODEL\n")
cat("================================================================\n\n")

cutoff_ptau <- 0.8  # Plasma p-tau181 positivity threshold (pg/mL)

cat(sprintf("Training SILA model (dt=0.25, val0=%.1f, maxi=200)...\n", cutoff_ptau))

res_sila <- tryCatch(
      sila(sila_df, dt = 0.25, val0 = cutoff_ptau, maxi = 200),
      error = function(e) {
            cat(sprintf("ERROR in sila(): %s\n", e$message))
            cat("Trying with maxi=500...\n")
            tryCatch(
                  sila(sila_df, dt = 0.25, val0 = cutoff_ptau, maxi = 500),
                  error = function(e2) {
                        cat(sprintf("ERROR with maxi=500: %s\n", e2$message))
                        NULL
                  }
            )
      }
)

if (is.null(res_sila)) {
      cat("\nSILA FAILED TO CONVERGE.\n")
      stop("SILA failed for plasma PTAU181. Check data or adjust parameters.")
}

cat(sprintf("  tsila: %d rows (trajectory curve points)\n", nrow(res_sila$tsila)))
cat(sprintf("  adtime range: [%.2f, %.2f] years from threshold\n",
            min(res_sila$tsila$adtime), max(res_sila$tsila$adtime)))
cat(sprintf("  val range: [%.3f, %.3f]\n",
            min(res_sila$tsila$val), max(res_sila$tsila$val)))
cat(sprintf("  nsubs range: [%d, %d] (subjects contributing per value)\n",
            min(res_sila$tsila$nsubs), max(res_sila$tsila$nsubs)))


###############################################################################
## SECTION 6: Estimate Individual EAOA
###############################################################################

cat("\n================================================================\n")
cat("SECTION 6: ESTIMATE INDIVIDUAL EAOA\n")
cat("================================================================\n\n")

# Build estimation input from ALL valid subjects (>= 1 measurement)
sila_df_all <- tibble(
      subid = as.numeric(ptau_valid$SUBJECT_ID),
      age   = ptau_valid$blood_age,
      val   = ptau_valid$PTAU181
)
sila_df_all <- sila_df_all[!duplicated(sila_df_all[, c("subid", "age")]), ]
cat(sprintf("Estimating for ALL valid subjects: %d obs from %d subjects (trained on %d)\n",
            nrow(sila_df_all), length(unique(sila_df_all$subid)),
            length(unique(sila_df$subid))))

cat("Running sila_estimate(align_event='last')...\n")
resfit <- sila_estimate(res_sila$tsila, sila_df_all, align_event = "last")
cat(sprintf("  sila_estimate output: %d rows\n", nrow(resfit)))

# Extract one row per subject: last observation
resfit_last <- do.call(rbind, lapply(split(resfit, resfit$subid), function(x) {
      x <- x[x$age == max(x$age), ]
      x[1, ]
}))
cat(sprintf("  Unique subjects: %d\n", nrow(resfit_last)))

# Identify onset age column
onset_col <- NULL
if ("estage0" %in% names(resfit_last)) {
      onset_col <- "estage0"
} else if ("estaget0" %in% names(resfit_last)) {
      onset_col <- "estaget0"
} else {
      age_cols <- grep("est.*age|age.*est", names(resfit_last), value = TRUE)
      if (length(age_cols) > 0) onset_col <- age_cols[1]
}

if (!is.null(onset_col)) {
      resfit_last$EAOA_plasma <- resfit_last[[onset_col]]
      n_valid <- sum(!is.na(resfit_last$EAOA_plasma))
      cat(sprintf("  Valid EAOA (col=%s): %d (%.1f%%)\n", onset_col,
                  n_valid, n_valid / nrow(resfit_last) * 100))

      eaoa_vals <- resfit_last$EAOA_plasma[!is.na(resfit_last$EAOA_plasma)]
      if (length(eaoa_vals) > 0) {
            cat(sprintf("  EAOA distribution: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
                        mean(eaoa_vals), sd(eaoa_vals),
                        min(eaoa_vals), max(eaoa_vals)))
      }
      if ("estpos" %in% names(resfit_last)) {
            n_pos <- sum(resfit_last$estpos == TRUE, na.rm = TRUE)
            cat(sprintf("  SILA-estimated positive: %d (%.1f%%)\n",
                        n_pos, n_pos / nrow(resfit_last) * 100))
      }
} else {
      cat("WARNING: Could not identify onset age column.\n")
      cat("  Available columns:", paste(names(resfit_last), collapse = ", "), "\n")
}

resfit_last$SUBJECT_ID <- as.integer(resfit_last$subid)


###############################################################################
## SECTION 7: Merge with Survival Data and Covariates
###############################################################################

cat("\n================================================================\n")
cat("SECTION 7: MERGE WITH SURVIVAL DATA\n")
cat("================================================================\n\n")

# Merge EAOA with survival + covariates
plasma_analysis <- merge(
      resfit_last[, c("SUBJECT_ID", "EAOA_plasma", "estpos")],
      dat_surv[, c("SUBJECT_ID", "LETTERCODE", "baseline.age", "onset.age",
                    "d", "censor.age", "Sex_F", "EDUC", "EDUC_z")],
      by = "SUBJECT_ID", all.x = TRUE
)
plasma_analysis <- merge(plasma_analysis, dat_apoe, by = "SUBJECT_ID", all.x = TRUE)

# Add JHU baseline age
plasma_analysis <- merge(plasma_analysis, jhu_first, by = "SUBJECT_ID", all.x = TRUE)

cat(sprintf("Merged analysis dataset: %d subjects\n", nrow(plasma_analysis)))
cat(sprintf("  Events: %d\n", sum(plasma_analysis$d == 1, na.rm = TRUE)))
cat(sprintf("  Censored: %d\n", sum(plasma_analysis$d == 0, na.rm = TRUE)))
cat(sprintf("  Missing survival data: %d\n", sum(is.na(plasma_analysis$d))))

# Degeneracy check
if (!is.null(onset_col)) {
      pos_at_entry <- !is.na(plasma_analysis$EAOA_plasma) &
            plasma_analysis$EAOA_plasma <= plasma_analysis$jhu_baseline_age
      cat(sprintf("\n  Positive at JHU entry: %d (%.1f%%)\n",
                  sum(pos_at_entry, na.rm = TRUE),
                  mean(pos_at_entry, na.rm = TRUE) * 100))
}


###############################################################################
## SECTION 8: Save Intermediate Results
###############################################################################

cat("\n================================================================\n")
cat("SECTION 8: SAVE INTERMEDIATE\n")
cat("================================================================\n\n")

# Save SILA objects and merged analysis data
res_sila_plasma <- res_sila
resfit_plasma <- resfit
resfit_last_plasma <- resfit_last
sila_df_plasma <- sila_df

save(res_sila_plasma, resfit_plasma, resfit_last_plasma, sila_df_plasma,
     plasma_analysis, dat_surv,
     file = file.path(data_dir, "BIOCARD_plasma_SILA_intermediate.rda"))
cat("Saved: data/BIOCARD_plasma_SILA_intermediate.rda\n")
cat("  Contains: SILA model, per-subject estimates, merged analysis data\n")

# Final summary
cat("\n================================================================\n")
cat("SUMMARY\n")
cat("================================================================\n\n")

cat("--- Plasma p-tau181 (NTK/Simoa, JHU phase, val0=0.8) ---\n")
cat(sprintf("Subjects in SILA: %d\n", nrow(resfit_last)))
if (!is.null(onset_col)) {
      cat(sprintf("Valid EAOA: %d\n", sum(!is.na(resfit_last$EAOA_plasma))))
}
if ("estpos" %in% names(resfit_last)) {
      cat(sprintf("Estimated positive: %d (%.1f%%)\n",
                  sum(resfit_last$estpos == TRUE, na.rm = TRUE),
                  mean(resfit_last$estpos == TRUE, na.rm = TRUE) * 100))
}
cat(sprintf("Analysis-ready (with survival data): %d subjects, %d events\n",
            sum(!is.na(plasma_analysis$d)), sum(plasma_analysis$d == 1, na.rm = TRUE)))

# Provenance metadata
plasma_sila_provenance <- data.frame(
      script = "BIOCARD_plasma_sila.R",
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version = paste0(R.version$major, ".", R.version$minor),
      stringsAsFactors = FALSE
)
write.csv(plasma_sila_provenance,
          file = file.path(data_dir, "BIOCARD_plasma_SILA_provenance.csv"),
          row.names = FALSE)

cat("\n=== BIOCARD Plasma SILA complete ===\n")
