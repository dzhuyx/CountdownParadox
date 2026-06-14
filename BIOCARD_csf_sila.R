###############################################################################
## BIOCARD CSF SILA Trajectory Fitting — AB42/AB40 & p-tau181
##
## Purpose: Run the silaR package on BIOCARD CSF longitudinal data to estimate
##   trajectory-based age at onset (EAOA) for two biomarkers:
##   1. AB42/AB40 ratio (negated so increasing = abnormal; val0 = 0.085)
##   2. p-tau181 (already increasing; val0 = 35)
##
## Uses SILA for cross-cohort consistency (ADNI PET and ADNI plasma also
## use SILA via the silaR package).
##
## Input:  Same raw CSF file as BIOCARD_data_extraction.R
## Output: data/BIOCARD_CSF_SILA_intermediate.rda
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
## SECTION 1: Load Data (replicating BIOCARD_data_extraction.R Parts 1-3)
###############################################################################

cat("================================================================\n")
cat("SECTION 1: LOAD DATA AND APPLY EXCLUSIONS\n")
cat("================================================================\n\n")

# --- 1a: Demographics & Diagnosis (for DOB, survival) --- #
dat_demo <- read_excel(path = file.path(biocard_dir,
                                         "BIOCARD_Demographics_2024.08.07.xlsx"))

dat_dx <- read_excel(path = file.path(biocard_dir,
                                       "BIOCARD_DiagnosisData_2024.09.08.xlsx"))
dat_dx$DIAGDATE <- do.call("c", lapply(dat_dx$DIAGDATE, function(x) {
      as.Date(strsplit(as.character(x), " UTC")[[1]], "%Y-%m-%d")
}))
dat_dx$DOB <- do.call("c", lapply(dat_dx$DOB, function(x) {
      as.Date(strsplit(as.character(x), " UTC")[[1]], "%Y-%m-%d")
}))
dat_dob <- dat_dx[!duplicated(dat_dx$SUBJECT_ID), c("SUBJECT_ID", "DOB")]

cat(sprintf("Demographics: %d subjects\n", nrow(dat_demo)))

# --- 1b: CSF Lumipulse biomarkers --- #
dat_csf <- read_excel(path = file.path(biocard_dir,
                                        "BIOCARD_CSF_Lumipulse_NFL_GFAP_Data_2025.02.24.xlsx"))
dat_csf$CSF_DATE <- do.call("c", lapply(dat_csf$CSF_DATE, function(x) {
      as.Date(strsplit(as.character(x), " UTC")[[1]], "%Y-%m-%d")
}))
cat(sprintf("CSF Lumipulse: %d rows, %d subjects\n",
            nrow(dat_csf), length(unique(dat_csf$SUBJECT_ID))))

# --- 1c: Compute CSF_age --- #
dat_csf <- merge(dat_csf, dat_dob, by = "SUBJECT_ID")
dat_csf$CSF_age <- as.numeric(dat_csf$CSF_DATE - dat_csf$DOB) / 365.25


###############################################################################
## SECTION 2: Apply Exclusion Criteria (same as BIOCARD_data_extraction.R)
###############################################################################

cat("\n================================================================\n")
cat("SECTION 2: APPLY EXCLUSION CRITERIA\n")
cat("================================================================\n\n")

# List A: Withdrawn subjects
listA <- read.csv(file = file.path(biocard_dir, "list_A_122021.csv"))
cat(sprintf("List A (withdrawn): %d subjects\n", nrow(listA)))

# List B: Impaired at baseline
listB <- read_excel(path = file.path(biocard_dir,
                                      "LIST_B_IMPAIRED_AT_BASELINE.09.22.2015.xlsx"))
cat(sprintf("List B (impaired at baseline): %d subjects\n", nrow(listB)))

# New BIOCARD enrollees (SUBJECT_ID >= 400)
list_new <- dat_demo$SUBJECT_ID[which(dat_demo$SUBJECT_ID >= 400)]
cat(sprintf("New enrollees (ID >= 400): %d subjects\n", length(list_new)))

# Combine all exclusion IDs
# The onset-before-CSF-baseline exclusion is intentionally not applied here.
# SILA estimates biomarker trajectory crossing ages independently of MCI onset.
# Subjects whose onset predated their first CSF measurement still contribute
# valid biomarker trajectory data. The onset-based filter was unjustified and
# reduced the SILA fitting population, biasing the population trajectory curve.
ID_exclude <- unique(c(
      listA$ID,
      listB$STUDY_ID,
      list_new
))
# Note: We do NOT exclude single-CSF-measure subjects at this stage.
# SILA fits a population trajectory and can use all subjects.
# However, sila_estimate() needs multiple observations for reliable
# individual estimates, so subjects with 1 visit will get an estimate
# but it may be less reliable.

# Apply exclusions
dat_csf_orig <- dat_csf  # preserve for accounting
dat_csf <- dat_csf[dat_csf$SUBJECT_ID %in%
                         setdiff(unique(dat_csf$SUBJECT_ID), ID_exclude), ]
cat(sprintf("Included subjects after exclusions: %d\n",
            length(unique(dat_csf$SUBJECT_ID))))


###############################################################################
## SECTION 3: Prepare SILA Input — AB42/AB40
###############################################################################

cat("\n================================================================\n")
cat("SECTION 3: PREPARE silaR INPUT — AB42/AB40\n")
cat("================================================================\n\n")

# Negate AB42/AB40 so that INCREASING values = more abnormal
# Original AB42/AB40 decreases with pathology; EABA-derived cutoff = 0.085
# After negation: val0 = -0.085 (on negated scale)
dat_csf$AB42AB40_neg <- -dat_csf$AB42AB40

# Filter to valid AB42/AB40 measurements
ab_valid <- dat_csf[!is.na(dat_csf$AB42AB40_neg) & !is.na(dat_csf$CSF_age), ]

# Identify subjects with >=2 visits for AB42/AB40
n_per_ab <- tapply(ab_valid$SUBJECT_ID, ab_valid$SUBJECT_ID, length)
ids_ab_ge2 <- as.integer(names(n_per_ab[n_per_ab >= 2]))
ab_multi <- ab_valid[ab_valid$SUBJECT_ID %in% ids_ab_ge2, ]

cat(sprintf("Subjects with valid AB42/AB40: %d total, %d with >=2 visits\n",
            length(unique(ab_valid$SUBJECT_ID)), length(ids_ab_ge2)))
cat(sprintf("Total measurements (>=2 visits): %d\n", nrow(ab_multi)))

# Visit counts
time_spans_ab <- do.call(rbind, lapply(split(ab_multi, ab_multi$SUBJECT_ID), function(x) {
      data.frame(SUBJECT_ID = x$SUBJECT_ID[1],
                 n_visits = nrow(x),
                 span_yr = max(x$CSF_age) - min(x$CSF_age))
}))
cat(sprintf("Visit counts: %d with 2, %d with 3, %d with 4+\n",
            sum(time_spans_ab$n_visits == 2),
            sum(time_spans_ab$n_visits == 3),
            sum(time_spans_ab$n_visits >= 4)))
cat(sprintf("Time span: mean = %.2f yr, median = %.2f yr, range = [%.2f, %.2f]\n",
            mean(time_spans_ab$span_yr), median(time_spans_ab$span_yr),
            min(time_spans_ab$span_yr), max(time_spans_ab$span_yr)))

# Create silaR input
sila_df_ab <- tibble(
      subid = as.numeric(ab_multi$SUBJECT_ID),
      age   = ab_multi$CSF_age,
      val   = ab_multi$AB42AB40_neg
)
# Remove duplicate ages within subject
sila_df_ab <- sila_df_ab[!duplicated(sila_df_ab[, c("subid", "age")]), ]

cat(sprintf("\nsilaR input: %d observations from %d subjects\n",
            nrow(sila_df_ab), length(unique(sila_df_ab$subid))))
cat(sprintf("  val (negated AB42/AB40): mean = %.4f, SD = %.4f, range = [%.4f, %.4f]\n",
            mean(sila_df_ab$val), sd(sila_df_ab$val),
            min(sila_df_ab$val), max(sila_df_ab$val)))
cat(sprintf("  Age: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(sila_df_ab$age), sd(sila_df_ab$age),
            min(sila_df_ab$age), max(sila_df_ab$age)))


###############################################################################
## SECTION 4: Train SILA — AB42/AB40
###############################################################################

cat("\n================================================================\n")
cat("SECTION 4: TRAIN SILA — AB42/AB40\n")
cat("================================================================\n\n")

cutoff_ab <- -0.085  # val0: threshold on negated scale (AB42/AB40 < 0.085 = positive)

cat(sprintf("Training SILA model (dt=0.25, val0=%.3f, maxi=200)...\n", cutoff_ab))

res_sila_ab <- tryCatch(
      sila(sila_df_ab, dt = 0.25, val0 = cutoff_ab, maxi = 200),
      error = function(e) {
            cat(sprintf("ERROR in sila(): %s\n", e$message))
            cat("Trying with maxi=500...\n")
            tryCatch(
                  sila(sila_df_ab, dt = 0.25, val0 = cutoff_ab, maxi = 500),
                  error = function(e2) {
                        cat(sprintf("ERROR with maxi=500: %s\n", e2$message))
                        NULL
                  }
            )
      }
)

if (is.null(res_sila_ab)) {
      cat("\nAB42/AB40 SILA FAILED TO CONVERGE.\n")
      stop("SILA failed for AB42/AB40. Check data or adjust parameters.")
}

cat(sprintf("  tsila: %d rows (trajectory curve points)\n", nrow(res_sila_ab$tsila)))
cat(sprintf("  adtime range: [%.2f, %.2f] years from threshold\n",
            min(res_sila_ab$tsila$adtime), max(res_sila_ab$tsila$adtime)))
cat(sprintf("  val range: [%.4f, %.4f]\n",
            min(res_sila_ab$tsila$val), max(res_sila_ab$tsila$val)))
cat(sprintf("  nsubs range: [%d, %d] (subjects contributing per value)\n",
            min(res_sila_ab$tsila$nsubs), max(res_sila_ab$tsila$nsubs)))


###############################################################################
## SECTION 5: Estimate Individual EAOA — AB42/AB40
###############################################################################

cat("\n================================================================\n")
cat("SECTION 5: ESTIMATE INDIVIDUAL EAOA — AB42/AB40\n")
cat("================================================================\n\n")

# Build estimation input from ALL valid subjects (>= 1 measurement),
# not just the >= 2 training set. SILA was trained on >= 2, but
# sila_estimate() only uses each subject's last observation for alignment.
sila_df_ab_all <- tibble(
      subid = as.numeric(ab_valid$SUBJECT_ID),
      age   = ab_valid$CSF_age,
      val   = ab_valid$AB42AB40_neg
)
sila_df_ab_all <- sila_df_ab_all[!duplicated(sila_df_ab_all[, c("subid", "age")]), ]
cat(sprintf("Estimating for ALL valid subjects: %d obs from %d subjects (trained on %d)\n",
            nrow(sila_df_ab_all), length(unique(sila_df_ab_all$subid)),
            length(unique(sila_df_ab$subid))))

cat("Running sila_estimate(align_event='last')...\n")
resfit_ab <- sila_estimate(res_sila_ab$tsila, sila_df_ab_all, align_event = "last")
cat(sprintf("  sila_estimate output: %d rows\n", nrow(resfit_ab)))

# Extract one row per subject: use the last observation's estimate
resfit_last_ab <- do.call(rbind, lapply(split(resfit_ab, resfit_ab$subid), function(x) {
      x <- x[x$age == max(x$age), ]
      x[1, ]
}))
cat(sprintf("  Unique subjects: %d\n", nrow(resfit_last_ab)))

# Identify onset age column
onset_col_ab <- NULL
if ("estage0" %in% names(resfit_last_ab)) {
      onset_col_ab <- "estage0"
} else if ("estaget0" %in% names(resfit_last_ab)) {
      onset_col_ab <- "estaget0"
} else {
      age_cols <- grep("est.*age|age.*est", names(resfit_last_ab), value = TRUE)
      if (length(age_cols) > 0) onset_col_ab <- age_cols[1]
}

if (!is.null(onset_col_ab)) {
      resfit_last_ab$EAOA_AB <- resfit_last_ab[[onset_col_ab]]
      n_valid_ab <- sum(!is.na(resfit_last_ab$EAOA_AB))
      cat(sprintf("  Valid EAOA (col=%s): %d (%.1f%%)\n", onset_col_ab,
                  n_valid_ab, n_valid_ab / nrow(resfit_last_ab) * 100))

      eaoa_vals_ab <- resfit_last_ab$EAOA_AB[!is.na(resfit_last_ab$EAOA_AB)]
      if (length(eaoa_vals_ab) > 0) {
            cat(sprintf("  EAOA distribution: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
                        mean(eaoa_vals_ab), sd(eaoa_vals_ab),
                        min(eaoa_vals_ab), max(eaoa_vals_ab)))
      }
      if ("estpos" %in% names(resfit_last_ab)) {
            n_pos_ab <- sum(resfit_last_ab$estpos == TRUE, na.rm = TRUE)
            cat(sprintf("  SILA-estimated positive: %d (%.1f%%)\n",
                        n_pos_ab, n_pos_ab / nrow(resfit_last_ab) * 100))
      }
} else {
      cat("WARNING: Could not identify onset age column in sila_estimate output.\n")
      cat("  Available columns:", paste(names(resfit_last_ab), collapse = ", "), "\n")
}

resfit_last_ab$SUBJECT_ID <- as.integer(resfit_last_ab$subid)


###############################################################################
## SECTION 6: Prepare SILA Input — p-tau181
###############################################################################

cat("\n================================================================\n")
cat("SECTION 6: PREPARE silaR INPUT — p-tau181\n")
cat("================================================================\n\n")

# p-tau181 is already increasing — no negation needed
ptau_valid <- dat_csf[!is.na(dat_csf$PTAU181) & !is.na(dat_csf$CSF_age), ]

n_per_ptau <- tapply(ptau_valid$SUBJECT_ID, ptau_valid$SUBJECT_ID, length)
ids_ptau_ge2 <- as.integer(names(n_per_ptau[n_per_ptau >= 2]))
ptau_multi <- ptau_valid[ptau_valid$SUBJECT_ID %in% ids_ptau_ge2, ]

cat(sprintf("Subjects with valid p-tau181: %d total, %d with >=2 visits\n",
            length(unique(ptau_valid$SUBJECT_ID)), length(ids_ptau_ge2)))
cat(sprintf("Total measurements (>=2 visits): %d\n", nrow(ptau_multi)))

time_spans_ptau <- do.call(rbind, lapply(split(ptau_multi, ptau_multi$SUBJECT_ID), function(x) {
      data.frame(SUBJECT_ID = x$SUBJECT_ID[1],
                 n_visits = nrow(x),
                 span_yr = max(x$CSF_age) - min(x$CSF_age))
}))
cat(sprintf("Visit counts: %d with 2, %d with 3, %d with 4+\n",
            sum(time_spans_ptau$n_visits == 2),
            sum(time_spans_ptau$n_visits == 3),
            sum(time_spans_ptau$n_visits >= 4)))
cat(sprintf("Time span: mean = %.2f yr, median = %.2f yr, range = [%.2f, %.2f]\n",
            mean(time_spans_ptau$span_yr), median(time_spans_ptau$span_yr),
            min(time_spans_ptau$span_yr), max(time_spans_ptau$span_yr)))

sila_df_ptau <- tibble(
      subid = as.numeric(ptau_multi$SUBJECT_ID),
      age   = ptau_multi$CSF_age,
      val   = ptau_multi$PTAU181
)
sila_df_ptau <- sila_df_ptau[!duplicated(sila_df_ptau[, c("subid", "age")]), ]

cat(sprintf("\nsilaR input: %d observations from %d subjects\n",
            nrow(sila_df_ptau), length(unique(sila_df_ptau$subid))))
cat(sprintf("  val (PTAU181): mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(sila_df_ptau$val), sd(sila_df_ptau$val),
            min(sila_df_ptau$val), max(sila_df_ptau$val)))


###############################################################################
## SECTION 7: Train SILA — p-tau181
###############################################################################

cat("\n================================================================\n")
cat("SECTION 7: TRAIN SILA — p-tau181\n")
cat("================================================================\n\n")

cutoff_ptau <- 35  # p-tau181 positivity threshold

cat(sprintf("Training SILA model (dt=0.25, val0=%.0f, maxi=200)...\n", cutoff_ptau))

res_sila_ptau <- tryCatch(
      sila(sila_df_ptau, dt = 0.25, val0 = cutoff_ptau, maxi = 200),
      error = function(e) {
            cat(sprintf("ERROR in sila(): %s\n", e$message))
            cat("Trying with maxi=500...\n")
            tryCatch(
                  sila(sila_df_ptau, dt = 0.25, val0 = cutoff_ptau, maxi = 500),
                  error = function(e2) {
                        cat(sprintf("ERROR with maxi=500: %s\n", e2$message))
                        NULL
                  }
            )
      }
)

if (is.null(res_sila_ptau)) {
      cat("\np-tau181 SILA FAILED TO CONVERGE.\n")
      stop("SILA failed for p-tau181. Check data or adjust parameters.")
}

cat(sprintf("  tsila: %d rows (trajectory curve points)\n", nrow(res_sila_ptau$tsila)))
cat(sprintf("  adtime range: [%.2f, %.2f] years from threshold\n",
            min(res_sila_ptau$tsila$adtime), max(res_sila_ptau$tsila$adtime)))
cat(sprintf("  val range: [%.1f, %.1f]\n",
            min(res_sila_ptau$tsila$val), max(res_sila_ptau$tsila$val)))
cat(sprintf("  nsubs range: [%d, %d] (subjects contributing per value)\n",
            min(res_sila_ptau$tsila$nsubs), max(res_sila_ptau$tsila$nsubs)))


###############################################################################
## SECTION 8: Estimate Individual EAOA — p-tau181
###############################################################################

cat("\n================================================================\n")
cat("SECTION 8: ESTIMATE INDIVIDUAL EAOA — p-tau181\n")
cat("================================================================\n\n")

# Build estimation input from ALL valid subjects (>= 1 measurement)
sila_df_ptau_all <- tibble(
      subid = as.numeric(ptau_valid$SUBJECT_ID),
      age   = ptau_valid$CSF_age,
      val   = ptau_valid$PTAU181
)
sila_df_ptau_all <- sila_df_ptau_all[!duplicated(sila_df_ptau_all[, c("subid", "age")]), ]
cat(sprintf("Estimating for ALL valid subjects: %d obs from %d subjects (trained on %d)\n",
            nrow(sila_df_ptau_all), length(unique(sila_df_ptau_all$subid)),
            length(unique(sila_df_ptau$subid))))

cat("Running sila_estimate(align_event='last')...\n")
resfit_ptau <- sila_estimate(res_sila_ptau$tsila, sila_df_ptau_all, align_event = "last")
cat(sprintf("  sila_estimate output: %d rows\n", nrow(resfit_ptau)))

resfit_last_ptau <- do.call(rbind, lapply(split(resfit_ptau, resfit_ptau$subid), function(x) {
      x <- x[x$age == max(x$age), ]
      x[1, ]
}))
cat(sprintf("  Unique subjects: %d\n", nrow(resfit_last_ptau)))

onset_col_ptau <- NULL
if ("estage0" %in% names(resfit_last_ptau)) {
      onset_col_ptau <- "estage0"
} else if ("estaget0" %in% names(resfit_last_ptau)) {
      onset_col_ptau <- "estaget0"
} else {
      age_cols <- grep("est.*age|age.*est", names(resfit_last_ptau), value = TRUE)
      if (length(age_cols) > 0) onset_col_ptau <- age_cols[1]
}

if (!is.null(onset_col_ptau)) {
      resfit_last_ptau$EAOA_PTAU <- resfit_last_ptau[[onset_col_ptau]]
      n_valid_ptau <- sum(!is.na(resfit_last_ptau$EAOA_PTAU))
      cat(sprintf("  Valid EAOA (col=%s): %d (%.1f%%)\n", onset_col_ptau,
                  n_valid_ptau, n_valid_ptau / nrow(resfit_last_ptau) * 100))

      eaoa_vals_ptau <- resfit_last_ptau$EAOA_PTAU[!is.na(resfit_last_ptau$EAOA_PTAU)]
      if (length(eaoa_vals_ptau) > 0) {
            cat(sprintf("  EAOA distribution: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
                        mean(eaoa_vals_ptau), sd(eaoa_vals_ptau),
                        min(eaoa_vals_ptau), max(eaoa_vals_ptau)))
      }
      if ("estpos" %in% names(resfit_last_ptau)) {
            n_pos_ptau <- sum(resfit_last_ptau$estpos == TRUE, na.rm = TRUE)
            cat(sprintf("  SILA-estimated positive: %d (%.1f%%)\n",
                        n_pos_ptau, n_pos_ptau / nrow(resfit_last_ptau) * 100))
      }
} else {
      cat("WARNING: Could not identify onset age column.\n")
      cat("  Available columns:", paste(names(resfit_last_ptau), collapse = ", "), "\n")
}

resfit_last_ptau$SUBJECT_ID <- as.integer(resfit_last_ptau$subid)


###############################################################################
## SECTION 9: Save SILA Intermediate Results
###############################################################################

cat("\n================================================================\n")
cat("SECTION 10: SAVE INTERMEDIATE\n")
cat("================================================================\n\n")

# Save all SILA objects for both biomarkers
save(res_sila_ab, resfit_ab, resfit_last_ab, sila_df_ab,
     res_sila_ptau, resfit_ptau, resfit_last_ptau, sila_df_ptau,
     file = file.path(data_dir, "BIOCARD_CSF_SILA_intermediate.rda"))
cat("Saved: data/BIOCARD_CSF_SILA_intermediate.rda\n")
cat("  Contains: AB42/AB40 SILA + p-tau181 SILA results\n")

# Final summary
cat("\n================================================================\n")
cat("SUMMARY\n")
cat("================================================================\n\n")

cat("--- AB42/AB40 (negated, val0=-0.085) ---\n")
cat(sprintf("Subjects in SILA: %d\n", nrow(resfit_last_ab)))
if (!is.null(onset_col_ab)) {
      cat(sprintf("Valid EAOA: %d\n", sum(!is.na(resfit_last_ab$EAOA_AB))))
}
if ("estpos" %in% names(resfit_last_ab)) {
      cat(sprintf("Estimated positive: %d\n",
                  sum(resfit_last_ab$estpos == TRUE, na.rm = TRUE)))
}

cat("\n--- p-tau181 (val0=35) ---\n")
cat(sprintf("Subjects in SILA: %d\n", nrow(resfit_last_ptau)))
if (!is.null(onset_col_ptau)) {
      cat(sprintf("Valid EAOA: %d\n", sum(!is.na(resfit_last_ptau$EAOA_PTAU))))
}
if ("estpos" %in% names(resfit_last_ptau)) {
      cat(sprintf("Estimated positive: %d\n",
                  sum(resfit_last_ptau$estpos == TRUE, na.rm = TRUE)))
}

# Provenance metadata
csf_sila_provenance <- data.frame(
      script = "BIOCARD_csf_sila.R",
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version = paste0(R.version$major, ".", R.version$minor),
      stringsAsFactors = FALSE
)
write.csv(csf_sila_provenance,
          file = file.path(data_dir, "BIOCARD_CSF_SILA_provenance.csv"),
          row.names = FALSE)

cat("\n=== BIOCARD CSF SILA complete ===\n")
