###############################################################################
## ADNI Plasma SILA Trajectory Fitting — Fujirebio pT217_F & C2N %p-tau217
##
## Purpose: Run the silaR package on plasma p-tau217 data to estimate
##   trajectory-based age at positivity (EAOA) for two biomarkers:
##   1. Fujirebio pT217_F (val0 = 0.300, ~479 subjects with ≥2 visits)
##   2. C2N %p-tau217 (val0 = 4.06, ~143 subjects with ≥2 visits)
##
## Input:  data/ADNI_Biospecimen_28Feb2026/UPENN_PLASMA_FUJIREBIO_QUANTERIX_28Feb2026.csv
##         data/C2N_Precivity_28Feb2026/C2N_PRECIVITYAD2_PLASMA_28Feb2026.csv
##         ADNI_2026_data/PTDEMOG_11Feb2026.csv (for demographics/DOB)
## Output: data/ADNI_plasma_SILA_intermediate.rda
##         (contains Fujirebio and C2N SILA results)
##
## Usage:  source() from reproducibility/ directory or via run_all.R
##
## Author: Yuxin Zhu
## Date: March 2026 (reproducibility copy April 2026)
###############################################################################

rm(list = ls())
library(survival)
library(silaR)
library(tibble)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT", "/Users/daisyzhu/Documents/Research Projects/CountdownParadox_BiomarkerPositivity/CountdownParadox_Analysis")
data_dir     <- file.path(project_root, "data")
bio_dir      <- file.path(data_dir, "ADNI_Biospecimen_28Feb2026")
c2n_dir      <- file.path(data_dir, "C2N_Precivity_28Feb2026")
data_2026    <- file.path(project_root, "ADNI_2026_data")


###############################################################################
## SECTION 1: Load Data
###############################################################################

cat("================================================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("================================================================\n\n")

# --- 1a: Demographics (for DOB -> age) --- #
demo_file <- file.path(data_2026, "PTDEMOG_11Feb2026.csv")
dat_demo <- read.csv(demo_file, na.strings = -4)
dat_demo$DOB <- as.Date(paste0("01/", dat_demo$PTDOB), "%d/%m/%Y")
dat_demo$update_stamp <- as.Date(dat_demo$update_stamp, "%Y-%m-%d")
dat_demo <- do.call(rbind, lapply(split(dat_demo, dat_demo$RID), function(x) {
      x[x$update_stamp == max(x$update_stamp, na.rm = TRUE), ][1, ]
}))
cat(sprintf("Demographics: %d subjects\n", nrow(dat_demo)))

# --- 1b: Fujirebio plasma --- #
plasma_file <- file.path(bio_dir, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_28Feb2026.csv")
plasma <- read.csv(plasma_file, stringsAsFactors = FALSE, na.strings = c("", "NA", "-4"))

# Recode ADNI negative missing codes
for (col in c("pT217_F", "AB42_F", "AB40_F", "AB42_AB40_F",
              "pT217_AB42_F", "NfL_Q", "GFAP_Q", "NfL_F", "GFAP_F")) {
      if (col %in% names(plasma)) {
            bad <- !is.na(plasma[[col]]) & plasma[[col]] < 0
            plasma[[col]][bad] <- NA
      }
}

plasma$EXAMDATE <- as.Date(plasma$EXAMDATE)
plasma <- merge(plasma, dat_demo[, c("RID", "DOB")], by = "RID", all.x = TRUE)
plasma$age_at_exam <- as.numeric(plasma$EXAMDATE - plasma$DOB) / 365.25
cat(sprintf("Fujirebio plasma: %d measurements, %d subjects\n",
            nrow(plasma), length(unique(plasma$RID))))

# Exclude 78 participants removed by ADNI (data quality notice, March 2026)
adni_exclude_rids <- unique(dat_demo$RID[grepl("^381_S_", dat_demo$PTID)])
n_excl_plasma <- length(unique(plasma$RID[plasma$RID %in% adni_exclude_rids]))
plasma <- plasma[!plasma$RID %in% adni_exclude_rids, ]
cat(sprintf("Excluded %d subjects (%d flagged 381_S) from plasma data\n",
            n_excl_plasma, length(adni_exclude_rids)))


###############################################################################
## SECTION 2: Filter to Multi-Visit Subjects
###############################################################################

cat("\n================================================================\n")
cat("SECTION 2: FILTER TO MULTI-VISIT SUBJECTS\n")
cat("================================================================\n\n")

# Filter to valid pT217_F measurements
plasma_valid <- plasma[!is.na(plasma$pT217_F) & !is.na(plasma$age_at_exam), ]

# Identify subjects with >=2 visits
n_per_subj <- tapply(plasma_valid$RID, plasma_valid$RID, length)
rids_ge2 <- as.integer(names(n_per_subj[n_per_subj >= 2]))
plasma_multi <- plasma_valid[plasma_valid$RID %in% rids_ge2, ]

cat(sprintf("Subjects with >=2 plasma visits: %d\n", length(rids_ge2)))
cat(sprintf("Total measurements: %d\n", nrow(plasma_multi)))

# Time span for multi-visit subjects
time_spans <- do.call(rbind, lapply(split(plasma_multi, plasma_multi$RID), function(x) {
      data.frame(RID = x$RID[1],
                 n_visits = nrow(x),
                 span_yr = max(x$age_at_exam) - min(x$age_at_exam))
}))
cat(sprintf("Visit counts: %d with 2, %d with 3, %d with 4+\n",
            sum(time_spans$n_visits == 2),
            sum(time_spans$n_visits == 3),
            sum(time_spans$n_visits >= 4)))
cat(sprintf("Time span: mean = %.2f yr, median = %.2f yr, range = [%.2f, %.2f]\n",
            mean(time_spans$span_yr), median(time_spans$span_yr),
            min(time_spans$span_yr), max(time_spans$span_yr)))


###############################################################################
## SECTION 3: Prepare silaR Input
###############################################################################

cat("\n================================================================\n")
cat("SECTION 3: PREPARE silaR INPUT\n")
cat("================================================================\n\n")

# Create silaR input: subid (numeric), age, val
sila_df <- tibble(
      subid = as.numeric(plasma_multi$RID),
      age   = plasma_multi$age_at_exam,
      val   = plasma_multi$pT217_F
)

# Remove duplicate ages within subject (silaR requirement)
sila_df <- sila_df[!duplicated(sila_df[, c("subid", "age")]), ]

cat(sprintf("silaR input: %d observations from %d subjects\n",
            nrow(sila_df), length(unique(sila_df$subid))))
cat(sprintf("  pT217_F: mean = %.3f, SD = %.3f, range = [%.3f, %.3f]\n",
            mean(sila_df$val), sd(sila_df$val),
            min(sila_df$val), max(sila_df$val)))
cat(sprintf("  Age: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(sila_df$age), sd(sila_df$age),
            min(sila_df$age), max(sila_df$age)))


###############################################################################
## SECTION 4: Train SILA Model
###############################################################################

cat("\n================================================================\n")
cat("SECTION 4: TRAIN SILA MODEL\n")
cat("================================================================\n\n")

cutoff_fuj <- 0.300  # Fujirebio positivity threshold

cat(sprintf("Training SILA model (dt=0.25, val0=%.3f, maxi=200)...\n", cutoff_fuj))

res_sila <- tryCatch(
      sila(sila_df, dt = 0.25, val0 = cutoff_fuj, maxi = 200),
      error = function(e) {
            cat(sprintf("ERROR in sila(): %s\n", e$message))
            cat("Trying with maxi=500...\n")
            tryCatch(
                  sila(sila_df, dt = 0.25, val0 = cutoff_fuj, maxi = 500),
                  error = function(e2) {
                        cat(sprintf("ERROR with maxi=500: %s\n", e2$message))
                        NULL
                  }
            )
      }
)

if (is.null(res_sila)) {
      cat("\nSILA FAILED TO CONVERGE — saving empty intermediate.\n")
      resfit <- NULL
      resfit_last_plasma <- NULL

      save(res_sila, resfit, resfit_last_plasma, sila_df,
           file = file.path(data_dir, "ADNI_plasma_SILA_intermediate.rda"))
      cat("Saved: data/ADNI_plasma_SILA_intermediate.rda (empty — SILA failed)\n")
      cat("\n=== ADNI Plasma SILA complete (with failure) ===\n")
      stop("SILA failed to converge. Check data or adjust parameters.")
}

# Check convergence
cat(sprintf("  tsila: %d rows (trajectory curve points)\n", nrow(res_sila$tsila)))
cat(sprintf("  adtime range: [%.2f, %.2f] years from threshold\n",
            min(res_sila$tsila$adtime), max(res_sila$tsila$adtime)))
cat(sprintf("  val range: [%.3f, %.3f]\n",
            min(res_sila$tsila$val), max(res_sila$tsila$val)))
cat(sprintf("  nsubs range: [%d, %d] (subjects contributing per value)\n",
            min(res_sila$tsila$nsubs), max(res_sila$tsila$nsubs)))


###############################################################################
## SECTION 5: Estimate Individual Onset Ages
###############################################################################

cat("\n================================================================\n")
cat("SECTION 5: ESTIMATE INDIVIDUAL ONSET AGES (EAOA)\n")
cat("================================================================\n\n")

# Build estimation input from ALL valid subjects (>= 1 measurement)
sila_df_all <- tibble(
      subid = as.numeric(plasma_valid$RID),
      age   = plasma_valid$age_at_exam,
      val   = plasma_valid$pT217_F
)
sila_df_all <- sila_df_all[!duplicated(sila_df_all[, c("subid", "age")]), ]
cat(sprintf("Estimating for ALL valid subjects: %d obs from %d subjects (trained on %d)\n",
            nrow(sila_df_all), length(unique(sila_df_all$subid)),
            length(unique(sila_df$subid))))

cat("Running sila_estimate(align_event='last')...\n")
resfit <- sila_estimate(res_sila$tsila, sila_df_all, align_event = "last")

cat(sprintf("  sila_estimate output: %d rows\n", nrow(resfit)))

# Extract one row per subject: use the last observation's estaget0
resfit_last_plasma <- do.call(rbind, lapply(split(resfit, resfit$subid), function(x) {
      x <- x[x$age == max(x$age), ]
      x[1, ]
}))

cat(sprintf("  Unique subjects: %d\n", nrow(resfit_last_plasma)))

# Check which column name SILA uses for estimated onset age
# (estaget0 in PET code, but sila_estimate() docs say estage0)
onset_col <- NULL
if ("estaget0" %in% names(resfit_last_plasma)) {
      onset_col <- "estaget0"
} else if ("estage0" %in% names(resfit_last_plasma)) {
      onset_col <- "estage0"
} else {
      cat("WARNING: Neither estaget0 nor estage0 found in sila_estimate output!\n")
      cat("  Available columns: ", paste(names(resfit_last_plasma), collapse = ", "), "\n")
      # Try to find any column with "est" and "age" in the name
      age_cols <- grep("est.*age|age.*est", names(resfit_last_plasma), value = TRUE)
      if (length(age_cols) > 0) {
            onset_col <- age_cols[1]
            cat(sprintf("  Using fallback column: %s\n", onset_col))
      }
}

if (!is.null(onset_col)) {
      n_valid_eaoa <- sum(!is.na(resfit_last_plasma[[onset_col]]))
      cat(sprintf("  Valid EAOA (non-NA, col=%s): %d (%.1f%%)\n", onset_col,
                  n_valid_eaoa, n_valid_eaoa / nrow(resfit_last_plasma) * 100))

      eaoa_vals <- resfit_last_plasma[[onset_col]][!is.na(resfit_last_plasma[[onset_col]])]
      if (length(eaoa_vals) > 0) {
            cat(sprintf("  EAOA distribution: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
                        mean(eaoa_vals), sd(eaoa_vals), min(eaoa_vals), max(eaoa_vals)))
      }

      # How many are positive (estpos)?
      if ("estpos" %in% names(resfit_last_plasma)) {
            n_pos <- sum(resfit_last_plasma$estpos == TRUE, na.rm = TRUE)
            cat(sprintf("  SILA-estimated positive (estpos = TRUE): %d (%.1f%%)\n",
                        n_pos, n_pos / nrow(resfit_last_plasma) * 100))
      }

      # Standardize onset column name for downstream use
      resfit_last_plasma$EAOA_plasma <- resfit_last_plasma[[onset_col]]
}


###############################################################################
## SECTION 6: Compare with Threshold-Based Age at First Positive (Fujirebio)
###############################################################################

cat("\n================================================================\n")
cat("SECTION 6: SILA vs THRESHOLD COMPARISON (Fujirebio)\n")
cat("================================================================\n\n")

# Compute first-positive-scan age
first_pos_plasma <- do.call(rbind, lapply(split(plasma_valid, plasma_valid$RID), function(x) {
      x <- x[order(x$EXAMDATE), ]
      pos_idx <- which(x$pT217_F >= cutoff_fuj)
      if (length(pos_idx) == 0) return(NULL)
      data.frame(RID = x$RID[1], first_pos_age = x$age_at_exam[pos_idx[1]])
}))

resfit_last_plasma$RID <- as.integer(resfit_last_plasma$subid)
resfit_last_plasma <- merge(resfit_last_plasma, first_pos_plasma, by = "RID", all.x = TRUE)

if (!is.null(onset_col)) {
      both_valid <- !is.na(resfit_last_plasma$EAOA_plasma) &
            !is.na(resfit_last_plasma$first_pos_age)
      if (sum(both_valid) > 5) {
            cor_val <- cor(resfit_last_plasma$EAOA_plasma[both_valid],
                           resfit_last_plasma$first_pos_age[both_valid])
            mean_diff <- mean(resfit_last_plasma$EAOA_plasma[both_valid] -
                                    resfit_last_plasma$first_pos_age[both_valid])
            cat(sprintf("Correlation SILA EAOA vs first-positive-plasma age: r = %.3f (n = %d)\n",
                        cor_val, sum(both_valid)))
            cat(sprintf("Mean EAOA - first_pos: %.1f years\n", mean_diff))
      } else {
            cat(sprintf("Too few subjects with both SILA EAOA and first-positive (n = %d)\n",
                        sum(both_valid)))
      }
}

# Store Fujirebio results
res_sila_fuj <- res_sila
resfit_fuj <- resfit
resfit_last_fuj <- resfit_last_plasma
sila_df_fuj <- sila_df


###############################################################################
## SECTION 7: C2N %p-tau217 — Load and Prepare
###############################################################################

cat("\n================================================================\n")
cat("SECTION 7: C2N %p-tau217 — LOAD AND PREPARE\n")
cat("================================================================\n\n")

# --- 7a: Load C2N plasma data --- #
c2n_file <- file.path(c2n_dir, "C2N_PRECIVITYAD2_PLASMA_28Feb2026.csv")
c2n <- read.csv(c2n_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))

# Recode ADNI negative missing codes
for (col in c("pT217_C2N", "npT217_C2N", "AB42_C2N", "AB40_C2N",
              "AB42_AB40_C2N", "pT217_npT217_C2N", "APS2_C2N")) {
      if (col %in% names(c2n)) {
            bad <- !is.na(c2n[[col]]) & c2n[[col]] < 0
            c2n[[col]][bad] <- NA
      }
}

c2n$EXAMDATE <- as.Date(c2n$EXAMDATE)
c2n <- merge(c2n, dat_demo[, c("RID", "DOB")], by = "RID", all.x = TRUE)
c2n$age_at_exam <- as.numeric(c2n$EXAMDATE - c2n$DOB) / 365.25
cat(sprintf("C2N plasma: %d measurements, %d subjects\n",
            nrow(c2n), length(unique(c2n$RID))))

# Exclude 381_S subjects (same exclusion as Fujirebio above)
n_excl_c2n <- length(unique(c2n$RID[c2n$RID %in% adni_exclude_rids]))
c2n <- c2n[!c2n$RID %in% adni_exclude_rids, ]
cat(sprintf("Excluded %d C2N subjects (381_S data quality)\n", n_excl_c2n))

# --- 7b: Filter to valid %p-tau217 and multi-visit --- #
c2n_valid <- c2n[!is.na(c2n$pT217_npT217_C2N) & !is.na(c2n$age_at_exam), ]

n_per_c2n <- tapply(c2n_valid$RID, c2n_valid$RID, length)
rids_c2n_ge2 <- as.integer(names(n_per_c2n[n_per_c2n >= 2]))
c2n_multi <- c2n_valid[c2n_valid$RID %in% rids_c2n_ge2, ]

cat(sprintf("C2N subjects with >=2 visits: %d\n", length(rids_c2n_ge2)))
cat(sprintf("Total measurements: %d\n", nrow(c2n_multi)))

time_spans_c2n <- do.call(rbind, lapply(split(c2n_multi, c2n_multi$RID), function(x) {
      data.frame(RID = x$RID[1],
                 n_visits = nrow(x),
                 span_yr = max(x$age_at_exam) - min(x$age_at_exam))
}))
cat(sprintf("Visit counts: %d with 2, %d with 3, %d with 4+\n",
            sum(time_spans_c2n$n_visits == 2),
            sum(time_spans_c2n$n_visits == 3),
            sum(time_spans_c2n$n_visits >= 4)))
cat(sprintf("Time span: mean = %.2f yr, median = %.2f yr, range = [%.2f, %.2f]\n",
            mean(time_spans_c2n$span_yr), median(time_spans_c2n$span_yr),
            min(time_spans_c2n$span_yr), max(time_spans_c2n$span_yr)))


###############################################################################
## SECTION 8: C2N — Prepare silaR Input and Train SILA
###############################################################################

cat("\n================================================================\n")
cat("SECTION 8: C2N SILA MODEL\n")
cat("================================================================\n\n")

sila_df_c2n <- tibble(
      subid = as.numeric(c2n_multi$RID),
      age   = c2n_multi$age_at_exam,
      val   = c2n_multi$pT217_npT217_C2N
)
sila_df_c2n <- sila_df_c2n[!duplicated(sila_df_c2n[, c("subid", "age")]), ]

cat(sprintf("C2N silaR input: %d observations from %d subjects\n",
            nrow(sila_df_c2n), length(unique(sila_df_c2n$subid))))
cat(sprintf("  pT217_npT217_C2N: mean = %.2f, SD = %.2f, range = [%.2f, %.2f]\n",
            mean(sila_df_c2n$val), sd(sila_df_c2n$val),
            min(sila_df_c2n$val), max(sila_df_c2n$val)))

cutoff_c2n <- 4.06  # C2N %p-tau217 positivity threshold (Petersen 2026)
cat(sprintf("Training C2N SILA model (dt=0.25, val0=%.2f, maxi=200)...\n", cutoff_c2n))

res_sila_c2n <- tryCatch(
      sila(sila_df_c2n, dt = 0.25, val0 = cutoff_c2n, maxi = 200),
      error = function(e) {
            cat(sprintf("ERROR in sila() for C2N: %s\n", e$message))
            cat("Trying with maxi=500...\n")
            tryCatch(
                  sila(sila_df_c2n, dt = 0.25, val0 = cutoff_c2n, maxi = 500),
                  error = function(e2) {
                        cat(sprintf("ERROR with maxi=500: %s\n", e2$message))
                        NULL
                  }
            )
      }
)

resfit_c2n <- NULL
resfit_last_c2n <- NULL

if (is.null(res_sila_c2n)) {
      cat("\nC2N SILA FAILED TO CONVERGE.\n")
      cat("NOTE: n=143 subjects is small for SILA population trajectory fitting.\n")
      cat("TIRA (LME-based) may handle sparse data better.\n")
} else {
      cat(sprintf("  tsila: %d rows\n", nrow(res_sila_c2n$tsila)))
      cat(sprintf("  adtime range: [%.2f, %.2f]\n",
                  min(res_sila_c2n$tsila$adtime), max(res_sila_c2n$tsila$adtime)))
      cat(sprintf("  val range: [%.2f, %.2f]\n",
                  min(res_sila_c2n$tsila$val), max(res_sila_c2n$tsila$val)))

      # Build estimation input from ALL valid C2N subjects (>= 1 measurement)
      sila_df_c2n_all <- tibble(
            subid = as.numeric(c2n_valid$RID),
            age   = c2n_valid$age_at_exam,
            val   = c2n_valid$pT217_npT217_C2N
      )
      sila_df_c2n_all <- sila_df_c2n_all[!duplicated(sila_df_c2n_all[, c("subid", "age")]), ]
      cat(sprintf("Estimating for ALL valid C2N subjects: %d obs from %d subjects (trained on %d)\n",
                  nrow(sila_df_c2n_all), length(unique(sila_df_c2n_all$subid)),
                  length(unique(sila_df_c2n$subid))))

      cat("Running sila_estimate(align_event='last') for C2N...\n")
      resfit_c2n <- sila_estimate(res_sila_c2n$tsila, sila_df_c2n_all, align_event = "last")
      cat(sprintf("  sila_estimate output: %d rows\n", nrow(resfit_c2n)))

      resfit_last_c2n <- do.call(rbind, lapply(split(resfit_c2n, resfit_c2n$subid), function(x) {
            x <- x[x$age == max(x$age), ]
            x[1, ]
      }))
      cat(sprintf("  Unique subjects: %d\n", nrow(resfit_last_c2n)))

      # Check onset column name
      onset_col_c2n <- NULL
      if ("estaget0" %in% names(resfit_last_c2n)) {
            onset_col_c2n <- "estaget0"
      } else if ("estage0" %in% names(resfit_last_c2n)) {
            onset_col_c2n <- "estage0"
      } else {
            age_cols <- grep("est.*age|age.*est", names(resfit_last_c2n), value = TRUE)
            if (length(age_cols) > 0) onset_col_c2n <- age_cols[1]
      }

      if (!is.null(onset_col_c2n)) {
            resfit_last_c2n$EAOA_plasma <- resfit_last_c2n[[onset_col_c2n]]
            n_valid_c2n <- sum(!is.na(resfit_last_c2n$EAOA_plasma))
            cat(sprintf("  Valid EAOA (col=%s): %d (%.1f%%)\n", onset_col_c2n,
                        n_valid_c2n, n_valid_c2n / nrow(resfit_last_c2n) * 100))
            eaoa_vals_c2n <- resfit_last_c2n$EAOA_plasma[!is.na(resfit_last_c2n$EAOA_plasma)]
            if (length(eaoa_vals_c2n) > 0) {
                  cat(sprintf("  EAOA distribution: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
                              mean(eaoa_vals_c2n), sd(eaoa_vals_c2n),
                              min(eaoa_vals_c2n), max(eaoa_vals_c2n)))
            }
            if ("estpos" %in% names(resfit_last_c2n)) {
                  n_pos_c2n <- sum(resfit_last_c2n$estpos == TRUE, na.rm = TRUE)
                  cat(sprintf("  SILA-estimated positive: %d (%.1f%%)\n",
                              n_pos_c2n, n_pos_c2n / nrow(resfit_last_c2n) * 100))
            }
      }

      resfit_last_c2n$RID <- as.integer(resfit_last_c2n$subid)
}


###############################################################################
## SECTION 9: Save All Intermediate Results
###############################################################################

cat("\n================================================================\n")
cat("SECTION 9: SAVE INTERMEDIATE\n")
cat("================================================================\n\n")

# Save both Fujirebio and C2N results
# Fujirebio: res_sila_fuj, resfit_fuj, resfit_last_fuj, sila_df_fuj
# C2N:       res_sila_c2n, resfit_c2n, resfit_last_c2n, sila_df_c2n
save(res_sila_fuj, resfit_fuj, resfit_last_fuj, sila_df_fuj,
     res_sila_c2n, resfit_c2n, resfit_last_c2n, sila_df_c2n,
     file = file.path(data_dir, "ADNI_plasma_SILA_intermediate.rda"))
cat("Saved: data/ADNI_plasma_SILA_intermediate.rda\n")
cat("  Contains: Fujirebio SILA + C2N SILA results\n")

# Final summary
cat("\n================================================================\n")
cat("SUMMARY\n")
cat("================================================================\n\n")

cat("--- Fujirebio pT217_F ---\n")
cat(sprintf("Total subjects: %d\n", nrow(resfit_last_fuj)))
if (!is.null(onset_col)) {
      cat(sprintf("Valid EAOA: %d\n", sum(!is.na(resfit_last_fuj$EAOA_plasma))))
}
if ("estpos" %in% names(resfit_last_fuj)) {
      cat(sprintf("Estimated positive: %d\n",
                  sum(resfit_last_fuj$estpos == TRUE, na.rm = TRUE)))
}

cat("\n--- C2N %%p-tau217 ---\n")
if (!is.null(resfit_last_c2n)) {
      cat(sprintf("Total subjects: %d\n", nrow(resfit_last_c2n)))
      cat(sprintf("Valid EAOA: %d\n", sum(!is.na(resfit_last_c2n$EAOA_plasma))))
      if ("estpos" %in% names(resfit_last_c2n)) {
            cat(sprintf("Estimated positive: %d\n",
                        sum(resfit_last_c2n$estpos == TRUE, na.rm = TRUE)))
      }
} else {
      cat("C2N SILA failed to converge — no results available.\n")
}

cat("\n=== ADNI Plasma SILA complete ===\n")
