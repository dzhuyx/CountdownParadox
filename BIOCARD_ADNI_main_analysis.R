###############################################################################
## Component C: Main Analysis — 5 Biomarkers × 3 Models
##
## Purpose: Run P1 (countdown), P2 (TVC z_only), P3 (TVC interaction) for
##   all five biomarker–cohort combinations, compiling a unified results table.
##
##   C1: BIOCARD CSF AB42/AB40 (SILA EAOA from Component A)
##   C2: BIOCARD CSF p-tau181  (SILA EAOA from Component A)
##   C3: BIOCARD plasma p-tau181 (SILA EAOA from Component B)
##   C4: ADNI amyloid PET (carry forward existing results)
##   C5: ADNI plasma p-tau217 Fujirebio (carry forward existing results)
##
## Input:
##   data/BIOCARD_CSF_SILA_intermediate.rda   (Component A)
##   data/BIOCARD_plasma_SILA_intermediate.rda (Component B)
##   data/analysis_data_merged.rda             (survival + covariates)
##   results/ADNI_countdown_vs_tvc_2026.csv    (C4 existing)
##   results/ADNI_plasma_all_models.csv        (C5 existing)
##
## Output:
##   results/main_results_all_biomarkers.csv
##   results/degeneracy_all_biomarkers.csv
##   results/sample_sizes_all_biomarkers.csv
##
## Author: Yuxin Zhu
## Date: March 2026
###############################################################################

rm(list = ls())
library(survival)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT", "/Users/daisyzhu/Documents/Research Projects/CountdownParadox_BiomarkerPositivity/CountdownParadox_Analysis")
out_dir      <- project_root
data_dir     <- file.path(project_root, "data")
results_dir  <- file.path(project_root, "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

covariates <- c("Sex_F", "EDUC_z", "apoe4")


###############################################################################
## Helper functions (copied from BIOCARD_reanalysis_v2.R)
###############################################################################

reportcox <- function(fit) {
      fit_s <- summary(fit)
      result <- round(cbind(fit_s$conf.int[, c(1, 3, 4)],
                            fit_s$coefficients[, 5]), 2)
      result <- as.data.frame(result)
      result[which(result[, 4] == 0), 4] <- "< 0.01"
      names(result) <- c("hazard ratio", "lower .95", "upper .95", "p-value")
      result <- rbind(result, rep(" ", ncol(result)))
      result[nrow(result), 1] <- paste0("n = ", fit_s$n)
      result[nrow(result), 2] <- paste0("nevent = ", fit_s$nevent)
      return(result)
}

extract_hr <- function(fit, var_name, n_override = NULL) {
      # n_override: when fit is a counting-process Cox model, fit_s$n returns
      # the number of (tstart, tstop) interval rows, not unique subjects.
      # Pass length(unique(data$id)) to report subject count instead.
      fit_s <- summary(fit)
      n_val <- if (!is.null(n_override)) n_override else fit_s$n
      idx <- which(rownames(fit_s$coefficients) == var_name)
      if (length(idx) == 0) return(data.frame(
            hr = NA, lower95 = NA, upper95 = NA, pvalue = NA,
            n = n_val, nevent = fit_s$nevent
      ))
      data.frame(
            hr = fit_s$conf.int[idx, 1],
            lower95 = fit_s$conf.int[idx, 3],
            upper95 = fit_s$conf.int[idx, 4],
            pvalue = fit_s$coefficients[idx, 5],
            n = n_val,
            nevent = fit_s$nevent
      )
}

create_tvc_data <- function(dat, A_col, entry_col, exit_col, event_col,
                            estpos_col = NULL) {
      # Use estpos indicator to classify never-positive subjects.
      # When estpos_col is provided (SILA analyses):
      #   estpos=0 → never_positive; estpos=1 → use A for timing; NA A → skip (safety net)
      # When estpos_col is NULL (e.g., first_pos_age):
      #   NA A → never_positive (backward compatible)
      A         <- dat[[A_col]]
      entry_age <- dat[[entry_col]]
      exit_age  <- dat[[exit_col]]
      event     <- dat[[event_col]]
      estpos    <- if (!is.null(estpos_col)) dat[[estpos_col]] else NULL

      rows <- list()
      for (i in seq_along(A)) {
            if (!is.null(estpos) && !is.na(estpos[i]) && !estpos[i]) {
                  # estpos=0: SILA assessed as never reaching positivity
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 0, A = 0,
                        status = "never_positive"
                  )
            } else if (is.na(A[i])) {
                  if (!is.null(estpos_col)) {
                        # SILA mode: estpos=1 but EAOA=NA → treat as never-positive
                        rows[[length(rows) + 1]] <- data.frame(
                              id = i, tstart = entry_age[i], tstop = exit_age[i],
                              event = event[i], Z_tv = 0, A = 0,
                              status = "never_positive"
                        )
                        next
                  }
                  # Non-SILA mode: NA means never-positive (e.g., first_pos_age)
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 0, A = 0,
                        status = "never_positive"
                  )
            } else if (A[i] <= entry_age[i]) {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 1, A = A[i],
                        status = "already_positive"
                  )
            } else if (A[i] < exit_age[i]) {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = A[i],
                        event = 0, Z_tv = 0, A = A[i],
                        status = "transition_pre"
                  )
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = A[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 1, A = A[i],
                        status = "transition_post"
                  )
            } else {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 0, A = 0,
                        status = "never_positive"
                  )
            }
      }
      result <- do.call(rbind, rows)
      result <- result[result$tstop > result$tstart, ]
      result
}


###############################################################################
## Generic function: run P1/P2/P3 for one biomarker
###############################################################################

run_three_models <- function(dat, eaoa_col, entry_col, exit_col, event_col,
                             biomarker_label, cohort_label,
                             censor_col = NULL, estpos_col = NULL) {
      #
      # dat: data frame with EAOA, survival times, covariates
      # eaoa_col: column name for EAOA
      # estpos_col: column with SILA estpos indicator (1=positive, 0=never positive)
      # exit_col: column for event/censor time in survival models (onset.age for BIOCARD)
      # censor_col: column for last-visit age (censor.age). If NULL, defaults to exit_col.
      #             Used only for follow-up computation.
      # biomarker_label: e.g. "CSF_AB42_AB40"
      # cohort_label: e.g. "BIOCARD"
      #
      cat(sprintf("\n================================================================\n"))
      cat(sprintf("  %s — %s\n", cohort_label, biomarker_label))
      cat(sprintf("================================================================\n\n"))

      dat_c <- dat[complete.cases(dat[, covariates]), ]
      # Exclude subjects with zero follow-up (entry >= exit)
      n_zero_fu <- sum(dat_c[[entry_col]] >= dat_c[[exit_col]])
      if (n_zero_fu > 0) {
            cat(sprintf("Excluding %d subject(s) with zero follow-up\n", n_zero_fu))
            dat_c <- dat_c[which(dat_c[[entry_col]] < dat_c[[exit_col]]), ]
      }
      n_total <- nrow(dat_c)
      n_pos   <- if (!is.null(estpos_col)) sum(dat_c[[estpos_col]] == 1, na.rm = TRUE) else sum(!is.na(dat_c[[eaoa_col]]))
      n_events <- sum(dat_c[[event_col]])
      cat(sprintf("Sample: %d total, %d positive, %d events\n",
                  n_total, n_pos, n_events))

      results <- data.frame()

      # ---- P1: Countdown ---- #
      cat("\n--- P1: Countdown (T-Z ~ Z) ---\n")

      dat_c$time_to_dx <- dat_c[[exit_col]] - dat_c[[eaoa_col]]
      pos_filter <- if (!is.null(estpos_col)) which(dat_c[[estpos_col]] == 1) else which(!is.na(dat_c[[eaoa_col]]))
      dat_p1 <- dat_c[intersect(pos_filter, which(dat_c$time_to_dx > 0)), ]

      if (nrow(dat_p1) >= 10 && sum(dat_p1[[event_col]]) >= 5) {
            Z_mean <- mean(dat_p1[[eaoa_col]])
            Z_sd   <- sd(dat_p1[[eaoa_col]])
            dat_p1$Z_std <- (dat_p1[[eaoa_col]] - Z_mean) / Z_sd
            cat(sprintf("Z standardization: mean = %.1f, SD = %.1f\n", Z_mean, Z_sd))
            cat(sprintf("Countdown sample: n = %d, events = %d\n",
                        nrow(dat_p1), sum(dat_p1[[event_col]])))

            fml_p1 <- as.formula(paste0("Surv(time_to_dx, ", event_col,
                                         ") ~ Sex_F + EDUC_z + apoe4 + Z_std"))
            fit_p1 <- coxph(fml_p1, data = dat_p1)
            print(reportcox(fit_p1))

            hr_p1 <- extract_hr(fit_p1, "Z_std")
      } else {
            cat(sprintf("SKIPPED: insufficient (n=%d, events=%d)\n",
                        nrow(dat_p1), sum(dat_p1[[event_col]])))
            hr_p1 <- data.frame(hr = NA, lower95 = NA, upper95 = NA, pvalue = NA,
                                n = nrow(dat_p1), nevent = sum(dat_p1[[event_col]]))
      }
      hr_p1$biomarker <- biomarker_label
      hr_p1$cohort    <- cohort_label
      hr_p1$model     <- "P1_countdown"
      hr_p1$parameter <- "Z_std"
      results <- rbind(results, hr_p1)

      # ---- P2: TVC z_only ---- #
      cat("\n--- P2: TVC z_only ---\n")

      tvc_long <- create_tvc_data(dat_c, eaoa_col, entry_col, exit_col, event_col,
                                   estpos_col)

      # Merge covariates
      covar_df <- data.frame(
            id = seq_len(nrow(dat_c)),
            Sex_F = dat_c$Sex_F,
            EDUC_z = dat_c$EDUC_z,
            apoe4 = dat_c$apoe4
      )
      tvc_long <- merge(tvc_long, covar_df, by = "id", all.x = TRUE)
      tvc_long <- tvc_long[complete.cases(tvc_long[, covariates]), ]

      n_pos_rows <- sum(tvc_long$Z_tv == 1)
      n_neg_rows <- sum(tvc_long$Z_tv == 0)
      n_tvc_subj <- length(unique(tvc_long$id))
      cat(sprintf("TVC data: %d rows (%d Z_tv=1, %d Z_tv=0), %d subjects\n",
                  nrow(tvc_long), n_pos_rows, n_neg_rows, n_tvc_subj))

      # ---- Tier 1 assertions ---- #
      # Check 1.2: TVC subject count = analysis N (no silent drops)
      stopifnot("TVC subject count must equal analysis N" =
                n_tvc_subj == nrow(dat_c))
      # Check 1.3: estpos=0 subjects must all have Z_tv=0
      if (!is.null(estpos_col)) {
            estpos_0_ids <- which(dat_c[[estpos_col]] == 0)
            if (length(estpos_0_ids) > 0) {
                  stopifnot("estpos=0 subjects must have Z_tv=0" =
                            all(tvc_long$Z_tv[tvc_long$id %in% estpos_0_ids] == 0))
                  # Check 1.5: estpos=0 subjects must have A=0 (no extrapolated EAOA used)
                  stopifnot("estpos=0 subjects must have A=0 in TVC" =
                            all(tvc_long$A[tvc_long$id %in% estpos_0_ids] == 0))
            }
      }
      cat("  Assertions passed: subject count, estpos↔Z_tv concordance\n")

      fit_p2 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + Sex_F + EDUC_z + apoe4,
                       data = tvc_long)
      print(reportcox(fit_p2))

      hr_p2 <- extract_hr(fit_p2, "Z_tv", n_override = n_tvc_subj)
      hr_p2$biomarker <- biomarker_label
      hr_p2$cohort    <- cohort_label
      hr_p2$model     <- "P2_tvc_z_only"
      hr_p2$parameter <- "Z_tv"
      results <- rbind(results, hr_p2)

      # ---- P3: TVC interaction ---- #
      cat("\n--- P3: TVC interaction ---\n")

      # Standardize A among positive subjects — use unique subjects, not TVC rows
      # (transition subjects have 2 rows; using all rows double-counts them)
      pos_unique <- unique(tvc_long[tvc_long$A > 0, c("id", "A")])$A
      A_mean <- mean(pos_unique)
      A_sd   <- sd(pos_unique)
      tvc_long$A_z <- ifelse(tvc_long$A > 0, (tvc_long$A - A_mean) / A_sd, 0)
      cat(sprintf("A standardization: mean = %.1f, SD = %.1f (n_unique = %d)\n",
                  A_mean, A_sd, length(pos_unique)))

      # Check degeneracy
      all_positive <- all(tvc_long$Z_tv == 1)
      cat(sprintf("All Z_tv=1? %s\n", all_positive))

      if (all_positive) {
            # B8: When all subjects are already positive at entry, Z_tv has no
            # variation. The P3 model reduces to ~ A_z + covariates (no Z_tv term).
            # The gamma (A_z) coefficient is still interpretable as the modulation
            # of hazard by EAOA among positive subjects. Manuscript methods should
            # note this reduced specification for degenerate biomarkers.
            cat("DEGENERATE: all Z_tv = 1. Fitting reduced model ~ A_z only.\n")
            fit_p3 <- tryCatch({
                  coxph(Surv(tstart, tstop, event) ~ A_z + Sex_F + EDUC_z + apoe4,
                        data = tvc_long)
            }, error = function(e) {
                  cat(sprintf("  Model error: %s\n", e$message))
                  NULL
            })

            if (!is.null(fit_p3)) {
                  print(reportcox(fit_p3))

                  hr_gamma <- extract_hr(fit_p3, "A_z", n_override = n_tvc_subj)
                  hr_gamma$biomarker <- biomarker_label
                  hr_gamma$cohort    <- cohort_label
                  hr_gamma$model     <- "P3_tvc_interaction"
                  hr_gamma$parameter <- "gamma_A"
                  results <- rbind(results, hr_gamma)

                  results <- rbind(results, data.frame(
                        hr = NA, lower95 = NA, upper95 = NA, pvalue = NA,
                        n = n_tvc_subj, nevent = summary(fit_p3)$nevent,
                        biomarker = biomarker_label, cohort = cohort_label,
                        model = "P3_tvc_interaction", parameter = "beta_Z_tv"
                  ))
            }
      } else {
            fit_p3 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + A_z:Z_tv +
                                   Sex_F + EDUC_z + apoe4,
                             data = tvc_long)
            print(reportcox(fit_p3))

            # Beta (Z_tv)
            hr_beta <- extract_hr(fit_p3, "Z_tv", n_override = n_tvc_subj)
            hr_beta$biomarker <- biomarker_label
            hr_beta$cohort    <- cohort_label
            hr_beta$model     <- "P3_tvc_interaction"
            hr_beta$parameter <- "beta_Z_tv"
            results <- rbind(results, hr_beta)

            # Gamma (A_z:Z_tv interaction)
            gamma_name <- setdiff(rownames(summary(fit_p3)$coefficients),
                                  c("Z_tv", "Sex_F", "EDUC_z", "apoe4"))
            if (length(gamma_name) == 1) {
                  hr_gamma <- extract_hr(fit_p3, gamma_name, n_override = n_tvc_subj)
            } else {
                  hr_gamma <- data.frame(hr = NA, lower95 = NA, upper95 = NA,
                                          pvalue = NA, n = n_tvc_subj,
                                          nevent = summary(fit_p3)$nevent)
            }
            hr_gamma$biomarker <- biomarker_label
            hr_gamma$cohort    <- cohort_label
            hr_gamma$model     <- "P3_tvc_interaction"
            hr_gamma$parameter <- "gamma_A"
            results <- rbind(results, hr_gamma)

            # Correlation between beta and gamma
            vc <- vcov(fit_p3)
            z_idx <- which(rownames(vc) == "Z_tv")
            g_idx <- which(rownames(vc) == gamma_name)
            if (length(z_idx) == 1 && length(g_idx) == 1) {
                  cor_bg <- cov2cor(vc)[z_idx, g_idx]
                  cat(sprintf("  vcov correlation (beta, gamma) = %.3f\n", cor_bg))
            }
      }

      # ---- Degeneracy summary ---- #
      A_vals    <- dat_c[[eaoa_col]]
      entry     <- dat_c[[entry_col]]
      exit      <- dat_c[[exit_col]]

      if (!is.null(estpos_col)) {
            pos_idx       <- which(dat_c[[estpos_col]] == 1)
            n_valid       <- length(pos_idx)
            n_already_pos <- sum(A_vals[pos_idx] <= entry[pos_idx])
            n_transition  <- sum(A_vals[pos_idx] > entry[pos_idx] & A_vals[pos_idx] < exit[pos_idx])
            n_never_pos   <- n_total - n_already_pos - n_transition
      } else {
            valid_idx     <- which(!is.na(A_vals))
            n_valid       <- length(valid_idx)
            n_already_pos <- sum(A_vals[valid_idx] <= entry[valid_idx])
            n_transition  <- sum(A_vals[valid_idx] > entry[valid_idx] & A_vals[valid_idx] < exit[valid_idx])
            n_never_pos   <- sum(A_vals[valid_idx] >= exit[valid_idx]) + (n_total - n_valid)
      }
      pct_already   <- if (n_valid > 0) n_already_pos / n_valid * 100 else NA

      cat(sprintf("\nDegeneracy: valid=%d, already_pos=%d (%.0f%%), transition=%d, never_pos=%d\n",
                  n_valid, n_already_pos, pct_already, n_transition, n_never_pos))
      # Check 1.4: Degeneracy arithmetic
      stopifnot("Degeneracy must sum to n_total" =
                n_already_pos + n_transition + n_never_pos == n_total)

      # Follow-up = censor.age - baseline.age (total observation)
      # Person-time = onset.age - baseline.age (time at risk)
      censor_age <- if (!is.null(censor_col)) dat_c[[censor_col]] else exit
      degen <- data.frame(
            biomarker = biomarker_label, cohort = cohort_label,
            n_total = n_total, n_events = n_events,
            n_valid_eaoa = n_valid, n_already_pos = n_already_pos,
            pct_already_pos = round(pct_already, 1),
            n_transition = n_transition, n_never_pos = n_never_pos,
            followup_median = round(median(censor_age - entry), 1),
            persontime_median = round(median(exit - entry), 1),
            stringsAsFactors = FALSE
      )

      return(list(results = results, degeneracy = degen))
}


###############################################################################
## SECTION 1: Load BIOCARD data
###############################################################################

cat("================================================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("================================================================\n\n")

# --- BIOCARD survival + covariates --- #
load(file.path(data_dir, "analysis_data_merged.rda"))
cat(sprintf("analysis_data: %d subjects, %d events\n",
            nrow(analysis_data), sum(analysis_data$d)))

# --- Safeguard assertions (prevent censor.age/onset.age confusion) --- #
# For BIOCARD events: onset.age < censor.age (onset is earlier than last visit)
stopifnot("onset.age must be <= censor.age for all subjects" =
          all(analysis_data$onset.age <= analysis_data$censor.age))
stopifnot("onset.age must be < censor.age for at least some events" =
          any(analysis_data$onset.age[analysis_data$d == 1] <
              analysis_data$censor.age[analysis_data$d == 1]))
# baseline.age must be positive and reasonable
stopifnot("baseline.age must be > 0" = all(analysis_data$baseline.age > 0))
stopifnot("onset.age must be > baseline.age" =
          all(analysis_data$onset.age > analysis_data$baseline.age))
cat("  Assertions passed: onset.age/censor.age/baseline.age verified.\n")

# --- BIOCARD CSF SILA (Component A) --- #
load(file.path(data_dir, "BIOCARD_CSF_SILA_intermediate.rda"))
cat(sprintf("CSF SILA AB: %d subjects, %d positive\n",
            nrow(resfit_last_ab), sum(resfit_last_ab$estpos)))
cat(sprintf("CSF SILA PTAU: %d subjects, %d positive\n",
            nrow(resfit_last_ptau), sum(resfit_last_ptau$estpos)))

# Merge SILA EAOA into analysis_data
# Keep EAOA values for all subjects (including estpos=0).
# estpos indicator is passed to create_tvc_data to classify never-positive.
ab_eaoa <- resfit_last_ab[, c("SUBJECT_ID", "EAOA_AB", "estpos")]
names(ab_eaoa)[3] <- "estpos_ab"

ptau_eaoa <- resfit_last_ptau[, c("SUBJECT_ID", "EAOA_PTAU", "estpos")]
names(ptau_eaoa)[3] <- "estpos_ptau"

dat_biocard_csf <- merge(analysis_data, ab_eaoa[, c("SUBJECT_ID", "EAOA_AB", "estpos_ab")],
                         by = "SUBJECT_ID", all.x = TRUE)
dat_biocard_csf <- merge(dat_biocard_csf, ptau_eaoa[, c("SUBJECT_ID", "EAOA_PTAU", "estpos_ptau")],
                         by = "SUBJECT_ID", all.x = TRUE)

cat(sprintf("After merge: %d subjects, AB positive=%d, PTAU positive=%d\n",
            nrow(dat_biocard_csf),
            sum(dat_biocard_csf$estpos_ab == 1, na.rm = TRUE),
            sum(dat_biocard_csf$estpos_ptau == 1, na.rm = TRUE)))
# Check: estpos must be defined for all SILA subjects
stopifnot("estpos_ab must be defined for all merged subjects" =
          all(!is.na(dat_biocard_csf$estpos_ab)))
stopifnot("estpos_ptau must be defined for all merged subjects" =
          all(!is.na(dat_biocard_csf$estpos_ptau)))

# --- BIOCARD plasma SILA (Component B) --- #
load(file.path(data_dir, "BIOCARD_plasma_SILA_intermediate.rda"))

# Keep EAOA values for all subjects; estpos used by create_tvc_data
cat(sprintf("Plasma SILA: %d subjects, %d positive (estpos=1), %d events\n",
            nrow(plasma_analysis),
            sum(plasma_analysis$estpos == 1, na.rm = TRUE),
            sum(plasma_analysis$d)))


###############################################################################
## SECTION 2: C1 — BIOCARD CSF AB42/AB40
###############################################################################

all_results <- data.frame()
all_degen   <- data.frame()

out_c1 <- run_three_models(
      dat = dat_biocard_csf,
      eaoa_col = "EAOA_AB",
      entry_col = "baseline.age",
      exit_col = "onset.age",
      event_col = "d",
      biomarker_label = "CSF_AB42_AB40",
      cohort_label = "BIOCARD",
      censor_col = "censor.age",
      estpos_col = "estpos_ab"
)
all_results <- rbind(all_results, out_c1$results)
all_degen   <- rbind(all_degen, out_c1$degeneracy)


###############################################################################
## SECTION 3: C2 — BIOCARD CSF p-tau181
###############################################################################

out_c2 <- run_three_models(
      dat = dat_biocard_csf,
      eaoa_col = "EAOA_PTAU",
      entry_col = "baseline.age",
      exit_col = "onset.age",
      event_col = "d",
      biomarker_label = "CSF_pTau181",
      cohort_label = "BIOCARD",
      censor_col = "censor.age",
      estpos_col = "estpos_ptau"
)
all_results <- rbind(all_results, out_c2$results)
all_degen   <- rbind(all_degen, out_c2$degeneracy)


###############################################################################
## SECTION 4: C3 — BIOCARD plasma p-tau181
###############################################################################

out_c3 <- run_three_models(
      dat = plasma_analysis,
      eaoa_col = "EAOA_plasma",
      entry_col = "baseline.age",
      exit_col = "onset.age",
      event_col = "d",
      biomarker_label = "Plasma_pTau181",
      cohort_label = "BIOCARD",
      censor_col = "censor.age",
      estpos_col = "estpos"
)
all_results <- rbind(all_results, out_c3$results)
all_degen   <- rbind(all_degen, out_c3$degeneracy)


###############################################################################
## SECTION 5: C4 — ADNI amyloid PET (carry forward)
###############################################################################

cat("\n================================================================\n")
cat("  C4: ADNI Amyloid PET — Carry Forward Existing Results\n")
cat("================================================================\n\n")

adni_pet_file <- file.path(results_dir, "ADNI_countdown_vs_tvc_2026.csv")
if (file.exists(adni_pet_file)) {
      adni_pet <- read.csv(adni_pet_file, stringsAsFactors = FALSE)
      # Keep only SILA method (not first_positive sensitivity)
      adni_pet <- adni_pet[adni_pet$method == "SILA", ]

      # Map to unified format
      for (i in seq_len(nrow(adni_pet))) {
            row <- adni_pet[i, ]
            model_map <- c("A1_countdown" = "P1_countdown",
                           "A2_tvc_z_only" = "P2_tvc_z_only",
                           "A3_tvc_interaction" = "P3_tvc_interaction")
            param_map <- c("EAOA_z" = "Z_std")

            unified_model <- model_map[row$model]
            unified_param <- ifelse(row$parameter %in% names(param_map),
                                    param_map[row$parameter], row$parameter)

            all_results <- rbind(all_results, data.frame(
                  hr = row$hr, lower95 = row$lower95, upper95 = row$upper95,
                  pvalue = row$pvalue, n = row$n, nevent = row$nevent,
                  biomarker = "Amyloid_PET_FBP", cohort = "ADNI",
                  model = unified_model, parameter = unified_param,
                  stringsAsFactors = FALSE
            ))
      }
      cat(sprintf("Carried forward %d ADNI PET results (SILA method)\n", nrow(adni_pet)))
} else {
      cat("WARNING: ADNI PET results file not found!\n")
}

# ADNI PET degeneracy
adni_pet_degen_file <- file.path(results_dir, "ADNI_degeneracy_2026.csv")
if (file.exists(adni_pet_degen_file)) {
      adni_degen <- read.csv(adni_pet_degen_file, stringsAsFactors = FALSE)
      # Reformat to match our unified format
      # Get n_events from carried-forward results (use P2 for full-cohort event count)
      adni_pet_nevent <- all_results$nevent[all_results$cohort == "ADNI" &
            all_results$biomarker == "Amyloid_PET_FBP" & all_results$model == "P2_tvc_z_only"]
      adni_pet_nevent <- if (length(adni_pet_nevent) > 0) adni_pet_nevent[1] else NA

      all_degen <- rbind(all_degen, data.frame(
            biomarker = "Amyloid_PET_FBP", cohort = "ADNI",
            n_total = adni_degen$n_total,
            n_events = adni_pet_nevent,
            n_valid_eaoa = adni_degen$n_valid,
            n_already_pos = adni_degen$n_already_pos,
            pct_already_pos = adni_degen$pct_already_pos,
            n_transition = adni_degen$n_transition,
            n_never_pos = adni_degen$n_never_pos,
            followup_median = NA,  # populated by the descriptives script
            persontime_median = NA,
            stringsAsFactors = FALSE
      ))
      cat("Carried forward ADNI PET degeneracy data\n")
}


###############################################################################
## SECTION 6: C5 — ADNI plasma p-tau217 Fujirebio (carry forward)
###############################################################################

cat("\n================================================================\n")
cat("  C5: ADNI Plasma p-tau217 Fujirebio — Carry Forward\n")
cat("================================================================\n\n")

adni_plasma_file <- file.path(results_dir, "ADNI_plasma_all_models.csv")
if (file.exists(adni_plasma_file)) {
      adni_plasma <- read.csv(adni_plasma_file, stringsAsFactors = FALSE)
      # Keep only Fujirebio SILA, P1/P2/P3 (no landmark)
      fuj_sila <- adni_plasma[adni_plasma$config == "Fujirebio_SILA" &
                              adni_plasma$model %in% c("P1_countdown", "P2_tvc_z_only",
                                                        "P3_tvc_interaction"), ]

      for (i in seq_len(nrow(fuj_sila))) {
            row <- fuj_sila[i, ]
            all_results <- rbind(all_results, data.frame(
                  hr = row$hr, lower95 = row$lower95, upper95 = row$upper95,
                  pvalue = row$pvalue, n = row$n, nevent = row$nevent,
                  biomarker = "Plasma_pTau217_Fuji", cohort = "ADNI",
                  model = row$model, parameter = row$parameter,
                  stringsAsFactors = FALSE
            ))
      }
      cat(sprintf("Carried forward %d ADNI plasma Fujirebio SILA results\n",
                  nrow(fuj_sila)))
} else {
      cat("WARNING: ADNI plasma results file not found!\n")
}

# ADNI plasma degeneracy
adni_plasma_degen_file <- file.path(results_dir, "ADNI_plasma_degeneracy.csv")
if (file.exists(adni_plasma_degen_file)) {
      ap_degen <- read.csv(adni_plasma_degen_file, stringsAsFactors = FALSE)
      fuj_degen <- ap_degen[ap_degen$config == "Fujirebio_SILA", ]
      if (nrow(fuj_degen) > 0) {
            # Get n_events from carried-forward results (use P2 for full-cohort event count)
            adni_plasma_nevent <- all_results$nevent[all_results$cohort == "ADNI" &
                  all_results$biomarker == "Plasma_pTau217_Fuji" &
                  all_results$model == "P2_tvc_z_only"]
            adni_plasma_nevent <- if (length(adni_plasma_nevent) > 0) adni_plasma_nevent[1] else NA

            all_degen <- rbind(all_degen, data.frame(
                  biomarker = "Plasma_pTau217_Fuji", cohort = "ADNI",
                  n_total = fuj_degen$n_total,
                  n_events = adni_plasma_nevent,
                  n_valid_eaoa = fuj_degen$n_valid,
                  n_already_pos = fuj_degen$n_already_pos,
                  pct_already_pos = fuj_degen$pct_already_pos,
                  n_transition = fuj_degen$n_transition,
                  n_never_pos = fuj_degen$n_never_pos,
                  followup_median = NA,  # populated by the descriptives script
                  persontime_median = NA,
                  stringsAsFactors = FALSE
            ))
            cat("Carried forward ADNI plasma degeneracy data\n")
      }
}


###############################################################################
## SECTION 7: Compile and Save
###############################################################################

cat("\n================================================================\n")
cat("SECTION 7: COMPILE AND SAVE\n")
cat("================================================================\n\n")

# Format results
all_results$HR_95CI <- ifelse(is.na(all_results$hr), "N/A",
                               sprintf("%.2f (%.2f, %.2f)", all_results$hr,
                                       all_results$lower95, all_results$upper95))
all_results$p_formatted <- ifelse(is.na(all_results$pvalue), "N/A",
                                   ifelse(all_results$pvalue < 0.001, "< 0.001",
                                          ifelse(all_results$pvalue < 0.01, "< 0.01",
                                                 sprintf("%.3f", all_results$pvalue))))

write.csv(all_results, file.path(results_dir, "main_results_all_biomarkers.csv"),
          row.names = FALSE)
cat("Saved: results/main_results_all_biomarkers.csv\n")

write.csv(all_degen, file.path(results_dir, "degeneracy_all_biomarkers.csv"),
          row.names = FALSE)
cat("Saved: results/degeneracy_all_biomarkers.csv\n")

# Sample sizes per biomarker-cohort
sample_sizes <- all_degen[, c("biomarker", "cohort", "n_total", "n_events",
                               "n_valid_eaoa", "n_already_pos")]
write.csv(sample_sizes, file.path(results_dir, "sample_sizes_all_biomarkers.csv"),
          row.names = FALSE)
cat("Saved: results/sample_sizes_all_biomarkers.csv\n")


###############################################################################
## SECTION 8: Summary Table
###############################################################################

cat("\n================================================================\n")
cat("SUMMARY: MAIN RESULTS TABLE\n")
cat("================================================================\n\n")

# Print in a structured format
biomarker_order <- c("CSF_AB42_AB40", "CSF_pTau181", "Plasma_pTau181",
                     "Amyloid_PET_FBP", "Plasma_pTau217_Fuji")

cat(sprintf("%-22s  %-8s  %-22s  %-22s  %-22s  %-22s\n",
            "Biomarker", "Cohort", "P1 Countdown",
            "P2 TVC z_only", "P3 beta (Z_tv)", "P3 gamma (A)"))
cat(strrep("-", 130), "\n")

for (bm in biomarker_order) {
      r <- all_results[all_results$biomarker == bm, ]
      if (nrow(r) == 0) next

      cohort <- r$cohort[1]

      p1 <- r[r$model == "P1_countdown", ]
      p1_str <- if (nrow(p1) > 0) p1$HR_95CI[1] else "N/A"

      p2 <- r[r$model == "P2_tvc_z_only", ]
      p2_str <- if (nrow(p2) > 0) p2$HR_95CI[1] else "N/A"

      p3_beta <- r[r$model == "P3_tvc_interaction" & r$parameter == "beta_Z_tv", ]
      p3b_str <- if (nrow(p3_beta) > 0) p3_beta$HR_95CI[1] else "N/A"

      p3_gamma <- r[r$model == "P3_tvc_interaction" & r$parameter == "gamma_A", ]
      p3g_str <- if (nrow(p3_gamma) > 0) p3_gamma$HR_95CI[1] else "N/A"

      cat(sprintf("%-22s  %-8s  %-22s  %-22s  %-22s  %-22s\n",
                  bm, cohort, p1_str, p2_str, p3b_str, p3g_str))
}

cat("\n\nDegeneracy Summary:\n")
cat(sprintf("%-22s  %-8s  %5s  %5s  %5s  %5s  %5s  %5s\n",
            "Biomarker", "Cohort", "N", "Event", "Valid", "AlrPos", "%Pos", "Trans"))
cat(strrep("-", 80), "\n")
for (i in seq_len(nrow(all_degen))) {
      d <- all_degen[i, ]
      cat(sprintf("%-22s  %-8s  %5d  %5s  %5d  %5d  %5.1f  %5d\n",
                  d$biomarker, d$cohort, d$n_total,
                  ifelse(is.na(d$n_events), "?", as.character(d$n_events)),
                  d$n_valid_eaoa, d$n_already_pos,
                  d$pct_already_pos, d$n_transition))
}

cat("\n=== Component C: Main Analysis Complete ===\n")
