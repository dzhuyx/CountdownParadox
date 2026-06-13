###############################################################################
## ADNI SILA Re-Analysis v2: Countdown Paradox — TVC-First (2026 Data)
##
## Purpose: Apply countdown and TVC analyses to ADNI amyloid PET data.
##   Model A1: Countdown — Surv(T-Z, d) ~ EAOA_z + covariates
##   Model A2: TVC z_only (secondary/supportive) — Surv(tstart, tstop, event) ~ Z_tv + covariates
##   Model A3: TVC interaction (primary) — Surv(tstart, tstop, event) ~ Z_tv + A:Z_tv + covariates
##   Sensitivity: first-positive-scan as alternative Z
##
## Input:  data/ADNI_SILA_intermediate_2026.rda (saved from prior SILA run)
## Output: results/ADNI_countdown_vs_tvc_2026.csv
##         results/ADNI_degeneracy_2026.csv
##
## Usage:  source() from reproducibility/ directory or via run_all.R
##
## Author: Yuxin Zhu
## Date: February 2026 (reproducibility copy April 2026)
###############################################################################

rm(list = ls())
library(survival)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT", "/Users/daisyzhu/Documents/Research Projects/CountdownParadox_BiomarkerPositivity/CountdownParadox_Analysis")
data_dir     <- file.path(project_root, "data")
res_dir      <- file.path(project_root, "results")
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)


###############################################################################
## Helper functions
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
      # n_override: for counting-process fits, fit_s$n is interval rows, not subjects.
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

# TVC data creation (from BIOCARD_ADNI_tvc_degeneracy_check.R)
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
      # Remove rows with zero or negative interval (e.g., subjects with no follow-up)
      result <- result[result$tstop > result$tstart, ]
      result
}


###############################################################################
## SECTION 1: Load Intermediate Data
###############################################################################

cat("================================================================\n")
cat("SECTION 1: LOAD INTERMEDIATE DATA\n")
cat("================================================================\n\n")

load(file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"))

# Merge EAOA into CN survival dataset — restrict to subjects with SILA estimates
# (SILA trained on >=2 PET scans, applied to all subjects with >=1 scan)
eaoa_df <- resfit_last[, c("RID", "estaget0", "estpos", "first_pos_age")]
names(eaoa_df)[2] <- "EAOA"
dat_analysis <- merge(dat_surv_cn, eaoa_df, by = "RID", all.x = FALSE)
cat(sprintf("Merged SILA estimates: %d of %d CN subjects\n",
            nrow(dat_analysis), nrow(dat_surv_cn)))
stopifnot("All subjects must have SILA EAOA after restriction" =
          all(!is.na(dat_analysis$EAOA)))

# Exclude 78 participants removed by ADNI (data quality notice, March 2026)
data_2026 <- file.path(project_root, "ADNI_2026_data")
adni_exclude_demo <- read.csv(file.path(data_2026, "PTDEMOG_11Feb2026.csv"))
adni_exclude_rids <- unique(adni_exclude_demo$RID[grepl("^381_S_", adni_exclude_demo$PTID)])
n_excluded <- sum(dat_analysis$RID %in% adni_exclude_rids)
dat_analysis <- dat_analysis[!dat_analysis$RID %in% adni_exclude_rids, ]
cat(sprintf("Excluded %d of %d flagged 381_S subjects\n", n_excluded, length(adni_exclude_rids)))
rm(adni_exclude_demo)

# Add demographics
dat_analysis <- merge(dat_analysis, dat_demo[, c("RID", "PTGENDER", "PTEDUCAT")],
                      by = "RID", all.x = TRUE)
dat_analysis$Sex_F <- as.integer(dat_analysis$PTGENDER == 2)
dat_analysis$EDUC_z <- scale(dat_analysis$PTEDUCAT)[, 1]

# Add APOE4
dat_analysis <- merge(dat_analysis, dat_apoe[, c("RID", "apoe4")],
                      by = "RID", all.x = TRUE)

# Compute countdown time
dat_analysis$time_EAOA2DX <- dat_analysis$onset.age - dat_analysis$EAOA
dat_analysis$time_fps2DX  <- dat_analysis$onset.age - dat_analysis$first_pos_age

# Standardize EAOA among estpos=1 subjects (SILA gives EAOA to all subjects, but
# sub-threshold subjects have extrapolated onset ages >> 100 that inflate the SD)
eaoa_pos_pool <- which(!is.na(dat_analysis$EAOA) & dat_analysis$estpos == 1)
eaoa_mean_pos <- mean(dat_analysis$EAOA[eaoa_pos_pool])
eaoa_sd_pos   <- sd(dat_analysis$EAOA[eaoa_pos_pool])
dat_analysis$EAOA_z <- (dat_analysis$EAOA - eaoa_mean_pos) / eaoa_sd_pos

# First-positive-scan scaling
valid_fps <- !is.na(dat_analysis$first_pos_age)
dat_analysis$fps_z <- NA
dat_analysis$fps_z[valid_fps] <- scale(dat_analysis$first_pos_age[valid_fps])[, 1]

covariates <- c("Sex_F", "EDUC_z", "apoe4")
dat_complete <- dat_analysis[complete.cases(dat_analysis[, covariates]), ]

cat(sprintf("CN subjects: %d total, %d with complete covariates\n",
            nrow(dat_analysis), nrow(dat_complete)))
cat(sprintf("  With SILA EAOA: %d (estpos=1: %d)\n",
            sum(!is.na(dat_complete$EAOA)),
            sum(!is.na(dat_complete$EAOA) & dat_complete$estpos == 1)))
cat(sprintf("  EAOA_z scaling pool (estpos=1): mean = %.1f, SD = %.1f\n",
            eaoa_mean_pos, eaoa_sd_pos))
cat(sprintf("  Events: %d (%.1f%%)\n\n", sum(dat_complete$d == 1),
            mean(dat_complete$d == 1) * 100))

# Save the PET analysis cohort (n = 575) for manuscript descriptives.
# This is the exact cohort the TV models below are fit on.
save(dat_complete, file = file.path(data_dir, "ADNI_pet_analysis_cohort.rda"))
cat(sprintf("Saved PET analysis cohort (%d subjects) -> data/ADNI_pet_analysis_cohort.rda\n",
            nrow(dat_complete)))

hr_results <- data.frame()


###############################################################################
## SECTION 2: Model A1 — Countdown Analysis
###############################################################################

cat("================================================================\n")
cat("MODEL A1: COUNTDOWN ANALYSIS (SILA EAOA)\n")
cat("================================================================\n\n")

dat_a1 <- dat_complete[which(!is.na(dat_complete$EAOA) &
                              dat_complete$estpos == 1 &
                              dat_complete$time_EAOA2DX > 0), ]

# Re-standardize EAOA within P1 subset (consistent with BIOCARD_ADNI_main_analysis.R)
eaoa_mean_a1 <- mean(dat_a1$EAOA)
eaoa_sd_a1   <- sd(dat_a1$EAOA)
dat_a1$EAOA_z <- (dat_a1$EAOA - eaoa_mean_a1) / eaoa_sd_a1
cat(sprintf("A1 sample: n = %d, events = %d\n", nrow(dat_a1), sum(dat_a1$d == 1)))
cat(sprintf("A1 Z standardization: mean = %.1f, SD = %.1f\n", eaoa_mean_a1, eaoa_sd_a1))

fit_a1 <- coxph(Surv(time_EAOA2DX, d) ~ Sex_F + EDUC_z + apoe4 + EAOA_z,
                 data = dat_a1)
cat("Cox model:\n")
print(reportcox(fit_a1))

hr_a1 <- extract_hr(fit_a1, "EAOA_z")
hr_a1$model <- "A1_countdown"
hr_a1$method <- "SILA"
hr_a1$parameter <- "EAOA_z"
hr_a1$tracer <- "FBP"
hr_results <- rbind(hr_results, hr_a1)

cat("\n")


###############################################################################
## SECTION 3: Model A2 — TVC z_only Analysis (Primary)
###############################################################################

cat("================================================================\n")
cat("MODEL A2: TVC z_only ANALYSIS (PRIMARY)\n")
cat("================================================================\n\n")

# Create TVC data for all CN subjects
# estpos used to classify never-positive (not NA in EAOA)
tvc_long <- create_tvc_data(dat_complete, "EAOA",
                            "baseline.age", "onset.age", "d",
                            estpos_col = "estpos")

# Merge covariates
covar_df <- data.frame(
      id = seq_len(nrow(dat_complete)),
      Sex_F = dat_complete$Sex_F,
      EDUC_z = dat_complete$EDUC_z,
      apoe4 = dat_complete$apoe4
)
tvc_long <- merge(tvc_long, covar_df, by = "id", all.x = TRUE)

n_pos_rows <- sum(tvc_long$Z_tv == 1)
n_neg_rows <- sum(tvc_long$Z_tv == 0)
n_tvc_subj <- length(unique(tvc_long$id))
cat(sprintf("TVC data: %d rows (%d Z_tv=1, %d Z_tv=0), %d unique subjects\n",
            nrow(tvc_long), n_pos_rows, n_neg_rows, n_tvc_subj))

# ---- Tier 1 assertions ---- #
stopifnot("TVC subject count must equal analysis N" =
          n_tvc_subj == nrow(dat_complete))
estpos_0_ids <- which(dat_complete$estpos == 0)
if (length(estpos_0_ids) > 0) {
      stopifnot("estpos=0 subjects must have Z_tv=0" =
                all(tvc_long$Z_tv[tvc_long$id %in% estpos_0_ids] == 0))
      stopifnot("estpos=0 subjects must have A=0 in TVC" =
                all(tvc_long$A[tvc_long$id %in% estpos_0_ids] == 0))
}
cat("  Assertions passed: subject count, estpos↔Z_tv concordance\n")

fit_a2 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + Sex_F + EDUC_z + apoe4,
                 data = tvc_long)
cat("Cox model:\n")
print(reportcox(fit_a2))

hr_a2 <- extract_hr(fit_a2, "Z_tv", n_override = n_tvc_subj)
hr_a2$model <- "A2_tvc_z_only"
hr_a2$method <- "SILA"
hr_a2$parameter <- "Z_tv"
hr_a2$tracer <- "FBP"
hr_results <- rbind(hr_results, hr_a2)

cat("\n")


###############################################################################
## SECTION 4: Model A3 — TVC Interaction Analysis (Secondary)
###############################################################################

cat("================================================================\n")
cat("MODEL A3: TVC INTERACTION ANALYSIS (SECONDARY)\n")
cat("================================================================\n\n")

# Standardize A among positive subjects — use unique subjects, not TVC rows
# Transition subjects have 2 rows with A > 0; counting them twice biases mean/SD
pos_A <- unique(tvc_long[tvc_long$A > 0, c("id", "A")])$A
A_mean <- mean(pos_A)
A_sd   <- sd(pos_A)
tvc_long$A_z <- ifelse(tvc_long$A > 0, (tvc_long$A - A_mean) / A_sd, 0)
cat(sprintf("A standardization (unique subjects): mean = %.1f, SD = %.1f\n", A_mean, A_sd))

# Check degeneracy
all_positive <- all(tvc_long$Z_tv == 1)
cat(sprintf("All Z_tv=1? %s\n", all_positive))

if (all_positive) {
      cat("DEGENERATE: all Z_tv = 1. Fitting reduced model ~ A_z only.\n")
      fit_a3 <- coxph(Surv(tstart, tstop, event) ~ A_z + Sex_F + EDUC_z + apoe4,
                       data = tvc_long)
      cat("Cox model (reduced):\n")
      print(reportcox(fit_a3))

      hr_gamma <- extract_hr(fit_a3, "A_z", n_override = n_tvc_subj)
      hr_gamma$model <- "A3_tvc_interaction"
      hr_gamma$method <- "SILA"
      hr_gamma$parameter <- "gamma_A"
      hr_gamma$tracer <- "FBP"
      hr_results <- rbind(hr_results, hr_gamma)

      hr_results <- rbind(hr_results, data.frame(
            hr = NA, lower95 = NA, upper95 = NA, pvalue = NA,
            n = n_tvc_subj, nevent = summary(fit_a3)$nevent,
            model = "A3_tvc_interaction", method = "SILA",
            parameter = "beta_Z_tv", tracer = "FBP"
      ))
} else {
      fit_a3 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + A_z:Z_tv +
                             Sex_F + EDUC_z + apoe4,
                       data = tvc_long)
      cat("Cox model:\n")
      print(reportcox(fit_a3))

      # Beta (Z_tv)
      hr_beta <- extract_hr(fit_a3, "Z_tv", n_override = n_tvc_subj)
      hr_beta$model <- "A3_tvc_interaction"
      hr_beta$method <- "SILA"
      hr_beta$parameter <- "beta_Z_tv"
      hr_beta$tracer <- "FBP"
      hr_results <- rbind(hr_results, hr_beta)

      # Gamma (A:Z_tv)
      gamma_name <- setdiff(rownames(summary(fit_a3)$coefficients),
                            c("Z_tv", "Sex_F", "EDUC_z", "apoe4"))
      if (length(gamma_name) == 1) {
            hr_gamma <- extract_hr(fit_a3, gamma_name, n_override = n_tvc_subj)
      } else {
            hr_gamma <- data.frame(hr = NA, lower95 = NA, upper95 = NA,
                                    pvalue = NA, n = n_tvc_subj,
                                    nevent = summary(fit_a3)$nevent)
      }
      hr_gamma$model <- "A3_tvc_interaction"
      hr_gamma$method <- "SILA"
      hr_gamma$parameter <- "gamma_A"
      hr_gamma$tracer <- "FBP"
      hr_results <- rbind(hr_results, hr_gamma)

      # Correlation
      vc <- vcov(fit_a3)
      z_idx <- which(rownames(vc) == "Z_tv")
      g_idx <- which(rownames(vc) == gamma_name)
      if (length(z_idx) == 1 && length(g_idx) == 1) {
            cor_bg <- cov2cor(vc)[z_idx, g_idx]
            cat(sprintf("  vcov correlation (beta, gamma) = %.3f\n", cor_bg))
      }
}

cat("\n")


###############################################################################
## SECTION 5: Sensitivity — First Positive Scan as Z
###############################################################################

cat("================================================================\n")
cat("SENSITIVITY: FIRST POSITIVE SCAN AS Z\n")
cat("================================================================\n\n")

# A1-sensitivity: Countdown with first-positive-scan age
dat_fps_a1 <- dat_complete[!is.na(dat_complete$first_pos_age) &
                                  dat_complete$time_fps2DX > 0, ]

cat(sprintf("First-positive countdown: n = %d, events = %d\n",
            nrow(dat_fps_a1), sum(dat_fps_a1$d == 1)))

if (nrow(dat_fps_a1) >= 10 && sum(dat_fps_a1$d == 1) >= 5) {
      fit_fps_a1 <- coxph(Surv(time_fps2DX, d) ~ Sex_F + EDUC_z + apoe4 + fps_z,
                           data = dat_fps_a1)
      cat("Cox model (countdown, first-positive-scan):\n")
      print(reportcox(fit_fps_a1))

      hr_fps <- extract_hr(fit_fps_a1, "fps_z")
      hr_fps$model <- "A1_countdown"
      hr_fps$method <- "first_positive"
      hr_fps$parameter <- "fps_z"
      hr_fps$tracer <- "FBP"
      hr_results <- rbind(hr_results, hr_fps)
}

# TVC with first-positive-scan as Z
cat("\n--- TVC z_only with first-positive-scan ---\n")
tvc_fps_long <- create_tvc_data(dat_complete, "first_pos_age",
                                "baseline.age", "onset.age", "d")
covar_df <- data.frame(
      id = seq_len(nrow(dat_complete)),
      Sex_F = dat_complete$Sex_F,
      EDUC_z = dat_complete$EDUC_z,
      apoe4 = dat_complete$apoe4
)
tvc_fps_long <- merge(tvc_fps_long, covar_df, by = "id", all.x = TRUE)

n_fps_tvc_subj <- length(unique(tvc_fps_long$id))
cat(sprintf("TVC (first-positive-scan): %d rows (%d Z_tv=1, %d Z_tv=0), %d subjects\n",
            nrow(tvc_fps_long), sum(tvc_fps_long$Z_tv == 1),
            sum(tvc_fps_long$Z_tv == 0), n_fps_tvc_subj))

fit_fps_tvc <- coxph(Surv(tstart, tstop, event) ~ Z_tv + Sex_F + EDUC_z + apoe4,
                      data = tvc_fps_long)
cat("Cox model:\n")
print(reportcox(fit_fps_tvc))

hr_fps_tvc <- extract_hr(fit_fps_tvc, "Z_tv", n_override = n_fps_tvc_subj)
hr_fps_tvc$model <- "A2_tvc_z_only"
hr_fps_tvc$method <- "first_positive"
hr_fps_tvc$parameter <- "Z_tv"
hr_fps_tvc$tracer <- "FBP"
hr_results <- rbind(hr_results, hr_fps_tvc)

cat("\n")


###############################################################################
## SECTION 6: Degeneracy Summary
###############################################################################

cat("================================================================\n")
cat("DEGENERACY SUMMARY\n")
cat("================================================================\n\n")

A_adni     <- dat_complete$EAOA
entry_adni <- dat_complete$baseline.age
exit_adni  <- dat_complete$onset.age

# Use estpos to classify positivity, not NA in EAOA
pos_idx       <- which(dat_complete$estpos == 1)
n_valid       <- length(pos_idx)
n_already_pos <- sum(A_adni[pos_idx] <= entry_adni[pos_idx])
pct_already   <- if (n_valid > 0) n_already_pos / n_valid * 100 else NA
n_transition  <- sum(A_adni[pos_idx] > entry_adni[pos_idx] & A_adni[pos_idx] < exit_adni[pos_idx])
n_never_pos   <- nrow(dat_complete) - n_already_pos - n_transition

stopifnot("Degeneracy must sum to n_total" =
          n_already_pos + n_transition + n_never_pos == nrow(dat_complete))

cat(sprintf("  ADNI Amyloid (FBP):\n"))
cat(sprintf("    N total: %d\n", nrow(dat_complete)))
cat(sprintf("    N with valid EAOA: %d\n", n_valid))
cat(sprintf("    Already positive at entry: %d (%.1f%%)\n", n_already_pos, pct_already))
cat(sprintf("    Transition during follow-up: %d (%.1f%%)\n",
            n_transition, n_transition / nrow(dat_complete) * 100))
cat(sprintf("    Never positive: %d (%.1f%%)\n",
            n_never_pos, n_never_pos / nrow(dat_complete) * 100))

degen_adni <- data.frame(
      biomarker = "Amyloid_FBP", n_total = nrow(dat_complete),
      n_valid = n_valid, n_already_pos = n_already_pos,
      pct_already_pos = round(pct_already, 1),
      n_transition = n_transition, n_never_pos = n_never_pos,
      stringsAsFactors = FALSE
)

write.csv(degen_adni, file.path(res_dir, "ADNI_degeneracy_2026.csv"),
          row.names = FALSE)
cat("\nSaved: results/ADNI_degeneracy_2026.csv\n\n")


###############################################################################
## SECTION 7: Output Tables
###############################################################################

cat("================================================================\n")
cat("SAVING OUTPUT TABLES\n")
cat("================================================================\n\n")

# Format results
hr_out <- hr_results
hr_out$HR_95CI <- ifelse(is.na(hr_out$hr), "N/A",
                          sprintf("%.2f (%.2f, %.2f)", hr_out$hr,
                                  hr_out$lower95, hr_out$upper95))
hr_out$p_formatted <- ifelse(is.na(hr_out$pvalue), "N/A",
                              ifelse(hr_out$pvalue < 0.001, "< 0.001",
                                     ifelse(hr_out$pvalue < 0.01, "< 0.01",
                                            sprintf("%.3f", hr_out$pvalue))))

write.csv(hr_out, file.path(res_dir, "ADNI_countdown_vs_tvc_2026.csv"),
          row.names = FALSE)
cat("Saved: results/ADNI_countdown_vs_tvc_2026.csv\n")


###############################################################################
## SUMMARY PRINT
###############################################################################

cat("\n================================================================\n")
cat("SUMMARY: KEY RESULTS\n")
cat("================================================================\n\n")

cat("HR comparison (amyloid PET, FBP tracer):\n\n")
for (i in seq_len(nrow(hr_out))) {
      cat(sprintf("  %-45s  HR = %-22s  p = %-8s  (n=%s, events=%s)\n",
                  paste0("[", hr_out$method[i], "] ",
                         hr_out$model[i], " (", hr_out$parameter[i], ")"),
                  hr_out$HR_95CI[i],
                  hr_out$p_formatted[i],
                  hr_out$n[i],
                  hr_out$nevent[i]))
}

cat("\n=== ADNI SILA re-analysis v2 (2026 data) complete ===\n")
