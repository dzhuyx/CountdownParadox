###############################################################################
## ADNI Plasma Re-Analysis: Countdown Paradox — Trajectory-Only
##
## Purpose: Apply countdown, TVC, and landmark analyses to ADNI plasma p-tau217
##   data using trajectory-estimated age-at-positivity (Z) only.
##
##   Model P1: Countdown — Surv(T-Z, d) ~ Z_std + covariates
##   Model P2: TVC z_only — Surv(tstart, tstop, event) ~ Z_tv + covariates
##   Model P3: TVC interaction — Surv(tstart, tstop, event) ~ Z_tv + A:Z_tv + covariates
##   Model P4: Landmark continuous — Surv(T-L, event_L) ~ Z_std + covariates
##   Model P5: Landmark binary — Surv(T-L, event_L) ~ positive_by_L + covariates
##
##   4 Configs (trajectory-only):
##     1. Fujirebio_SILA (val0=0.300)
##     2. Fujirebio_TIRA (tip_point=0.300)
##     3. C2N_SILA (val0=4.06)
##     4. C2N_TIRA (tip_point=4.06)
##
## Input:  data/ADNI_SILA_intermediate_2026.rda  (dat_surv_cn, dat_demo, dat_apoe)
##         data/ADNI_plasma_SILA_intermediate.rda (SILA trajectory Z: Fuji + C2N)
##         data/ADNI_plasma_TIRA_intermediate.rda (TIRA trajectory Z: Fuji + C2N)
## Output: results/ADNI_plasma_all_models.csv
##         results/ADNI_plasma_degeneracy.csv
##         results/ADNI_plasma_positivity_summary.csv
##
## Usage:  source() from reproducibility/ directory or via run_all.R
##
## Author: Yuxin Zhu
## Date: March 2026 (reproducibility copy April 2026)
###############################################################################

rm(list = ls())
library(survival)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set it to <project>/CountdownParadox_Analysis.")
data_dir     <- file.path(project_root, "data")
data_2026    <- file.path(project_root, "ADNI_2026_data")
res_dir      <- file.path(project_root, "results")
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# Landmark times
LANDMARK_TIMES <- c(70, 75, 80)


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
                        # (SILA predicted positivity but couldn't compute valid EAOA)
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
## Run all five models for a given trajectory config
###############################################################################

run_all_models <- function(dat, config_label, method_label) {
      #
      # dat must have columns:
      #   age_at_positivity  = Z (NA if never positive)
      #   baseline.age, censor.age, d  = survival
      #   Sex_F, EDUC_z, apoe4  = covariates
      #
      cat(sprintf("\n================================================================\n"))
      cat(sprintf("CONFIG: %s | METHOD: %s\n", config_label, method_label))
      cat(sprintf("================================================================\n\n"))

      covariates <- c("Sex_F", "EDUC_z", "apoe4")
      dat_c <- dat[complete.cases(dat[, c(covariates, "baseline.age",
                                          "censor.age", "d")]), ]
      # Exclude subjects with zero follow-up (baseline.age >= censor.age)
      n_zero_fu <- sum(dat_c$baseline.age >= dat_c$censor.age)
      if (n_zero_fu > 0) {
            cat(sprintf("Excluding %d subject(s) with zero follow-up\n", n_zero_fu))
            dat_c <- dat_c[which(dat_c$baseline.age < dat_c$censor.age), ]
      }

      n_total <- nrow(dat_c)
      # Use ever_positive for positivity classification
      n_pos   <- sum(dat_c$ever_positive == 1, na.rm = TRUE)
      n_neg   <- sum(dat_c$ever_positive == 0, na.rm = TRUE)

      cat(sprintf("Sample: %d total (%d positive, %d never positive)\n",
                  n_total, n_pos, n_neg))
      cat(sprintf("Events: %d (%.1f%%)\n", sum(dat_c$d == 1),
                  mean(dat_c$d == 1) * 100))

      results <- data.frame()

      # ---- Model P1: Countdown ---- #
      cat("\n--- Model P1: Countdown (T-Z ~ Z) ---\n")

      dat_c$time_to_dx <- dat_c$censor.age - dat_c$age_at_positivity
      # Use ever_positive instead of !is.na(age_at_positivity)
      dat_p1 <- dat_c[which(dat_c$ever_positive == 1 &
                             dat_c$time_to_dx > 0), ]

      if (nrow(dat_p1) >= 10 && sum(dat_p1$d == 1) >= 5) {
            Z_mean <- mean(dat_p1$age_at_positivity)
            Z_sd   <- sd(dat_p1$age_at_positivity)
            dat_p1$Z_std <- (dat_p1$age_at_positivity - Z_mean) / Z_sd
            cat(sprintf("Z standardization: mean = %.1f, SD = %.1f\n", Z_mean, Z_sd))
            cat(sprintf("Countdown sample: n = %d, events = %d\n",
                        nrow(dat_p1), sum(dat_p1$d == 1)))

            fit_p1 <- coxph(Surv(time_to_dx, d) ~ Sex_F + EDUC_z + apoe4 + Z_std,
                            data = dat_p1)
            print(reportcox(fit_p1))

            hr_p1 <- extract_hr(fit_p1, "Z_std")
            hr_p1$model <- "P1_countdown"
            hr_p1$parameter <- "Z_std"
            hr_p1$config <- config_label
            hr_p1$method <- method_label
            hr_p1$landmark_time <- NA
            results <- rbind(results, hr_p1)
      } else {
            cat(sprintf("SKIPPED: insufficient sample (n=%d, events=%d)\n",
                        nrow(dat_p1), sum(dat_p1$d == 1)))
            results <- rbind(results, data.frame(
                  hr = NA, lower95 = NA, upper95 = NA, pvalue = NA,
                  n = nrow(dat_p1), nevent = sum(dat_p1$d == 1),
                  model = "P1_countdown", parameter = "Z_std",
                  config = config_label, method = method_label,
                  landmark_time = NA
            ))
      }

      # ---- Model P2: TVC z_only ---- #
      cat("\n--- Model P2: TVC z_only ---\n")

      tvc_long <- create_tvc_data(dat_c, "age_at_positivity",
                                  "baseline.age", "censor.age", "d",
                                  estpos_col = "ever_positive")
      covar_df <- data.frame(
            id = seq_len(nrow(dat_c)),
            Sex_F = dat_c$Sex_F,
            EDUC_z = dat_c$EDUC_z,
            apoe4 = dat_c$apoe4
      )
      tvc_long <- merge(tvc_long, covar_df, by = "id", all.x = TRUE)

      n_pos_rows <- sum(tvc_long$Z_tv == 1)
      n_neg_rows <- sum(tvc_long$Z_tv == 0)
      cat(sprintf("TVC data: %d rows (%d Z_tv=1, %d Z_tv=0), %d subjects\n",
                  nrow(tvc_long), n_pos_rows, n_neg_rows,
                  length(unique(tvc_long$id))))

      # ---- Tier 1 assertions ---- #
      n_tvc_subj <- length(unique(tvc_long$id))
      stopifnot("TVC subject count must equal analysis N" =
                n_tvc_subj == nrow(dat_c))
      estpos_0_ids <- which(dat_c$ever_positive == 0)
      if (length(estpos_0_ids) > 0) {
            stopifnot("estpos=0 subjects must have Z_tv=0" =
                      all(tvc_long$Z_tv[tvc_long$id %in% estpos_0_ids] == 0))
            stopifnot("estpos=0 subjects must have A=0 in TVC" =
                      all(tvc_long$A[tvc_long$id %in% estpos_0_ids] == 0))
      }
      cat("  Assertions passed: subject count, estpos↔Z_tv concordance\n")

      fit_p2 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + Sex_F + EDUC_z + apoe4,
                       data = tvc_long)
      print(reportcox(fit_p2))

      hr_p2 <- extract_hr(fit_p2, "Z_tv", n_override = n_tvc_subj)
      hr_p2$model <- "P2_tvc_z_only"
      hr_p2$parameter <- "Z_tv"
      hr_p2$config <- config_label
      hr_p2$method <- method_label
      hr_p2$landmark_time <- NA
      results <- rbind(results, hr_p2)

      # ---- Model P3: TVC interaction ---- #
      cat("\n--- Model P3: TVC interaction ---\n")

      # Use unique subjects, not TVC rows — transition subjects have 2 rows
      pos_A <- unique(tvc_long[tvc_long$A > 0, c("id", "A")])$A
      A_mean <- mean(pos_A)
      A_sd   <- sd(pos_A)
      tvc_long$A_z <- ifelse(tvc_long$A > 0, (tvc_long$A - A_mean) / A_sd, 0)
      cat(sprintf("A standardization (unique subjects): mean = %.1f, SD = %.1f\n", A_mean, A_sd))

      all_positive <- all(tvc_long$Z_tv == 1)
      cat(sprintf("All Z_tv=1? %s\n", all_positive))

      if (all_positive) {
            cat("DEGENERATE: all Z_tv = 1. Fitting reduced model ~ A_z only.\n")
            fit_p3 <- coxph(Surv(tstart, tstop, event) ~ A_z + Sex_F + EDUC_z + apoe4,
                            data = tvc_long)
            print(reportcox(fit_p3))

            hr_gamma <- extract_hr(fit_p3, "A_z", n_override = n_tvc_subj)
            hr_gamma$model <- "P3_tvc_interaction"
            hr_gamma$parameter <- "gamma_A"
            hr_gamma$config <- config_label
            hr_gamma$method <- method_label
            hr_gamma$landmark_time <- NA
            results <- rbind(results, hr_gamma)

            results <- rbind(results, data.frame(
                  hr = NA, lower95 = NA, upper95 = NA, pvalue = NA,
                  n = n_tvc_subj, nevent = summary(fit_p3)$nevent,
                  model = "P3_tvc_interaction", parameter = "beta_Z_tv",
                  config = config_label, method = method_label,
                  landmark_time = NA
            ))
      } else {
            fit_p3 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + A_z:Z_tv +
                                   Sex_F + EDUC_z + apoe4,
                            data = tvc_long)
            print(reportcox(fit_p3))

            hr_beta <- extract_hr(fit_p3, "Z_tv", n_override = n_tvc_subj)
            hr_beta$model <- "P3_tvc_interaction"
            hr_beta$parameter <- "beta_Z_tv"
            hr_beta$config <- config_label
            hr_beta$method <- method_label
            hr_beta$landmark_time <- NA
            results <- rbind(results, hr_beta)

            gamma_name <- setdiff(rownames(summary(fit_p3)$coefficients),
                                  c("Z_tv", "Sex_F", "EDUC_z", "apoe4"))
            if (length(gamma_name) == 1) {
                  hr_gamma <- extract_hr(fit_p3, gamma_name, n_override = n_tvc_subj)
            } else {
                  hr_gamma <- data.frame(hr = NA, lower95 = NA, upper95 = NA,
                                          pvalue = NA, n = n_tvc_subj,
                                          nevent = summary(fit_p3)$nevent)
            }
            hr_gamma$model <- "P3_tvc_interaction"
            hr_gamma$parameter <- "gamma_A"
            hr_gamma$config <- config_label
            hr_gamma$method <- method_label
            hr_gamma$landmark_time <- NA
            results <- rbind(results, hr_gamma)

            vc <- vcov(fit_p3)
            z_idx <- which(rownames(vc) == "Z_tv")
            g_idx <- which(rownames(vc) == gamma_name)
            if (length(z_idx) == 1 && length(g_idx) == 1) {
                  cor_bg <- cov2cor(vc)[z_idx, g_idx]
                  cat(sprintf("  vcov correlation (beta, gamma) = %.3f\n", cor_bg))
            }
      }

      # ---- Model P4: Landmark continuous (at L = 70, 75, 80) ---- #
      for (L in LANDMARK_TIMES) {
            cat(sprintf("\n--- Model P4: Landmark continuous (L=%d) ---\n", L))

            # Restrict to: positive by L (Z < L), event-free at L
            # Use ever_positive instead of !is.na(age_at_positivity)
            lm_dat <- dat_c[which(dat_c$ever_positive == 1 &
                                  dat_c$age_at_positivity < L &
                                  pmin(dat_c$censor.age,
                                       ifelse(dat_c$d == 1, dat_c$censor.age, Inf)) > L), ]
            lm_dat$time_from_L <- lm_dat$censor.age - L
            lm_dat$event_L <- as.integer(lm_dat$d == 1 & lm_dat$censor.age > L)
            lm_dat <- lm_dat[lm_dat$time_from_L > 0, ]

            if (nrow(lm_dat) >= 10 && sum(lm_dat$event_L) >= 5) {
                  Z_mean_L <- mean(lm_dat$age_at_positivity)
                  Z_sd_L   <- sd(lm_dat$age_at_positivity)
                  lm_dat$Z_std <- (lm_dat$age_at_positivity - Z_mean_L) / Z_sd_L
                  cat(sprintf("  Sample: n = %d, events = %d\n",
                              nrow(lm_dat), sum(lm_dat$event_L)))

                  fit_p4 <- tryCatch(
                        coxph(Surv(time_from_L, event_L) ~ Z_std + Sex_F + EDUC_z + apoe4,
                              data = lm_dat),
                        error = function(e) NULL
                  )

                  if (!is.null(fit_p4)) {
                        print(reportcox(fit_p4))
                        hr_p4 <- extract_hr(fit_p4, "Z_std")
                  } else {
                        hr_p4 <- data.frame(hr = NA, lower95 = NA, upper95 = NA,
                                            pvalue = NA, n = nrow(lm_dat),
                                            nevent = sum(lm_dat$event_L))
                  }
            } else {
                  cat(sprintf("  SKIPPED: insufficient (n=%d, events=%d)\n",
                              nrow(lm_dat), sum(lm_dat$event_L)))
                  hr_p4 <- data.frame(hr = NA, lower95 = NA, upper95 = NA,
                                      pvalue = NA, n = nrow(lm_dat),
                                      nevent = sum(lm_dat$event_L))
            }
            hr_p4$model <- "P4_landmark_continuous"
            hr_p4$parameter <- "Z_std"
            hr_p4$config <- config_label
            hr_p4$method <- method_label
            hr_p4$landmark_time <- L
            results <- rbind(results, hr_p4)
      }

      # ---- Model P5: Landmark binary (at L = 70, 75, 80) ---- #
      for (L in LANDMARK_TIMES) {
            cat(sprintf("\n--- Model P5: Landmark binary (L=%d) ---\n", L))

            # All event-free at L; binary indicator: positive_by_L = (Z < L)
            lm_dat <- dat_c[pmin(dat_c$censor.age,
                                 ifelse(dat_c$d == 1, dat_c$censor.age, Inf)) > L, ]
            # Use ever_positive instead of !is.na(age_at_positivity)
            lm_dat$positive_by_L <- as.integer(!is.na(lm_dat$ever_positive) &
                                                      lm_dat$ever_positive == 1 &
                                                      lm_dat$age_at_positivity < L)
            lm_dat$time_from_L <- lm_dat$censor.age - L
            lm_dat$event_L <- as.integer(lm_dat$d == 1 & lm_dat$censor.age > L)
            lm_dat <- lm_dat[lm_dat$time_from_L > 0, ]

            n_pos_L <- sum(lm_dat$positive_by_L == 1)
            n_neg_L <- sum(lm_dat$positive_by_L == 0)
            cat(sprintf("  Sample: n = %d (pos=%d, neg=%d), events = %d\n",
                        nrow(lm_dat), n_pos_L, n_neg_L, sum(lm_dat$event_L)))

            if (nrow(lm_dat) >= 10 && n_pos_L >= 3 && n_neg_L >= 3 &&
                sum(lm_dat$event_L) >= 5) {

                  fit_p5 <- tryCatch(
                        coxph(Surv(time_from_L, event_L) ~ positive_by_L +
                                    Sex_F + EDUC_z + apoe4,
                              data = lm_dat),
                        error = function(e) NULL
                  )

                  if (!is.null(fit_p5)) {
                        print(reportcox(fit_p5))
                        hr_p5 <- extract_hr(fit_p5, "positive_by_L")
                  } else {
                        hr_p5 <- data.frame(hr = NA, lower95 = NA, upper95 = NA,
                                            pvalue = NA, n = nrow(lm_dat),
                                            nevent = sum(lm_dat$event_L))
                  }
            } else {
                  cat(sprintf("  SKIPPED: insufficient\n"))
                  hr_p5 <- data.frame(hr = NA, lower95 = NA, upper95 = NA,
                                      pvalue = NA, n = nrow(lm_dat),
                                      nevent = sum(lm_dat$event_L))
            }
            hr_p5$model <- "P5_landmark_binary"
            hr_p5$parameter <- "positive_by_L"
            hr_p5$config <- config_label
            hr_p5$method <- method_label
            hr_p5$landmark_time <- L
            results <- rbind(results, hr_p5)
      }

      # ---- Degeneracy summary ---- #
      A_vals    <- dat_c$age_at_positivity
      entry     <- dat_c$baseline.age
      exit      <- dat_c$censor.age

      # Use ever_positive for positivity classification
      pos_idx       <- which(dat_c$ever_positive == 1)
      n_valid       <- length(pos_idx)
      n_already_pos <- sum(A_vals[pos_idx] <= entry[pos_idx])
      pct_already   <- if (n_valid > 0) n_already_pos / n_valid * 100 else NA
      n_transition  <- sum(A_vals[pos_idx] > entry[pos_idx] & A_vals[pos_idx] < exit[pos_idx])
      n_never_pos   <- n_total - n_already_pos - n_transition

      stopifnot("Degeneracy must sum to n_total" =
                n_already_pos + n_transition + n_never_pos == n_total)

      cat(sprintf("\nDegeneracy: valid=%d, already_pos=%d (%.0f%%), transition=%d, never_pos=%d\n",
                  n_valid, n_already_pos, pct_already, n_transition, n_never_pos))

      degen <- data.frame(
            config = config_label, method = method_label,
            n_total = n_total, n_valid = n_valid,
            n_already_pos = n_already_pos,
            pct_already_pos = round(pct_already, 1),
            n_transition = n_transition, n_never_pos = n_never_pos,
            stringsAsFactors = FALSE
      )

      return(list(results = results, degeneracy = degen))
}


###############################################################################
## SECTION 1: Load Data
###############################################################################

cat("================================================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("================================================================\n\n")

# --- 1a: Survival outcome and demographics from SILA intermediate --- #
load(file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"))
cat(sprintf("Loaded SILA intermediate: dat_surv_cn (%d subjects)\n", nrow(dat_surv_cn)))

# --- 1b: Demographics --- #
demo_file <- file.path(data_2026, "PTDEMOG_11Feb2026.csv")
dat_demo <- read.csv(demo_file, na.strings = -4)
dat_demo$DOB <- as.Date(paste0("01/", dat_demo$PTDOB), "%d/%m/%Y")
dat_demo$update_stamp <- as.Date(dat_demo$update_stamp, "%Y-%m-%d")
dat_demo <- do.call(rbind, lapply(split(dat_demo, dat_demo$RID), function(x) {
      x[x$update_stamp == max(x$update_stamp, na.rm = TRUE), ][1, ]
}))
cat(sprintf("Demographics: %d subjects\n", nrow(dat_demo)))

# --- 1c: APOE --- #
apoe_file <- file.path(data_2026, "All_Subjects_APOERES_15Feb2026.csv")
dat_apoe <- read.csv(apoe_file, na.strings = -4)
dat_apoe$apoe4 <- as.integer(grepl("4", dat_apoe$GENOTYPE))
dat_apoe <- do.call(rbind, lapply(split(dat_apoe, dat_apoe$RID), function(x) x[1, ]))
cat(sprintf("APOE: %d subjects\n", nrow(dat_apoe)))

# --- 1d: SILA trajectory-based Z (Fujirebio + C2N from ADNI_plasma_sila.R) --- #
sila_plasma_file <- file.path(data_dir, "ADNI_plasma_SILA_intermediate.rda")
sila_fuj_z_df <- NULL
sila_c2n_z_df <- NULL
if (file.exists(sila_plasma_file)) {
      sila_env <- new.env()
      load(sila_plasma_file, envir = sila_env)

      # Fujirebio SILA
      if (exists("resfit_last_fuj", envir = sila_env)) {
            rflp <- get("resfit_last_fuj", envir = sila_env)
            eaoa_col <- if ("EAOA_plasma" %in% names(rflp)) "EAOA_plasma"
                        else if ("estaget0" %in% names(rflp)) "estaget0"
                        else if ("estage0" %in% names(rflp)) "estage0"
                        else NULL
            estpos_col <- if ("estpos" %in% names(rflp)) "estpos" else NULL
            if (!is.null(eaoa_col)) {
                  sila_fuj_z_df <- data.frame(
                        RID = rflp$RID,
                        age_at_positivity = rflp[[eaoa_col]],
                        ever_positive = if (!is.null(estpos_col)) rflp[[estpos_col]] else !is.na(rflp[[eaoa_col]])
                  )
                  # Keep age_at_positivity for all subjects; ever_positive used by create_tvc_data
                  cat(sprintf("SILA Fujirebio Z: %d subjects (%d positive)\n",
                              nrow(sila_fuj_z_df), sum(sila_fuj_z_df$ever_positive, na.rm = TRUE)))
            }
      }

      # C2N SILA
      if (exists("resfit_last_c2n", envir = sila_env)) {
            rlc <- get("resfit_last_c2n", envir = sila_env)
            if (!is.null(rlc)) {
                  eaoa_col_c <- if ("EAOA_plasma" %in% names(rlc)) "EAOA_plasma"
                                else if ("estaget0" %in% names(rlc)) "estaget0"
                                else if ("estage0" %in% names(rlc)) "estage0"
                                else NULL
                  estpos_col_c <- if ("estpos" %in% names(rlc)) "estpos" else NULL
                  if (!is.null(eaoa_col_c)) {
                        sila_c2n_z_df <- data.frame(
                              RID = rlc$RID,
                              age_at_positivity = rlc[[eaoa_col_c]],
                              ever_positive = if (!is.null(estpos_col_c)) rlc[[estpos_col_c]] else !is.na(rlc[[eaoa_col_c]])
                        )
                        # Keep age_at_positivity for all subjects
                        cat(sprintf("SILA C2N Z: %d subjects (%d positive)\n",
                                    nrow(sila_c2n_z_df), sum(sila_c2n_z_df$ever_positive, na.rm = TRUE)))
                  } else {
                        cat("WARNING: C2N SILA intermediate has no EAOA column.\n")
                  }
            } else {
                  cat("C2N SILA result is NULL (convergence failure).\n")
            }
      } else {
            cat("C2N SILA results not found in intermediate (may not have been run).\n")
      }

      # Backward compatibility: old format had resfit_last_plasma (Fuji only)
      if (is.null(sila_fuj_z_df) && exists("resfit_last_plasma", envir = sila_env)) {
            rflp <- get("resfit_last_plasma", envir = sila_env)
            eaoa_col <- if ("EAOA_plasma" %in% names(rflp)) "EAOA_plasma"
                        else if ("estaget0" %in% names(rflp)) "estaget0"
                        else if ("estage0" %in% names(rflp)) "estage0"
                        else NULL
            estpos_col <- if ("estpos" %in% names(rflp)) "estpos" else NULL
            if (!is.null(eaoa_col)) {
                  sila_fuj_z_df <- data.frame(
                        RID = rflp$RID,
                        age_at_positivity = rflp[[eaoa_col]],
                        ever_positive = if (!is.null(estpos_col)) rflp[[estpos_col]] else !is.na(rflp[[eaoa_col]])
                  )
                  # Keep age_at_positivity for all subjects
                  cat(sprintf("SILA Fujirebio Z (legacy format): %d subjects (%d positive)\n",
                              nrow(sila_fuj_z_df), sum(sila_fuj_z_df$ever_positive, na.rm = TRUE)))
            }
      }

      rm(sila_env)
} else {
      cat("SILA plasma intermediate not found — SILA configs will be skipped.\n")
}

# --- 1e: TIRA trajectory-based Z (from ADNI_plasma_tira.R) --- #
tira_plasma_file <- file.path(data_dir, "ADNI_plasma_TIRA_intermediate.rda")
tira_fuj_z_df <- NULL
tira_c2n_z_df <- NULL
if (file.exists(tira_plasma_file)) {
      tira_env <- new.env()
      load(tira_plasma_file, envir = tira_env)

      # Fujirebio TIRA
      if (exists("tira_fuj_result", envir = tira_env)) {
            tfr <- get("tira_fuj_result", envir = tira_env)
            if (!is.null(tfr) && "est_onset_age" %in% names(tfr)) {
                  tira_fuj_z_df <- data.frame(
                        RID = tfr$RID,
                        age_at_positivity = tfr$est_onset_age,
                        ever_positive = !is.na(tfr$est_onset_age)
                  )
                  cat(sprintf("TIRA Fujirebio Z: %d subjects (%d with est_onset_age)\n",
                              nrow(tira_fuj_z_df),
                              sum(!is.na(tira_fuj_z_df$age_at_positivity))))
            }
      }

      # C2N TIRA
      if (exists("tira_c2n_result", envir = tira_env)) {
            tcr <- get("tira_c2n_result", envir = tira_env)
            if (!is.null(tcr) && "est_onset_age" %in% names(tcr)) {
                  tira_c2n_z_df <- data.frame(
                        RID = tcr$RID,
                        age_at_positivity = tcr$est_onset_age,
                        ever_positive = !is.na(tcr$est_onset_age)
                  )
                  cat(sprintf("TIRA C2N Z: %d subjects (%d with est_onset_age)\n",
                              nrow(tira_c2n_z_df),
                              sum(!is.na(tira_c2n_z_df$age_at_positivity))))
            }
      }
      rm(tira_env)
} else {
      cat("TIRA plasma intermediate not found — TIRA configs will be skipped.\n")
}


###############################################################################
## SECTION 2: Positivity Summary
###############################################################################

cat("\n================================================================\n")
cat("SECTION 2: POSITIVITY SUMMARY\n")
cat("================================================================\n\n")

pos_summary <- data.frame()

add_pos <- function(df, label, method_name) {
      if (is.null(df)) return(NULL)
      data.frame(
            config = label, method = method_name,
            n_total = nrow(df),
            n_positive = sum(df$ever_positive, na.rm = TRUE),
            pct_positive = round(mean(df$ever_positive, na.rm = TRUE) * 100, 1),
            stringsAsFactors = FALSE
      )
}

pos_summary <- rbind(pos_summary,
      add_pos(sila_fuj_z_df, "Fujirebio_SILA", "SILA"),
      add_pos(tira_fuj_z_df, "Fujirebio_TIRA", "TIRA"),
      add_pos(sila_c2n_z_df, "C2N_SILA", "SILA"),
      add_pos(tira_c2n_z_df, "C2N_TIRA", "TIRA")
)

cat("Positivity summary:\n")
print(pos_summary)

write.csv(pos_summary, file.path(res_dir, "ADNI_plasma_positivity_summary.csv"),
          row.names = FALSE)
cat("Saved: results/ADNI_plasma_positivity_summary.csv\n")


###############################################################################
## SECTION 3: Merge Z with Survival Data
###############################################################################

cat("\n================================================================\n")
cat("SECTION 3: MERGE Z WITH SURVIVAL DATA\n")
cat("================================================================\n\n")

# Base survival data from dat_surv_cn (CN at baseline)
# dat_surv_cn already has censor.age (identical to onset.age); drop onset.age to avoid duplication
dat_base <- dat_surv_cn
dat_base$onset.age <- NULL

# Add demographics
dat_base <- merge(dat_base, dat_demo[, c("RID", "PTGENDER", "PTEDUCAT")],
                  by = "RID", all.x = TRUE)
dat_base$Sex_F <- as.integer(dat_base$PTGENDER == 2)
dat_base$EDUC_z <- scale(dat_base$PTEDUCAT)[, 1]

# Add APOE4
dat_base <- merge(dat_base, dat_apoe[, c("RID", "apoe4")],
                  by = "RID", all.x = TRUE)

covariates <- c("Sex_F", "EDUC_z", "apoe4")
dat_base <- dat_base[complete.cases(dat_base[, covariates]), ]

# Exclude 78 participants removed by ADNI (data quality notice, March 2026)
adni_exclude_rids <- unique(dat_demo$RID[grepl("^381_S_", dat_demo$PTID)])
n_excluded <- sum(dat_base$RID %in% adni_exclude_rids)
dat_base <- dat_base[!dat_base$RID %in% adni_exclude_rids, ]
cat(sprintf("Excluded %d of %d flagged 381_S subjects\n", n_excluded, length(adni_exclude_rids)))

cat(sprintf("CN subjects with complete covariates: %d\n", nrow(dat_base)))
cat(sprintf("  Events: %d (%.1f%%)\n\n", sum(dat_base$d == 1),
            mean(dat_base$d == 1) * 100))

# Function to merge Z into base dataset — restrict to trajectory-eligible subjects
# Include ever_positive for estpos-based classification in create_tvc_data
merge_Z <- function(base, z_df) {
      merged <- merge(base, z_df[, c("RID", "age_at_positivity", "ever_positive")],
                      by = "RID", all.x = FALSE)
      merged
}

# Save the plasma analysis cohort (n = 801, Fujirebio SILA) for manuscript
# descriptives. Reconstructs exactly the cohort that
# run_all_models uses internally for the primary Fujirebio_SILA config:
# dat_base -> trajectory-eligible (merge_Z) -> complete covariates + nonzero
# follow-up. This is a distinct subject set from the PET cohort.
if (exists("sila_fuj_z_df") && !is.null(sila_fuj_z_df)) {
      dat_plasma_cohort <- merge_Z(dat_base, sila_fuj_z_df)
      dat_plasma_cohort <- dat_plasma_cohort[complete.cases(dat_plasma_cohort[, c("Sex_F", "EDUC_z", "apoe4", "baseline.age", "censor.age", "d")]), ]
      dat_plasma_cohort <- dat_plasma_cohort[which(dat_plasma_cohort$baseline.age < dat_plasma_cohort$censor.age), ]
      dat_plasma_cohort$onset.age <- dat_plasma_cohort$censor.age  # ADNI: event age = censor age
      save(dat_plasma_cohort, file = file.path(data_dir, "ADNI_plasma_analysis_cohort.rda"))
      cat(sprintf("Saved plasma analysis cohort (%d subjects) -> data/ADNI_plasma_analysis_cohort.rda\n",
                  nrow(dat_plasma_cohort)))
}


###############################################################################
## SECTION 4: Run Analyses — 4 Trajectory Configs × 5 Models
###############################################################################

cat("================================================================\n")
cat("SECTION 4: RUN ANALYSES\n")
cat("================================================================\n\n")

all_results  <- data.frame()
all_degen    <- data.frame()

# Define 4 trajectory-only configurations
configs <- list()

if (!is.null(sila_fuj_z_df)) {
      configs[[length(configs) + 1]] <- list(
            z_df = sila_fuj_z_df,
            label = "Fujirebio_SILA",
            method = "SILA"
      )
      cat("Config 1: Fujirebio SILA\n")
}
if (!is.null(tira_fuj_z_df)) {
      configs[[length(configs) + 1]] <- list(
            z_df = tira_fuj_z_df,
            label = "Fujirebio_TIRA",
            method = "TIRA"
      )
      cat("Config 2: Fujirebio TIRA\n")
}
if (!is.null(sila_c2n_z_df)) {
      configs[[length(configs) + 1]] <- list(
            z_df = sila_c2n_z_df,
            label = "C2N_SILA",
            method = "SILA"
      )
      cat("Config 3: C2N SILA\n")
}
if (!is.null(tira_c2n_z_df)) {
      configs[[length(configs) + 1]] <- list(
            z_df = tira_c2n_z_df,
            label = "C2N_TIRA",
            method = "TIRA"
      )
      cat("Config 4: C2N TIRA\n")
}

cat(sprintf("\nTotal configs: %d\n", length(configs)))

for (cfg in configs) {
      dat_merged <- merge_Z(dat_base, cfg$z_df)

      n_positive_in_cohort <- sum(dat_merged$ever_positive == 1, na.rm = TRUE)
      n_never_pos <- sum(dat_merged$ever_positive == 0, na.rm = TRUE)
      cat(sprintf("\n--- %s ---\n", cfg$label))
      cat(sprintf("Trajectory-eligible subjects: %d (restricted from %d CN)\n",
                  nrow(dat_merged), nrow(dat_base)))
      stopifnot("All subjects must have ever_positive defined after restriction" =
                all(!is.na(dat_merged$ever_positive)))
      cat(sprintf("Positive: %d, Never positive: %d\n", n_positive_in_cohort, n_never_pos))

      if (n_positive_in_cohort >= 10) {
            out <- run_all_models(dat_merged, cfg$label, cfg$method)
            all_results <- rbind(all_results, out$results)
            all_degen   <- rbind(all_degen, out$degeneracy)
      } else {
            cat(sprintf("SKIPPED: only %d positive subjects in CN cohort\n",
                        n_positive_in_cohort))
      }
}


###############################################################################
## SECTION 5: Output Tables
###############################################################################

cat("\n================================================================\n")
cat("SECTION 5: SAVING OUTPUT TABLES\n")
cat("================================================================\n\n")

# Format results
if (nrow(all_results) > 0) {
      all_results$HR_95CI <- ifelse(is.na(all_results$hr), "N/A",
                                     sprintf("%.2f (%.2f, %.2f)", all_results$hr,
                                             all_results$lower95, all_results$upper95))
      all_results$p_formatted <- ifelse(is.na(all_results$pvalue), "N/A",
                                         ifelse(all_results$pvalue < 0.001, "< 0.001",
                                                ifelse(all_results$pvalue < 0.01, "< 0.01",
                                                       sprintf("%.3f", all_results$pvalue))))

      write.csv(all_results, file.path(res_dir, "ADNI_plasma_all_models.csv"),
                row.names = FALSE)
      cat("Saved: results/ADNI_plasma_all_models.csv\n")
}

if (nrow(all_degen) > 0) {
      write.csv(all_degen, file.path(res_dir, "ADNI_plasma_degeneracy.csv"),
                row.names = FALSE)
      cat("Saved: results/ADNI_plasma_degeneracy.csv\n")
}


###############################################################################
## SUMMARY
###############################################################################

cat("\n================================================================\n")
cat("SUMMARY: KEY RESULTS\n")
cat("================================================================\n\n")

if (nrow(all_results) > 0) {
      for (i in seq_len(nrow(all_results))) {
            lm_str <- if (!is.na(all_results$landmark_time[i]))
                  sprintf(" [L=%d]", all_results$landmark_time[i]) else ""
            cat(sprintf("  %-18s %-15s %-22s  HR = %-22s  p = %-8s  (n=%s, events=%s)%s\n",
                        all_results$config[i],
                        all_results$method[i],
                        paste0(all_results$model[i], " (", all_results$parameter[i], ")"),
                        all_results$HR_95CI[i],
                        all_results$p_formatted[i],
                        all_results$n[i],
                        all_results$nevent[i],
                        lm_str))
      }
}

cat(sprintf("\nTotal model fits: %d\n", nrow(all_results)))
cat(sprintf("Configs: %d\n", length(configs)))
cat("\n=== ADNI Plasma Re-Analysis complete ===\n")
