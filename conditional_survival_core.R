###############################################################################
## conditional_survival_core.R
##
## SHARED CORE (Option A). In this repository it is sourced by the analyses added
## in the A&D revision (R1_S6_beyond_curve.R, R1_S9_cause_specific.R), so they start
## from the SAME audited fits. Repository copy of CountdownParadox_Analysis/
## conditional_survival_core.R: only project_root is changed (read from CP_PROJECT_ROOT).
##
## Provides, in the caller's environment:
##   - paths (project_root, data_dir, results_dir), `covariates`, `COEF_TOL`
##   - helpers: extract_all_hr, create_tvc_data, coefn, inter_coef,
##              run_all_models_capture, make_stepH, evalH,
##              cond_F, solve_median, eta_std, eta_tvbc, eta_tvaabc, sh_of
##   - `fits`  : per-biomarker captured Cox fits (P1/P2/P3), standardization
##               constants, aabc_pool, AND the subject-level analysis frame
##               (`$dat`, `$cols`) used to fit them (needed for calibration cohorts)
##   - validation gate G0 (refit coefficients vs saved CSV) — HALTS on mismatch
##
## This file does NOT rm(list=ls()) and does NOT save anything or make figures.
## The 15 fits are IDENTICAL to extract_full_coefficients.R.
##
## Author: Yuxin Zhu
###############################################################################

suppressPackageStartupMessages(library(survival))

project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set CP_PROJECT_ROOT to the analysis root (the folder that holds data/ and results/).")
data_dir     <- file.path(project_root, "data")
results_dir  <- file.path(project_root, "results")

covariates <- c("Sex_F", "EDUC_z", "apoe4")
COEF_TOL   <- 1e-5        # tolerance for refit-vs-saved coefficient validation

cat("[core] loading shared fitting core...\n")


###############################################################################
## Helper functions (verbatim from extract_full_coefficients.R)
###############################################################################

extract_all_hr <- function(fit, biomarker, cohort, model, n_override = NULL) {
      fit_s <- summary(fit)
      n_coef <- nrow(fit_s$coefficients)
      n_val <- if (!is.null(n_override)) n_override else fit_s$n
      data.frame(
            parameter = rownames(fit_s$coefficients),
            hr        = fit_s$conf.int[1:n_coef, 1],
            lower95   = fit_s$conf.int[1:n_coef, 3],
            upper95   = fit_s$conf.int[1:n_coef, 4],
            pvalue    = fit_s$coefficients[1:n_coef, 5],
            n         = n_val,
            nevent    = fit_s$nevent,
            biomarker = biomarker,
            cohort    = cohort,
            model     = model,
            stringsAsFactors = FALSE
      )
}

create_tvc_data <- function(dat, A_col, entry_col, exit_col, event_col,
                            estpos_col = NULL) {
      A         <- dat[[A_col]]
      entry_age <- dat[[entry_col]]
      exit_age  <- dat[[exit_col]]
      event     <- dat[[event_col]]
      estpos    <- if (!is.null(estpos_col)) dat[[estpos_col]] else NULL

      rows <- list()
      for (i in seq_along(A)) {
            if (!is.null(estpos) && !is.na(estpos[i]) && !estpos[i]) {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 0, A = 0,
                        status = "never_positive")
            } else if (is.na(A[i])) {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 0, A = 0,
                        status = "never_positive")
            } else if (A[i] <= entry_age[i]) {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 1, A = A[i],
                        status = "already_positive")
            } else if (A[i] < exit_age[i]) {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = A[i],
                        event = 0, Z_tv = 0, A = A[i],
                        status = "transition_pre")
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = A[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 1, A = A[i],
                        status = "transition_post")
            } else {
                  rows[[length(rows) + 1]] <- data.frame(
                        id = i, tstart = entry_age[i], tstop = exit_age[i],
                        event = event[i], Z_tv = 0, A = 0,
                        status = "never_positive")
            }
      }
      result <- do.call(rbind, rows)
      result <- result[result$tstop > result$tstart, ]
      result
}

coefn <- function(fit, nm) {
      cf <- coef(fit)
      if (nm %in% names(cf)) cf[[nm]] else NA_real_
}
inter_coef <- function(fit) {              # the single Z_tv x A_z interaction term
      cf <- coef(fit)
      ix <- grep(":", names(cf))
      if (length(ix) == 1) cf[[ix]] else NA_real_
}


###############################################################################
## run_all_models_capture: identical fits to extract_full_coefficients.R,
##   returns fitted coxph objects + standardization constants + positive-AABC
##   pool + the subject-level analysis frame (`$dat`, `$cols`) for calibration.
###############################################################################

run_all_models_capture <- function(dat, eaoa_col, entry_col, exit_col, event_col,
                                    biomarker_label, cohort_label,
                                    estpos_col = NULL) {
      cat(sprintf("\n=== %s — %s ===\n", cohort_label, biomarker_label))

      dat_c <- dat[complete.cases(dat[, covariates]), ]
      dat_c <- dat_c[which(dat_c[[entry_col]] < dat_c[[exit_col]]), ]

      all_coefs <- data.frame()
      out <- list(biomarker = biomarker_label, cohort = cohort_label,
                  fit_p1 = NULL, fit_p2 = NULL, fit_p3 = NULL, fit_p3_bh = NULL,
                  p1 = NULL, p3 = NULL, aabc_pool = NULL, p3_degenerate = NA)

      # ---- P1: Countdown / Standard ---- #
      pos_filter <- if (!is.null(estpos_col)) which(dat_c[[estpos_col]] == 1) else
            which(!is.na(dat_c[[eaoa_col]]))
      dat_c$time_to_dx <- dat_c[[exit_col]] - dat_c[[eaoa_col]]
      dat_p1 <- dat_c[intersect(pos_filter, which(dat_c$time_to_dx > 0)), ]

      if (nrow(dat_p1) >= 10 && sum(dat_p1[[event_col]]) >= 5) {
            Z_mean <- mean(dat_p1[[eaoa_col]])
            Z_sd   <- sd(dat_p1[[eaoa_col]])
            dat_p1$Z_std <- (dat_p1[[eaoa_col]] - Z_mean) / Z_sd
            fml_p1 <- as.formula(paste0("Surv(time_to_dx, ", event_col,
                                        ") ~ Sex_F + EDUC_z + apoe4 + Z_std"))
            fit_p1 <- coxph(fml_p1, data = dat_p1)
            all_coefs <- rbind(all_coefs,
                  extract_all_hr(fit_p1, biomarker_label, cohort_label, "P1_countdown"))
            out$fit_p1 <- fit_p1
            out$p1 <- list(Z_mean = Z_mean, Z_sd = Z_sd)
            cat(sprintf("  P1 fit: n=%d, events=%d, muZ=%.2f sdZ=%.2f\n",
                        nrow(dat_p1), sum(dat_p1[[event_col]]), Z_mean, Z_sd))
      } else {
            cat("  P1 SKIPPED (insufficient n/events)\n")
      }

      # ---- P2: TVC z_only (TV-BC) ---- #
      tvc_long <- create_tvc_data(dat_c, eaoa_col, entry_col, exit_col,
                                  event_col, estpos_col)
      covar_df <- data.frame(id = seq_len(nrow(dat_c)),
                             Sex_F = dat_c$Sex_F, EDUC_z = dat_c$EDUC_z,
                             apoe4 = dat_c$apoe4)
      tvc_long <- merge(tvc_long, covar_df, by = "id", all.x = TRUE)
      tvc_long <- tvc_long[complete.cases(tvc_long[, covariates]), ]
      n_tvc_subj <- length(unique(tvc_long$id))

      fit_p2 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + Sex_F + EDUC_z + apoe4,
                      data = tvc_long)
      all_coefs <- rbind(all_coefs,
            extract_all_hr(fit_p2, biomarker_label, cohort_label,
                           "P2_tvc_z_only", n_override = n_tvc_subj))
      out$fit_p2 <- fit_p2

      # ---- P3: TVC interaction (TV-AABC) ---- #
      pos_unique <- unique(tvc_long[tvc_long$A > 0, c("id", "A")])$A
      A_mean <- mean(pos_unique)
      A_sd   <- sd(pos_unique)
      tvc_long$A_z <- ifelse(tvc_long$A > 0, (tvc_long$A - A_mean) / A_sd, 0)
      all_positive <- all(tvc_long$Z_tv == 1)

      if (all_positive) {
            fit_p3 <- coxph(Surv(tstart, tstop, event) ~ A_z + Sex_F + EDUC_z + apoe4,
                            data = tvc_long)
            all_coefs <- rbind(all_coefs,
                  extract_all_hr(fit_p3, biomarker_label, cohort_label,
                                 "P3_tvc_interaction", n_override = n_tvc_subj))
      } else {
            fit_p3 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + A_z:Z_tv +
                                  Sex_F + EDUC_z + apoe4, data = tvc_long)
            all_coefs <- rbind(all_coefs,
                  extract_all_hr(fit_p3, biomarker_label, cohort_label,
                                 "P3_tvc_interaction", n_override = n_tvc_subj))
            # basehaz()/survfit() refuse an interaction without its lower-order
            # term; refit with an explicit product (numerically identical design)
            # solely to extract the baseline cumulative hazard for TV-AABC.
            tvc_long$AzZtv <- tvc_long$A_z * tvc_long$Z_tv
            fit_p3_bh <- coxph(Surv(tstart, tstop, event) ~ Z_tv + AzZtv +
                                     Sex_F + EDUC_z + apoe4, data = tvc_long)
            if (max(abs(inter_coef(fit_p3) - coefn(fit_p3_bh, "AzZtv")),
                    abs(coefn(fit_p3, "Z_tv") - coefn(fit_p3_bh, "Z_tv"))) > 1e-8)
                  stop("P3 reparametrization changed the fit; aborting.")
            out$fit_p3_bh <- fit_p3_bh
      }
      out$fit_p3 <- fit_p3
      out$p3 <- list(A_mean = A_mean, A_sd = A_sd)
      out$p3_degenerate <- all_positive
      out$aabc_pool <- pos_unique
      out$coefs <- all_coefs

      # subject-level frame + column map (for calibration cohorts; does not affect fits)
      out$dat  <- dat_c
      out$cols <- list(eaoa = eaoa_col, entry = entry_col, exit = exit_col,
                       event = event_col, estpos = estpos_col)

      cat(sprintf("  P2/P3 fit: %d subjects, %d positive, muA=%.2f sdA=%.2f, degenerate=%s\n",
                  n_tvc_subj, length(pos_unique), A_mean, A_sd, all_positive))
      out
}


###############################################################################
## Baseline hazard step functions (no extrapolation) + conditional helpers
###############################################################################

make_stepH <- function(fit) {
      bh <- basehaz(fit, centered = FALSE)
      bh <- bh[order(bh$time), ]
      list(time = bh$time, hazard = bh$hazard,
           tmin = min(bh$time), tmax = max(bh$time))
}

# Right-continuous step; 0 below first event time; NA strictly above tmax.
evalH <- function(sh, x) {
      idx <- findInterval(x, sh$time)
      H <- ifelse(idx == 0, 0, sh$hazard[pmax(idx, 1)])
      H[x > sh$tmax + 1e-9] <- NA_real_
      H
}

# Conditional cumulative incidence F(t) = 1 - S_c from `anchor` over width t.
cond_F <- function(sh, anchor, eta, t) {
      inc <- evalH(sh, anchor + t) - evalH(sh, anchor)
      1 - exp(-inc * exp(eta))
}

# Median residual time: smallest t with S_c(t) = 0.5; linear interp in Delta.
solve_median <- function(sh, anchor, eta) {
      thr <- log(2) / exp(eta)
      H0 <- evalH(sh, anchor)
      if (is.na(H0)) return(list(t_med = NA_real_, reached = FALSE, lower_bound = NA_real_))
      ev <- sh$time[sh$time > anchor + 1e-9 & sh$time <= sh$tmax + 1e-9]
      if (length(ev) == 0)
            return(list(t_med = NA_real_, reached = FALSE,
                        lower_bound = max(0, sh$tmax - anchor)))
      Delta <- evalH(sh, ev) - H0
      tt <- ev - anchor
      if (Delta[length(Delta)] < thr)
            return(list(t_med = NA_real_, reached = FALSE, lower_bound = sh$tmax - anchor))
      j <- which(Delta >= thr)[1]
      if (j == 1L) { t_lo <- 0; D_lo <- 0 } else { t_lo <- tt[j - 1]; D_lo <- Delta[j - 1] }
      t_hi <- tt[j]; D_hi <- Delta[j]
      t_med <- if (D_hi > D_lo) t_lo + (thr - D_lo) / (D_hi - D_lo) * (t_hi - t_lo) else t_hi
      list(t_med = t_med, reached = TRUE, lower_bound = NA_real_)
}

# Per-method linear predictors at Sex_F = 0, EDUC_z = 0 (mirror Section 7).
eta_std <- function(f, aabc, ap)
      coefn(f$fit_p1, "apoe4") * ap +
      coefn(f$fit_p1, "Z_std") * (aabc - f$p1$Z_mean) / f$p1$Z_sd
eta_tvbc <- function(f, ap)
      coefn(f$fit_p2, "Z_tv") + coefn(f$fit_p2, "apoe4") * ap
eta_tvaabc <- function(f, aabc, ap) {
      if (f$p3_degenerate)
            coefn(f$fit_p3, "A_z") * (aabc - f$p3$A_mean) / f$p3$A_sd +
                  coefn(f$fit_p3, "apoe4") * ap
      else
            coefn(f$fit_p3_bh, "Z_tv") +
                  coefn(f$fit_p3_bh, "AzZtv") * (aabc - f$p3$A_mean) / f$p3$A_sd +
                  coefn(f$fit_p3, "apoe4") * ap
}

sh_of <- function(f, method) {
      if (method == "Standard") make_stepH(f$fit_p1)
      else if (method == "TV-BC") make_stepH(f$fit_p2)
      else make_stepH(if (!is.null(f$fit_p3_bh)) f$fit_p3_bh else f$fit_p3)
}


###############################################################################
## Data loading (Sections 1-3; verbatim structure from extract_full_coefficients.R)
###############################################################################

cat("\n[core] Loading BIOCARD data...\n")
load(file.path(data_dir, "analysis_data_merged.rda"))
stopifnot(all(analysis_data$onset.age <= analysis_data$censor.age))

load(file.path(data_dir, "BIOCARD_CSF_SILA_intermediate.rda"))
ab_eaoa <- resfit_last_ab[, c("SUBJECT_ID", "EAOA_AB", "estpos")]
names(ab_eaoa)[3] <- "estpos_ab"
ptau_eaoa <- resfit_last_ptau[, c("SUBJECT_ID", "EAOA_PTAU", "estpos")]
names(ptau_eaoa)[3] <- "estpos_ptau"

dat_biocard_csf <- merge(analysis_data,
      ab_eaoa[, c("SUBJECT_ID", "EAOA_AB", "estpos_ab")], by = "SUBJECT_ID", all.x = TRUE)
dat_biocard_csf <- merge(dat_biocard_csf,
      ptau_eaoa[, c("SUBJECT_ID", "EAOA_PTAU", "estpos_ptau")], by = "SUBJECT_ID", all.x = TRUE)

load(file.path(data_dir, "BIOCARD_plasma_SILA_intermediate.rda"))
plasma_analysis <- plasma_analysis[which(plasma_analysis$onset.age > plasma_analysis$jhu_baseline_age), ]

cat("[core] Loading ADNI PET data...\n")
load(file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"))
eaoa_df <- resfit_last[, c("RID", "estaget0", "estpos", "first_pos_age")]
names(eaoa_df)[2] <- "EAOA"
dat_adni_pet <- merge(dat_surv_cn, eaoa_df, by = "RID", all.x = FALSE)

data_2026 <- file.path(project_root, "ADNI_2026_data")
if (file.exists(file.path(data_2026, "PTDEMOG_11Feb2026.csv"))) {
      adni_exclude_demo <- read.csv(file.path(data_2026, "PTDEMOG_11Feb2026.csv"))
      adni_exclude_rids <- unique(adni_exclude_demo$RID[grepl("^381_S_", adni_exclude_demo$PTID)])
      dat_adni_pet <- dat_adni_pet[!dat_adni_pet$RID %in% adni_exclude_rids, ]
      rm(adni_exclude_demo)
}
dat_adni_pet <- merge(dat_adni_pet, dat_demo[, c("RID", "PTGENDER", "PTEDUCAT")],
                      by = "RID", all.x = TRUE)
dat_adni_pet$Sex_F <- as.integer(dat_adni_pet$PTGENDER == 2)
dat_adni_pet$EDUC_z <- scale(dat_adni_pet$PTEDUCAT)[, 1]
dat_adni_pet <- merge(dat_adni_pet, dat_apoe[, c("RID", "apoe4")], by = "RID", all.x = TRUE)

cat("[core] Loading ADNI plasma data...\n")
sila_plasma_file <- file.path(data_dir, "ADNI_plasma_SILA_intermediate.rda")
has_plasma <- FALSE
if (file.exists(sila_plasma_file)) {
      sila_env <- new.env()
      load(sila_plasma_file, envir = sila_env)
      rflp <- NULL
      if (exists("resfit_last_fuj", envir = sila_env)) {
            rflp <- get("resfit_last_fuj", envir = sila_env)
      } else if (exists("resfit_last_plasma", envir = sila_env)) {
            rflp <- get("resfit_last_plasma", envir = sila_env)
      }
      if (!is.null(rflp)) {
            eaoa_col_p <- if ("EAOA_plasma" %in% names(rflp)) "EAOA_plasma" else
                  if ("estaget0" %in% names(rflp)) "estaget0" else NULL
            estpos_col_p <- if ("estpos" %in% names(rflp)) "estpos" else NULL
            if (!is.null(eaoa_col_p)) {
                  sila_fuj_z_df <- data.frame(
                        RID = rflp$RID,
                        age_at_positivity = rflp[[eaoa_col_p]],
                        ever_positive = if (!is.null(estpos_col_p)) rflp[[estpos_col_p]]
                                        else !is.na(rflp[[eaoa_col_p]]))
                  dat_base_plasma <- dat_surv_cn
                  dat_base_plasma <- merge(dat_base_plasma,
                        dat_demo[, c("RID", "PTGENDER", "PTEDUCAT")], by = "RID", all.x = TRUE)
                  dat_base_plasma$Sex_F <- as.integer(dat_base_plasma$PTGENDER == 2)
                  dat_base_plasma$EDUC_z <- scale(dat_base_plasma$PTEDUCAT)[, 1]
                  dat_base_plasma <- merge(dat_base_plasma,
                        dat_apoe[, c("RID", "apoe4")], by = "RID", all.x = TRUE)
                  dat_base_plasma <- dat_base_plasma[complete.cases(dat_base_plasma[, covariates]), ]
                  if (exists("adni_exclude_rids")) {
                        dat_base_plasma <- dat_base_plasma[!dat_base_plasma$RID %in% adni_exclude_rids, ]
                  }
                  dat_adni_plasma <- merge(dat_base_plasma,
                        sila_fuj_z_df[, c("RID", "age_at_positivity", "ever_positive")],
                        by = "RID", all.x = FALSE)
                  has_plasma <- TRUE
            }
      }
      rm(sila_env)
}
stopifnot(has_plasma)   # rigorous: do not silently proceed without the 5th biomarker


###############################################################################
## SECTION 4: Run all 5 biomarkers, capturing fits
###############################################################################

fits <- list()
fits[["CSF_AB42_AB40"]] <- run_all_models_capture(
      dat_biocard_csf, "EAOA_AB", "baseline.age", "onset.age", "d",
      "CSF_AB42_AB40", "BIOCARD", "estpos_ab")
fits[["CSF_pTau181"]] <- run_all_models_capture(
      dat_biocard_csf, "EAOA_PTAU", "baseline.age", "onset.age", "d",
      "CSF_pTau181", "BIOCARD", "estpos_ptau")
fits[["Plasma_pTau181"]] <- run_all_models_capture(
      plasma_analysis, "EAOA_plasma", "jhu_baseline_age", "onset.age", "d",
      "Plasma_pTau181", "BIOCARD", "estpos")
fits[["Amyloid_PET_FBP"]] <- run_all_models_capture(
      dat_adni_pet, "EAOA", "baseline.age", "onset.age", "d",
      "Amyloid_PET_FBP", "ADNI", "estpos")
fits[["Plasma_pTau217_Fuji"]] <- run_all_models_capture(
      dat_adni_plasma, "age_at_positivity", "baseline.age", "censor.age", "d",
      "Plasma_pTau217_Fuji", "ADNI", "ever_positive")


###############################################################################
## SECTION 5: VALIDATION GATE G0 — refit coefficients must match the saved CSV
###############################################################################

cat("\n[core] --- Validation gate G0: refit vs saved coefficients ---\n")
saved <- read.csv(file.path(results_dir, "full_coefficients_all_models.csv"),
                  stringsAsFactors = FALSE)
refit <- do.call(rbind, lapply(fits, function(x) x$coefs))

key <- function(d) paste(d$biomarker, d$model, d$parameter, sep = "||")
saved$key <- key(saved); refit$key <- key(refit)
mrg <- merge(refit, saved, by = "key", suffixes = c("_new", "_saved"))

stopifnot(nrow(mrg) == nrow(refit))   # every refit row must be matched in saved
max_hr_diff <- max(abs(mrg$hr_new - mrg$hr_saved))
max_p_diff  <- max(abs(mrg$pvalue_new - mrg$pvalue_saved))
cat(sprintf("[core]   matched %d coefficients; max|dHR|=%.3e, max|dP|=%.3e\n",
            nrow(mrg), max_hr_diff, max_p_diff))
if (max_hr_diff > COEF_TOL || max_p_diff > COEF_TOL) {
      bad <- mrg[abs(mrg$hr_new - mrg$hr_saved) > COEF_TOL, ]
      print(bad[, c("biomarker_new", "model_new", "parameter_new", "hr_new", "hr_saved")])
      stop("VALIDATION FAILED: refit coefficients differ from saved CSV. Halting.")
}
cat("[core]   G0 PASS: refit reproduces the audited models exactly.\n")
