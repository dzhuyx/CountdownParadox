###############################################################################
## Extract Full Covariate Coefficients for Supplementary Table 8
##
## Purpose: Re-run all 15 Cox PH model fits (5 biomarkers x 3 models) and
##   extract ALL coefficients (not just target effects) for the full-model
##   supplementary table.
##
## Input:  Same .rda files as BIOCARD_ADNI_main_analysis.R and
##         ADNI_SILA_reanalysis_v2.R / ADNI_plasma_reanalysis.R
##
## Output: results/full_coefficients_all_models.csv
##
## Author: Yuxin Zhu
## Date: April 2026
###############################################################################

rm(list = ls())
library(survival)

project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set CP_PROJECT_ROOT to the analysis root (the folder that holds data/ and results/).")
data_dir     <- file.path(project_root, "data")
results_dir  <- file.path(project_root, "results")

covariates <- c("Sex_F", "EDUC_z", "apoe4")


###############################################################################
## Helper functions (from BIOCARD_ADNI_main_analysis.R)
###############################################################################

extract_all_hr <- function(fit, biomarker, cohort, model, n_override = NULL) {
      # n_override: for counting-process Cox fits, fit_s$n returns interval rows,
      # not unique subjects. Pass length(unique(data$id)) to report subjects.
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
                        status = "never_positive"
                  )
            } else if (is.na(A[i])) {
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


run_all_models <- function(dat, eaoa_col, entry_col, exit_col, event_col,
                           biomarker_label, cohort_label,
                           estpos_col = NULL) {
      cat(sprintf("\n=== %s — %s ===\n", cohort_label, biomarker_label))

      dat_c <- dat[complete.cases(dat[, covariates]), ]
      dat_c <- dat_c[which(dat_c[[entry_col]] < dat_c[[exit_col]]), ]

      all_coefs <- data.frame()

      # ---- P1: Countdown ---- #
      pos_filter <- if (!is.null(estpos_col)) {
            which(dat_c[[estpos_col]] == 1)
      } else {
            which(!is.na(dat_c[[eaoa_col]]))
      }
      dat_c$time_to_dx <- dat_c[[exit_col]] - dat_c[[eaoa_col]]
      dat_p1 <- dat_c[intersect(pos_filter, which(dat_c$time_to_dx > 0)), ]

      if (nrow(dat_p1) >= 10 && sum(dat_p1[[event_col]]) >= 5) {
            Z_mean <- mean(dat_p1[[eaoa_col]])
            Z_sd   <- sd(dat_p1[[eaoa_col]])
            dat_p1$Z_std <- (dat_p1[[eaoa_col]] - Z_mean) / Z_sd

            fml_p1 <- as.formula(paste0("Surv(time_to_dx, ", event_col,
                                         ") ~ Sex_F + EDUC_z + apoe4 + Z_std"))
            fit_p1 <- coxph(fml_p1, data = dat_p1)
            coefs_p1 <- extract_all_hr(fit_p1, biomarker_label, cohort_label,
                                        "P1_countdown")
            all_coefs <- rbind(all_coefs, coefs_p1)
            cat(sprintf("  P1: %d coefficients extracted\n", nrow(coefs_p1)))
      }

      # ---- P2: TVC z_only ---- #
      tvc_long <- create_tvc_data(dat_c, eaoa_col, entry_col, exit_col,
                                   event_col, estpos_col)
      covar_df <- data.frame(
            id = seq_len(nrow(dat_c)),
            Sex_F = dat_c$Sex_F,
            EDUC_z = dat_c$EDUC_z,
            apoe4 = dat_c$apoe4
      )
      tvc_long <- merge(tvc_long, covar_df, by = "id", all.x = TRUE)
      tvc_long <- tvc_long[complete.cases(tvc_long[, covariates]), ]

      n_tvc_subj <- length(unique(tvc_long$id))

      fit_p2 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + Sex_F + EDUC_z + apoe4,
                       data = tvc_long)
      coefs_p2 <- extract_all_hr(fit_p2, biomarker_label, cohort_label,
                                  "P2_tvc_z_only", n_override = n_tvc_subj)
      all_coefs <- rbind(all_coefs, coefs_p2)
      cat(sprintf("  P2: %d coefficients extracted\n", nrow(coefs_p2)))

      # ---- P3: TVC interaction ---- #
      pos_unique <- unique(tvc_long[tvc_long$A > 0, c("id", "A")])$A
      A_mean <- mean(pos_unique)
      A_sd   <- sd(pos_unique)
      tvc_long$A_z <- ifelse(tvc_long$A > 0, (tvc_long$A - A_mean) / A_sd, 0)

      all_positive <- all(tvc_long$Z_tv == 1)

      if (all_positive) {
            fit_p3 <- coxph(Surv(tstart, tstop, event) ~ A_z + Sex_F + EDUC_z + apoe4,
                             data = tvc_long)
            coefs_p3 <- extract_all_hr(fit_p3, biomarker_label, cohort_label,
                                        "P3_tvc_interaction", n_override = n_tvc_subj)
            all_coefs <- rbind(all_coefs, coefs_p3)
            cat(sprintf("  P3 (degenerate): %d coefficients extracted\n", nrow(coefs_p3)))
      } else {
            fit_p3 <- coxph(Surv(tstart, tstop, event) ~ Z_tv + A_z:Z_tv +
                                   Sex_F + EDUC_z + apoe4,
                             data = tvc_long)
            coefs_p3 <- extract_all_hr(fit_p3, biomarker_label, cohort_label,
                                        "P3_tvc_interaction", n_override = n_tvc_subj)
            all_coefs <- rbind(all_coefs, coefs_p3)
            cat(sprintf("  P3: %d coefficients extracted\n", nrow(coefs_p3)))
      }

      # Metadata needed for HR-vs-AABC figure (vcov, A standardization, range)
      fit_meta <- list(
            fit_p2 = fit_p2,
            fit_p3 = fit_p3,
            A_mean = A_mean,
            A_sd   = A_sd,
            aabc_range = range(pos_unique),
            all_positive = all_positive,
            biomarker = biomarker_label,
            cohort    = cohort_label
      )

      return(list(coefs = all_coefs, fit_meta = fit_meta))
}


###############################################################################
## SECTION 1: BIOCARD data loading (from BIOCARD_ADNI_main_analysis.R)
###############################################################################

cat("Loading BIOCARD data...\n")

load(file.path(data_dir, "analysis_data_merged.rda"))
stopifnot(all(analysis_data$onset.age <= analysis_data$censor.age))

load(file.path(data_dir, "BIOCARD_CSF_SILA_intermediate.rda"))
ab_eaoa <- resfit_last_ab[, c("SUBJECT_ID", "EAOA_AB", "estpos")]
names(ab_eaoa)[3] <- "estpos_ab"
ptau_eaoa <- resfit_last_ptau[, c("SUBJECT_ID", "EAOA_PTAU", "estpos")]
names(ptau_eaoa)[3] <- "estpos_ptau"

dat_biocard_csf <- merge(analysis_data,
                         ab_eaoa[, c("SUBJECT_ID", "EAOA_AB", "estpos_ab")],
                         by = "SUBJECT_ID", all.x = TRUE)
dat_biocard_csf <- merge(dat_biocard_csf,
                         ptau_eaoa[, c("SUBJECT_ID", "EAOA_PTAU", "estpos_ptau")],
                         by = "SUBJECT_ID", all.x = TRUE)

load(file.path(data_dir, "BIOCARD_plasma_SILA_intermediate.rda"))

# Plasma sub-study: left-truncate at the JHU-phase baseline; require CU there. Drop
# subjects with no follow-up after the JHU baseline (consistent with the main analysis).
plasma_analysis <- plasma_analysis[which(plasma_analysis$onset.age > plasma_analysis$jhu_baseline_age), ]


###############################################################################
## SECTION 2: ADNI PET data loading (from ADNI_SILA_reanalysis_v2.R)
###############################################################################

cat("Loading ADNI PET data...\n")

load(file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"))

eaoa_df <- resfit_last[, c("RID", "estaget0", "estpos", "first_pos_age")]
names(eaoa_df)[2] <- "EAOA"
dat_adni_pet <- merge(dat_surv_cn, eaoa_df, by = "RID", all.x = FALSE)

# Exclude 381_S subjects (ADNI March 2026 data quality notice)
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
dat_adni_pet <- merge(dat_adni_pet, dat_apoe[, c("RID", "apoe4")],
                      by = "RID", all.x = TRUE)

# Standardize EAOA among estpos=1 (same as original script)
eaoa_pos_pool <- which(!is.na(dat_adni_pet$EAOA) & dat_adni_pet$estpos == 1)
eaoa_mean_pos <- mean(dat_adni_pet$EAOA[eaoa_pos_pool])
eaoa_sd_pos   <- sd(dat_adni_pet$EAOA[eaoa_pos_pool])


###############################################################################
## SECTION 3: ADNI plasma data loading (from ADNI_plasma_reanalysis.R)
##
## The plasma script uses different conventions from the PET script:
##   - exit_col = "censor.age" (not "onset.age")
##   - eaoa_col = "age_at_positivity"
##   - estpos_col = "ever_positive"
##   - Object is resfit_last_fuj (not resfit_last_plasma)
###############################################################################

cat("Loading ADNI plasma data...\n")

sila_plasma_file <- file.path(data_dir, "ADNI_plasma_SILA_intermediate.rda")
has_plasma <- FALSE

if (file.exists(sila_plasma_file)) {
      sila_env <- new.env()
      load(sila_plasma_file, envir = sila_env)

      # Try resfit_last_fuj first (current format), then resfit_last_plasma (legacy)
      rflp <- NULL
      if (exists("resfit_last_fuj", envir = sila_env)) {
            rflp <- get("resfit_last_fuj", envir = sila_env)
      } else if (exists("resfit_last_plasma", envir = sila_env)) {
            rflp <- get("resfit_last_plasma", envir = sila_env)
      }

      if (!is.null(rflp)) {
            # Find EAOA column
            eaoa_col <- if ("EAOA_plasma" %in% names(rflp)) "EAOA_plasma"
                        else if ("estaget0" %in% names(rflp)) "estaget0"
                        else NULL
            estpos_col <- if ("estpos" %in% names(rflp)) "estpos" else NULL

            if (!is.null(eaoa_col)) {
                  sila_fuj_z_df <- data.frame(
                        RID = rflp$RID,
                        age_at_positivity = rflp[[eaoa_col]],
                        ever_positive = if (!is.null(estpos_col)) rflp[[estpos_col]]
                                        else !is.na(rflp[[eaoa_col]])
                  )
                  cat(sprintf("SILA Fujirebio Z: %d subjects (%d positive)\n",
                              nrow(sila_fuj_z_df),
                              sum(sila_fuj_z_df$ever_positive, na.rm = TRUE)))

                  # Build dat_base the same way as ADNI_plasma_reanalysis.R:
                  # uses censor.age (not onset.age) as exit time
                  dat_base_plasma <- dat_surv_cn
                  dat_base_plasma <- merge(dat_base_plasma,
                        dat_demo[, c("RID", "PTGENDER", "PTEDUCAT")],
                        by = "RID", all.x = TRUE)
                  dat_base_plasma$Sex_F <- as.integer(dat_base_plasma$PTGENDER == 2)
                  dat_base_plasma$EDUC_z <- scale(dat_base_plasma$PTEDUCAT)[, 1]
                  dat_base_plasma <- merge(dat_base_plasma,
                        dat_apoe[, c("RID", "apoe4")],
                        by = "RID", all.x = TRUE)
                  dat_base_plasma <- dat_base_plasma[
                        complete.cases(dat_base_plasma[, covariates]), ]

                  # Exclude 381_S subjects
                  if (exists("adni_exclude_rids")) {
                        dat_base_plasma <- dat_base_plasma[
                              !dat_base_plasma$RID %in% adni_exclude_rids, ]
                  }

                  # Merge Z
                  dat_adni_plasma <- merge(dat_base_plasma,
                        sila_fuj_z_df[, c("RID", "age_at_positivity", "ever_positive")],
                        by = "RID", all.x = FALSE)
                  has_plasma <- TRUE
                  cat(sprintf("ADNI plasma analysis dataset: %d subjects\n",
                              nrow(dat_adni_plasma)))
            }
      }
      rm(sila_env)
} else {
      cat("WARNING: ADNI plasma SILA file not found\n")
}


###############################################################################
## SECTION 4: Run all models
###############################################################################

all_coefs <- data.frame()
all_fits  <- list()

# C1: BIOCARD CSF AB42/AB40
out <- run_all_models(
      dat = dat_biocard_csf, eaoa_col = "EAOA_AB",
      entry_col = "baseline.age", exit_col = "onset.age", event_col = "d",
      biomarker_label = "CSF_AB42_AB40", cohort_label = "BIOCARD",
      estpos_col = "estpos_ab"
)
all_coefs <- rbind(all_coefs, out$coefs)
all_fits[["BIOCARD_CSF_AB42_AB40"]] <- out$fit_meta

# C2: BIOCARD CSF p-tau181
out <- run_all_models(
      dat = dat_biocard_csf, eaoa_col = "EAOA_PTAU",
      entry_col = "baseline.age", exit_col = "onset.age", event_col = "d",
      biomarker_label = "CSF_pTau181", cohort_label = "BIOCARD",
      estpos_col = "estpos_ptau"
)
all_coefs <- rbind(all_coefs, out$coefs)
all_fits[["BIOCARD_CSF_pTau181"]] <- out$fit_meta

# C3: BIOCARD plasma p-tau181
out <- run_all_models(
      dat = plasma_analysis, eaoa_col = "EAOA_plasma",
      entry_col = "jhu_baseline_age", exit_col = "onset.age", event_col = "d",
      biomarker_label = "Plasma_pTau181", cohort_label = "BIOCARD",
      estpos_col = "estpos"
)
all_coefs <- rbind(all_coefs, out$coefs)
all_fits[["BIOCARD_Plasma_pTau181"]] <- out$fit_meta

# C4: ADNI amyloid PET
out <- run_all_models(
      dat = dat_adni_pet, eaoa_col = "EAOA",
      entry_col = "baseline.age", exit_col = "onset.age", event_col = "d",
      biomarker_label = "Amyloid_PET_FBP", cohort_label = "ADNI",
      estpos_col = "estpos"
)
all_coefs <- rbind(all_coefs, out$coefs)
all_fits[["ADNI_Amyloid_PET_FBP"]] <- out$fit_meta

# C5: ADNI plasma p-tau217 Fujirebio
# Uses censor.age (not onset.age) and age_at_positivity/ever_positive columns
# (matching ADNI_plasma_reanalysis.R conventions)
if (has_plasma) {
      out <- run_all_models(
            dat = dat_adni_plasma, eaoa_col = "age_at_positivity",
            entry_col = "baseline.age", exit_col = "censor.age", event_col = "d",
            biomarker_label = "Plasma_pTau217_Fuji", cohort_label = "ADNI",
            estpos_col = "ever_positive"
      )
      all_coefs <- rbind(all_coefs, out$coefs)
      all_fits[["ADNI_Plasma_pTau217_Fuji"]] <- out$fit_meta
}


###############################################################################
## SECTION 5: Format and save
###############################################################################

# Add formatted columns
all_coefs$HR_95CI <- sprintf("%.2f (%.2f\u2013%.2f)",
                              all_coefs$hr, all_coefs$lower95, all_coefs$upper95)
all_coefs$p_formatted <- ifelse(
      all_coefs$pvalue < 0.001, "< 0.001",
      ifelse(all_coefs$pvalue < 0.01, "< 0.01",
             ifelse(all_coefs$pvalue < 0.05, "< 0.05",
                    sprintf("= %.2f", all_coefs$pvalue))))

out_file <- file.path(results_dir, "full_coefficients_all_models.csv")
write.csv(all_coefs, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))

fits_file <- file.path(results_dir, "fits_p2_p3_all.rds")
saveRDS(all_fits, fits_file)
cat(sprintf("Saved: %s (%d biomarker-cohorts)\n", fits_file, length(all_fits)))
cat(sprintf("Total coefficients: %d rows across %d unique models\n",
            nrow(all_coefs),
            length(unique(paste(all_coefs$biomarker, all_coefs$model)))))

# Print summary
cat("\n=== Summary ===\n")
for (bm in unique(all_coefs$biomarker)) {
      for (mod in unique(all_coefs$model[all_coefs$biomarker == bm])) {
            sub <- all_coefs[all_coefs$biomarker == bm & all_coefs$model == mod, ]
            cat(sprintf("\n%s — %s — %s (n=%d, events=%d):\n",
                        sub$cohort[1], bm, mod, sub$n[1], sub$nevent[1]))
            for (j in seq_len(nrow(sub))) {
                  cat(sprintf("  %-20s  %s; P %s\n",
                              sub$parameter[j], sub$HR_95CI[j], sub$p_formatted[j]))
            }
      }
}

cat("\nDone.\n")
