###############################################################################
## compute_manuscript_descriptives.R
##
## Stage 3b: Compute descriptive statistics for manuscript tables and figures.
## All outputs are CSV files consumed by create_natmed_tables.R and
## create_natmed_figures.R — no hardcoded values in downstream scripts.
##
## Input:
##   data/analysis_data_merged.rda             (BIOCARD survival + covariates)
##   data/BIOCARD_CSF_SILA_intermediate.rda    (CSF SILA EAOA)
##   data/BIOCARD_plasma_SILA_intermediate.rda (Plasma SILA EAOA)
##   data/ADNI_SILA_intermediate_2026.rda      (ADNI PET survival + SILA)
##   data/ADNI_plasma_SILA_intermediate.rda    (ADNI Plasma SILA)
##
## Output:
##   results/table1_demographics.csv
##   results/person_years_by_ztv.csv
##   results/eaoa_summary.csv
##
## Author: Yuxin Zhu
## Date: April 2026
###############################################################################

rm(list = ls())
library(survival)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT", "/Users/daisyzhu/Documents/Research Projects/CountdownParadox_BiomarkerPositivity/CountdownParadox_Analysis")
out_dir      <- project_root
data_dir     <- file.path(project_root, "data")
results_dir  <- file.path(project_root, "results")
data_2026    <- file.path(project_root, "ADNI_2026_data")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

covariates <- c("Sex_F", "EDUC_z", "apoe4")


###############################################################################
## Helper: create_tvc_data (same as in BIOCARD_ADNI_main_analysis.R)
###############################################################################

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
## Helper: compute person-years and events by Z_tv state
###############################################################################

compute_person_years <- function(tvc_long) {
      py_z0 <- sum(tvc_long$tstop[tvc_long$Z_tv == 0] -
                    tvc_long$tstart[tvc_long$Z_tv == 0])
      py_z1 <- sum(tvc_long$tstop[tvc_long$Z_tv == 1] -
                    tvc_long$tstart[tvc_long$Z_tv == 1])
      ev_z0 <- sum(tvc_long$event[tvc_long$Z_tv == 0])
      ev_z1 <- sum(tvc_long$event[tvc_long$Z_tv == 1])
      data.frame(py_z0 = round(py_z0), py_z1 = round(py_z1),
                 ev_z0 = ev_z0, ev_z1 = ev_z1)
}


###############################################################################
## Helper: compute demographics for a dataset
###############################################################################

compute_demographics <- function(dat, entry_col, exit_col, event_col,
                                 censor_col = NULL, educ_raw_col = NULL) {
      n <- nrow(dat)
      n_events <- sum(dat[[event_col]])

      # Follow-up = censor.age - baseline.age (total observation, B2 convention)
      censor_age <- if (!is.null(censor_col)) dat[[censor_col]] else dat[[exit_col]]
      followup <- censor_age - dat[[entry_col]]

      # Entry age
      entry_age <- dat[[entry_col]]

      # Education (raw, not z-scored)
      educ <- if (!is.null(educ_raw_col)) dat[[educ_raw_col]] else NULL

      data.frame(
            n = n,
            n_events = n_events,
            pct_events = round(n_events / n * 100, 1),
            pct_female = round(mean(dat$Sex_F, na.rm = TRUE) * 100, 1),
            educ_mean = if (!is.null(educ)) round(mean(educ, na.rm = TRUE), 1) else NA,
            educ_sd = if (!is.null(educ)) round(sd(educ, na.rm = TRUE), 1) else NA,
            n_apoe4 = sum(dat$apoe4 == 1, na.rm = TRUE),
            pct_apoe4 = round(mean(dat$apoe4 == 1, na.rm = TRUE) * 100, 1),
            entry_age_mean = round(mean(entry_age), 1),
            entry_age_sd = round(sd(entry_age), 1),
            entry_age_min = round(min(entry_age)),
            entry_age_max = round(max(entry_age)),
            followup_mean = round(mean(followup), 1),
            followup_sd = round(sd(followup), 1),
            followup_median = round(median(followup), 1),
            stringsAsFactors = FALSE
      )
}


###############################################################################
## SECTION 1: Load BIOCARD data
###############################################################################

cat("================================================================\n")
cat("SECTION 1: LOAD BIOCARD DATA\n")
cat("================================================================\n\n")

# Survival + covariates
load(file.path(data_dir, "analysis_data_merged.rda"))
cat(sprintf("analysis_data: %d subjects, %d events\n",
            nrow(analysis_data), sum(analysis_data$d)))

# CSF SILA
load(file.path(data_dir, "BIOCARD_CSF_SILA_intermediate.rda"))

# SILA diagnostic: capture estpos breakdown BEFORE filtering
sila_diag <- list()
sila_diag[["BC_AB"]] <- data.frame(
      biomarker = "CSF_AB42_AB40", cohort = "BIOCARD",
      n_sila_total = nrow(resfit_last_ab),
      n_estpos_true = sum(resfit_last_ab$estpos == TRUE, na.rm = TRUE),
      n_estpos_false = sum(resfit_last_ab$estpos == FALSE, na.rm = TRUE),
      mean_eaoa_all = round(mean(resfit_last_ab$EAOA_AB, na.rm = TRUE), 1),
      sd_eaoa_all = round(sd(resfit_last_ab$EAOA_AB, na.rm = TRUE), 1),
      mean_eaoa_pos = round(mean(resfit_last_ab$EAOA_AB[which(resfit_last_ab$estpos == TRUE)], na.rm = TRUE), 1),
      sd_eaoa_pos = round(sd(resfit_last_ab$EAOA_AB[which(resfit_last_ab$estpos == TRUE)], na.rm = TRUE), 1),
      mean_eaoa_neg = round(mean(resfit_last_ab$EAOA_AB[which(resfit_last_ab$estpos == FALSE)], na.rm = TRUE), 1),
      sd_eaoa_neg = round(sd(resfit_last_ab$EAOA_AB[which(resfit_last_ab$estpos == FALSE)], na.rm = TRUE), 1)
)
sila_diag[["BC_PTAU"]] <- data.frame(
      biomarker = "CSF_pTau181", cohort = "BIOCARD",
      n_sila_total = nrow(resfit_last_ptau),
      n_estpos_true = sum(resfit_last_ptau$estpos == TRUE, na.rm = TRUE),
      n_estpos_false = sum(resfit_last_ptau$estpos == FALSE, na.rm = TRUE),
      mean_eaoa_all = round(mean(resfit_last_ptau$EAOA_PTAU, na.rm = TRUE), 1),
      sd_eaoa_all = round(sd(resfit_last_ptau$EAOA_PTAU, na.rm = TRUE), 1),
      mean_eaoa_pos = round(mean(resfit_last_ptau$EAOA_PTAU[which(resfit_last_ptau$estpos == TRUE)], na.rm = TRUE), 1),
      sd_eaoa_pos = round(sd(resfit_last_ptau$EAOA_PTAU[which(resfit_last_ptau$estpos == TRUE)], na.rm = TRUE), 1),
      mean_eaoa_neg = round(mean(resfit_last_ptau$EAOA_PTAU[which(resfit_last_ptau$estpos == FALSE)], na.rm = TRUE), 1),
      sd_eaoa_neg = round(sd(resfit_last_ptau$EAOA_PTAU[which(resfit_last_ptau$estpos == FALSE)], na.rm = TRUE), 1)
)

# Keep EAOA values for all subjects; estpos used by create_tvc_data
ab_eaoa <- resfit_last_ab[, c("SUBJECT_ID", "EAOA_AB", "estpos")]
names(ab_eaoa)[3] <- "estpos_ab"

ptau_eaoa <- resfit_last_ptau[, c("SUBJECT_ID", "EAOA_PTAU", "estpos")]
names(ptau_eaoa)[3] <- "estpos_ptau"

dat_bc_csf <- merge(analysis_data, ab_eaoa[, c("SUBJECT_ID", "EAOA_AB", "estpos_ab")],
                    by = "SUBJECT_ID", all.x = TRUE)
dat_bc_csf <- merge(dat_bc_csf, ptau_eaoa[, c("SUBJECT_ID", "EAOA_PTAU", "estpos_ptau")],
                    by = "SUBJECT_ID", all.x = TRUE)

# Keep only subjects with complete covariates (matching main analysis)
dat_bc_csf <- dat_bc_csf[complete.cases(dat_bc_csf[, covariates]), ]
cat(sprintf("BIOCARD CSF: %d subjects (complete covariates), AB positive=%d, PTAU positive=%d\n",
            nrow(dat_bc_csf),
            sum(dat_bc_csf$estpos_ab == 1, na.rm = TRUE),
            sum(dat_bc_csf$estpos_ptau == 1, na.rm = TRUE)))

# Plasma SILA
load(file.path(data_dir, "BIOCARD_plasma_SILA_intermediate.rda"))

sila_diag[["BC_Plasma"]] <- data.frame(
      biomarker = "Plasma_pTau181", cohort = "BIOCARD",
      n_sila_total = nrow(plasma_analysis),
      n_estpos_true = sum(plasma_analysis$estpos == TRUE, na.rm = TRUE),
      n_estpos_false = sum(plasma_analysis$estpos == FALSE, na.rm = TRUE),
      mean_eaoa_all = round(mean(plasma_analysis$EAOA_plasma, na.rm = TRUE), 1),
      sd_eaoa_all = round(sd(plasma_analysis$EAOA_plasma, na.rm = TRUE), 1),
      mean_eaoa_pos = round(mean(plasma_analysis$EAOA_plasma[which(plasma_analysis$estpos == TRUE)], na.rm = TRUE), 1),
      sd_eaoa_pos = round(sd(plasma_analysis$EAOA_plasma[which(plasma_analysis$estpos == TRUE)], na.rm = TRUE), 1),
      mean_eaoa_neg = round(mean(plasma_analysis$EAOA_plasma[which(plasma_analysis$estpos == FALSE)], na.rm = TRUE), 1),
      sd_eaoa_neg = round(sd(plasma_analysis$EAOA_plasma[which(plasma_analysis$estpos == FALSE)], na.rm = TRUE), 1)
)

# Keep EAOA values; estpos used by create_tvc_data
dat_bc_plasma <- plasma_analysis[complete.cases(plasma_analysis[, covariates]), ]
cat(sprintf("BIOCARD Plasma: %d subjects (complete covariates), %d positive (estpos=1)\n",
            nrow(dat_bc_plasma), sum(dat_bc_plasma$estpos == 1, na.rm = TRUE)))


###############################################################################
## SECTION 2: Load ADNI data
###############################################################################

cat("\n================================================================\n")
cat("SECTION 2: LOAD ADNI DATA\n")
cat("================================================================\n\n")

# ADNI PET SILA intermediate
load(file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"))

# Merge EAOA
eaoa_df <- resfit_last[, c("RID", "estaget0", "estpos", "first_pos_age")]
names(eaoa_df)[2] <- "EAOA"

sila_diag[["ADNI_PET"]] <- data.frame(
      biomarker = "Amyloid_PET_FBP", cohort = "ADNI",
      n_sila_total = nrow(resfit_last),
      n_estpos_true = sum(resfit_last$estpos == TRUE, na.rm = TRUE),
      n_estpos_false = sum(resfit_last$estpos == FALSE, na.rm = TRUE),
      mean_eaoa_all = round(mean(eaoa_df$EAOA, na.rm = TRUE), 1),
      sd_eaoa_all = round(sd(eaoa_df$EAOA, na.rm = TRUE), 1),
      mean_eaoa_pos = round(mean(eaoa_df$EAOA[which(eaoa_df$estpos == TRUE)], na.rm = TRUE), 1),
      sd_eaoa_pos = round(sd(eaoa_df$EAOA[which(eaoa_df$estpos == TRUE)], na.rm = TRUE), 1),
      mean_eaoa_neg = round(mean(eaoa_df$EAOA[which(eaoa_df$estpos == FALSE)], na.rm = TRUE), 1),
      sd_eaoa_neg = round(sd(eaoa_df$EAOA[which(eaoa_df$estpos == FALSE)], na.rm = TRUE), 1)
)

dat_adni <- merge(dat_surv_cn, eaoa_df, by = "RID", all.x = TRUE)

# Exclude 381_S subjects (data quality notice, March 2026)
adni_exclude_demo <- read.csv(file.path(data_2026, "PTDEMOG_11Feb2026.csv"))
adni_exclude_rids <- unique(adni_exclude_demo$RID[grepl("^381_S_", adni_exclude_demo$PTID)])
dat_adni <- dat_adni[!dat_adni$RID %in% adni_exclude_rids, ]
rm(adni_exclude_demo)

# Add demographics
dat_adni <- merge(dat_adni, dat_demo[, c("RID", "PTGENDER", "PTEDUCAT")],
                  by = "RID", all.x = TRUE)
dat_adni$Sex_F <- as.integer(dat_adni$PTGENDER == 2)
dat_adni$EDUC_z <- scale(dat_adni$PTEDUCAT)[, 1]

# Add APOE4
dat_adni <- merge(dat_adni, dat_apoe[, c("RID", "apoe4")],
                  by = "RID", all.x = TRUE)

# Keep EAOA values; estpos used by create_tvc_data

# Complete covariates
dat_adni_complete <- dat_adni[complete.cases(dat_adni[, covariates]), ]
cat(sprintf("ADNI PET: %d subjects (complete covariates), %d positive (estpos=1), %d events\n",
            nrow(dat_adni_complete),
            sum(dat_adni_complete$estpos == 1, na.rm = TRUE),
            sum(dat_adni_complete$d)))

# ADNI Plasma SILA
sila_plasma_file <- file.path(data_dir, "ADNI_plasma_SILA_intermediate.rda")
dat_adni_plasma <- NULL
if (file.exists(sila_plasma_file)) {
      sila_env <- new.env()
      load(sila_plasma_file, envir = sila_env)

      if (exists("resfit_last_fuj", envir = sila_env)) {
            rflp <- get("resfit_last_fuj", envir = sila_env)
            eaoa_col <- if ("EAOA_plasma" %in% names(rflp)) "EAOA_plasma"
                        else if ("estaget0" %in% names(rflp)) "estaget0"
                        else NULL
            if (!is.null(eaoa_col)) {
                  plasma_z_df <- data.frame(
                        RID = rflp$RID,
                        EAOA_plasma = rflp[[eaoa_col]],
                        estpos_plasma = if ("estpos" %in% names(rflp)) rflp[["estpos"]] else !is.na(rflp[[eaoa_col]])
                  )
                  sila_diag[["ADNI_Plasma"]] <- data.frame(
                        biomarker = "Plasma_pTau217_Fuji", cohort = "ADNI",
                        n_sila_total = nrow(rflp),
                        n_estpos_true = sum(plasma_z_df$estpos_plasma == TRUE, na.rm = TRUE),
                        n_estpos_false = sum(plasma_z_df$estpos_plasma == FALSE, na.rm = TRUE),
                        mean_eaoa_all = round(mean(plasma_z_df$EAOA_plasma, na.rm = TRUE), 1),
                        sd_eaoa_all = round(sd(plasma_z_df$EAOA_plasma, na.rm = TRUE), 1),
                        mean_eaoa_pos = round(mean(plasma_z_df$EAOA_plasma[which(plasma_z_df$estpos_plasma == TRUE)], na.rm = TRUE), 1),
                        sd_eaoa_pos = round(sd(plasma_z_df$EAOA_plasma[which(plasma_z_df$estpos_plasma == TRUE)], na.rm = TRUE), 1),
                        mean_eaoa_neg = round(mean(plasma_z_df$EAOA_plasma[which(plasma_z_df$estpos_plasma == FALSE)], na.rm = TRUE), 1),
                        sd_eaoa_neg = round(sd(plasma_z_df$EAOA_plasma[which(plasma_z_df$estpos_plasma == FALSE)], na.rm = TRUE), 1)
                  )

                  # Keep EAOA values; estpos_plasma used by create_tvc_data

                  # Merge into ADNI base (same cohort, different biomarker)
                  dat_adni_plasma <- merge(dat_adni_complete,
                                           plasma_z_df[, c("RID", "EAOA_plasma", "estpos_plasma")],
                                           by = "RID", all.x = TRUE)
                  cat(sprintf("ADNI Plasma: %d subjects, %d positive (estpos=1)\n",
                              nrow(dat_adni_plasma),
                              sum(dat_adni_plasma$estpos_plasma == 1, na.rm = TRUE)))
            }
      }
}


###############################################################################
## SECTION 3: Table 1 Demographics
###############################################################################

cat("\n================================================================\n")
cat("SECTION 3: TABLE 1 DEMOGRAPHICS\n")
cat("================================================================\n\n")

demo_rows <- list()

# BIOCARD CSF
demo_rows[["BIOCARD_CSF"]] <- compute_demographics(
      dat_bc_csf, "baseline.age", "onset.age", "d",
      censor_col = "censor.age", educ_raw_col = "EDUC"
)

# BIOCARD Plasma
demo_rows[["BIOCARD_Plasma"]] <- compute_demographics(
      dat_bc_plasma, "baseline.age", "onset.age", "d",
      censor_col = "censor.age", educ_raw_col = "EDUC"
)

# ADNI PET and ADNI Plasma are DISTINCT analysis cohorts (n = 575 and n = 801).
# Each cohort is loaded from its own saved analysis-cohort file. We load
# the exact analysis cohorts saved by the reanalysis scripts, so the demographics
# reproduce the manuscript Table 1.
pet_cohort_file    <- file.path(data_dir, "ADNI_pet_analysis_cohort.rda")
plasma_cohort_file <- file.path(data_dir, "ADNI_plasma_analysis_cohort.rda")
stopifnot("ADNI PET analysis cohort not found — run ADNI_SILA_reanalysis_v2.R first" =
          file.exists(pet_cohort_file))
stopifnot("ADNI plasma analysis cohort not found — run ADNI_plasma_reanalysis.R first" =
          file.exists(plasma_cohort_file))

load(pet_cohort_file)     # loads dat_complete (PET analysis cohort, n = 575)
demo_rows[["ADNI_PET"]] <- compute_demographics(
      dat_complete, "baseline.age", "onset.age", "d",
      censor_col = NULL, educ_raw_col = "PTEDUCAT"
)

load(plasma_cohort_file)  # loads dat_plasma_cohort (plasma analysis cohort, n = 801)
demo_rows[["ADNI_Plasma"]] <- compute_demographics(
      dat_plasma_cohort, "baseline.age", "onset.age", "d",
      censor_col = NULL, educ_raw_col = "PTEDUCAT"
)

table1 <- do.call(rbind, demo_rows)
table1$cohort_subset <- names(demo_rows)
table1 <- table1[, c("cohort_subset", names(table1)[names(table1) != "cohort_subset"])]

write.csv(table1, file.path(results_dir, "table1_demographics.csv"), row.names = FALSE)
cat("Saved: results/table1_demographics.csv\n")
print(table1)


###############################################################################
## SECTION 4: Person-Years by Z_tv State
###############################################################################

cat("\n================================================================\n")
cat("SECTION 4: PERSON-YEARS BY Z_tv STATE\n")
cat("================================================================\n\n")

py_rows <- list()

# BIOCARD CSF AB42/40
tvc_ab <- create_tvc_data(dat_bc_csf, "EAOA_AB", "baseline.age", "onset.age", "d",
                          estpos_col = "estpos_ab")
py_ab <- compute_person_years(tvc_ab)
py_ab$biomarker <- "CSF_AB42_AB40"
py_ab$cohort <- "BIOCARD"
py_rows[["BC_AB"]] <- py_ab

# BIOCARD CSF p-tau181
tvc_ptau <- create_tvc_data(dat_bc_csf, "EAOA_PTAU", "baseline.age", "onset.age", "d",
                            estpos_col = "estpos_ptau")
py_ptau <- compute_person_years(tvc_ptau)
py_ptau$biomarker <- "CSF_pTau181"
py_ptau$cohort <- "BIOCARD"
py_rows[["BC_PTAU"]] <- py_ptau

# BIOCARD Plasma p-tau181
tvc_plasma <- create_tvc_data(dat_bc_plasma, "EAOA_plasma", "baseline.age", "onset.age", "d",
                              estpos_col = "estpos")
py_plasma <- compute_person_years(tvc_plasma)
py_plasma$biomarker <- "Plasma_pTau181"
py_plasma$cohort <- "BIOCARD"
py_rows[["BC_Plasma"]] <- py_plasma

# ADNI Amyloid PET — use the 575-subject analysis cohort
# (the cohort the TV models are fit on).
tvc_adni_pet <- create_tvc_data(dat_complete, "EAOA", "baseline.age", "onset.age", "d",
                                estpos_col = "estpos")
py_adni_pet <- compute_person_years(tvc_adni_pet)
py_adni_pet$biomarker <- "Amyloid_PET_FBP"
py_adni_pet$cohort <- "ADNI"
py_rows[["ADNI_PET"]] <- py_adni_pet

# ADNI Plasma p-tau217 — use the 801-subject plasma analysis cohort.
if (exists("dat_plasma_cohort")) {
      tvc_adni_plasma <- create_tvc_data(dat_plasma_cohort, "age_at_positivity",
                                          "baseline.age", "onset.age", "d",
                                          estpos_col = "ever_positive")
      py_adni_plasma <- compute_person_years(tvc_adni_plasma)
      py_adni_plasma$biomarker <- "Plasma_pTau217_Fuji"
      py_adni_plasma$cohort <- "ADNI"
      py_rows[["ADNI_Plasma"]] <- py_adni_plasma
}

person_years <- do.call(rbind, py_rows)
person_years <- person_years[, c("biomarker", "cohort", "py_z0", "py_z1", "ev_z0", "ev_z1")]

write.csv(person_years, file.path(results_dir, "person_years_by_ztv.csv"), row.names = FALSE)
cat("Saved: results/person_years_by_ztv.csv\n")
print(person_years)


###############################################################################
## SECTION 5: EAOA Summary Statistics
###############################################################################

cat("\n================================================================\n")
cat("SECTION 5: EAOA SUMMARY\n")
cat("================================================================\n\n")

eaoa_rows <- list()

# BIOCARD CSF AB42/40 — use estpos to identify positive subjects
ab_pos <- dat_bc_csf$EAOA_AB[which(dat_bc_csf$estpos_ab == 1)]
eaoa_rows[["BC_AB"]] <- data.frame(
      biomarker = "CSF_AB42_AB40", cohort = "BIOCARD",
      n_pos = length(ab_pos),
      mean_eaoa = round(mean(ab_pos), 1),
      sd_eaoa = round(sd(ab_pos), 1),
      min_eaoa = round(min(ab_pos), 1),
      q25_eaoa = round(quantile(ab_pos, 0.25), 1),
      median_eaoa = round(median(ab_pos), 1),
      q75_eaoa = round(quantile(ab_pos, 0.75), 1),
      max_eaoa = round(max(ab_pos), 1)
)

# BIOCARD CSF p-tau181
ptau_pos <- dat_bc_csf$EAOA_PTAU[which(dat_bc_csf$estpos_ptau == 1)]
eaoa_rows[["BC_PTAU"]] <- data.frame(
      biomarker = "CSF_pTau181", cohort = "BIOCARD",
      n_pos = length(ptau_pos),
      mean_eaoa = round(mean(ptau_pos), 1),
      sd_eaoa = round(sd(ptau_pos), 1),
      min_eaoa = round(min(ptau_pos), 1),
      q25_eaoa = round(quantile(ptau_pos, 0.25), 1),
      median_eaoa = round(median(ptau_pos), 1),
      q75_eaoa = round(quantile(ptau_pos, 0.75), 1),
      max_eaoa = round(max(ptau_pos), 1)
)

# BIOCARD Plasma p-tau181
plasma_pos <- dat_bc_plasma$EAOA_plasma[which(dat_bc_plasma$estpos == 1)]
eaoa_rows[["BC_Plasma"]] <- data.frame(
      biomarker = "Plasma_pTau181", cohort = "BIOCARD",
      n_pos = length(plasma_pos),
      mean_eaoa = round(mean(plasma_pos), 1),
      sd_eaoa = round(sd(plasma_pos), 1),
      min_eaoa = round(min(plasma_pos), 1),
      q25_eaoa = round(quantile(plasma_pos, 0.25), 1),
      median_eaoa = round(median(plasma_pos), 1),
      q75_eaoa = round(quantile(plasma_pos, 0.75), 1),
      max_eaoa = round(max(plasma_pos), 1)
)

# ADNI Amyloid PET — use the 575-subject analysis cohort.
# estpos=1 = biomarker-positive.
adni_pet_pos <- dat_complete$EAOA[which(dat_complete$estpos == 1)]
stopifnot("ADNI PET: all estpos=1 subjects must have valid EAOA" =
          all(!is.na(adni_pet_pos)))
eaoa_rows[["ADNI_PET"]] <- data.frame(
      biomarker = "Amyloid_PET_FBP", cohort = "ADNI",
      n_pos = length(adni_pet_pos),
      mean_eaoa = round(mean(adni_pet_pos), 1),
      sd_eaoa = round(sd(adni_pet_pos), 1),
      min_eaoa = round(min(adni_pet_pos), 1),
      q25_eaoa = round(quantile(adni_pet_pos, 0.25), 1),
      median_eaoa = round(median(adni_pet_pos), 1),
      q75_eaoa = round(quantile(adni_pet_pos, 0.75), 1),
      max_eaoa = round(max(adni_pet_pos), 1)
)

# ADNI Plasma p-tau217 — use the 801-subject plasma analysis cohort
# (the cohort the TV models are fit on).
if (exists("dat_plasma_cohort")) {
      adni_plasma_pos <- dat_plasma_cohort$age_at_positivity[which(dat_plasma_cohort$ever_positive == 1)]
      if (length(adni_plasma_pos) > 0) {
            eaoa_rows[["ADNI_Plasma"]] <- data.frame(
                  biomarker = "Plasma_pTau217_Fuji", cohort = "ADNI",
                  n_pos = length(adni_plasma_pos),
                  mean_eaoa = round(mean(adni_plasma_pos), 1),
                  sd_eaoa = round(sd(adni_plasma_pos), 1),
                  min_eaoa = round(min(adni_plasma_pos), 1),
                  q25_eaoa = round(quantile(adni_plasma_pos, 0.25), 1),
                  median_eaoa = round(median(adni_plasma_pos), 1),
                  q75_eaoa = round(quantile(adni_plasma_pos, 0.75), 1),
                  max_eaoa = round(max(adni_plasma_pos), 1)
            )
      }
}

eaoa_summary <- do.call(rbind, eaoa_rows)
rownames(eaoa_summary) <- NULL

write.csv(eaoa_summary, file.path(results_dir, "eaoa_summary.csv"), row.names = FALSE)
cat("Saved: results/eaoa_summary.csv\n")
print(eaoa_summary)

# SILA estpos diagnostic: full breakdown of SILA estimation vs filtering
cat("\n--- SILA estpos diagnostic ---\n")
cat("SILA estimates a crossing age (EAOA) for every subject, but subjects whose\n")
cat("biomarker has not yet crossed the positivity threshold (estpos=FALSE) have\n")
cat("extrapolated values that can be far in the future. Only estpos=TRUE subjects\n")
cat("(biomarker estimated above threshold at last observation) are used for EAOA\n")
cat("summary statistics and downstream analyses.\n\n")

sila_diag_df <- do.call(rbind, sila_diag)
rownames(sila_diag_df) <- NULL
write.csv(sila_diag_df, file.path(results_dir, "eaoa_sila_diagnostic.csv"), row.names = FALSE)
cat("Saved: results/eaoa_sila_diagnostic.csv\n")
print(sila_diag_df)


###############################################################################
## Summary
###############################################################################

cat("\n================================================================\n")
cat("ALL DESCRIPTIVES COMPLETE\n")
cat("================================================================\n\n")

cat("Output files:\n")
cat("  results/table1_demographics.csv\n")
cat("  results/person_years_by_ztv.csv\n")
cat("  results/eaoa_summary.csv\n")
cat("  results/eaoa_sila_diagnostic.csv\n")
