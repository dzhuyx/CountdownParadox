###############################################################################
## ADNI PET SILA Trajectory Fitting — Amyloid PET (FBP tracer)
##
## Purpose: Apply SILA algorithm to ADNI amyloid PET data (2026 download) to
##   estimate amyloid onset age (EAOA = estaget0) for each subject.
##
## SILA trajectory fitting for ADNI amyloid PET. Countdown and landmark
## analyses are performed separately in ADNI_SILA_reanalysis_v2.R.
##
## Input:  ADNI_2026_data/*.csv (raw ADNI data files)
## Output: data/ADNI_SILA_intermediate_2026.rda
##         (contains: res_sila, resfit, resfit_last, sila_df,
##          dat_surv_cn, dat_av45, dat_demo, dat_apoe, first_pos_scan)
##
## Author: Yuxin Zhu
## Date: February 2026 (reproducibility copy April 2026)
###############################################################################

rm(list = ls())
library(survival)
library(silaR)
library(tibble)

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set CP_PROJECT_ROOT to the analysis root (the folder that holds data/ and results/).")
data_2026    <- file.path(project_root, "ADNI_2026_data")
data_dir     <- file.path(project_root, "data")
res_dir      <- file.path(project_root, "results")


###############################################################################
## Helper Functions
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

extract_hr <- function(fit, var_name) {
      fit_s <- summary(fit)
      idx <- which(rownames(fit_s$coefficients) == var_name)
      if (length(idx) == 0) return(data.frame(
            hr = NA, lower95 = NA, upper95 = NA, pvalue = NA,
            n = fit_s$n, nevent = fit_s$nevent
      ))
      data.frame(
            hr = fit_s$conf.int[idx, 1],
            lower95 = fit_s$conf.int[idx, 3],
            upper95 = fit_s$conf.int[idx, 4],
            pvalue = fit_s$coefficients[idx, 5],
            n = fit_s$n,
            nevent = fit_s$nevent
      )
}


###############################################################################
## SECTION 1: Setup and Data Loading (2026 Data)
###############################################################################

cat("================================================================\n")
cat("SECTION 1: SETUP AND DATA LOADING (2026 DATA)\n")
cat("================================================================\n\n")

# --- Amyloid PET (2026: UCBERKELEY_AMY_6MM, multi-tracer) --- #
dat_av45 <- read.csv(file.path(data_2026, "All_Subjects_UCBERKELEY_AMY_6MM_15Feb2026.csv"),
                      na.strings = -4)
dat_av45$SCANDATE <- as.Date(dat_av45$SCANDATE, "%Y-%m-%d")
cat(sprintf("Amyloid PET (all tracers): %d rows, %d subjects\n",
            nrow(dat_av45), length(unique(dat_av45$RID))))

# Report tracer distribution before filtering
cat("  Tracer distribution:\n")
print(table(dat_av45$TRACER, useNA = "ifany"))

# Report QC distribution before filtering
cat("  QC flag distribution:\n")
print(table(dat_av45$qc_flag, useNA = "ifany"))

# Filter to FBP tracer (primary analysis — consistent with old FBP-only file)
dat_av45 <- dat_av45[dat_av45$TRACER == "FBP", ]
cat(sprintf("After FBP filter: %d rows, %d subjects\n",
            nrow(dat_av45), length(unique(dat_av45$RID))))

# QC filter: Pass (2) or Partial pass (1)
dat_av45 <- dat_av45[!is.na(dat_av45$qc_flag) & dat_av45$qc_flag >= 1, ]
cat(sprintf("After QC filter (qc_flag >= 1): %d rows, %d subjects\n",
            nrow(dat_av45), length(unique(dat_av45$RID))))

# --- Demographics (2026: PTDOB replaces PTDOBYY+PTDOBMM) --- #
dat_demo <- read.csv(file.path(data_2026, "PTDEMOG_11Feb2026.csv"), na.strings = -4)
dat_demo$DOB <- as.Date(paste0("01/", dat_demo$PTDOB), "%d/%m/%Y")
dat_demo$update_stamp <- as.Date(dat_demo$update_stamp, "%Y-%m-%d")
dat_demo <- do.call(rbind, lapply(split(dat_demo, dat_demo$RID), function(x) {
      x[x$update_stamp == max(x$update_stamp, na.rm = TRUE), ][1, ]
}))
cat(sprintf("Demographics: %d unique subjects\n", nrow(dat_demo)))

# --- Diagnosis (2026: DIAGNOSIS pre-mapped, replaces DXCURREN+DXCHANGE) --- #
dat_dx <- read.csv(file.path(data_2026, "DXSUM_11Feb2026.csv"), na.strings = -4)
dat_dx$EXAMDATE <- as.Date(dat_dx$EXAMDATE, "%Y-%m-%d")
dat_dx$dx <- dat_dx$DIAGNOSIS
cat(sprintf("Diagnosis: %d rows\n", nrow(dat_dx)))

# --- ApoE (2026: GENOTYPE text field replaces APGEN1+APGEN2) --- #
dat_apoe <- read.csv(file.path(data_2026, "All_Subjects_APOERES_15Feb2026.csv"), na.strings = -4)
dat_apoe$apoe4 <- as.integer(grepl("4", dat_apoe$GENOTYPE))
# Already 1 row per subject in 2026 file, but deduplicate just in case
dat_apoe <- do.call(rbind, lapply(split(dat_apoe, dat_apoe$RID), function(x) x[1, ]))
cat(sprintf("ApoE: %d unique subjects\n", nrow(dat_apoe)))

# --- ARM (enrollment diagnosis) --- #
dat_arm <- read.csv(file.path(data_2026, "ARM_11Feb2026.csv"), na.strings = -4)
dat_arm$dxarm <- NA
dat_arm$dxarm[dat_arm$ARM %in% c(1, 4, 7)] <- 1  # NL
dat_arm$dxarm[dat_arm$ARM %in% c(2, 5, 8)] <- 2  # MCI
dat_arm$dxarm[dat_arm$ARM %in% c(3, 6, 9)] <- 3  # AD
# Note: ADNI4 ARM codes 10-11 will have dxarm = NA (filtered downstream)
cat(sprintf("ARM: %d rows\n", nrow(dat_arm)))

# --- Add age to amyloid PET (using SCANDATE) --- #
dat_av45 <- merge(dat_av45, dat_demo[, c("RID", "DOB")], by = "RID", all.x = TRUE)
dat_av45$age_at_scan <- as.numeric(dat_av45$SCANDATE - dat_av45$DOB) / 365.25

# --- Construct survival data --- #
dat_dx <- merge(dat_dx, dat_demo[, c("RID", "DOB")], by = "RID", all.x = TRUE)
dat_dx$age <- as.numeric(dat_dx$EXAMDATE - dat_dx$DOB) / 365.25

dat_surv <- do.call(rbind, lapply(split(dat_dx, dat_dx$RID), function(x) {
      x <- x[complete.cases(x[, c("dx", "age")]), , drop = FALSE]
      if (nrow(x) == 0) return(NULL)
      x <- x[order(x$EXAMDATE), , drop = FALSE]
      baseline.age <- x$age[1]
      baseline.dx  <- x$dx[1]
      if (rev(x$dx)[1] %in% c(2, 3)) {
            d <- 1
            onset.pos <- which(x$dx %in% c(2, 3))[1]
            onset.age <- x$age[onset.pos]
      } else {
            d <- 0
            onset.age <- rev(x$age)[1]
      }
      data.frame(RID = x$RID[1],
                 baseline.age = baseline.age,
                 baseline.dx = baseline.dx,
                 onset.age = onset.age,
                 time = onset.age - baseline.age,
                 d = d)
}))

# Deduplicate ARM and merge enrollment dx
dat_arm_dedup <- do.call(rbind, lapply(split(dat_arm, dat_arm$RID), function(x) {
      x_valid <- x[!is.na(x$dxarm), ]
      if (nrow(x_valid) > 0) return(x_valid[1, ]) else return(x[1, ])
}))
dat_surv <- merge(dat_surv, dat_arm_dedup[, c("RID", "dxarm")], by = "RID", all.x = TRUE)

# Restrict to CN at baseline (first dx visit == 1)
dat_surv_cn <- dat_surv[dat_surv$baseline.dx == 1, ]
# Also use censor.age = onset.age for consistency with BIOCARD code
dat_surv_cn$censor.age <- dat_surv_cn$onset.age

cat(sprintf("\nSurvival data: %d total, %d CN at baseline\n", nrow(dat_surv), nrow(dat_surv_cn)))
cat(sprintf("  Events: %d (%.1f%%)\n", sum(dat_surv_cn$d == 1),
            mean(dat_surv_cn$d == 1) * 100))
cat(sprintf("  Age range: [%.1f, %.1f]\n", min(dat_surv_cn$onset.age),
            max(dat_surv_cn$onset.age)))

# --- Identify subjects with >=2 amyloid PET scans --- #
n_av45_per_subj <- tapply(dat_av45$RID, dat_av45$RID, length)
rids_ge2 <- as.integer(names(n_av45_per_subj[n_av45_per_subj >= 2]))
rids_cn_ge2 <- intersect(rids_ge2, dat_surv_cn$RID)

cat(sprintf("  Amyloid PET >= 2 scans: %d subjects\n", length(rids_ge2)))
cat(sprintf("  Amyloid PET >= 2 scans + CN at baseline: %d subjects\n", length(rids_cn_ge2)))

cat("\n")


###############################################################################
## SECTION 2: Prepare silaR Input
###############################################################################

cat("================================================================\n")
cat("SECTION 2: PREPARE silaR INPUT\n")
cat("================================================================\n\n")

suvr_var <- "SUMMARY_SUVR"

# Filter to subjects with >=2 scans (all subjects, not just CN — SILA trains on whole population)
dat_sila_input <- dat_av45[dat_av45$RID %in% rids_ge2 &
                                 !is.na(dat_av45[[suvr_var]]) &
                                 !is.na(dat_av45$age_at_scan), ]

# Create silaR input: subid (numeric), age, val
sila_df <- tibble(
      subid = as.numeric(dat_sila_input$RID),
      age   = dat_sila_input$age_at_scan,
      val   = dat_sila_input[[suvr_var]]
)

# Remove duplicate ages within subject (silaR requirement)
sila_df <- sila_df[!duplicated(sila_df[, c("subid", "age")]), ]

cat(sprintf("silaR input: %d observations from %d subjects\n",
            nrow(sila_df), length(unique(sila_df$subid))))
cat(sprintf("  SUVR: mean = %.3f, SD = %.3f, range = [%.3f, %.3f]\n",
            mean(sila_df$val), sd(sila_df$val),
            min(sila_df$val), max(sila_df$val)))
cat(sprintf("  Age: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(sila_df$age), sd(sila_df$age),
            min(sila_df$age), max(sila_df$age)))

cat("\n")


###############################################################################
## SECTION 3: Train SILA Model
###############################################################################

cat("================================================================\n")
cat("SECTION 3: TRAIN SILA MODEL\n")
cat("================================================================\n\n")

cutoff_av45 <- 1.11  # FBP positivity threshold (val0 in SILA)

cat("Training SILA model (dt=0.25, val0=1.11, maxi=200)...\n")
res_sila <- sila(sila_df, dt = 0.25, val0 = cutoff_av45, maxi = 200)

# Check convergence
cat(sprintf("  tsila: %d rows (trajectory curve points)\n", nrow(res_sila$tsila)))
cat(sprintf("  adtime range: [%.2f, %.2f] years from threshold\n",
            min(res_sila$tsila$adtime), max(res_sila$tsila$adtime)))
cat(sprintf("  val range: [%.3f, %.3f]\n",
            min(res_sila$tsila$val), max(res_sila$tsila$val)))
cat(sprintf("  nsubs range: [%d, %d] (subjects contributing per value)\n",
            min(res_sila$tsila$nsubs), max(res_sila$tsila$nsubs)))

cat("\n")


###############################################################################
## SECTION 4: Estimate Individual Onset Ages
###############################################################################

cat("================================================================\n")
cat("SECTION 4: ESTIMATE INDIVIDUAL ONSET AGES (EAOA)\n")
cat("================================================================\n\n")

# Build estimation input from ALL subjects with valid PET data (>= 1 scan),
# not just the >= 2 training set. SILA was trained on >= 2, but
# sila_estimate() only uses each subject's last observation for alignment.
dat_sila_all <- dat_av45[!is.na(dat_av45[[suvr_var]]) &
                         !is.na(dat_av45$age_at_scan), ]
sila_df_all <- tibble(
      subid = as.numeric(dat_sila_all$RID),
      age   = dat_sila_all$age_at_scan,
      val   = dat_sila_all[[suvr_var]]
)
sila_df_all <- sila_df_all[!duplicated(sila_df_all[, c("subid", "age")]), ]
cat(sprintf("Estimating for ALL valid subjects: %d obs from %d subjects (trained on %d)\n",
            nrow(sila_df_all), length(unique(sila_df_all$subid)),
            length(unique(sila_df$subid))))

cat("Running sila_estimate(align_event='last')...\n")
resfit <- sila_estimate(res_sila$tsila, sila_df_all, align_event = "last")

cat(sprintf("  sila_estimate output: %d rows\n", nrow(resfit)))

# Extract one row per subject: use the last observation's estaget0
resfit_last <- do.call(rbind, lapply(split(resfit, resfit$subid), function(x) {
      x <- x[x$age == max(x$age), ]
      x[1, ]
}))

cat(sprintf("  Unique subjects: %d\n", nrow(resfit_last)))

# estaget0 = SILA-estimated age at threshold crossing (our Z = EAOA)
n_valid_eaoa <- sum(!is.na(resfit_last$estaget0))
cat(sprintf("  Valid EAOA (non-NA): %d (%.1f%%)\n", n_valid_eaoa,
            n_valid_eaoa / nrow(resfit_last) * 100))

eaoa_vals <- resfit_last$estaget0[!is.na(resfit_last$estaget0)]
cat(sprintf("  EAOA distribution: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(eaoa_vals), sd(eaoa_vals), min(eaoa_vals), max(eaoa_vals)))

# How many are positive (estpos)?
n_pos <- sum(resfit_last$estpos == TRUE, na.rm = TRUE)
cat(sprintf("  SILA-estimated positive (estpos = TRUE): %d (%.1f%%)\n",
            n_pos, n_pos / nrow(resfit_last) * 100))

# Compare with first-positive-scan age
first_pos_scan <- do.call(rbind, lapply(split(dat_av45, dat_av45$RID), function(x) {
      x <- x[complete.cases(x[, c(suvr_var, "age_at_scan")]), ]
      x <- x[order(x$SCANDATE), ]
      pos_idx <- which(x[[suvr_var]] >= cutoff_av45)
      if (length(pos_idx) == 0) return(NULL)
      data.frame(RID = x$RID[1], first_pos_age = x$age_at_scan[pos_idx[1]])
}))

resfit_last$RID <- as.integer(resfit_last$subid)
resfit_last <- merge(resfit_last, first_pos_scan, by = "RID", all.x = TRUE)

both_valid <- !is.na(resfit_last$estaget0) & !is.na(resfit_last$first_pos_age)
if (sum(both_valid) > 5) {
      cor_val <- cor(resfit_last$estaget0[both_valid],
                     resfit_last$first_pos_age[both_valid])
      cat(sprintf("  Correlation SILA EAOA vs first-positive-scan age: r = %.3f (n = %d)\n",
                  cor_val, sum(both_valid)))
      cat(sprintf("  Mean EAOA - first_pos: %.1f years (SILA estimates earlier onset)\n",
                  mean(resfit_last$estaget0[both_valid] - resfit_last$first_pos_age[both_valid])))
}

cat("\n")

# Save intermediate results
save(res_sila, resfit, resfit_last, sila_df, dat_surv_cn, dat_av45, dat_demo, dat_apoe,
     first_pos_scan,
     file = file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"))
cat("Saved: data/ADNI_SILA_intermediate_2026.rda\n")

cat("\n=== ADNI PET SILA complete ===\n")
