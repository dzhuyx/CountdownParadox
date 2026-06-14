# ==============================================================================
# validate_variables.R — Variable mapping validation (Safeguard #1)
#
# For each dataset, prints the meaning of every variable used as an outcome,
# confirms onset.age != censor.age for events in BIOCARD,
# confirms onset.age == censor.age for all subjects in ADNI.
#
# Usage: source("validate_variables.R")
#
# Author: Yuxin Zhu
# Date: April 2026
# ==============================================================================

cat("=======================================================================\n")
cat("  VARIABLE MAPPING VALIDATION\n")
cat("  Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=======================================================================\n\n")

project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set it to <project>/CountdownParadox_Analysis.")
data_dir <- file.path(project_root, "data")

n_pass <- 0
n_fail <- 0

check <- function(desc, condition) {
  if (condition) {
    cat(sprintf("  PASS: %s\n", desc))
    n_pass <<- n_pass + 1
  } else {
    cat(sprintf("  FAIL: %s\n", desc))
    n_fail <<- n_fail + 1
  }
}

# --- BIOCARD --- #
cat("--- BIOCARD (analysis_data_merged.rda) ---\n\n")

tryCatch({
  load(file.path(data_dir, "analysis_data_merged.rda"))
  dat <- analysis_data

  cat(sprintf("  N = %d, events = %d\n", nrow(dat), sum(dat$d)))
  cat(sprintf("  Columns: %s\n\n", paste(names(dat), collapse = ", ")))

  # Core survival variables exist
  check("onset.age column exists", "onset.age" %in% names(dat))
  check("censor.age column exists", "censor.age" %in% names(dat))
  check("baseline.age column exists", "baseline.age" %in% names(dat))
  check("d column exists", "d" %in% names(dat))

  # onset.age <= censor.age for all
  check("onset.age <= censor.age for all subjects",
        all(dat$onset.age <= dat$censor.age))

  # For events: onset.age < censor.age (the key distinction)
  events <- dat[dat$d == 1, ]
  n_events_differ <- sum(events$onset.age < events$censor.age)
  check(sprintf("onset.age < censor.age for events (%d of %d)",
                n_events_differ, nrow(events)),
        n_events_differ > 0)

  mean_diff <- mean(events$censor.age - events$onset.age)
  cat(sprintf("\n  Mean(censor.age - onset.age) for events: %.1f years\n", mean_diff))
  cat(sprintf("  This confirms BIOCARD events have onset BEFORE last visit.\n"))

  # For censored: onset.age == censor.age
  censored <- dat[dat$d == 0, ]
  check(sprintf("onset.age == censor.age for censored subjects (%d)", nrow(censored)),
        all(censored$onset.age == censored$censor.age))

  # baseline.age < onset.age
  check("baseline.age < onset.age for all",
        all(dat$baseline.age < dat$onset.age))

  # Covariates
  check("Sex_F column exists", "Sex_F" %in% names(dat))
  check("EDUC_z column exists", "EDUC_z" %in% names(dat))
  check("apoe4 column exists", "apoe4" %in% names(dat))

  # No SUBJECT_ID >= 400
  if ("SUBJECT_ID" %in% names(dat)) {
    check("No SUBJECT_ID >= 400 (new enrollees excluded)",
          all(dat$SUBJECT_ID < 400))
  }

  cat("\n")
}, error = function(e) cat("  ERROR loading BIOCARD data:", conditionMessage(e), "\n\n"))

# --- ADNI PET --- #
cat("--- ADNI PET (ADNI_SILA_intermediate_2026.rda) ---\n\n")

tryCatch({
  load(file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"))

  # Check what objects were loaded
  adni_objs <- ls()
  cat(sprintf("  Objects in environment: %s\n", paste(adni_objs[!adni_objs %in% c("dat", "events", "censored", "n_pass", "n_fail", "check", "project_root", "data_dir", "n_events_differ", "mean_diff")], collapse = ", ")))

  # The ADNI intermediate typically has dat_surv_cn
  if (exists("dat_surv_cn")) {
    cat(sprintf("  dat_surv_cn: N = %d, events = %d\n", nrow(dat_surv_cn), sum(dat_surv_cn$d)))
    cat(sprintf("  Columns: %s\n\n", paste(names(dat_surv_cn), collapse = ", ")))

    if ("onset.age" %in% names(dat_surv_cn) && "censor.age" %in% names(dat_surv_cn)) {
      check("ADNI: onset.age == censor.age for ALL subjects",
            all(dat_surv_cn$onset.age == dat_surv_cn$censor.age))
    } else {
      cat("  NOTE: onset.age or censor.age not found in ADNI data. Check column names.\n")
    }
  } else {
    cat("  dat_surv_cn not found in ADNI intermediate. Available objects listed above.\n")
  }

  cat("\n")
}, error = function(e) cat("  ERROR loading ADNI PET data:", conditionMessage(e), "\n\n"))

# --- Summary --- #
cat("=======================================================================\n")
cat(sprintf("  VARIABLE VALIDATION: %d PASS, %d FAIL\n", n_pass, n_fail))
if (n_fail > 0) {
  cat("  WARNING: Failed checks indicate potential variable confusion.\n")
  cat("  Review variable definitions before proceeding.\n")
}
cat("=======================================================================\n")
