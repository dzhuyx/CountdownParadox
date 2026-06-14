# ==============================================================================
# run_all.R — Master reproducibility pipeline
#
# Runs all real-data analyses for the manuscript end-to-end.
# Does NOT duplicate code — sources scripts from this reproducibility folder.
#
# Usage: Rscript run_all.R   (from any directory)
#
# Phases:
#   1. Data extraction & SILA fitting (BIOCARD + ADNI, from raw data)
#   1b. Variable validation (safeguard)
#   2. Main analysis (canonical 5-biomarker × 3-model, carries forward ADNI results)
#   3b. Manuscript descriptives (demographics, person-years, EAOA summary)
#   5. Tables and figures
#
# Note: Simulation summarization (Phase 4) is run separately on cluster.
#       See run_simulations.R for documentation.
#
# Author: Yuxin Zhu
# Date: April 2026
# ==============================================================================

cat("=== Countdown Paradox: Manuscript Reproducibility Pipeline ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -- Paths -------------------------------------------------------------------
# All analysis scripts live in this folder. Locate it portably (works with
# `Rscript run_all.R` or source()), derive the project root from it, and
# expose both to the child scripts via environment variables.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("^--file=", "", fa[length(fa)]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(getwd())
}
script_dir   <- Sys.getenv("CP_SCRIPT_DIR")
if (script_dir == "") script_dir <- get_script_dir()
project_dir  <- dirname(dirname(script_dir))
project_root <- file.path(project_dir, "CountdownParadox_Analysis")
Sys.setenv(CP_PROJECT_DIR = project_dir, CP_PROJECT_ROOT = project_root)

cat("Script directory:", script_dir, "\n\n")

# Helper: source a script in an isolated environment so rm(list=ls())
# inside the script doesn't clear our variables here.
run_script <- function(script_name) {
  path <- file.path(script_dir, script_name)
  if (!file.exists(path)) {
    cat("  ERROR: Script not found:", path, "\n")
    return(invisible(FALSE))
  }
  tryCatch({
    env <- new.env(parent = globalenv())
    sys.source(path, envir = env)
    cat("  Done.\n")
    invisible(TRUE)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    invisible(FALSE)
  })
}

# ==============================================================================
# Phase 1: Data Extraction & SILA (BIOCARD + ADNI)
# ==============================================================================
cat("=== Phase 1: Data Extraction & SILA ===\n")

# 1a. BIOCARD data extraction → data/analysis_data_merged.rda
cat("[1/7] BIOCARD_data_extraction.R...\n")
run_script("BIOCARD_data_extraction.R")

# 1b. BIOCARD CSF SILA → data/BIOCARD_CSF_SILA_intermediate.rda
cat("[2/7] BIOCARD_csf_sila.R...\n")
run_script("BIOCARD_csf_sila.R")

# 1c. BIOCARD Plasma SILA → data/BIOCARD_plasma_SILA_intermediate.rda
cat("[3/7] BIOCARD_plasma_sila.R...\n")
run_script("BIOCARD_plasma_sila.R")

# 1d. ADNI PET SILA → data/ADNI_SILA_intermediate_2026.rda
cat("[4/7] ADNI_pet_sila.R...\n")
run_script("ADNI_pet_sila.R")

# 1e. ADNI Plasma SILA → data/ADNI_plasma_SILA_intermediate.rda
cat("[5/7] ADNI_plasma_sila.R...\n")
run_script("ADNI_plasma_sila.R")

# 1f. ADNI PET reanalysis → results/ADNI_countdown_vs_tvc_2026.csv, results/ADNI_degeneracy_2026.csv
cat("[6/7] ADNI_SILA_reanalysis_v2.R...\n")
run_script("ADNI_SILA_reanalysis_v2.R")

# 1g. ADNI Plasma reanalysis → results/ADNI_plasma_all_models.csv, results/ADNI_plasma_degeneracy.csv
cat("[7/7] ADNI_plasma_reanalysis.R...\n")
run_script("ADNI_plasma_reanalysis.R")

cat("\n")

# ==============================================================================
# Phase 1b: Variable Validation (Safeguard)
# ==============================================================================
cat("=== Phase 1b: Variable Validation ===\n")
cat("[1/1] validate_variables.R...\n")
run_script("validate_variables.R")
cat("\n")

# ==============================================================================
# Phase 2: Main Analysis (Canonical)
# ==============================================================================
cat("=== Phase 2: Main Analysis ===\n")

# BIOCARD_ADNI_main_analysis.R produces:
#   results/main_results_all_biomarkers.csv
#   results/degeneracy_all_biomarkers.csv
#   results/sample_sizes_all_biomarkers.csv
cat("[1/2] BIOCARD_ADNI_main_analysis.R...\n")
run_script("BIOCARD_ADNI_main_analysis.R")

# extract_full_coefficients.R produces:
#   results/full_coefficients_all_models.csv  (Supp Table 8 source)
cat("[2/2] extract_full_coefficients.R...\n")
run_script("extract_full_coefficients.R")
cat("\n")

# ==============================================================================
# Phase 3b: Manuscript Descriptives
# ==============================================================================
cat("=== Phase 3b: Manuscript Descriptives ===\n")

# compute_manuscript_descriptives.R produces:
#   results/table1_demographics.csv
#   results/person_years_by_ztv.csv
#   results/eaoa_summary.csv
cat("[1/2] compute_manuscript_descriptives.R...\n")
run_script("compute_manuscript_descriptives.R")

# adni_event_types.R produces:
#   results/adni_event_types.csv  (progressor-without-MCI counts: PET 7/106, plasma 7/162)
cat("[2/2] adni_event_types.R...\n")
run_script("adni_event_types.R")
cat("\n")

# ==============================================================================
# Phase 4: Simulation Summarization — SKIPPED
# ==============================================================================
cat("=== Phase 4: Simulation Summarization — SKIPPED ===\n")
cat("  Simulations are re-run on cluster. See run_simulations.R.\n")
cat("  Pre-existing summary CSVs in CountdownParadox_Manuscript_Simulations/results/ are used.\n\n")

# ==============================================================================
# Phase 5: Tables & Figures
# ==============================================================================
cat("=== Phase 5: Tables & Figures ===\n")

cat("[1/3] create_manuscript_tables.R...\n")
run_script("create_manuscript_tables.R")

cat("[2/3] create_manuscript_figures.R...\n")
run_script("create_manuscript_figures.R")

# Manuscript Figure 3 (HR as a function of AABC, 5-panel).
cat("[3/3] figure_hr_aabc_panel.R...\n")
run_script("figure_hr_aabc_panel.R")

cat("\n=== Pipeline complete ===\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
