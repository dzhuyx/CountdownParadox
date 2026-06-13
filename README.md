# Reproducibility Pipeline

Reproduces all analyses, tables, and figures for the Nature Medicine manuscript on the countdown paradox.

## Quick Start

```r
source("run_all.R")
```

This runs the full pipeline end-to-end (~30 min for real-data analyses). Simulations are not re-run; only summarized from existing per-replicate results.

## Pipeline Overview

`run_all.R` orchestrates 6 phases in dependency order:

| Phase | Scripts | Output | Runtime |
|-------|---------|--------|---------|
| 1. Data extraction & SILA | BIOCARD_data_extraction.R, BIOCARD_csf_sila.R, BIOCARD_plasma_sila.R, ADNI_pet_sila.R, ADNI_plasma_sila.R, ADNI_SILA_reanalysis_v2.R, ADNI_plasma_reanalysis.R | data/*.rda, results/ADNI_*.csv | ~5 min |
| 1b. Variable validation | validate_variables.R | console (13 checks) | ~1 sec |
| 2. Main analysis | BIOCARD_ADNI_main_analysis.R, extract_full_coefficients.R | main_results_all_biomarkers.csv, degeneracy_all_biomarkers.csv, sample_sizes_all_biomarkers.csv, full_coefficients_all_models.csv | ~5 min |
| 3b. Manuscript descriptives | compute_manuscript_descriptives.R, adni_event_types.R | table1_demographics.csv, person_years_by_ztv.csv, eaoa_summary.csv, adni_event_types.csv | ~1 min |
| 4. Simulation summary | Skipped (run on cluster separately) | results/study{1,2}/summary_results.csv | N/A |
| 5. Tables & figures | create_natmed_tables.R, create_natmed_figures.R | results/manuscript_tables/, results/manuscript_figures/ | ~5 min |

## Analysis Scripts

### Phase 1: Data Extraction & SILA

| Script | Location | Input | Output | Description |
|--------|----------|-------|--------|-------------|
| `BIOCARD_data_extraction.R` | reproducibility/ | Raw BIOCARD Excel files | data/analysis_data_merged.rda | Extracts BIOCARD survival data (age at entry, onset/censor, covariates) from raw clinical data |
| `BIOCARD_csf_sila.R` | reproducibility/ | CSF biomarker data + analysis_data_merged.rda | data/BIOCARD_CSF_SILA_intermediate.rda | Fits Sampled Iterative Local Approximation (SILA) to estimate age at CSF biomarker positivity (AB42/40 and p-tau181) |
| `BIOCARD_plasma_sila.R` | reproducibility/ | Plasma biomarker data + analysis_data_merged.rda | data/BIOCARD_plasma_SILA_intermediate.rda | Fits SILA for plasma p-tau181 |
| `ADNI_pet_sila.R` | reproducibility/ | ADNI PET data files | data/ADNI_SILA_intermediate_2026.rda | Fits SILA for ADNI amyloid PET |
| `ADNI_plasma_sila.R` | reproducibility/ | ADNI plasma data files | data/ADNI_plasma_SILA_intermediate.rda | Fits SILA for ADNI plasma p-tau217 |
| `ADNI_SILA_reanalysis_v2.R` | reproducibility/ | data/ADNI_SILA_intermediate_2026.rda | results/ADNI_countdown_vs_tvc_2026.csv, results/ADNI_degeneracy_2026.csv | ADNI PET countdown + TVC analysis |
| `ADNI_plasma_reanalysis.R` | reproducibility/ | data/ADNI_plasma_SILA_intermediate.rda | results/ADNI_plasma_all_models.csv, results/ADNI_plasma_degeneracy.csv | ADNI plasma countdown + TVC analysis |

### Phase 2: Main Analysis

| Script | Location | Input | Output | Description |
|--------|----------|-------|--------|-------------|
| `BIOCARD_ADNI_main_analysis.R` | reproducibility/ | Phase 1 outputs + ADNI pre-computed CSVs | main_results_all_biomarkers.csv, degeneracy_all_biomarkers.csv, sample_sizes_all_biomarkers.csv | **Canonical analysis script**. Runs standard countdown (P1), TV-BC (P2), and TV-AABC (P3) models for all 5 biomarker-cohort combinations. BIOCARD results computed live; ADNI results carried forward from pre-computed CSVs. Includes runtime assertions for onset.age/censor.age. |
| `extract_full_coefficients.R` | reproducibility/ | Phase 1 .rda files | full_coefficients_all_models.csv | Re-runs all 15 Cox PH models and extracts ALL coefficients (not just target effects) for Supp Table S8. |

### Phase 1b: Variable Validation (Safeguard)

| Script | Location | Input | Output | Description |
|--------|----------|-------|--------|-------------|
| `validate_variables.R` | reproducibility/ | data/*.rda | console (13 PASS/FAIL checks) | Confirms onset.age != censor.age for BIOCARD events, onset.age == censor.age for ADNI, and all required columns exist |

### Phase 3b: Manuscript Descriptives

| Script | Location | Input | Output | Description |
|--------|----------|-------|--------|-------------|
| `compute_manuscript_descriptives.R` | reproducibility/ | Phase 1-2 outputs | table1_demographics.csv, person_years_by_ztv.csv, eaoa_summary.csv | Demographics, person-years by Z_tv state, EAOA distribution summaries |
| `adni_event_types.R` | reproducibility/ | DXSUM + ADNI analysis-cohort .rda | results/adni_event_types.csv | Tabulates ADNI progressors with a direct CN→dementia transition (no documented MCI) |

### Phase 4: Simulations

| Script | Location | Input | Output | Description |
|--------|----------|-------|--------|-------------|
| `study1_simulation.R` | CountdownParadox_Manuscript_Simulations/ | None (self-contained) | results/study1/all_results.rds | Study 1: 16 null-effect scenarios (Z independent of T). 1000 replicates. Standard + TVC. |
| `study2_simulation.R` | CountdownParadox_Manuscript_Simulations/ | None (self-contained) | results/study2/all_results.rds | Study 2: 42+ scenarios with LMM trajectories. Null + alternative effects. |
| `study1_summarize.R` | CountdownParadox_Manuscript_Simulations/ | results/study1/all_results.rds | results/study1/summary_results.csv | Computes Type I error, mean HR, coverage from per-replicate results |
| `study2_summarize.R` | CountdownParadox_Manuscript_Simulations/ | results/study2/all_results.rds | results/study2/summary_results.csv | Computes power, direction accuracy, bias from per-replicate results |
| `run_simulations.R` | reproducibility/ | N/A | results/study{1,2}/ | Runnable local runner: executes both simulation scripts (1000 replicates each) and the summarizers (~11–12 h on a single core). The cluster workflow is documented below. |

### Phase 5: Tables & Figures

| Script | Location | Input | Output | Description |
|--------|----------|-------|--------|-------------|
| `create_natmed_tables.R` | reproducibility/ | main_results_all_biomarkers.csv, degeneracy_all_biomarkers.csv, table1_demographics.csv | CountdownParadox_Analysis/results/manuscript_tables/*.csv | Generates Tables 1-3 for the manuscript. All values read dynamically from upstream CSVs. |
| `create_natmed_figures.R` | reproducibility/ | Phase 2-3b outputs + simulation summaries | results/manuscript_figures/*.pdf,*.png | Generates all main + supplementary figures. Uses Okabe-Ito colorblind-safe palette. Person-years and EAOA values read dynamically. |

### Validation Scripts

| Script | Description |
|--------|-------------|
| `validate_variables.R` | Confirms onset.age/censor.age semantics (13 checks). Run as Phase 1b. |

## Simulation Workflow on HPC Cluster

The full simulations (Study 1 + Study 2) take ~6 hours total and were run on the JHPCE cluster at Johns Hopkins.

### How simulations were run

1. **Job submission**: Each study submitted as a single-core batch job
   ```bash
   sbatch --mem=8G --time=6:00:00 --wrap="cd $SIM_DIR && Rscript study1_simulation.R"
   sbatch --mem=16G --time=8:00:00 --wrap="cd $SIM_DIR && Rscript study2_simulation.R"
   ```

2. **Seed management**: Deterministic seeds (Study 1: 12345, Study 2: 54321). Results are reproducible given the same R version (4.5.0).

3. **Output**: Per-replicate results saved to `.rds` files. These are the archival format — all summary statistics are computed post-hoc.

4. **Re-summarization**: After any change to summary statistics, run the summarizer scripts (~5 seconds each). No need to re-run the full simulations.

### Why simulations are not re-run in the pipeline

- They are independent of the real-data analyses
- Full runs take ~6 hours
- Per-replicate results are saved and can be re-summarized at any time

## Data Dependencies

| File | Source | Notes |
|------|--------|-------|
| `data/analysis_data_merged.rda` | BIOCARD_data_extraction.R | Regenerated by pipeline |
| `data/BIOCARD_CSF_SILA_intermediate.rda` | BIOCARD_csf_sila.R | Regenerated by pipeline |
| `data/BIOCARD_plasma_SILA_intermediate.rda` | BIOCARD_plasma_sila.R | Regenerated by pipeline |
| `data/ADNI_SILA_intermediate_2026.rda` | ADNI_pet_sila.R | Regenerated by pipeline |
| `data/ADNI_plasma_SILA_intermediate.rda` | ADNI_plasma_sila.R | Regenerated by pipeline |
| `results/ADNI_countdown_vs_tvc_2026.csv` | ADNI_SILA_reanalysis_v2.R | Carried forward by canonical script |
| `results/ADNI_plasma_all_models.csv` | ADNI_plasma_reanalysis.R | Carried forward by canonical script |

## Requirements

- **R version**: 4.5.0
- **Key packages**: survival, dplyr, tidyr, ggplot2, patchwork, scales, lme4, MASS
- **Directory structure**: Scripts expect the standard project layout with `CountdownParadox_Analysis/` and `CountdownParadox_Manuscript_Simulations/` as siblings under the project root.
