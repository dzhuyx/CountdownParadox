# Reproducibility Pipeline

Reproduces all analyses, tables, and figures for the manuscript on the countdown paradox.

## Quick Start

```r
# 1. Set CP_BIOCARD_DIR to the folder with the raw BIOCARD data files.
# 2. Place the raw ADNI data files in <project_root>/ADNI_2026_data/
#    (both are listed under Data Dependencies). Then run:
Sys.setenv(CP_BIOCARD_DIR = "/path/to/raw/BIOCARD/data")
source("run_all.R")
```

This runs the full pipeline end-to-end (\~30 min for real-data analyses). `run_all.R` locates itself and resolves all project paths automatically: in the manuscript authors' directory tree it uses the existing `CountdownParadox_Analysis/` folder, and from a standalone clone it falls back to a self-contained layout **under the repository** (creating `CountdownParadox_Analysis/data/` and `.../results/`). Two raw-data inputs must be in place first: (1) set `CP_BIOCARD_DIR` to the raw BIOCARD data folder, and (2) put the raw ADNI data files in `<project_root>/ADNI_2026_data/`. Any location can be redirected via the environment variables listed under Requirements. Simulations are not re-run; only summarized from existing per-replicate results.

## Pipeline Overview

`run_all.R` orchestrates 6 phases in dependency order:

| Phase                       | Scripts                                                                                                                                                                       | Output                                                                                                                                       | Runtime |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- | ------- |
| 1. Data extraction & SILA   | BIOCARD\_data\_extraction.R, BIOCARD\_csf\_sila.R, BIOCARD\_plasma\_sila.R, ADNI\_pet\_sila.R, ADNI\_plasma\_sila.R, ADNI\_SILA\_reanalysis\_v2.R, ADNI\_plasma\_reanalysis.R | data/*.rda, results/ADNI\_*.csv                                                                                                              | \~5 min |
| 1b. Variable validation     | validate\_variables.R                                                                                                                                                         | console (13 checks)                                                                                                                          | \~1 sec |
| 2. Main analysis            | BIOCARD\_ADNI\_main\_analysis.R, extract\_full\_coefficients.R                                                                                                                | main\_results\_all\_biomarkers.csv, degeneracy\_all\_biomarkers.csv, sample\_sizes\_all\_biomarkers.csv, full\_coefficients\_all\_models.csv | \~5 min |
| 3b. Manuscript descriptives | compute\_manuscript\_descriptives.R, adni\_event\_types.R                                                                                                                     | table1\_demographics.csv, person\_years\_by\_ztv.csv, eaoa\_summary.csv, adni\_event\_types.csv                                              | \~1 min |
| 4. Simulation summary       | Skipped (run on cluster separately)                                                                                                                                           | results/study{1,2}/summary\_results.csv                                                                                                      | N/A     |
| 5. Tables & figures         | create\_manuscript\_tables.R, create\_manuscript\_figures.R                                                                                                                           | results/manuscript\_tables/, results/manuscript\_figures/                                                                                    | \~5 min |

## Analysis Scripts

### Phase 1: Data Extraction & SILA

| Script | Input | Output | Description |
| --------------------------- | -------------------------------------------------- | ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------- |
| `BIOCARD_data_extraction.R` | Raw BIOCARD Excel files | data/analysis\_data\_merged.rda | Extracts BIOCARD survival data (age at entry, onset/censor, covariates) from raw clinical data |
| `BIOCARD_csf_sila.R` | CSF biomarker data + analysis\_data\_merged.rda | data/BIOCARD\_CSF\_SILA\_intermediate.rda | Fits Sampled Iterative Local Approximation (SILA) to estimate age at CSF biomarker positivity (AB42/40 and p-tau181) |
| `BIOCARD_plasma_sila.R` | Plasma biomarker data + analysis\_data\_merged.rda | data/BIOCARD\_plasma\_SILA\_intermediate.rda | Fits SILA for plasma p-tau181 |
| `ADNI_pet_sila.R` | ADNI PET data files | data/ADNI\_SILA\_intermediate\_2026.rda | Fits SILA for ADNI amyloid PET |
| `ADNI_plasma_sila.R` | ADNI plasma data files | data/ADNI\_plasma\_SILA\_intermediate.rda | Fits SILA for ADNI plasma p-tau217 |
| `ADNI_SILA_reanalysis_v2.R` | data/ADNI\_SILA\_intermediate\_2026.rda | results/ADNI\_countdown\_vs\_tvc\_2026.csv, results/ADNI\_degeneracy\_2026.csv | ADNI PET countdown + TVC analysis |
| `ADNI_plasma_reanalysis.R` | data/ADNI\_plasma\_SILA\_intermediate.rda | results/ADNI\_plasma\_all\_models.csv, results/ADNI\_plasma\_degeneracy.csv | ADNI plasma countdown + TVC analysis |

### Phase 2: Main Analysis

| Script | Input | Output | Description |
| ------------------------------ | ---------------------------------------- | ------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `BIOCARD_ADNI_main_analysis.R` | Phase 1 outputs + ADNI pre-computed CSVs | main\_results\_all\_biomarkers.csv, degeneracy\_all\_biomarkers.csv, sample\_sizes\_all\_biomarkers.csv | **Canonical analysis script**. Runs standard countdown (P1), TV-BC (P2), and TV-AABC (P3) models for all 5 biomarker-cohort combinations. BIOCARD results computed live; ADNI results carried forward from pre-computed CSVs. Includes runtime assertions for onset.age/censor.age. |
| `extract_full_coefficients.R` | Phase 1 .rda files | full\_coefficients\_all\_models.csv | Re-runs all 15 Cox PH models and extracts ALL coefficients (not just target effects) for Supp Table S8. |

### Phase 1b: Variable Validation (Safeguard)

| Script | Input | Output | Description |
| ---------------------- | ----------- | ----------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `validate_variables.R` | data/\*.rda | console (13 PASS/FAIL checks) | Confirms onset.age != censor.age for BIOCARD events, onset.age == censor.age for ADNI, and all required columns exist |

### Phase 3b: Manuscript Descriptives

| Script | Input | Output | Description |
| ----------------------------------- | --------------------------------- | ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `compute_manuscript_descriptives.R` | Phase 1-2 outputs | table1\_demographics.csv, person\_years\_by\_ztv.csv, eaoa\_summary.csv | Demographics, person-years by Z\_tv state, EAOA distribution summaries |
| `adni_event_types.R` | DXSUM + ADNI analysis-cohort .rda | results/adni\_event\_types.csv | Tabulates ADNI progressors with a direct CN→dementia transition (no documented MCI) |

### Phase 4: Simulations

The simulation scripts (`study1_simulation.R`, `study2_simulation.R`, `study1_summarize.R`, `study2_summarize.R`) reside in the sibling `CountdownParadox_Manuscript_Simulations/` directory, not in this repository; `run_simulations.R` (in this repository) is the runner that invokes them. The main pipeline only summarizes their saved per-replicate results.

| Script | Input | Output | Description |
| --------------------- | ------------------------------- | ----------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `study1_simulation.R` | None (self-contained) | results/study1/all\_results.rds | Study 1: 16 null-effect scenarios (Z independent of T). 1000 replicates. Standard + TVC. |
| `study2_simulation.R` | None (self-contained) | results/study2/all\_results.rds | Study 2: 42+ scenarios with LMM trajectories. Null + alternative effects. |
| `study1_summarize.R` | results/study1/all\_results.rds | results/study1/summary\_results.csv | Computes Type I error, mean HR, coverage from per-replicate results |
| `study2_summarize.R` | results/study2/all\_results.rds | results/study2/summary\_results.csv | Computes power, direction accuracy, bias from per-replicate results |
| `run_simulations.R` | N/A | results/study{1,2}/ | Runnable local runner: executes both simulation scripts (1000 replicates each) and the summarizers (\~11–12 h on a single core). The cluster workflow is documented below. |

### Phase 5: Tables & Figures

| Script | Input | Output | Description |
| ------------------------- | --------------------------------------------------------------------------------------------- | ------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------- |
| `create_manuscript_tables.R` | main\_results\_all\_biomarkers.csv, degeneracy\_all\_biomarkers.csv, table1\_demographics.csv | results/manuscript\_tables/\*.csv | Generates Tables 1-3 for the manuscript. All values read dynamically from upstream CSVs. |
| `create_manuscript_figures.R` | Phase 2-3b outputs + simulation summaries | results/manuscript\_figures/*.pdf,*.png | Generates all main + supplementary figures. Uses Okabe-Ito colorblind-safe palette. Person-years and EAOA values read dynamically. |

### Validation Scripts

| Script                 | Description                                                           |
| ---------------------- | --------------------------------------------------------------------- |
| `validate_variables.R` | Confirms onset.age/censor.age semantics (13 checks). Run as Phase 1b. |

## Simulation Workflow on HPC Cluster

The full simulations (Study 1 + Study 2) take \~6 hours total and were run on the JHPCE cluster at Johns Hopkins.

### How simulations were run

1. **Job submission**: Each study submitted as a single-core batch job

   ```bash
   sbatch --mem=8G --time=6:00:00 --wrap="cd $SIM_DIR && Rscript study1_simulation.R"
   sbatch --mem=16G --time=8:00:00 --wrap="cd $SIM_DIR && Rscript study2_simulation.R"
   ```

2. **Seed management**: Deterministic seeds (Study 1: 12345, Study 2: 54321). Results are reproducible given the same R version (4.5.0).

3. **Output**: Per-replicate results saved to `.rds` files. These are the archival format — all summary statistics are computed post-hoc.

4. **Re-summarization**: After any change to summary statistics, run the summarizer scripts (\~5 seconds each).

## Data Dependencies

| File                                        | Source                       | Notes                               |
| ------------------------------------------- | ---------------------------- | ----------------------------------- |
| Raw BIOCARD data (Excel files) | User-provided | Folder pointed to by `CP_BIOCARD_DIR` |
| Raw ADNI data (CSV files, e.g. `DXSUM_*.csv`, `PTDEMOG_*.csv`) | User-provided | Place in `<project_root>/ADNI_2026_data/` (`ADNI_2026_data` is the required folder name) |
| `data/analysis_data_merged.rda`             | BIOCARD\_data\_extraction.R  | Regenerated by pipeline             |
| `data/BIOCARD_CSF_SILA_intermediate.rda`    | BIOCARD\_csf\_sila.R         | Regenerated by pipeline             |
| `data/BIOCARD_plasma_SILA_intermediate.rda` | BIOCARD\_plasma\_sila.R      | Regenerated by pipeline             |
| `data/ADNI_SILA_intermediate_2026.rda`      | ADNI\_pet\_sila.R            | Regenerated by pipeline             |
| `data/ADNI_plasma_SILA_intermediate.rda`    | ADNI\_plasma\_sila.R         | Regenerated by pipeline             |
| `results/ADNI_countdown_vs_tvc_2026.csv`    | ADNI\_SILA\_reanalysis\_v2.R | Carried forward by canonical script |
| `results/ADNI_plasma_all_models.csv`        | ADNI\_plasma\_reanalysis.R   | Carried forward by canonical script |

## Requirements

- **R version**: 4.5.0

- **Key packages** (real-data analysis, tables, figures): survival, dplyr, tibble, tidyr, ggplot2, patchwork, scales, readxl

- **SILA estimation**: `silaR` (v0.0.0.9000) — an R port (by M. Bilgel) of the SILA algorithm of Betthauser et al. Not on CRAN; install from GitHub:

  ```r
  # install.packages("remotes")
  remotes::install_github("Betthauser-Neuro-Lab/silaR")
  ```

- **Simulations** (scripts in `CountdownParadox_Manuscript_Simulations/`): additionally require lme4 and MASS

- **Directory structure**: `run_all.R` resolves the project root automatically — it uses the authors' layout (`CountdownParadox_Analysis/` and `CountdownParadox_Manuscript_Simulations/` as siblings under a common project root) when present, and otherwise treats the repository folder itself as the project root (a self-contained clone). No paths are hardcoded; all I/O is relative to the resolved root.

- **Environment variables** — `run_all.R` respects any of these that are already set and derives the rest, so the pipeline runs from a checkout in any location:

  - `CP_BIOCARD_DIR` — **required**: folder with the raw BIOCARD data files.
  - `CP_PROJECT_ROOT` — optional: where analysis I/O lives (`data/`, `results/`, and the raw `ADNI_2026_data/`). Defaults to `CountdownParadox_Analysis` under the resolved project root (i.e. under the repository for a standalone clone).
  - `CP_PROJECT_DIR` — optional: the parent holding `CountdownParadox_Analysis/` and `CountdownParadox_Manuscript_Simulations/`; set it to relocate both.
  - `CP_SIM_DIR` — optional (used by `run_simulations.R`): the `CountdownParadox_Manuscript_Simulations` folder.
  - `CP_SCRIPT_DIR` — optional: overrides the detected location of this repository.

## License

The code in this repository is released under the MIT License (see `LICENSE`).
The `silaR` package and the other R packages listed above are distributed
separately under their own licenses.
