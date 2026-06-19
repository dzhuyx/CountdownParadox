# ==============================================================================
# run_simulations.R — Local simulation runner for external readers
#
# Runs both manuscript simulation scripts in full mode (1000 replicates each)
# and executes the summarizers. Reproduces the per-replicate coefficient
# estimates and summary statistics used in the manuscript.
#
# Audience:
#   External readers who have cloned the repository and installed R
#   dependencies, and who want to regenerate simulation outputs on a
#   single local machine without cluster access.
#
# Wall clock:
#   Approximately 11–12 hours on a single modern CPU core. Study 1 is fast
#   (~1 hour); Study 2 dominates (~10–11 hours) and is driven by the n=2000
#   PED twin scenario. There is no within-R parallelization in this runner.
#
# Manuscript authors:
#   The published results were produced on the JHPCE cluster as a 75-task
#   SLURM array distributed alongside the simulation scripts. This local runner
#   is deterministic given the seed and R version — it produces the same
#   per-replicate estimates as the cluster run.
#
# Usage:
#   # Interactive (will prompt before running):
#   Rscript run_simulations.R
#
#   # CI / non-interactive:
#   RUN_FULL_SIMS=1 Rscript run_simulations.R
#
# R version used for published results: 4.4.0
# Required packages: survival, dplyr, lme4, MASS
#
# Author: Yuxin Zhu
# Date: April 2026
# ==============================================================================

# -- Locate the simulations directory ----------------------------------------
# The simulation scripts (study1_simulation.R, study2_simulation.R, and the
# summarizers) are bundled in this repository. By default they are sourced from
# the current working directory (run this from the repository root); set
# CP_SIM_DIR to override, e.g. to point at a separate working copy.
sim_dir <- Sys.getenv("CP_SIM_DIR")
if (sim_dir == "" || !dir.exists(sim_dir)) {
  sim_dir <- getwd()
}
.sim_scripts <- c("study1_simulation.R", "study2_simulation.R",
                  "study1_summarize.R", "study2_summarize.R")
.missing <- .sim_scripts[!file.exists(file.path(sim_dir, .sim_scripts))]
if (length(.missing)) {
  stop(sprintf(paste0("Simulation scripts not found in '%s': %s.\n",
                      "Run this from the repository root, or set CP_SIM_DIR to ",
                      "the folder that holds them."),
               sim_dir, paste(.missing, collapse = ", ")))
}

# -- Banner + confirmation gate ----------------------------------------------
cat("============================================================\n")
cat("Countdown Paradox — local simulation runner\n")
cat("============================================================\n")
cat(sprintf("Simulations directory: %s\n", sim_dir))
cat("\nThis will run:\n")
cat("  1. study1_simulation.R  (16 configs x 1000 reps, ~1 hour)\n")
cat("  2. study2_simulation.R  (59 scenarios x 1000 reps, ~10-11 hours)\n")
cat("  3. study1_summarize.R   (post-hoc summary, ~5 sec)\n")
cat("  4. study2_summarize.R   (post-hoc summary, ~5 sec)\n")
cat("\nTotal wall clock: ~11-12 hours on a single core.\n")
cat("Outputs: results/study{1,2}/all_results.rds and summary_results.csv\n")
cat("============================================================\n\n")

confirm_env <- Sys.getenv("RUN_FULL_SIMS", unset = "")
if (!nzchar(confirm_env)) {
  if (interactive()) {
    ans <- readline("Proceed with full simulation run? (type 'yes' to continue): ")
    if (!identical(tolower(trimws(ans)), "yes")) {
      cat("Aborted by user.\n")
      quit(save = "no", status = 0)
    }
  } else {
    cat("Non-interactive run detected. Set RUN_FULL_SIMS=1 to skip this prompt.\n")
    cat("Aborting. Example:  RUN_FULL_SIMS=1 Rscript run_simulations.R\n")
    quit(save = "no", status = 0)
  }
}

# -- Set working directory so getwd()-relative results paths resolve ---------
old_wd <- getwd()
setwd(sim_dir)
on.exit(setwd(old_wd), add = TRUE)

t_master_start <- Sys.time()

# ==============================================================================
# Study 1
# ==============================================================================
cat("\n=== Study 1: study1_simulation.R ===\n")
t1 <- Sys.time()
source("study1_simulation.R", echo = FALSE)
t2 <- Sys.time()
cat(sprintf(">> Study 1 wall clock: %.1f min\n\n",
            as.numeric(difftime(t2, t1, units = "mins"))))

cat("=== Study 1 summary: study1_summarize.R ===\n")
source("study1_summarize.R", echo = FALSE)
cat("\n")

# ==============================================================================
# Study 2
# ==============================================================================
cat("=== Study 2: study2_simulation.R ===\n")
t3 <- Sys.time()
source("study2_simulation.R", echo = FALSE)
t4 <- Sys.time()
cat(sprintf(">> Study 2 wall clock: %.1f min\n\n",
            as.numeric(difftime(t4, t3, units = "mins"))))

cat("=== Study 2 summary: study2_summarize.R ===\n")
source("study2_summarize.R", echo = FALSE)

# ==============================================================================
t_master_end <- Sys.time()
cat("\n============================================================\n")
cat(sprintf("Total wall clock: %.2f hours\n",
            as.numeric(difftime(t_master_end, t_master_start, units = "hours"))))
cat(sprintf("Outputs in %s/results/study{1,2}/\n", sim_dir))
cat("  - all_results.rds      (per-replicate coefficient estimates)\n")
cat("  - summary_results.csv  (aggregate: Type I error, bias, coverage, power)\n")
cat("============================================================\n")
