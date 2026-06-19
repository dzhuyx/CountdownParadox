# =============================================================================
# Study 1: Summarize Results (standalone)
#
# Loads raw per-replicate results from all_results.rds and computes summary
# statistics. To add or modify summary statistics, edit this file and re-run
# — no need to re-run the simulation.
#
# Input:  results/study1/all_results.rds
# Output: results/study1/summary_results.csv
#
# Author: Yuxin Zhu
# Date: February 2026
# =============================================================================

library(dplyr)

# =============================================================================
# LOAD RAW RESULTS
# =============================================================================

results_dir <- file.path(getwd(), "results", "study1")
rds_path <- file.path(results_dir, "all_results.rds")

if (!file.exists(rds_path)) {
  stop("Raw results not found: ", rds_path,
       "\nRun study1_simulation.R first.")
}

cat("Loading raw results...\n")
results <- readRDS(rds_path)
cat(sprintf("  %d rows loaded (%d unique scenarios, %d methods)\n",
            nrow(results), length(unique(results$scenario)),
            length(unique(results$method))))

# =============================================================================
# SUMMARIZE
# =============================================================================

summarize_results <- function(results) {
  results %>%
    group_by(scenario, scenario_name, family, n_target, method) %>%
    summarize(
      # --- Counts ---
      n_sims           = n(),
      n_valid          = sum(!is.na(log_hr)),

      # --- Sample size ---
      mean_n_analyzed  = mean(n_analyzed, na.rm = TRUE),

      # --- Event rate ---
      mean_event_rate  = mean(event_rate, na.rm = TRUE),

      # --- Var(Z) ---
      mean_var_Z       = mean(var_Z, na.rm = TRUE),

      # --- Beta (primary coefficient) inference ---
      type1_error      = mean(p_value < 0.05, na.rm = TRUE),
      mean_p_value     = mean(p_value, na.rm = TRUE),
      median_p_value   = median(p_value, na.rm = TRUE),
      mean_hr          = mean(hr, na.rm = TRUE),
      median_hr        = median(hr, na.rm = TRUE),
      q025_hr          = quantile(hr, 0.025, na.rm = TRUE),
      q975_hr          = quantile(hr, 0.975, na.rm = TRUE),
      mean_log_hr      = mean(log_hr, na.rm = TRUE),
      bias             = mean(log_hr, na.rm = TRUE),
      empirical_se     = sd(log_hr, na.rm = TRUE),
      mean_model_se    = mean(se, na.rm = TRUE),
      se_ratio         = mean(se, na.rm = TRUE) / sd(log_hr, na.rm = TRUE),
      coverage         = mean(ci_lower <= 1 & ci_upper >= 1, na.rm = TRUE),

      # --- Gamma (interaction coefficient) inference ---
      n_valid_gamma      = sum(!is.na(log_hr_gamma)),
      type1_error_gamma  = mean(p_gamma < 0.05, na.rm = TRUE),
      mean_hr_gamma      = mean(hr_gamma, na.rm = TRUE),
      mean_log_hr_gamma  = mean(log_hr_gamma, na.rm = TRUE),
      mean_se_gamma      = mean(se_gamma, na.rm = TRUE),
      coverage_gamma     = mean(ci_lower_gamma <= 1 & ci_upper_gamma >= 1, na.rm = TRUE),

      # --- TVC-specific diagnostics ---
      mean_n_already_pos = mean(n_already_pos, na.rm = TRUE),
      mean_n_transition  = mean(n_transition, na.rm = TRUE),
      mean_n_never_pos   = mean(n_never_pos, na.rm = TRUE),
      convergence_rate   = mean(converged, na.rm = TRUE),
      mean_vcov_cor      = mean(vcov_cor_beta_gamma, na.rm = TRUE),

      # --- General diagnostics ---
      mean_cor_Z_TminusZ = mean(cor_Z_TminusZ, na.rm = TRUE),
      mean_cor_Z_T       = mean(cor_Z_T, na.rm = TRUE),
      mean_mean_Z        = mean(mean_Z, na.rm = TRUE),
      mean_sd_Z          = mean(sd_Z, na.rm = TRUE),
      mean_mean_T        = mean(mean_T, na.rm = TRUE),
      mean_sd_T          = mean(sd_T, na.rm = TRUE),
      mean_followup      = mean(mean_followup, na.rm = TRUE),
      mean_pct_excluded  = mean(pct_excluded, na.rm = TRUE),

      .groups = "drop"
    )
}

summary_all <- summarize_results(results)

# =============================================================================
# PRINT KEY RESULTS
# =============================================================================

cat("\n=============================================================\n")
cat("SUMMARY: Naive (Countdown) Analysis\n")
cat("=============================================================\n")
naive_summary <- summary_all %>%
  filter(method == "naive") %>%
  select(scenario, family, n_target, mean_var_Z, mean_event_rate,
         type1_error, mean_hr, bias, coverage,
         mean_p_value, mean_cor_Z_TminusZ, mean_followup, mean_pct_excluded)
print(as.data.frame(naive_summary), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("SUMMARY: TVC z_only Analysis\n")
cat("=============================================================\n")
tvc_z_summary <- summary_all %>%
  filter(method == "tvc_z_only") %>%
  select(scenario, family, n_target, mean_n_analyzed,
         type1_error, mean_hr, bias, coverage,
         mean_p_value, median_p_value,
         mean_n_already_pos, mean_n_transition, mean_n_never_pos,
         convergence_rate)
print(as.data.frame(tvc_z_summary), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("SUMMARY: TVC Interaction Analysis (beta)\n")
cat("=============================================================\n")
tvc_int_summary <- summary_all %>%
  filter(method == "tvc_interaction") %>%
  select(scenario, family, n_target,
         type1_error, mean_hr, bias, coverage,
         type1_error_gamma, mean_hr_gamma, coverage_gamma,
         mean_vcov_cor, convergence_rate)
print(as.data.frame(tvc_int_summary), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("COMPARISON: Type I Error — Naive vs TVC z_only vs TVC Interaction\n")
cat("=============================================================\n")
comparison <- summary_all %>%
  select(scenario, family, method, type1_error, mean_hr, coverage) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = c(type1_error, mean_hr, coverage),
    names_sep = "_"
  )
print(as.data.frame(comparison), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("DIAGNOSTIC: Cov(Z, T-Z) ~ -Var(Z) verification\n")
cat("=============================================================\n")
diag_check <- summary_all %>%
  filter(method == "naive") %>%
  select(scenario, mean_var_Z, mean_cor_Z_TminusZ, mean_cor_Z_T,
         mean_followup, mean_mean_T)
print(as.data.frame(diag_check), row.names = FALSE, digits = 4)

# =============================================================================
# SAVE
# =============================================================================

csv_path <- file.path(results_dir, "summary_results.csv")
write.csv(summary_all, csv_path, row.names = FALSE)

cat("\n=============================================================\n")
cat(sprintf("Summary saved: %s\n", csv_path))
cat(sprintf("  %d rows (%d scenarios x %d methods)\n",
            nrow(summary_all), length(unique(summary_all$scenario)),
            length(unique(summary_all$method))))
cat("=============================================================\n")
