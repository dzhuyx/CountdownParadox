# =============================================================================
# Study 2: Summarize Results (standalone)
#
# Loads raw per-replicate results from all_results.rds and computes summary
# statistics. To add or modify summary statistics, edit this file and re-run
# — no need to re-run the simulation.
#
# Input:  results/study2/all_results.rds
# Output: results/study2/summary_results.csv
#
# Author: Yuxin Zhu
# Date: February 2026
# =============================================================================

library(dplyr)

# =============================================================================
# LOAD RAW RESULTS
# =============================================================================

results_dir <- file.path(getwd(), "results", "study2")
rds_path <- file.path(results_dir, "all_results.rds")

if (!file.exists(rds_path)) {
  stop("Raw results not found: ", rds_path,
       "\nRun study2_simulation.R first.")
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
    group_by(scenario, scenario_name, family, n_target, method,
             beta_true, true_hr) %>%
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
      reject_rate      = mean(p_value < 0.05, na.rm = TRUE),
      mean_p_value     = mean(p_value, na.rm = TRUE),
      median_p_value   = median(p_value, na.rm = TRUE),
      mean_hr          = mean(hr, na.rm = TRUE),
      median_hr        = median(hr, na.rm = TRUE),
      q025_hr          = quantile(hr, 0.025, na.rm = TRUE),
      q975_hr          = quantile(hr, 0.975, na.rm = TRUE),
      mean_log_hr      = mean(log_hr, na.rm = TRUE),
      empirical_se     = sd(log_hr, na.rm = TRUE),
      mean_model_se    = mean(se, na.rm = TRUE),
      se_ratio         = mean(se, na.rm = TRUE) / sd(log_hr, na.rm = TRUE),

      # --- Coverage of HR=1 ---
      coverage_null    = mean(ci_lower <= 1 & ci_upper >= 1, na.rm = TRUE),

      # --- Direction accuracy ---
      pct_hr_above_1   = mean(hr > 1, na.rm = TRUE),
      pct_hr_below_1   = mean(hr < 1, na.rm = TRUE),

      # --- Gamma (interaction coefficient) inference ---
      n_valid_gamma      = sum(!is.na(log_hr_gamma)),
      reject_rate_gamma  = mean(p_gamma < 0.05, na.rm = TRUE),
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

      # --- Correlations ---
      mean_cor_Z_TminusZ = mean(cor_Z_TminusZ, na.rm = TRUE),
      mean_cor_Z_T       = mean(cor_Z_T, na.rm = TRUE),
      mean_cor_Z_Y       = mean(cor_Z_Y, na.rm = TRUE),
      mean_cor_Y_T       = mean(cor_Y_T, na.rm = TRUE),

      # --- Descriptive ---
      mean_mean_Z        = mean(mean_Z, na.rm = TRUE),
      mean_sd_Z          = mean(sd_Z, na.rm = TRUE),
      mean_mean_T        = mean(mean_T, na.rm = TRUE),
      mean_sd_T          = mean(sd_T, na.rm = TRUE),
      mean_followup      = mean(mean_followup, na.rm = TRUE),

      # --- LMM recovery ---
      mean_beta_0_hat    = mean(beta_0_hat, na.rm = TRUE),
      mean_beta_1_hat    = mean(beta_1_hat, na.rm = TRUE),
      mean_sigma_0_hat   = mean(sigma_0_hat, na.rm = TRUE),
      mean_sigma_1_hat   = mean(sigma_1_hat, na.rm = TRUE),
      mean_rho_hat       = mean(rho_hat, na.rm = TRUE),
      mean_sigma_eps_hat = mean(sigma_eps_hat, na.rm = TRUE),

      # --- Z estimation ---
      mean_n_Z_hat_NA    = mean(n_Z_hat_NA, na.rm = TRUE),
      mean_pct_Z_hat_NA  = mean(pct_Z_hat_NA, na.rm = TRUE),
      mean_visits        = mean(mean_visits, na.rm = TRUE),

      .groups = "drop"
    )
}

summary_all <- summarize_results(results)

# =============================================================================
# PRINT KEY RESULTS
# =============================================================================

cat("\n=============================================================\n")
cat("SUMMARY: Naive Analysis — Effect Size Sweep\n")
cat("=============================================================\n")
naive_effect <- summary_all %>%
  filter(method == "naive",
         grepl("^S3[a-g]-(ped|bio)$", scenario)) %>%
  select(scenario, family, beta_true, mean_event_rate,
         reject_rate, mean_hr, mean_log_hr, pct_hr_above_1)
print(as.data.frame(naive_effect), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("SUMMARY: TVC z_only — Effect Size Sweep\n")
cat("=============================================================\n")
tvc_z_effect <- summary_all %>%
  filter(method == "tvc_z_only",
         grepl("^S3[a-g]-(ped|bio)$", scenario)) %>%
  select(scenario, family, beta_true, mean_n_analyzed,
         reject_rate, mean_hr, mean_log_hr,
         coverage_null, pct_hr_above_1,
         mean_n_already_pos, mean_n_transition, mean_n_never_pos)
print(as.data.frame(tvc_z_effect), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("SUMMARY: TVC Interaction — Effect Size Sweep (beta)\n")
cat("=============================================================\n")
tvc_int_effect <- summary_all %>%
  filter(method == "tvc_interaction",
         grepl("^S3[a-g]-(ped|bio)$", scenario)) %>%
  select(scenario, family, beta_true,
         reject_rate, mean_hr, mean_log_hr,
         reject_rate_gamma, mean_hr_gamma,
         mean_vcov_cor, convergence_rate)
print(as.data.frame(tvc_int_effect), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("COMPARISON: Direction Accuracy — Naive vs TVC z_only\n")
cat("=============================================================\n")
direction <- summary_all %>%
  filter(method %in% c("naive", "tvc_z_only"),
         grepl("^S3[a-g]-(ped|bio)$", scenario)) %>%
  select(scenario, family, beta_true, method, mean_hr, pct_hr_above_1, pct_hr_below_1) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = c(mean_hr, pct_hr_above_1, pct_hr_below_1),
    names_sep = "_"
  )
print(as.data.frame(direction), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("SUMMARY: TVC z_only Power by Sample Size (beta=0.5)\n")
cat("=============================================================\n")
power_n <- summary_all %>%
  filter(method == "tvc_z_only",
         scenario %in% c("S3h-ped", "S3b-ped", "S3i-ped", "S3j-ped",
                          "S3h-bio", "S3b-bio", "S3i-bio")) %>%
  select(scenario, family, n_target, mean_n_analyzed,
         reject_rate, mean_hr, coverage_null)
print(as.data.frame(power_n), row.names = FALSE, digits = 3)

cat("\n=============================================================\n")
cat("DIAGNOSTIC: Correlation Structure by Effect Size\n")
cat("=============================================================\n")
cor_diag <- summary_all %>%
  filter(method == "naive",
         grepl("^S3[a-g]-(ped|bio)$", scenario)) %>%
  select(scenario, family, beta_true,
         mean_cor_Z_T, mean_cor_Z_Y, mean_cor_Y_T,
         mean_cor_Z_TminusZ)
print(as.data.frame(cor_diag), row.names = FALSE, digits = 4)

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
