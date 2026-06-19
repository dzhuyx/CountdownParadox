# =============================================================================
# Study 1: Manuscript Simulation — Direct Generation (Null Effect)
#
# Demonstrates that regressing T-Z on Z produces spurious associations when
# Z and T are independent, regardless of their marginal distributions.
# The time-varying covariate (TVC) analysis correctly maintains Type I error.
#
# 15 configurations: 9 base scenarios (S1-S8 pedagogical + S9 BIOCARD)
#   plus 6 sample-size variants (S1 and S9 at n=150, 200, 1000).
#
# Methods compared:
#   1. Naive/countdown: Surv(T-Z, event) ~ Z — always biased
#   2. TVC z_only: Surv(tstart, tstop, event) ~ Z_tv — valid
#   3. TVC interaction: Surv(tstart, tstop, event) ~ Z_tv + A:Z_tv — valid
#
# Source: Adapted from archived study1_simulation.R (landmark replaced with TVC)
# Specification: manuscript_simulation_plan.md
#
# Author: Yuxin Zhu
# Date: February 2026
# =============================================================================

library(survival)
library(dplyr)

# =============================================================================
# SECTION 1: MODE CONTROL
# =============================================================================

MODE <- "full"  # "prototype" or "full"
# prototype: n_sim=10, configs=c("S1","S9")  (~10 sec)
# full:      n_sim=1000, all 15 configs       (~60-70 min)

if (MODE == "prototype") {
  N_SIM <- 10
  CONFIG_SUBSET <- c("S1", "S9")
  cat("=== PROTOTYPE MODE: 10 reps, S1 + S9 only ===\n\n")
} else {
  N_SIM <- 1000
  CONFIG_SUBSET <- NULL  # Run all
  cat("=== FULL MODE: 1000 reps, all 15 configs ===\n\n")
}

SEED_BASE <- 12345
POOL_MULT <- 2  # Oversample multiplier for all Study 1 scenarios

# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

rtruncnorm <- function(n, mean, sd, lower, upper) {
  x <- rnorm(n, mean, sd)
  while (any(x < lower | x > upper)) {
    bad <- x < lower | x > upper
    x[bad] <- rnorm(sum(bad), mean, sd)
  }
  x
}

# =============================================================================
# SECTION 3: SCENARIO DEFINITIONS (15 configurations)
# =============================================================================

scenarios <- list(

  # --- Pedagogical Scenarios (S1-S8) ---

  S1 = list(
    name = "S1: Uniform-Uniform",
    generate_Z = function(n) runif(n, 50, 75),
    generate_T = function(n) runif(n, 60, 90),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  S2 = list(
    name = "S2: Normal-Normal",
    generate_Z = function(n) rtruncnorm(n, 62, 5, 45, 80),
    generate_T = function(n) rtruncnorm(n, 75, 7, 55, 95),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  S3 = list(
    name = "S3: Normal-Weibull",
    generate_Z = function(n) rtruncnorm(n, 62, 5, 45, 80),
    generate_T = function(n) {
      T_raw <- rweibull(n, shape = 5, scale = 78)
      T_raw[T_raw < 50] <- 50 + runif(sum(T_raw < 50), 0, 5)
      T_raw
    },
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  S4 = list(
    name = "S4: Gamma-Weibull",
    generate_Z = function(n) {
      z <- rgamma(n, shape = 20, rate = 20/22) + 40
      pmin(pmax(z, 45), 85)
    },
    generate_T = function(n) 55 + rweibull(n, shape = 3, scale = 20),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  S5 = list(
    name = "S5: Beta-Exponential",
    generate_Z = function(n) 45 + 35 * rbeta(n, 2, 5),
    generate_T = function(n) pmin(60 + rexp(n, rate = 0.05), 95),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  S6 = list(
    name = "S6: Light Censoring",
    generate_Z = function(n) rtruncnorm(n, 62, 5, 45, 80),
    generate_T = function(n) rtruncnorm(n, 75, 7, 55, 95),
    entry_min = 50, entry_max = 60, duration = 35,
    max_age = 95, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  S7 = list(
    name = "S7: Heavy Censoring",
    generate_Z = function(n) rtruncnorm(n, 62, 5, 45, 80),
    generate_T = function(n) rtruncnorm(n, 75, 7, 55, 95),
    entry_min = 60, entry_max = 75, duration = 12,
    max_age = 90, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  S8 = list(
    name = "S8: High Var(Z)",
    generate_Z = function(n) rtruncnorm(n, 62, 10, 35, 90),
    generate_T = function(n) rtruncnorm(n, 75, 7, 55, 95),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 500, family = "PED",
    pool_mult = POOL_MULT
  ),

  # --- BIOCARD-Calibrated Scenario (S9) ---

  S9 = list(
    name = "S9: BIOCARD-calibrated",
    generate_Z = function(n) rtruncnorm(n, 53, 10, 25, 82),
    generate_T = function(n) rtruncnorm(n, 73, 10, 45, 96),
    entry_min = 40, entry_max = 65, duration = 20,
    max_age = 95, n = 500, family = "BIO",
    pool_mult = POOL_MULT
  ),

  # --- Sample Size Variants ---

  S1_n150 = list(
    name = "S1_n150: Uniform-Uniform (n=150)",
    generate_Z = function(n) runif(n, 50, 75),
    generate_T = function(n) runif(n, 60, 90),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 150, family = "PED",
    pool_mult = POOL_MULT
  ),

  S1_n200 = list(
    name = "S1_n200: Uniform-Uniform (n=200)",
    generate_Z = function(n) runif(n, 50, 75),
    generate_T = function(n) runif(n, 60, 90),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 200, family = "PED",
    pool_mult = POOL_MULT
  ),

  S1_n1000 = list(
    name = "S1_n1000: Uniform-Uniform (n=1000)",
    generate_Z = function(n) runif(n, 50, 75),
    generate_T = function(n) runif(n, 60, 90),
    entry_min = 55, entry_max = 70, duration = 20,
    max_age = 90, n = 1000, family = "PED",
    pool_mult = POOL_MULT
  ),

  S9_n150 = list(
    name = "S9_n150: BIOCARD-calibrated (n=150)",
    generate_Z = function(n) rtruncnorm(n, 53, 10, 25, 82),
    generate_T = function(n) rtruncnorm(n, 73, 10, 45, 96),
    entry_min = 40, entry_max = 65, duration = 20,
    max_age = 95, n = 150, family = "BIO",
    pool_mult = POOL_MULT
  ),

  S9_n200 = list(
    name = "S9_n200: BIOCARD-calibrated (n=200)",
    generate_Z = function(n) rtruncnorm(n, 53, 10, 25, 82),
    generate_T = function(n) rtruncnorm(n, 73, 10, 45, 96),
    entry_min = 40, entry_max = 65, duration = 20,
    max_age = 95, n = 200, family = "BIO",
    pool_mult = POOL_MULT
  ),

  S9_n1000 = list(
    name = "S9_n1000: BIOCARD-calibrated (n=1000)",
    generate_Z = function(n) rtruncnorm(n, 53, 10, 25, 82),
    generate_T = function(n) rtruncnorm(n, 73, 10, 45, 96),
    entry_min = 40, entry_max = 65, duration = 20,
    max_age = 95, n = 1000, family = "BIO",
    pool_mult = POOL_MULT
  )
)

# =============================================================================
# SECTION 4: CORE FUNCTIONS
# =============================================================================

generate_data <- function(n, scenario, seed = NULL) {
  # Oversample-filter-subsample design (manuscript plan Section 2.3)
  if (!is.null(seed)) set.seed(seed)

  pool_mult <- scenario$pool_mult
  n_pool <- n * pool_mult

  # Step 1: Generate pool
  Z_pool     <- scenario$generate_Z(n_pool)
  T_pool     <- scenario$generate_T(n_pool)
  entry_pool <- runif(n_pool, scenario$entry_min, scenario$entry_max)

  # Step 2: Filter — exclude subjects with T <= entry_age
  eligible <- T_pool > entry_pool
  n_eligible <- sum(eligible)

  # Pool exhaustion handling
  current_mult <- pool_mult
  while (n_eligible < n) {
    current_mult <- current_mult + 1
    warning(sprintf("Pool exhausted (eligible=%d < n=%d). Increasing pool_mult to %d.",
                    n_eligible, n, current_mult))
    n_pool <- n * current_mult
    if (!is.null(seed)) set.seed(seed)
    Z_pool     <- scenario$generate_Z(n_pool)
    T_pool     <- scenario$generate_T(n_pool)
    entry_pool <- runif(n_pool, scenario$entry_min, scenario$entry_max)
    eligible   <- T_pool > entry_pool
    n_eligible <- sum(eligible)
  }

  # Step 3: Randomly sample n eligible subjects
  idx <- sample(which(eligible), n)
  Z               <- Z_pool[idx]
  T_event         <- T_pool[idx]
  study_entry_age <- entry_pool[idx]

  pct_excluded <- 1 - mean(eligible)

  # Step 4: Administrative censoring
  censor_age <- pmin(study_entry_age + scenario$duration, scenario$max_age)
  R <- T_event - Z
  observed_T <- pmin(T_event, censor_age)
  event <- as.integer(T_event <= censor_age)

  dat <- data.frame(
    id = 1:n,
    Z = Z,
    T_event = T_event,
    R = R,
    study_entry_age = study_entry_age,
    censor_age = censor_age,
    observed_T = observed_T,
    event = event
  )
  attr(dat, "pct_excluded") <- pct_excluded
  dat
}

run_naive_analysis <- function(dat) {
  dat_naive <- dat %>%
    filter(T_event > Z) %>%
    mutate(time_from_Z = observed_T - Z,
           event_naive = as.integer(event == 1 & T_event > Z)) %>%
    filter(time_from_Z > 0)

  if (nrow(dat_naive) < 10) {
    return(list(hr = NA, log_hr = NA, se = NA, p_value = NA,
                ci_lower = NA, ci_upper = NA, n_analyzed = nrow(dat_naive)))
  }

  fit <- tryCatch(coxph(Surv(time_from_Z, event_naive) ~ Z, data = dat_naive),
                  error = function(e) NULL)

  if (is.null(fit)) {
    return(list(hr = NA, log_hr = NA, se = NA, p_value = NA,
                ci_lower = NA, ci_upper = NA, n_analyzed = nrow(dat_naive)))
  }

  coef_summary <- summary(fit)$coefficients
  ci <- confint(fit)

  list(hr = exp(coef_summary["Z", "coef"]),
       log_hr = coef_summary["Z", "coef"],
       se = coef_summary["Z", "se(coef)"],
       p_value = coef_summary["Z", "Pr(>|z|)"],
       ci_lower = exp(ci["Z", 1]),
       ci_upper = exp(ci["Z", 2]),
       n_analyzed = nrow(dat_naive))
}

run_tvc_analysis <- function(dat, formulation, A_col, entry_col, exit_col, event_col) {
  # Build start-stop data and fit TVC Cox model.
  # Adapted from tvc_investigation.R run_tvc_analysis() with CI computation added.
  #
  # formulation: "z_only" or "interaction"
  # Column names: A_col = age at positivity, entry_col = study entry age,
  #               exit_col = observed exit age, event_col = event indicator

  A         <- dat[[A_col]]
  entry_age <- dat[[entry_col]]
  exit_age  <- dat[[exit_col]]
  event_raw <- dat[[event_col]]

  # Build start-stop rows
  rows <- list()
  n_already_pos <- 0
  n_transition  <- 0
  n_never_pos   <- 0

  for (i in seq_along(A)) {
    ai <- A[i]
    ei <- entry_age[i]
    xi <- exit_age[i]
    di <- event_raw[i]

    if (is.na(ai) || ai >= xi) {
      # Never positive during follow-up
      rows[[length(rows) + 1]] <- data.frame(
        id = i, tstart = ei, tstop = xi, event = di, Z_tv = 0, A = 0
      )
      n_never_pos <- n_never_pos + 1
    } else if (ai <= ei) {
      # Already positive at entry
      rows[[length(rows) + 1]] <- data.frame(
        id = i, tstart = ei, tstop = xi, event = di, Z_tv = 1, A = ai
      )
      n_already_pos <- n_already_pos + 1
    } else {
      # Transitions during follow-up (ei < ai < xi)
      rows[[length(rows) + 1]] <- data.frame(
        id = i, tstart = ei, tstop = ai, event = 0, Z_tv = 0, A = ai
      )
      rows[[length(rows) + 1]] <- data.frame(
        id = i, tstart = ai, tstop = xi, event = di, Z_tv = 1, A = ai
      )
      n_transition <- n_transition + 1
    }
  }

  dat_long <- do.call(rbind, rows)
  # Remove zero-length intervals
  dat_long <- dat_long[dat_long$tstop > dat_long$tstart, ]
  n_analyzed <- length(A)

  # Standardize A among positive subjects — unique subjects only (C1 fix)
  pos_A <- unique(dat_long[dat_long$A > 0, c("id", "A")])$A
  if (length(pos_A) > 1) {
    A_mean <- mean(pos_A)
    A_sd   <- sd(pos_A)
    dat_long$A_z <- ifelse(dat_long$A > 0, (dat_long$A - A_mean) / A_sd, 0)
  } else {
    dat_long$A_z <- dat_long$A
  }

  # Initialize result
  result <- list(
    hr = NA, log_hr = NA, se = NA, p_value = NA,
    ci_lower = NA, ci_upper = NA,
    hr_gamma = NA, log_hr_gamma = NA, se_gamma = NA, p_gamma = NA,
    ci_lower_gamma = NA, ci_upper_gamma = NA,
    n_analyzed = n_analyzed,
    n_already_pos = n_already_pos,
    n_transition = n_transition,
    n_never_pos = n_never_pos,
    converged = FALSE,
    vcov_cor_beta_gamma = NA
  )

  # Check degeneracy: all person-time is Z_tv = 1
  all_positive <- all(dat_long$Z_tv == 1)

  if (all_positive && formulation == "z_only") {
    # Z_tv is constant — cannot estimate beta
    return(result)
  }

  if (all_positive && formulation == "interaction") {
    # Z_tv = 1 for all, only A_z effect estimable
    tryCatch({
      fit <- coxph(Surv(tstart, tstop, event) ~ A_z, data = dat_long)
      fit_s <- summary(fit)
      ci <- confint(fit)
      result$converged     <- TRUE
      result$hr_gamma      <- exp(fit_s$coefficients["A_z", 1])
      result$log_hr_gamma  <- fit_s$coefficients["A_z", 1]
      result$se_gamma      <- fit_s$coefficients["A_z", 3]
      result$p_gamma       <- fit_s$coefficients["A_z", 5]
      result$ci_lower_gamma <- exp(ci["A_z", 1])
      result$ci_upper_gamma <- exp(ci["A_z", 2])
    }, error = function(e) NULL)
    return(result)
  }

  # Fit TVC model
  tryCatch({
    if (formulation == "interaction") {
      fit <- coxph(Surv(tstart, tstop, event) ~ Z_tv + A_z:Z_tv, data = dat_long)
    } else {
      fit <- coxph(Surv(tstart, tstop, event) ~ Z_tv, data = dat_long)
    }

    fit_s <- summary(fit)
    ci <- confint(fit)
    result$converged <- TRUE

    # Beta (Z_tv) coefficient
    result$hr      <- exp(fit_s$coefficients["Z_tv", 1])
    result$log_hr  <- fit_s$coefficients["Z_tv", 1]
    result$se      <- fit_s$coefficients["Z_tv", 3]
    result$p_value <- fit_s$coefficients["Z_tv", 5]
    result$ci_lower <- exp(ci["Z_tv", 1])
    result$ci_upper <- exp(ci["Z_tv", 2])

    # Gamma (interaction) coefficient
    if (formulation == "interaction") {
      coef_names <- rownames(fit_s$coefficients)
      gamma_row  <- setdiff(coef_names, "Z_tv")
      result$hr_gamma      <- exp(fit_s$coefficients[gamma_row, 1])
      result$log_hr_gamma  <- fit_s$coefficients[gamma_row, 1]
      result$se_gamma      <- fit_s$coefficients[gamma_row, 3]
      result$p_gamma       <- fit_s$coefficients[gamma_row, 5]
      result$ci_lower_gamma <- exp(ci[gamma_row, 1])
      result$ci_upper_gamma <- exp(ci[gamma_row, 2])
      # Vcov correlation between beta and gamma
      vc <- vcov(fit)
      if (nrow(vc) == 2) {
        result$vcov_cor_beta_gamma <- cov2cor(vc)[1, 2]
      }
    }
  }, error = function(e) NULL)

  result
}

run_single_sim <- function(n, scenario, seed) {
  dat <- generate_data(n = n, scenario = scenario, seed = seed)

  # Compute diagnostics
  var_Z         <- var(dat$Z)
  event_rate    <- mean(dat$event)
  cor_Z_TminusZ <- cor(dat$Z, dat$T_event - dat$Z)
  cor_Z_T       <- cor(dat$Z, dat$T_event)
  mean_Z        <- mean(dat$Z)
  sd_Z          <- sd(dat$Z)
  mean_T        <- mean(dat$T_event)
  sd_T          <- sd(dat$T_event)
  mean_followup <- mean(dat$observed_T - dat$study_entry_age)
  pct_excluded  <- attr(dat, "pct_excluded")

  diag <- list(
    var_Z = var_Z, event_rate = event_rate,
    cor_Z_TminusZ = cor_Z_TminusZ, cor_Z_T = cor_Z_T,
    mean_Z = mean_Z, sd_Z = sd_Z,
    mean_T = mean_T, sd_T = sd_T,
    mean_followup = mean_followup, pct_excluded = pct_excluded
  )

  # TVC-specific NA placeholders for naive method
  tvc_na <- list(
    hr_gamma = NA, log_hr_gamma = NA, se_gamma = NA, p_gamma = NA,
    ci_lower_gamma = NA, ci_upper_gamma = NA,
    n_already_pos = NA, n_transition = NA, n_never_pos = NA,
    converged = NA, vcov_cor_beta_gamma = NA
  )

  # --- Naive analysis ---
  naive_result         <- run_naive_analysis(dat)
  naive_result$method  <- "naive"
  naive_result <- c(naive_result, tvc_na, diag)

  # --- TVC z_only ---
  tvc_z         <- run_tvc_analysis(dat, "z_only", "Z", "study_entry_age", "observed_T", "event")
  tvc_z$method  <- "tvc_z_only"
  tvc_z <- c(tvc_z, diag)

  # --- TVC interaction ---
  tvc_int         <- run_tvc_analysis(dat, "interaction", "Z", "study_entry_age", "observed_T", "event")
  tvc_int$method  <- "tvc_interaction"
  tvc_int <- c(tvc_int, diag)

  all_results <- list(naive_result, tvc_z, tvc_int)
  bind_rows(lapply(all_results, as.data.frame))
}

run_scenario_simulation <- function(config_id, scenario, n_sim, seed_base,
                                    verbose = TRUE) {
  n <- scenario$n

  if (verbose) cat(sprintf("\nRunning %s (n=%d, %d reps)...\n",
                           scenario$name, n, n_sim))

  results_list <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    if (verbose && i %% 200 == 0) cat(sprintf("  Simulation %d / %d\n", i, n_sim))
    seed <- seed_base + i
    res <- run_single_sim(n = n, scenario = scenario, seed = seed)
    res$sim_id        <- i
    res$scenario      <- config_id
    res$scenario_name <- scenario$name
    res$family        <- scenario$family
    res$n_target      <- n
    results_list[[i]] <- res
  }

  bind_rows(results_list)
}

# NOTE: Summary statistics are computed by study1_summarize.R (standalone).
# This script only saves raw per-replicate results (all_results.rds).

# =============================================================================
# SECTION 5: MAIN EXECUTION
# =============================================================================
#
# When sourced by the Phase D cluster adapter (cluster/run_single_scenario.R),
# CLUSTER_MODE is set to TRUE before source(), and the block below is skipped.
# The adapter then picks one scenario by index and drives its own per-rep loop
# with incremental CSV writes. Local Rscript runs leave CLUSTER_MODE unset,
# so the guard passes and the block runs normally.

if (!exists("CLUSTER_MODE") || !isTRUE(CLUSTER_MODE)) {

config_ids <- names(scenarios)
if (!is.null(CONFIG_SUBSET)) {
  config_ids <- config_ids[config_ids %in% CONFIG_SUBSET]
}

cat("=============================================================\n")
cat("Study 1: Manuscript Simulation — Direct Generation (Null Effect)\n")
cat("=============================================================\n")
cat(sprintf("Mode: %s\n", MODE))
cat(sprintf("Replicates per scenario: %d\n", N_SIM))
cat(sprintf("Configurations to run: %d (%s)\n", length(config_ids),
            paste(config_ids, collapse = ", ")))
cat(sprintf("Seed base: %d\n", SEED_BASE))
cat(sprintf("Pool multiplier: %d\n", POOL_MULT))
cat(sprintf("Started at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n")

results_dir <- file.path(getwd(), "results", "study1")
figures_dir <- file.path(results_dir, "figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

all_results <- list()
t_start <- Sys.time()

for (idx in seq_along(config_ids)) {
  config_id <- config_ids[idx]
  scenario <- scenarios[[config_id]]

  scenario_seed_base <- SEED_BASE + (idx - 1) * 1000

  t1 <- Sys.time()
  results <- run_scenario_simulation(
    config_id  = config_id,
    scenario   = scenario,
    n_sim      = N_SIM,
    seed_base  = scenario_seed_base,
    verbose    = TRUE
  )
  t2 <- Sys.time()
  cat(sprintf("  Completed in %.1f seconds (%d/%d)\n",
              difftime(t2, t1, units = "secs"), idx, length(config_ids)))

  all_results[[config_id]] <- results
}

t_end <- Sys.time()
cat(sprintf("\nTotal runtime: %.1f minutes\n",
            difftime(t_end, t_start, units = "mins")))

combined_results <- bind_rows(all_results)

# =============================================================================
# SECTION 6: SAVE RAW RESULTS
# =============================================================================

saveRDS(combined_results, file.path(results_dir, "all_results.rds"))

# Provenance metadata (C2 fix)
provenance <- data.frame(
      script = "study1_simulation.R",
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version = paste0(R.version$major, ".", R.version$minor),
      n_sim = N_SIM,
      base_seed = SEED_BASE,
      n_rows = nrow(combined_results),
      stringsAsFactors = FALSE
)
write.csv(provenance, file.path(results_dir, "study1_provenance.csv"),
          row.names = FALSE)

cat("\n=============================================================\n")
cat(sprintf("Raw results saved to %s/\n", results_dir))
cat(sprintf("  all_results.rds: %d rows\n", nrow(combined_results)))
cat(sprintf("Completed at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("\nRun study1_summarize.R to compute summary statistics.\n")
cat("=============================================================\n")

}  # end if (!CLUSTER_MODE)
