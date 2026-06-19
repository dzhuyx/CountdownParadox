# =============================================================================
# Study 2: Manuscript Simulation — Trajectory-Based (True Effect)
#
# Demonstrates that when a true causal biomarker-disease relationship exists:
# 1. The naive analysis is still biased (always shows HR > 1)
# 2. The TVC analysis correctly identifies the true effect direction
# 3. The TVC analysis provides valid Type I error and power
#
# 59 unique scenarios: 36 pedagogical (PED) + 23 BIOCARD-calibrated (BIO)
# (42 originals + 17 Phase C null twins at off-base spec levels for Panel A)
#
# Methods compared:
#   1. Naive/countdown: Surv(T-Z, event) ~ Z_hat — always biased
#   2. TVC z_only: Surv(tstart, tstop, event) ~ Z_tv — valid
#   3. TVC interaction: Surv(tstart, tstop, event) ~ Z_tv + A:Z_tv — valid
#
# Source: Adapted from archived study2_simulation.R (landmark replaced with TVC)
# Specification: manuscript_simulation_plan.md Section 5
#
# Author: Yuxin Zhu
# Date: February 2026
# =============================================================================

library(lme4)
library(survival)
library(MASS)
library(dplyr)

# =============================================================================
# SECTION 1: MODE CONTROL
# =============================================================================

MODE <- "full"  # "prototype" or "full"
# prototype: n_sim=5, configs=c("S3b-ped","S3b-bio")  (~3 sec)
# full:      n_sim=1000, all 59 scenarios              (projected from pilot)

if (MODE == "prototype") {
  N_SIM <- 5
  CONFIG_SUBSET <- c("S3b-ped", "S3b-bio")
  cat("=== PROTOTYPE MODE: 5 reps, S3b-ped + S3b-bio only ===\n\n")
} else {
  N_SIM <- 1000
  CONFIG_SUBSET <- NULL
  cat("=== FULL MODE: 1000 reps, all 59 scenarios ===\n\n")
}

SEED_BASE <- 3026

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

create_params <- function(p) {
  p$Sigma <- matrix(
    c(p$sigma_0^2,
      p$rho * p$sigma_0 * p$sigma_1,
      p$rho * p$sigma_0 * p$sigma_1,
      p$sigma_1^2),
    nrow = 2
  )
  p$true_hr <- exp(p$beta_true)
  return(p)
}

# =============================================================================
# SECTION 3: SCENARIO DEFINITIONS (59 unique scenarios)
# =============================================================================

# --- Shared biomarker parameters ---

base_biomarker <- list(
  beta_0 = -1.5,
  beta_1 = 0.05,
  sigma_0 = 0.5,
  sigma_1 = 0.008,
  rho = -0.3,
  sigma_eps = 0.25,
  threshold = 1.5
)

# --- Pedagogical base ---
base_ped <- c(base_biomarker, list(
  entry_min = 50, entry_max = 65,
  visit_interval = 2, max_followup = 20, max_age = 85,
  lambda_0 = 0.03,
  n = 500, family = "PED"
))

# --- BIOCARD base ---
base_bio <- c(base_biomarker, list(
  entry_min = 40, entry_max = 70,
  visit_interval = 4, max_followup = 25, max_age = 95,
  lambda_0 = 0.02,
  n = 500, family = "BIO"
))

# --- Define all scenarios ---

scenarios <- list(

  # =========================================================================
  # Category 1: Effect Size Sweep
  # =========================================================================

  # PED effect size (7)
  `S3a-ped` = { p <- base_ped; p$name <- "S3a-ped: Null effect"; p$beta_true <- 0.0; create_params(p) },
  `S3b-ped` = { p <- base_ped; p$name <- "S3b-ped: Moderate harmful"; p$beta_true <- 0.5; create_params(p) },
  `S3c-ped` = { p <- base_ped; p$name <- "S3c-ped: Strong harmful"; p$beta_true <- 1.0; create_params(p) },
  `S3d-ped` = { p <- base_ped; p$name <- "S3d-ped: Moderate protective"; p$beta_true <- -0.5; create_params(p) },
  `S3e-ped` = { p <- base_ped; p$name <- "S3e-ped: Strong protective"; p$beta_true <- -1.0; create_params(p) },
  `S3f-ped` = { p <- base_ped; p$name <- "S3f-ped: Mild harmful"; p$beta_true <- 0.3; create_params(p) },
  `S3g-ped` = { p <- base_ped; p$name <- "S3g-ped: Mild protective"; p$beta_true <- -0.3; create_params(p) },

  # BIO effect size (7)
  `S3a-bio` = { p <- base_bio; p$name <- "S3a-bio: Null effect"; p$beta_true <- 0.0; create_params(p) },
  `S3b-bio` = { p <- base_bio; p$name <- "S3b-bio: Moderate harmful"; p$beta_true <- 0.5; create_params(p) },
  `S3c-bio` = { p <- base_bio; p$name <- "S3c-bio: Strong harmful"; p$beta_true <- 1.0; create_params(p) },
  `S3d-bio` = { p <- base_bio; p$name <- "S3d-bio: Moderate protective"; p$beta_true <- -0.5; create_params(p) },
  `S3e-bio` = { p <- base_bio; p$name <- "S3e-bio: Strong protective"; p$beta_true <- -1.0; create_params(p) },
  `S3f-bio` = { p <- base_bio; p$name <- "S3f-bio: Mild harmful"; p$beta_true <- 0.3; create_params(p) },
  `S3g-bio` = { p <- base_bio; p$name <- "S3g-bio: Mild protective"; p$beta_true <- -0.3; create_params(p) },

  # =========================================================================
  # Category 2: Sample Size Sweep (beta_true=0.5 and null twin at beta_true=0)
  # =========================================================================

  `S3h-ped`      = { p <- base_ped; p$name <- "S3h-ped: Small sample";                 p$beta_true <- 0.5; p$n <- 200;  create_params(p) },
  `S3h-ped-null` = { p <- base_ped; p$name <- "S3h-ped-null: Small sample (null)";     p$beta_true <- 0.0; p$n <- 200;  create_params(p) },
  `S3i-ped`      = { p <- base_ped; p$name <- "S3i-ped: Large sample";                 p$beta_true <- 0.5; p$n <- 1000; create_params(p) },
  `S3i-ped-null` = { p <- base_ped; p$name <- "S3i-ped-null: Large sample (null)";     p$beta_true <- 0.0; p$n <- 1000; create_params(p) },
  `S3j-ped`      = { p <- base_ped; p$name <- "S3j-ped: Very large sample";            p$beta_true <- 0.5; p$n <- 2000; create_params(p) },
  `S3j-ped-null` = { p <- base_ped; p$name <- "S3j-ped-null: Very large sample (null)"; p$beta_true <- 0.0; p$n <- 2000; create_params(p) },

  `S3h-bio`      = { p <- base_bio; p$name <- "S3h-bio: BIOCARD n=150";                p$beta_true <- 0.5; p$n <- 150;  create_params(p) },
  `S3h-bio-null` = { p <- base_bio; p$name <- "S3h-bio-null: BIOCARD n=150 (null)";    p$beta_true <- 0.0; p$n <- 150;  create_params(p) },
  `S3i-bio`      = { p <- base_bio; p$name <- "S3i-bio: BIOCARD large";                p$beta_true <- 0.5; p$n <- 1000; create_params(p) },
  `S3i-bio-null` = { p <- base_bio; p$name <- "S3i-bio-null: BIOCARD large (null)";    p$beta_true <- 0.0; p$n <- 1000; create_params(p) },

  # =========================================================================
  # Category 3: Trajectory Heterogeneity (beta_true=0.5 and null twin at beta_true=0, n=500)
  # =========================================================================

  `S3k-ped`      = { p <- base_ped; p$name <- "S3k-ped: Low heterogeneity";             p$beta_true <- 0.5; p$sigma_0 <- 0.25; p$sigma_1 <- 0.004; create_params(p) },
  `S3k-ped-null` = { p <- base_ped; p$name <- "S3k-ped-null: Low heterogeneity (null)"; p$beta_true <- 0.0; p$sigma_0 <- 0.25; p$sigma_1 <- 0.004; create_params(p) },
  `S3l-ped`      = { p <- base_ped; p$name <- "S3l-ped: High heterogeneity";            p$beta_true <- 0.5; p$sigma_0 <- 1.0;  p$sigma_1 <- 0.016; create_params(p) },
  `S3l-ped-null` = { p <- base_ped; p$name <- "S3l-ped-null: High heterogeneity (null)"; p$beta_true <- 0.0; p$sigma_0 <- 1.0; p$sigma_1 <- 0.016; create_params(p) },

  `S3k-bio`      = { p <- base_bio; p$name <- "S3k-bio: Low heterogeneity";             p$beta_true <- 0.5; p$sigma_0 <- 0.25; p$sigma_1 <- 0.004; create_params(p) },
  `S3k-bio-null` = { p <- base_bio; p$name <- "S3k-bio-null: Low heterogeneity (null)"; p$beta_true <- 0.0; p$sigma_0 <- 0.25; p$sigma_1 <- 0.004; create_params(p) },
  `S3l-bio`      = { p <- base_bio; p$name <- "S3l-bio: High heterogeneity";            p$beta_true <- 0.5; p$sigma_0 <- 1.0;  p$sigma_1 <- 0.016; create_params(p) },
  `S3l-bio-null` = { p <- base_bio; p$name <- "S3l-bio-null: High heterogeneity (null)"; p$beta_true <- 0.0; p$sigma_0 <- 1.0; p$sigma_1 <- 0.016; create_params(p) },

  # =========================================================================
  # Category 4: Disease Progression Rate (beta_true=0.5 and null twin at beta_true=0, n=500)
  # =========================================================================

  `S3m-ped`      = { p <- base_ped; p$name <- "S3m-ped: Slow disease";             p$beta_true <- 0.5; p$lambda_0 <- 0.02;  create_params(p) },
  `S3m-ped-null` = { p <- base_ped; p$name <- "S3m-ped-null: Slow disease (null)"; p$beta_true <- 0.0; p$lambda_0 <- 0.02;  create_params(p) },
  `S3n-ped`      = { p <- base_ped; p$name <- "S3n-ped: Fast disease";             p$beta_true <- 0.5; p$lambda_0 <- 0.05;  create_params(p) },
  `S3n-ped-null` = { p <- base_ped; p$name <- "S3n-ped-null: Fast disease (null)"; p$beta_true <- 0.0; p$lambda_0 <- 0.05;  create_params(p) },

  `S3m-bio`      = { p <- base_bio; p$name <- "S3m-bio: Slow disease";             p$beta_true <- 0.5; p$lambda_0 <- 0.015; create_params(p) },
  `S3m-bio-null` = { p <- base_bio; p$name <- "S3m-bio-null: Slow disease (null)"; p$beta_true <- 0.0; p$lambda_0 <- 0.015; create_params(p) },
  `S3n-bio`      = { p <- base_bio; p$name <- "S3n-bio: Fast disease";             p$beta_true <- 0.5; p$lambda_0 <- 0.03;  create_params(p) },
  `S3n-bio-null` = { p <- base_bio; p$name <- "S3n-bio-null: Fast disease (null)"; p$beta_true <- 0.0; p$lambda_0 <- 0.03;  create_params(p) },

  # =========================================================================
  # Category 5: Measurement Precision (beta_true=0.5 and null twin at beta_true=0, n=500)
  # =========================================================================

  `S3o-ped`      = { p <- base_ped; p$name <- "S3o-ped: Frequent visits";             p$beta_true <- 0.5; p$visit_interval <- 1; create_params(p) },
  `S3o-ped-null` = { p <- base_ped; p$name <- "S3o-ped-null: Frequent visits (null)"; p$beta_true <- 0.0; p$visit_interval <- 1; create_params(p) },
  `S3p-ped`      = { p <- base_ped; p$name <- "S3p-ped: Sparse visits";               p$beta_true <- 0.5; p$visit_interval <- 4; create_params(p) },
  `S3p-ped-null` = { p <- base_ped; p$name <- "S3p-ped-null: Sparse visits (null)";   p$beta_true <- 0.0; p$visit_interval <- 4; create_params(p) },
  `S3q-ped`      = { p <- base_ped; p$name <- "S3q-ped: High meas. error";            p$beta_true <- 0.5; p$sigma_eps <- 0.5;   create_params(p) },
  `S3q-ped-null` = { p <- base_ped; p$name <- "S3q-ped-null: High meas. error (null)"; p$beta_true <- 0.0; p$sigma_eps <- 0.5;  create_params(p) },

  `S3q-bio`      = { p <- base_bio; p$name <- "S3q-bio: High meas. error";            p$beta_true <- 0.5; p$sigma_eps <- 0.5;   create_params(p) },
  `S3q-bio-null` = { p <- base_bio; p$name <- "S3q-bio-null: High meas. error (null)"; p$beta_true <- 0.0; p$sigma_eps <- 0.5;  create_params(p) },

  # =========================================================================
  # Category 6: Cross-Product Scenarios
  # =========================================================================

  `S3r-ped` = { p <- base_ped; p$name <- "S3r-ped: Strong harmful + large n"; p$beta_true <- 1.0; p$n <- 1000; create_params(p) },
  `S3s-ped` = { p <- base_ped; p$name <- "S3s-ped: Strong harmful + very large n"; p$beta_true <- 1.0; p$n <- 2000; create_params(p) },
  `S3t-ped` = { p <- base_ped; p$name <- "S3t-ped: Strong protective + large n"; p$beta_true <- -1.0; p$n <- 1000; create_params(p) },
  `S3u-ped` = { p <- base_ped; p$name <- "S3u-ped: Strong harmful + high het"; p$beta_true <- 1.0; p$sigma_0 <- 1.0; p$sigma_1 <- 0.016; create_params(p) },
  `S3v-ped` = { p <- base_ped; p$name <- "S3v-ped: Strong protective + high het"; p$beta_true <- -1.0; p$sigma_0 <- 1.0; p$sigma_1 <- 0.016; create_params(p) },
  `S3w-ped` = { p <- base_ped; p$name <- "S3w-ped: Strong harmful + low het"; p$beta_true <- 1.0; p$sigma_0 <- 0.25; p$sigma_1 <- 0.004; create_params(p) },
  `S3x-ped` = { p <- base_ped; p$name <- "S3x-ped: Strong harmful + fast disease"; p$beta_true <- 1.0; p$lambda_0 <- 0.05; create_params(p) },
  `S3y-ped` = { p <- base_ped; p$name <- "S3y-ped: Strong harmful + high error"; p$beta_true <- 1.0; p$sigma_eps <- 0.5; create_params(p) },
  `S3z-ped` = { p <- base_ped; p$name <- "S3z-ped: Max power scenario"; p$beta_true <- 1.0; p$n <- 2000; p$lambda_0 <- 0.05; create_params(p) },

  `S3r-bio` = { p <- base_bio; p$name <- "S3r-bio: Strong harmful + large n"; p$beta_true <- 1.0; p$n <- 1000; create_params(p) },
  `S3s-bio` = { p <- base_bio; p$name <- "S3s-bio: Strong protective + large n"; p$beta_true <- -1.0; p$n <- 1000; create_params(p) }
)

cat(sprintf("Total scenarios defined: %d\n", length(scenarios)))

# =============================================================================
# SECTION 4: CORE FUNCTIONS
# =============================================================================

simulate_cohort_with_effect <- function(n, params) {
  # Generate random effects
  random_effects <- mvrnorm(n, mu = c(0, 0), Sigma = params$Sigma)
  colnames(random_effects) <- c("b0", "b1")

  # Study entry ages
  entry_age <- runif(n, params$entry_min, params$entry_max)

  # Individual trajectory parameters
  intercept_i <- params$beta_0 + random_effects[, "b0"]
  slope_i     <- params$beta_1 + random_effects[, "b1"]

  # Baseline biomarker (with measurement error)
  Y_baseline_true <- intercept_i + slope_i * entry_age
  Y_baseline      <- Y_baseline_true + rnorm(n, 0, params$sigma_eps)
  Y_baseline_mean <- mean(Y_baseline)
  Y_baseline_centered <- Y_baseline - Y_baseline_mean

  # Disease onset CONDITIONAL on Y_baseline
  rate_i       <- params$lambda_0 * exp(params$beta_true * Y_baseline_centered)
  waiting_time <- rexp(n, rate = rate_i)
  T_onset      <- entry_age + waiting_time

  # True Z
  true_Z <- (params$threshold - params$beta_0 - random_effects[, "b0"]) /
             (params$beta_1 + random_effects[, "b1"])
  true_Z[true_Z < 30]  <- 30
  true_Z[true_Z > 100] <- 100

  # Censoring
  censor_age    <- pmin(entry_age + params$max_followup, params$max_age)
  observed_time <- pmin(T_onset, censor_age)
  event         <- as.integer(T_onset <= censor_age)

  subject_data <- data.frame(
    subject_id = 1:n,
    entry_age = entry_age,
    T_onset = T_onset,
    true_Z = true_Z,
    censor_age = censor_age,
    observed_time = observed_time,
    event = event,
    b0 = random_effects[, "b0"],
    b1 = random_effects[, "b1"],
    Y_baseline = Y_baseline,
    Y_baseline_true = Y_baseline_true,
    Y_baseline_centered = Y_baseline_centered
  )

  # Generate longitudinal biomarker data
  biomarker_list <- vector("list", n)
  for (i in 1:n) {
    last_visit_age <- min(subject_data$observed_time[i], subject_data$censor_age[i])
    visit_ages <- seq(subject_data$entry_age[i], last_visit_age, by = params$visit_interval)
    if (length(visit_ages) < 2) {
      visit_ages <- c(subject_data$entry_age[i],
                      min(subject_data$observed_time[i], subject_data$censor_age[i]))
    }

    Y_true     <- intercept_i[i] + slope_i[i] * visit_ages
    Y_observed <- Y_true + rnorm(length(visit_ages), 0, params$sigma_eps)

    biomarker_list[[i]] <- data.frame(
      subject_id = i,
      visit = seq_along(visit_ages),
      age = visit_ages,
      Y = Y_observed,
      Y_true = Y_true
    )
  }

  biomarker_data <- do.call(rbind, biomarker_list)

  list(subject_data = subject_data, biomarker_data = biomarker_data,
       params = params, Y_baseline_mean = Y_baseline_mean)
}

estimate_Z_from_lmm <- function(cohort_data) {
  biomarker_data <- cohort_data$biomarker_data
  subject_data   <- cohort_data$subject_data
  params         <- cohort_data$params

  age_center <- 60
  biomarker_data$age_centered <- biomarker_data$age - age_center

  fit <- tryCatch(
    lmer(Y ~ age_centered + (1 + age_centered | subject_id),
         data = biomarker_data,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    result <- subject_data %>%
      mutate(Z_hat = NA, valid_Z = FALSE,
             beta_0_hat = NA, beta_1_hat = NA)
    return(list(
      data = result,
      model_info = list(
        beta_0_hat = NA, beta_1_hat = NA,
        sigma_0_hat = NA, sigma_1_hat = NA,
        rho_hat = NA, sigma_eps_hat = NA,
        n_valid = 0, n_invalid = nrow(subject_data)
      ),
      lmm_fit = NULL
    ))
  }

  beta_0_at_center <- fixef(fit)["(Intercept)"]
  beta_1_hat       <- fixef(fit)["age_centered"]
  beta_0_hat       <- beta_0_at_center - beta_1_hat * age_center

  re <- ranef(fit)$subject_id
  n <- nrow(subject_data)
  Z_hat   <- numeric(n)
  valid_Z <- logical(n)

  for (i in 1:n) {
    b0_hat <- re[as.character(i), "(Intercept)"]
    b1_hat <- re[as.character(i), "age_centered"]
    slope_i <- beta_1_hat + b1_hat

    if (is.na(slope_i) || slope_i <= 0) {
      Z_hat[i] <- NA
      valid_Z[i] <- FALSE
    } else {
      Z_hat[i] <- age_center + (params$threshold - beta_0_at_center - b0_hat) / slope_i
      valid_Z[i] <- TRUE
      Z_hat[i] <- pmin(pmax(Z_hat[i], 30), 100)
    }
  }

  result <- subject_data %>%
    mutate(Z_hat = Z_hat, valid_Z = valid_Z,
           beta_0_hat = beta_0_hat, beta_1_hat = beta_1_hat)

  model_info <- list(
    beta_0_hat  = beta_0_hat,
    beta_1_hat  = beta_1_hat,
    sigma_0_hat = as.numeric(attr(VarCorr(fit)$subject_id, "stddev")["(Intercept)"]),
    sigma_1_hat = as.numeric(attr(VarCorr(fit)$subject_id, "stddev")["age_centered"]),
    rho_hat     = as.numeric(attr(VarCorr(fit)$subject_id, "correlation")["(Intercept)", "age_centered"]),
    sigma_eps_hat = sigma(fit),
    n_valid   = sum(valid_Z),
    n_invalid = sum(!valid_Z)
  )

  list(data = result, model_info = model_info, lmm_fit = fit)
}

run_naive_analysis <- function(analysis_data) {
  dat <- analysis_data %>%
    filter(valid_Z, T_onset > Z_hat)

  if (nrow(dat) < 10) {
    return(list(hr = NA, log_hr = NA, se = NA, p_value = NA,
                ci_lower = NA, ci_upper = NA, n_analyzed = nrow(dat)))
  }

  dat$remaining_time     <- dat$T_onset - dat$Z_hat
  dat$observed_remaining <- pmin(dat$remaining_time, dat$censor_age - dat$Z_hat)
  dat$event_remaining    <- as.integer(dat$T_onset <= dat$censor_age)
  # Filter out subjects with non-positive remaining time (C3 fix — align with study1)
  dat <- dat[dat$observed_remaining > 0, ]

  tryCatch({
    fit <- coxph(Surv(observed_remaining, event_remaining) ~ Z_hat, data = dat)
    coef_summary <- summary(fit)$coefficients
    ci <- confint(fit)
    list(hr = exp(coef(fit)["Z_hat"]),
         log_hr = coef(fit)["Z_hat"],
         se = coef_summary["Z_hat", "se(coef)"],
         p_value = coef_summary["Z_hat", "Pr(>|z|)"],
         ci_lower = exp(ci["Z_hat", 1]),
         ci_upper = exp(ci["Z_hat", 2]),
         n_analyzed = nrow(dat))
  }, error = function(e) {
    list(hr = NA, log_hr = NA, se = NA, p_value = NA,
         ci_lower = NA, ci_upper = NA, n_analyzed = nrow(dat))
  })
}

run_tvc_analysis <- function(dat, formulation, A_col, entry_col, exit_col, event_col) {
  # Build start-stop data and fit TVC Cox model.
  # Adapted from tvc_investigation.R run_tvc_analysis() with CI computation added.
  #
  # formulation: "z_only" or "interaction"
  # Column names: A_col = age at positivity, entry_col = study entry age,
  #               exit_col = observed exit age, event_col = event indicator
  #
  # Subjects with NA in A_col are treated as never becoming positive.

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
    return(result)
  }

  if (all_positive && formulation == "interaction") {
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
      vc <- vcov(fit)
      if (nrow(vc) == 2) {
        result$vcov_cor_beta_gamma <- cov2cor(vc)[1, 2]
      }
    }
  }, error = function(e) NULL)

  result
}

run_single_simulation <- function(params, sim_id = 1) {
  # Simulate cohort with true effect
  cohort <- simulate_cohort_with_effect(params$n, params)

  # Estimate Z from LMM
  lmm_results  <- estimate_Z_from_lmm(cohort)
  analysis_data <- lmm_results$data
  model_info    <- lmm_results$model_info

  # Compute diagnostics
  valid_data <- analysis_data %>% filter(valid_Z)
  mean_Z     <- mean(valid_data$Z_hat, na.rm = TRUE)
  sd_Z       <- sd(valid_data$Z_hat, na.rm = TRUE)
  var_Z      <- var(valid_data$Z_hat, na.rm = TRUE)
  mean_T     <- mean(analysis_data$T_onset)
  sd_T       <- sd(analysis_data$T_onset)
  event_rate <- mean(analysis_data$event)
  mean_followup <- mean(analysis_data$observed_time - analysis_data$entry_age)
  n_Z_hat_NA    <- model_info$n_invalid
  pct_Z_hat_NA  <- model_info$n_invalid / nrow(analysis_data)
  mean_visits   <- nrow(cohort$biomarker_data) / params$n

  # Correlations
  cor_Z_TminusZ <- if (nrow(valid_data) > 2) cor(valid_data$Z_hat, valid_data$T_onset - valid_data$Z_hat) else NA
  cor_Z_T       <- if (nrow(valid_data) > 2) cor(valid_data$Z_hat, valid_data$T_onset) else NA
  cor_Z_Y       <- if (nrow(valid_data) > 2) cor(valid_data$Z_hat, valid_data$Y_baseline) else NA
  cor_Y_T       <- if (nrow(valid_data) > 2) cor(valid_data$Y_baseline, valid_data$T_onset) else NA

  # Shared diagnostics
  diag <- list(
    var_Z = var_Z, event_rate = event_rate,
    cor_Z_TminusZ = cor_Z_TminusZ, cor_Z_T = cor_Z_T,
    cor_Z_Y = cor_Z_Y, cor_Y_T = cor_Y_T,
    mean_Z = mean_Z, sd_Z = sd_Z,
    mean_T = mean_T, sd_T = sd_T,
    mean_followup = mean_followup,
    n_Z_hat_NA = n_Z_hat_NA, pct_Z_hat_NA = pct_Z_hat_NA,
    mean_visits = mean_visits,
    beta_0_hat = model_info$beta_0_hat,
    beta_1_hat = model_info$beta_1_hat,
    sigma_0_hat = model_info$sigma_0_hat,
    sigma_1_hat = model_info$sigma_1_hat,
    rho_hat = model_info$rho_hat,
    sigma_eps_hat = model_info$sigma_eps_hat,
    beta_true = params$beta_true,
    true_hr = params$true_hr
  )

  # TVC-specific NA placeholders for naive method
  tvc_na <- list(
    hr_gamma = NA, log_hr_gamma = NA, se_gamma = NA, p_gamma = NA,
    ci_lower_gamma = NA, ci_upper_gamma = NA,
    n_already_pos = NA, n_transition = NA, n_never_pos = NA,
    converged = NA, vcov_cor_beta_gamma = NA
  )

  # --- Naive analysis ---
  naive_result         <- run_naive_analysis(analysis_data)
  naive_result$method  <- "naive"
  naive_result <- c(naive_result, tvc_na, diag)

  # --- TVC z_only (uses ALL subjects, not just valid_Z) ---
  tvc_z         <- run_tvc_analysis(analysis_data, "z_only",
                                    "Z_hat", "entry_age", "observed_time", "event")
  tvc_z$method  <- "tvc_z_only"
  tvc_z <- c(tvc_z, diag)

  # --- TVC interaction ---
  tvc_int         <- run_tvc_analysis(analysis_data, "interaction",
                                      "Z_hat", "entry_age", "observed_time", "event")
  tvc_int$method  <- "tvc_interaction"
  tvc_int <- c(tvc_int, diag)

  all_results <- list(naive_result, tvc_z, tvc_int)
  bind_rows(lapply(all_results, as.data.frame))
}

run_scenario_simulation <- function(config_id, params, n_sim, seed_base,
                                    verbose = TRUE) {
  if (verbose) cat(sprintf("\nRunning %s (n=%d, beta=%.1f, %d reps)...\n",
                           params$name, params$n, params$beta_true, n_sim))

  results_list <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    if (verbose && i %% 50 == 0) cat(sprintf("  Simulation %d / %d\n", i, n_sim))
    set.seed(seed_base + i)

    tryCatch({
      res <- run_single_simulation(params, sim_id = i)
      res$sim_id        <- i
      res$scenario      <- config_id
      res$scenario_name <- params$name
      res$family        <- params$family
      res$n_target      <- params$n
      results_list[[i]] <- res
    }, error = function(e) {
      if (verbose) cat(sprintf("  Error in sim %d: %s\n", i, e$message))
      results_list[[i]] <- NULL
    })
  }

  results_list <- results_list[!sapply(results_list, is.null)]
  bind_rows(results_list)
}

# NOTE: Summary statistics are computed by study2_summarize.R (standalone).
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
cat("Study 2: Manuscript Simulation — Trajectory-Based (True Effect)\n")
cat("=============================================================\n")
cat(sprintf("Mode: %s\n", MODE))
cat(sprintf("Replicates per scenario: %d\n", N_SIM))
cat(sprintf("Configurations to run: %d (%s)\n", length(config_ids),
            paste(config_ids, collapse = ", ")))
cat(sprintf("Seed base: %d\n", SEED_BASE))
cat(sprintf("Started at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n")

results_dir <- file.path(getwd(), "results", "study2")
figures_dir <- file.path(results_dir, "figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

all_results <- list()
t_start <- Sys.time()

for (idx in seq_along(config_ids)) {
  config_id <- config_ids[idx]
  params <- scenarios[[config_id]]
  scenario_seed_base <- SEED_BASE + (idx - 1) * 1000

  t1 <- Sys.time()
  results <- run_scenario_simulation(
    config_id  = config_id,
    params     = params,
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
      script = "study2_simulation.R",
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version = paste0(R.version$major, ".", R.version$minor),
      n_sim = N_SIM,
      base_seed = SEED_BASE,
      n_rows = nrow(combined_results),
      stringsAsFactors = FALSE
)
write.csv(provenance, file.path(results_dir, "study2_provenance.csv"),
          row.names = FALSE)

cat("\n=============================================================\n")
cat(sprintf("Raw results saved to %s/\n", results_dir))
cat(sprintf("  all_results.rds: %d rows\n", nrow(combined_results)))
cat(sprintf("Completed at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("\nRun study2_summarize.R to compute summary statistics.\n")
cat("=============================================================\n")

}  # end if (!CLUSTER_MODE)
