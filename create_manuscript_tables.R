# ==============================================================================
# create_manuscript_tables.R
#
# Formatted CSV tables for the manuscript.
# Main tables + supplementary tables.
#
# Output: CSV files in <CP_PROJECT_ROOT>/results/manuscript_tables/
#   Table1.csv  — Demographics (4 cohort-subset columns)
#   Table2.csv  — TVC degeneracy assessment (5 biomarker-cohort rows)
#   Table3.csv  — Main results (standard + TV-BC + TV-AABC)
#   SuppTable_S1 — Study 1 complete results (15 scenarios × 4 methods)
#   SuppTable_S2 — Study 2 complete results (59 scenarios × 4 methods)
#   SuppTable_S3–S5 — Scenario specifications and degeneracy
#
# Columns match manuscript v1.6.1 supplementary tables:
#   Study 1: Type I Error, Mean HR, SD HR
#   Study 2: Rejection Rate, Mean HR, SD HR, % Correct Direction
# ==============================================================================

library(dplyr)
library(tidyr)

# -- Paths -------------------------------------------------------------------
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set it to the analysis I/O root.")
analysis_dir <- project_root
results_dir  <- file.path(project_root, "results")
# Simulation summary CSVs (study{1,2}/): under CP_SIM_DIR/results if set, else results/.
sim_dir      <- if (nzchar(Sys.getenv("CP_SIM_DIR"))) file.path(Sys.getenv("CP_SIM_DIR"), "results") else results_dir
out_dir      <- file.path(results_dir, "manuscript_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -- Helpers -----------------------------------------------------------------
fmt_pval <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("")
    if (x < 0.001) return("< 0.001")
    if (x < 0.01) return("< 0.01")
    if (x < 0.05) return("< 0.05")
    sprintf("= %.2f", x)
  })
}

fmt_hr <- function(x, d = 2) formatC(x, format = "f", digits = d)


# ==============================================================================
# TABLE 1: Demographics (4 cohort-subset columns)
# ==============================================================================

cat("\n===== TABLE 1: Demographics =====\n\n")

# Dynamic read from upstream CSVs (no hardcoded values)
demo <- read.csv(file.path(results_dir, "table1_demographics.csv"),
                 stringsAsFactors = FALSE)
degen <- read.csv(file.path(results_dir, "degeneracy_all_biomarkers.csv"),
                  stringsAsFactors = FALSE)

# Helper: format one cohort-subset column from demographics
fmt_demo <- function(d) {
  c(
    formatC(d$n, big.mark = ","),
    sprintf("%s (%s%%)", formatC(d$n_events, big.mark = ","), d$pct_events),
    sprintf("%d (%.1f%%)", round(d$n * d$pct_female / 100), d$pct_female),
    sprintf("%.1f (%.1f)", d$educ_mean, d$educ_sd),
    sprintf("%d (%.1f%%)", d$n_apoe4, d$pct_apoe4),
    sprintf("%.1f (%.1f)", d$entry_age_mean, d$entry_age_sd),
    sprintf("%d-%d", d$entry_age_min, d$entry_age_max),
    sprintf("%.1f (%.1f)", d$followup_mean, d$followup_sd),
    sprintf("%.1f", d$followup_median)
  )
}

char_labels_demo <- c(
  "N", "Events, n (%)", "Female, n (%)", "Education, mean (SD) yr",
  "APOE \u03b54 carrier, n (%)", "Entry age, mean (SD) yr", "Entry age, range",
  "Follow-up, mean (SD) yr", "Follow-up, median yr"
)

d_bc_csf  <- demo[demo$cohort_subset == "BIOCARD_CSF", ]
d_bc_pla  <- demo[demo$cohort_subset == "BIOCARD_Plasma", ]
d_adni_pet <- demo[demo$cohort_subset == "ADNI_PET", ]
d_adni_pla <- demo[demo$cohort_subset == "ADNI_Plasma", ]

tab1 <- data.frame(
  Characteristic = char_labels_demo,
  BIOCARD_CSF    = fmt_demo(d_bc_csf),
  BIOCARD_Plasma = fmt_demo(d_bc_pla),
  ADNI_PET       = fmt_demo(d_adni_pet),
  ADNI_Plasma    = fmt_demo(d_adni_pla),
  stringsAsFactors = FALSE
)

write.csv(tab1, file.path(out_dir, "Table1.csv"), row.names = FALSE)
cat("Table 1 saved.\n")
print(tab1, right = FALSE, row.names = FALSE)


# ==============================================================================
# TABLE 2: TVC Degeneracy Assessment (5 biomarker-cohort rows)
# ==============================================================================

cat("\n===== TABLE 2: TVC Degeneracy Assessment =====\n\n")

eaoa_summ <- read.csv(file.path(results_dir, "eaoa_summary.csv"),
                       stringsAsFactors = FALSE)

tab2 <- degen %>%
  left_join(eaoa_summ %>% select(biomarker, cohort, mean_eaoa, sd_eaoa),
            by = c("biomarker", "cohort")) %>%
  mutate(
    `Biomarker (Cohort)` = case_when(
      biomarker == "CSF_AB42_AB40"      ~ "CSF A\u03b242/A\u03b240 (BIOCARD)",
      biomarker == "CSF_pTau181"        ~ "CSF p-tau181 (BIOCARD)",
      biomarker == "Plasma_pTau181"     ~ "Plasma p-tau181 (BIOCARD)",
      biomarker == "Amyloid_PET_FBP"    ~ "Amyloid PET (ADNI)",
      biomarker == "Plasma_pTau217_Fuji"~ "Plasma p-tau217 (ADNI)"),
    `N with valid EAOA`             = as.character(n_valid_eaoa),
    `Already positive at entry, n (%)` = sprintf("%d (%.1f%%)", n_already_pos, pct_already_pos),
    `Transition during FU, n`       = as.character(n_transition),
    `Never positive, n`             = as.character(n_never_pos),
    `AABC, mean (SD) yr`            = sprintf("%.1f (%.1f)", mean_eaoa, sd_eaoa)
  ) %>%
  select(`Biomarker (Cohort)`, `N with valid EAOA`,
         `Already positive at entry, n (%)`,
         `Transition during FU, n`, `Never positive, n`,
         `AABC, mean (SD) yr`)

# Order: BIOCARD CSF AB, BIOCARD CSF ptau, BIOCARD Plasma, ADNI PET, ADNI Plasma
bio_order_t2 <- c("CSF A\u03b242/A\u03b240 (BIOCARD)", "CSF p-tau181 (BIOCARD)",
                  "Plasma p-tau181 (BIOCARD)", "Amyloid PET (ADNI)",
                  "Plasma p-tau217 (ADNI)")
tab2 <- tab2[match(bio_order_t2, tab2$`Biomarker (Cohort)`), ]
rownames(tab2) <- NULL

write.csv(tab2, file.path(out_dir, "Table2.csv"), row.names = FALSE)
cat("Table 2 saved.\n")
print(tab2, right = FALSE, row.names = FALSE)


# ==============================================================================
# TABLE 3: Main Results (merged Table 2 + Table 3)
#
# Layout: Biomarker (Cohort) | Positivity status | Standard Countdown cols |
#         TV-BC cols | TV-AABC cols
# - Combined "HR (95% CI); P" cells
# - Positivity status from degeneracy CSV
# - TV-BC n (events) shared with TV-AABC (footnote)
# ==============================================================================

cat("\n===== TABLE 3: Main Results =====\n\n")

main_res <- read.csv(file.path(results_dir, "main_results_all_biomarkers.csv"),
                     stringsAsFactors = FALSE)

# EN-DASH for CI ranges
EN <- "\u2013"

# Combined HR (95% CI); P formatter
fmt_hr_ci_p <- function(hr, lo, hi, p) {
  ci <- paste0(fmt_hr(hr), " (", fmt_hr(lo), EN, fmt_hr(hi), ")")
  paste0(ci, "; P ", fmt_pval(p))
}

# Biomarker display labels: biomarker first, cohort in parentheses
bio_labels <- c(
  "CSF_AB42_AB40.BIOCARD"     = "CSF A\u03b242/A\u03b240 (BIOCARD)",
  "CSF_pTau181.BIOCARD"       = "CSF p-tau181 (BIOCARD)",
  "Plasma_pTau181.BIOCARD"    = "Plasma p-tau181 (BIOCARD)",
  "Amyloid_PET_FBP.ADNI"      = "Amyloid PET (ADNI)",
  "Plasma_pTau217_Fuji.ADNI"  = "Plasma p-tau217 (ADNI)"
)

# Row order
bio_order <- names(bio_labels)

# Build one row per biomarker-cohort
tab3_rows <- lapply(bio_order, function(bk) {
  parts <- strsplit(bk, "\\.")[[1]]
  bm <- parts[1]; co <- parts[2]

  # Positivity status from degeneracy
  dg <- degen[degen$biomarker == bm & degen$cohort == co, ]
  pos_status <- sprintf("%d positive (%d at entry, %d transition)",
                        dg$n_valid_eaoa, dg$n_already_pos, dg$n_transition)

  # Extract model results
  get_row <- function(mod, par = NULL) {
    r <- main_res[main_res$biomarker == bm & main_res$cohort == co &
                  main_res$model == mod, ]
    if (!is.null(par)) r <- r[r$parameter == par, ]
    r
  }

  p1 <- get_row("P1_countdown")
  p2 <- get_row("P2_tvc_z_only")
  p3b <- get_row("P3_tvc_interaction", "beta_Z_tv")
  p3g <- get_row("P3_tvc_interaction", "gamma_A")

  data.frame(
    `Biomarker (Cohort)` = bio_labels[[bk]],
    `Positivity status`  = pos_status,
    `SC: n (events)`     = sprintf("%d (%d)", p1$n, p1$nevent),
    `SC: HR (95% CI); P` = fmt_hr_ci_p(p1$hr, p1$lower95, p1$upper95, p1$pvalue),
    `_1` = "",
    `TV-BC: n (events)*` = sprintf("%d (%d)", p2$n, p2$nevent),
    `TV-BC: HR (95% CI); P` = fmt_hr_ci_p(p2$hr, p2$lower95, p2$upper95, p2$pvalue),
    `_2` = "",
    `TV-AABC: HR-beta (95% CI); P`  = fmt_hr_ci_p(p3b$hr, p3b$lower95, p3b$upper95, p3b$pvalue),
    `TV-AABC: HR-gamma (95% CI); P` = fmt_hr_ci_p(p3g$hr, p3g$lower95, p3g$upper95, p3g$pvalue),
    stringsAsFactors = FALSE, check.names = FALSE
  )
})

tab3 <- do.call(rbind, tab3_rows)
rownames(tab3) <- NULL

write.csv(tab3, file.path(out_dir, "Table3.csv"), row.names = FALSE)
cat("Table 3 saved.\n")
cat("* n (events) under TV-BC same for TV-AABC\n")
cat("Spacer columns _1, _2 mark method group boundaries\n\n")
print(tab3[, !names(tab3) %in% c("_1", "_2")], right = FALSE, row.names = FALSE)


# ==============================================================================
# SUPPLEMENTARY TABLE S1: Study 1 Complete Simulation Results
# 15 scenarios x 4 methods
# Columns: AABC Variance, Event Rate, Type I Error, Mean HR, SD HR
# ==============================================================================

cat("\n===== Supp Table S1: Study 1 Complete Results =====\n\n")

s1 <- read.csv(file.path(sim_dir, "study1", "summary_results.csv"),
               stringsAsFactors = FALSE)

# Study 1 reports 15 scenarios (S1-S9, with S1 and S9 sample-size variants).
# study1_simulation.R emits these manuscript IDs directly (S8 = High Var(Z),
# S9 = BIOCARD-calibrated); no relabeling is needed here.

# Build 4 method series by splitting tvc_interaction into beta and gamma rows
# SD HR computed via delta method: SD(HR) ≈ mean_hr × empirical_se
s1_tab <- bind_rows(
  s1 %>% filter(method == "naive") %>%
    transmute(Scenario = scenario, Description = scenario_name, n = n_target,
              Method = "Standard countdown",
              AABC_Variance = mean_var_Z,
              Event_Rate_Pct = mean_event_rate * 100,
              Type_I_Error = type1_error, Mean_HR = mean_hr,
              SD_HR = mean_hr * empirical_se),
  s1 %>% filter(method == "tvc_z_only") %>%
    transmute(Scenario = scenario, Description = scenario_name, n = n_target,
              Method = "TV-BC",
              AABC_Variance = mean_var_Z,
              Event_Rate_Pct = mean_event_rate * 100,
              Type_I_Error = type1_error, Mean_HR = mean_hr,
              SD_HR = mean_hr * empirical_se),
  s1 %>% filter(method == "tvc_interaction") %>%
    transmute(Scenario = scenario, Description = scenario_name, n = n_target,
              Method = "TV-AABC \u03b2",
              AABC_Variance = mean_var_Z,
              Event_Rate_Pct = mean_event_rate * 100,
              Type_I_Error = type1_error, Mean_HR = mean_hr,
              SD_HR = mean_hr * empirical_se),
  s1 %>% filter(method == "tvc_interaction", !is.na(type1_error_gamma)) %>%
    transmute(Scenario = scenario, Description = scenario_name, n = n_target,
              Method = "TV-AABC \u03b3",
              AABC_Variance = mean_var_Z,
              Event_Rate_Pct = mean_event_rate * 100,
              Type_I_Error = type1_error_gamma, Mean_HR = mean_hr_gamma,
              SD_HR = mean_hr_gamma * mean_se_gamma)
) %>%
  arrange(Scenario, factor(Method, levels = c("Standard countdown", "TV-BC",
                                               "TV-AABC \u03b2", "TV-AABC \u03b3")))

s1_out <- s1_tab %>%
  mutate(Scenario = gsub("_n", ", n=", Scenario),
         Description = gsub("_n", ", n=", Description),
         AABC_Variance = sprintf("%.1f", AABC_Variance),
         Event_Rate_Pct = sprintf("%.1f%%", Event_Rate_Pct),
         Type_I_Error = sprintf("%.3f", Type_I_Error),
         Mean_HR = fmt_hr(Mean_HR),
         SD_HR = sprintf("%.3f", SD_HR))

write.csv(s1_out, file.path(out_dir, "SuppTable_S1_study1_complete.csv"),
          row.names = FALSE)
cat(sprintf("Supp Table S1 saved: %d rows (15 scenarios x 4 methods)\n", nrow(s1_out)))


# ==============================================================================
# SUPPLEMENTARY TABLE S2: Study 2 Complete Simulation Results
# 59 scenarios x 4 method series
# Columns: Reject Rate, Mean HR, SD HR, % Correct Direction (non-null only)
# ==============================================================================

cat("\n===== Supp Table S2: Study 2 Complete Results =====\n\n")

s2 <- read.csv(file.path(sim_dir, "study2", "summary_results.csv"),
               stringsAsFactors = FALSE)

# Pre-compute gamma direction metrics from per-replicate data
# (summary CSV has pct_hr_above_1/below_1 for beta only, not gamma)
cat("  Loading per-replicate data for gamma direction...\n")
s2_rds <- readRDS(file.path(sim_dir, "study2", "all_results.rds"))
gamma_dir <- s2_rds %>%
  filter(method == "tvc_interaction", !is.na(hr_gamma)) %>%
  group_by(scenario) %>%
  summarise(
    pct_hr_above_1_gamma = mean(hr_gamma > 1, na.rm = TRUE),
    pct_hr_below_1_gamma = mean(hr_gamma < 1, na.rm = TRUE),
    .groups = "drop")
rm(s2_rds)  # free memory

# Direction accuracy logic:
#   Standard: Z -> T-Z estimand. Structural bias always pushes HR > 1.
#     Correct direction when beta > 0 is HR < 1 (older Z = less exposure = lower risk).
#     Correct direction when beta < 0 is HR > 1 (older Z = less protection = higher risk).
#   TV-BC / TV-AABC beta: Z(t) -> hazard estimand.
#     Correct when beta > 0 is HR > 1; when beta < 0 is HR < 1.
#   TV-AABC gamma: A:Z(t) interaction. Direction is INVERTED relative to beta.
#     When beta > 0: earlier Z -> more exposure -> higher hazard -> gamma < 0 (HR < 1).
#     When beta < 0: earlier Z -> more protection -> lower hazard -> gamma > 0 (HR > 1).
#
# SD HR computed via delta method: SD(HR) ≈ mean_hr × empirical_se(log_hr)
# For gamma: SD(HR_γ) ≈ mean_hr_gamma × mean_se_gamma (model SE as proxy)

s2_tab <- bind_rows(
  # Standard countdown
  s2 %>% filter(method == "naive") %>%
    transmute(Scenario = scenario, Family = family, n = n_target,
              Beta_True = beta_true, True_HR = true_hr,
              Method = "Standard countdown",
              Reject_Rate = reject_rate, Mean_HR = mean_hr,
              SD_HR = mean_hr * empirical_se,
              Pct_Correct_Direction = case_when(
                beta_true == 0 ~ NA_real_,
                beta_true > 0  ~ pct_hr_below_1,
                beta_true < 0  ~ pct_hr_above_1)),
  # TV-BC
  s2 %>% filter(method == "tvc_z_only") %>%
    transmute(Scenario = scenario, Family = family, n = n_target,
              Beta_True = beta_true, True_HR = true_hr,
              Method = "TV-BC",
              Reject_Rate = reject_rate, Mean_HR = mean_hr,
              SD_HR = mean_hr * empirical_se,
              Pct_Correct_Direction = case_when(
                beta_true == 0 ~ NA_real_,
                beta_true > 0  ~ pct_hr_above_1,
                beta_true < 0  ~ pct_hr_below_1)),
  # TV-AABC beta
  s2 %>% filter(method == "tvc_interaction") %>%
    transmute(Scenario = scenario, Family = family, n = n_target,
              Beta_True = beta_true, True_HR = true_hr,
              Method = "TV-AABC \u03b2",
              Reject_Rate = reject_rate, Mean_HR = mean_hr,
              SD_HR = mean_hr * empirical_se,
              Pct_Correct_Direction = case_when(
                beta_true == 0 ~ NA_real_,
                beta_true > 0  ~ pct_hr_above_1,
                beta_true < 0  ~ pct_hr_below_1)),
  # TV-AABC gamma
  s2 %>% filter(method == "tvc_interaction") %>%
    left_join(gamma_dir, by = "scenario") %>%
    transmute(Scenario = scenario, Family = family, n = n_target,
              Beta_True = beta_true, True_HR = true_hr,
              Method = "TV-AABC \u03b3",
              Reject_Rate = reject_rate_gamma, Mean_HR = mean_hr_gamma,
              SD_HR = mean_hr_gamma * mean_se_gamma,
              Pct_Correct_Direction = case_when(
                beta_true == 0 ~ NA_real_,
                beta_true > 0  ~ pct_hr_below_1_gamma,
                beta_true < 0  ~ pct_hr_above_1_gamma))
) %>%
  arrange(Scenario, factor(Method, levels = c("Standard countdown", "TV-BC",
                                               "TV-AABC \u03b2", "TV-AABC \u03b3")))

s2_out <- s2_tab %>%
  mutate(Beta_True = sprintf("%.1f", Beta_True),
         True_HR = fmt_hr(True_HR),
         Reject_Rate = sprintf("%.3f", Reject_Rate),
         Mean_HR = fmt_hr(Mean_HR),
         SD_HR = sprintf("%.3f", SD_HR),
         Pct_Correct_Direction = ifelse(is.na(Pct_Correct_Direction), "",
                                        sprintf("%.1f%%", Pct_Correct_Direction * 100)))

write.csv(s2_out, file.path(out_dir, "SuppTable_S2_study2_complete.csv"),
          row.names = FALSE)
cat(sprintf("Supp Table S2 saved: %d rows (59 scenarios x 4 methods)\n", nrow(s2_out)))


# ==============================================================================
# SUPPLEMENTARY TABLE S3: Study 1 Scenario Specifications
# One row per scenario — DGP parameters, not results
# ==============================================================================

cat("\n===== Supp Table S3: Study 1 Scenario Specifications =====\n\n")

s1_specs <- data.frame(
  Scenario = c("S1", "S2", "S3", "S4", "S5",
               "S6", "S7", "S8", "S9",
               "S1, n=150", "S1, n=200", "S1, n=1000",
               "S9, n=150", "S9, n=200", "S9, n=1000"),
  n = c(rep(500, 9), 150, 200, 1000, 150, 200, 1000),
  Z_Distribution = c(
    "Uniform(50, 75)", "truncN(62, 5, 45, 80)", "truncN(62, 5, 45, 80)",
    "Gamma(20, 20/22) + 40", "Beta(2,5) scaled [45,80]",
    "truncN(62, 5, 45, 80)", "truncN(62, 5, 45, 80)",
    "truncN(62, 10, 35, 90)", "truncN(53, 10, 25, 82)",
    "Uniform(50, 75)", "Uniform(50, 75)", "Uniform(50, 75)",
    "truncN(53, 10, 25, 82)", "truncN(53, 10, 25, 82)", "truncN(53, 10, 25, 82)"),
  T_Distribution = c(
    "Uniform(60, 90)", "truncN(75, 7, 55, 95)", "Weibull(5, 78), floor 50",
    "Weibull(3, 20) + 55", "Exp(0.05) + 60, cap 95",
    "truncN(75, 7, 55, 95)", "truncN(75, 7, 55, 95)",
    "truncN(75, 7, 55, 95)", "truncN(73, 10, 45, 96)",
    "Uniform(60, 90)", "Uniform(60, 90)", "Uniform(60, 90)",
    "truncN(73, 10, 45, 96)", "truncN(73, 10, 45, 96)", "truncN(73, 10, 45, 96)"),
  Study_Entry = c(
    "U(55,70)", "U(55,70)", "U(55,70)", "U(55,70)", "U(55,70)",
    "U(50,60)", "U(60,75)",
    "U(55,70)", "U(40,65)",
    "U(55,70)", "U(55,70)", "U(55,70)",
    "U(40,65)", "U(40,65)", "U(40,65)"),
  Max_Followup = c(
    rep(20, 5), 35, 12, 20, 20,
    20, 20, 20, 20, 20, 20),
  Design_Note = c(
    "Base: Uniform-Uniform", "Base: Normal-Normal", "Normal Z, Weibull T",
    "Gamma Z, Weibull T", "Beta Z, Exponential T",
    "Light censoring (early entry, long FU)", "Heavy censoring (late entry, short FU)",
    "High Var(Z): SD=10 vs 5", "BIOCARD-calibrated parameters",
    "S1 at n=150", "S1 at n=200", "S1 at n=1000",
    "S9 at n=150", "S9 at n=200", "S9 at n=1000"),
  stringsAsFactors = FALSE
)

write.csv(s1_specs, file.path(out_dir, "SuppTable_S3_study1_specifications.csv"),
          row.names = FALSE)
cat(sprintf("Supp Table S3 saved: %d scenarios\n", nrow(s1_specs)))


# ==============================================================================
# SUPPLEMENTARY TABLE S4: Study 2 Scenario Specifications
# One row per scenario — DGP parameters, not results
# ==============================================================================

cat("\n===== Supp Table S4: Study 2 Scenario Specifications =====\n\n")

# Read scenario metadata from summary to extract n, beta_true, family
s2_meta <- s2 %>%
  filter(method == "naive") %>%
  select(Scenario = scenario, Family = family, n = n_target,
         Beta_True = beta_true, True_HR = true_hr)

# Shared biomarker defaults: beta_0=-1.5, beta_1=0.05, threshold=1.5
# Sigma and observation defaults differ by scenario

# Build specification columns from scenario IDs
s2_specs <- s2_meta %>%
  mutate(
    # sigma_0 (trajectory intercept SD)
    sigma_0 = case_when(
      grepl("S3k", Scenario)  ~ 0.25,   # low heterogeneity
      grepl("S3l", Scenario)  ~ 1.0,    # high heterogeneity
      grepl("S3u|S3v", Scenario) ~ 1.0, # cross-product high het
      grepl("S3w", Scenario)  ~ 0.25,   # cross-product low het
      TRUE ~ 0.5),                       # default
    # sigma_1 (trajectory slope SD)
    sigma_1 = case_when(
      grepl("S3k", Scenario)  ~ 0.004,
      grepl("S3l", Scenario)  ~ 0.016,
      grepl("S3u|S3v", Scenario) ~ 0.016,
      grepl("S3w", Scenario)  ~ 0.004,
      TRUE ~ 0.008),
    # lambda_0 (baseline hazard)
    lambda_0 = case_when(
      grepl("S3m", Scenario) & Family == "PED" ~ 0.02,
      grepl("S3m", Scenario) & Family == "BIO" ~ 0.015,
      grepl("S3n", Scenario) & Family == "PED" ~ 0.05,
      grepl("S3n", Scenario) & Family == "BIO" ~ 0.03,
      grepl("S3x|S3z", Scenario) ~ 0.05,  # cross-product high lambda
      Family == "PED" ~ 0.03,
      Family == "BIO" ~ 0.02),
    # visit_interval
    visit_interval = case_when(
      grepl("S3o", Scenario) ~ 1L,    # frequent visits
      grepl("S3p", Scenario) ~ 4L,    # sparse visits
      Family == "PED" ~ 2L,
      Family == "BIO" ~ 4L),
    # sigma_eps (measurement error)
    sigma_eps = case_when(
      grepl("S3q|S3y", Scenario) ~ 0.5,  # high measurement error
      TRUE ~ 0.25),
    # Category label
    Category = case_when(
      grepl("^S3[a-g]-", Scenario) ~ "Effect sweep",
      grepl("^S3[hi]-|^S3j-", Scenario) ~ "Sample size sweep",
      grepl("^S3[kl]-", Scenario) ~ "Heterogeneity sweep",
      grepl("^S3[mn]-", Scenario) ~ "Disease rate sweep",
      grepl("^S3[opq]-", Scenario) ~ "Measurement sweep",
      grepl("^S3[r-z]-", Scenario) ~ "Cross-product"),
    Is_Null_Twin = grepl("-null$", Scenario)
  ) %>%
  mutate(
    Beta_True = sprintf("%.1f", Beta_True),
    True_HR = fmt_hr(True_HR),
    sigma_0 = sprintf("%.2f", sigma_0),
    sigma_1 = sprintf("%.3f", sigma_1),
    lambda_0 = sprintf("%.3f", lambda_0),
    sigma_eps = sprintf("%.2f", sigma_eps)
  ) %>%
  select(Scenario, Family, Category, Is_Null_Twin, n,
         Beta_True, True_HR,
         sigma_0, sigma_1, lambda_0,
         visit_interval, sigma_eps) %>%
  arrange(Category, Scenario)

write.csv(s2_specs, file.path(out_dir, "SuppTable_S4_study2_specifications.csv"),
          row.names = FALSE)
cat(sprintf("Supp Table S4 saved: %d scenarios\n", nrow(s2_specs)))


# ==============================================================================
# SUPPLEMENTARY TABLE S5: Degeneracy Assessment (Real Data)
# ==============================================================================

cat("\n===== Supp Table S5: Degeneracy =====\n\n")

degen <- read.csv(file.path(results_dir, "degeneracy_all_biomarkers.csv"),
                   stringsAsFactors = FALSE)

degen_tab <- degen %>%
  mutate(
    Biomarker_Cohort = case_when(
      biomarker == "CSF_AB42_AB40" ~ "BIOCARD CSF AB42/40",
      biomarker == "CSF_pTau181" ~ "BIOCARD CSF p-tau181",
      biomarker == "Plasma_pTau181" ~ "BIOCARD Plasma p-tau181",
      biomarker == "Amyloid_PET_FBP" ~ "ADNI Amyloid PET",
      biomarker == "Plasma_pTau217_Fuji" ~ "ADNI Plasma p-tau217"),
    Pct_Already_Positive = sprintf("%.1f%%", pct_already_pos),
    Interpretation = case_when(
      pct_already_pos >= 50 ~ "Highly degenerate",
      pct_already_pos >= 20 ~ "Partially degenerate",
      TRUE ~ "Operational")
  ) %>%
  select(
    Biomarker_Cohort, N_Total = n_total, N_Events = n_events,
    N_Valid_EAOA = n_valid_eaoa, N_Already_Positive = n_already_pos,
    Pct_Already_Positive, N_Transition = n_transition,
    N_Never_Positive = n_never_pos, Interpretation
  )

write.csv(degen_tab, file.path(out_dir, "SuppTable_S5_degeneracy.csv"),
          row.names = FALSE)
cat("Supp Table S5 saved.\n")
print(degen_tab, right = FALSE, row.names = FALSE)


# ==============================================================================
# SUPPLEMENTARY TABLE S8: Full Model Coefficients
# All covariates from all 15 Cox PH models (5 biomarkers x 3 models)
# ==============================================================================

cat("\n===== Supp Table S8: Full Model Coefficients =====\n\n")

full_coefs_file <- file.path(results_dir, "full_coefficients_all_models.csv")
if (file.exists(full_coefs_file)) {
  full_coefs <- read.csv(full_coefs_file, stringsAsFactors = FALSE)

  # Biomarker display labels
  bm_labels <- c(
    "CSF_AB42_AB40"      = "CSF A\u03b242/A\u03b240",
    "CSF_pTau181"        = "CSF p-tau181",
    "Plasma_pTau181"     = "Plasma p-tau181",
    "Amyloid_PET_FBP"    = "Amyloid PET",
    "Plasma_pTau217_Fuji"= "Plasma p-tau217"
  )

  # Model display labels
  mod_labels <- c(
    "P1_countdown"       = "Standard countdown",
    "P2_tvc_z_only"      = "TV-BC",
    "P3_tvc_interaction" = "TV-AABC"
  )

  EN <- "\u2013"

  s8_tab <- full_coefs %>%
    mutate(
      Biomarker = bm_labels[biomarker],
      Cohort = cohort,
      Model = mod_labels[model],
      Parameter = parameter,
      `HR (95% CI)` = sprintf("%.2f (%.2f%s%.2f)", hr, lower95, EN, upper95),
      P = fmt_pval(pvalue),
      N = n,
      Events = nevent
    ) %>%
    select(Cohort, Biomarker, Model, Parameter, `HR (95% CI)`, P, N, Events)

  write.csv(s8_tab, file.path(out_dir, "SuppTable_S8_full_coefficients.csv"),
            row.names = FALSE)
  cat(sprintf("Supp Table S8 saved: %d rows across %d models\n",
              nrow(s8_tab),
              length(unique(paste(s8_tab$Cohort, s8_tab$Biomarker, s8_tab$Model)))))
} else {
  cat("  WARNING: full_coefficients_all_models.csv not found.\n")
  cat("  Run extract_full_coefficients.R first (Phase 2).\n")
}


# ==============================================================================
# MANUSCRIPT-EXACT SUPPLEMENT TABLES (docx Tables 1-8)
# Re-formats the validated results into the 8 tables exactly as they appear in
# the manuscript supplement, so the pipeline *produces* the manuscript tables:
#   * pretty headers; merged "Scenario: description" column
#   * Study 2 split into 6 grouped tables (Tables 2-7)
#   * LOGICAL method order throughout: Standard countdown, TV-BC, TV-AABC b, TV-AABC g
#   * Table 8: 15 model blocks (header rows), humanized parameters, merged HR+P
# Verified cell-identical to the manuscript supplement (v7).
# Output: SuppTable_1.csv ... SuppTable_8.csv
# ==============================================================================

cat("\n===== Manuscript-exact supplement tables (SuppTable_1..8) =====\n\n")

mt_B <- "β"; mt_G <- "γ"; mt_EN <- "–"; mt_EM <- "—"
mt_TIMES <- "×"; mt_EPS <- "ε"; mt_DASH <- "—"
mt_METH <- c("Standard countdown", "TV-BC", paste0("TV-AABC ", mt_B), paste0("TV-AABC ", mt_G))
mt_ord <- function(df) df[order(match(df$Method, mt_METH)), ]
mt_nm2 <- unique(read.csv(file.path(sim_dir, "study2", "summary_results.csv"),
                          stringsAsFactors = FALSE)[, c("scenario", "scenario_name")])

# ---- Table 1 (Study 1) ----
mt_T1_ORDER <- c("S1","S9","S9, n=1000","S9, n=150","S9, n=200",
                 "S1, n=1000","S1, n=150","S1, n=200","S2","S3","S4","S5","S6","S7","S8")
mt_t1 <- s1_out
mt_t1$Description <- gsub("High Var\\(Z\\)", "High AABC variance", mt_t1$Description)
mt_T1 <- do.call(rbind, lapply(mt_T1_ORDER, function(sc) mt_ord(mt_t1[mt_t1$Scenario == sc, ])))
mt_T1 <- mt_T1[, c("Description","n","Method","AABC_Variance",
                   "Event_Rate_Pct","Type_I_Error","Mean_HR","SD_HR")]
names(mt_T1) <- c("Scenario","N","Analytic approach","AABC variance",
                  "Event rate (%)","Type I error","Mean HR","SD HR")
write.csv(mt_T1, file.path(out_dir, "SuppTable_1.csv"), row.names = FALSE)

# ---- Tables 2-7 (Study 2) ----
mt_SPLIT <- list(
  `2` = c("S3e-bio","S3d-bio","S3g-bio","S3a-bio","S3f-bio","S3b-bio","S3c-bio",
          "S3e-ped","S3d-ped","S3g-ped","S3a-ped","S3f-ped","S3b-ped","S3c-ped"),
  `3` = c("S3h-bio","S3h-bio-null","S3h-ped","S3h-ped-null","S3i-bio","S3i-bio-null",
          "S3i-ped","S3i-ped-null","S3j-ped","S3j-ped-null"),
  `4` = c("S3k-bio","S3k-bio-null","S3k-ped","S3k-ped-null","S3l-bio","S3l-bio-null",
          "S3l-ped","S3l-ped-null"),
  `5` = c("S3m-bio","S3m-bio-null","S3m-ped","S3m-ped-null","S3n-bio","S3n-bio-null",
          "S3n-ped","S3n-ped-null"),
  `6` = c("S3o-ped","S3o-ped-null","S3p-ped","S3p-ped-null","S3q-bio","S3q-bio-null",
          "S3q-ped","S3q-ped-null"),
  `7` = c("S3r-bio","S3r-ped","S3s-bio","S3s-ped","S3t-ped","S3u-ped","S3v-ped",
          "S3w-ped","S3x-ped","S3y-ped","S3z-ped"))
mt_s2m <- merge(s2_out, mt_nm2, by.x = "Scenario", by.y = "scenario", all.x = TRUE)
for (mt_k in names(mt_SPLIT)) {
  mt_o <- do.call(rbind, lapply(mt_SPLIT[[mt_k]],
                                function(sc) mt_ord(mt_s2m[mt_s2m$Scenario == sc, ])))
  mt_pcd <- mt_o$Pct_Correct_Direction
  mt_pcd[is.na(mt_pcd) | mt_pcd == ""] <- mt_DASH
  mt_tb <- mt_o[, c("scenario_name","Family","Beta_True","True_HR","Method",
                    "Reject_Rate","Mean_HR","SD_HR")]
  mt_tb$pcd <- mt_pcd
  names(mt_tb) <- c("Scenario","Family", paste0("True effect size (", mt_B, ")"),
                    "True HR","Analytic approach","Rejection rate","Mean HR","SD HR",
                    "% Correct direction")
  write.csv(mt_tb, file.path(out_dir, sprintf("SuppTable_%s.csv", mt_k)), row.names = FALSE)
}

# ---- Table 8 (coefficients) ----
if (exists("s8_tab")) {
  mt_BIO <- c("CSF Aβ42/Aβ40" = "CSF Aβ42/Aβ40",
              "CSF p-tau181" = "CSF p-tau181", "Plasma p-tau181" = "plasma p-tau181",
              "Amyloid PET" = "amyloid PET", "Plasma p-tau217" = "plasma p-tau217")
  mt_PAR <- c("Sex_F" = "Sex (female)", "EDUC_z" = "Education (z)",
              "apoe4" = paste0("APOE", mt_EPS, "4"), "Z_std" = "AABC (z)",
              "Z_tv" = "Biomarker positivity",
              "Z_tv:A_z" = paste0("Biomarker positivity ", mt_TIMES, " AABC"))
  mt_COV <- c("Sex_F", "EDUC_z", "apoe4")
  mt_s8 <- s8_tab
  mt_s8$blk <- paste(mt_s8$Cohort, mt_s8$Biomarker, mt_s8$Model, sep = "||")
  mt_rows <- list(c("Parameter", "HR (95% CI; P)"))
  for (mt_b in unique(mt_s8$blk)) {
    mt_g <- mt_s8[mt_s8$blk == mt_b, ]
    mt_hdr <- sprintf("%s %s %s %s %s (n = %s; %s clinical onsets)",
                      mt_g$Cohort[1], mt_EN, mt_BIO[[mt_g$Biomarker[1]]], mt_EM,
                      mt_g$Model[1], mt_g$N[1], mt_g$Events[1])
    mt_rows[[length(mt_rows) + 1]] <- c(mt_hdr, "")
    mt_terms <- mt_g$Parameter[!mt_g$Parameter %in% mt_COV]
    for (mt_p in c(intersect(mt_COV, mt_g$Parameter), mt_terms)) {
      mt_r <- mt_g[mt_g$Parameter == mt_p, ]
      mt_rows[[length(mt_rows) + 1]] <- c(mt_PAR[[mt_p]],
                                          sprintf("%s; P %s", mt_r$`HR (95% CI)`, mt_r$P))
    }
  }
  write.table(do.call(rbind, mt_rows), file.path(out_dir, "SuppTable_8.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE)
}

cat("  Wrote SuppTable_1.csv ... SuppTable_8.csv (manuscript-exact, logical method order)\n")


# ==============================================================================
# Summary
# ==============================================================================

cat("\n===== All manuscript tables created =====\n")
cat("Output directory:", out_dir, "\n")
for (f in sort(list.files(out_dir, pattern = "Supp|Table"))) {
  cat("  ", f, "\n")
}
