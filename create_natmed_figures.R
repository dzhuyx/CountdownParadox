# ==============================================================================
# create_natmed_figures.R
#
# Publication-quality figures for the Nature Medicine manuscript.
# Output: PDF + PNG (300 dpi).
#
# Main figures:
#   Figure 1: TVC Schematic (external)
#   Figure 2: Study 1 — Type I error (3 panels: Standard, TV methods, sample size)
#   Figure 3: Study 2 — Type I error across spec sweeps + direction accuracy
#
# Supplementary:
#   Supp Figure 1: Type I error (15 Study 1 scenarios; S8 excluded, original S9/S10 IDs)
#   Supp Figure 2: HR boxplots under the null
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# -- Paths -------------------------------------------------------------------
project_dir <- Sys.getenv("CP_PROJECT_DIR", "/Users/daisyzhu/Documents/Research Projects/CountdownParadox_BiomarkerPositivity")
analysis_dir <- file.path(project_dir, "CountdownParadox_Analysis")
results_dir <- file.path(analysis_dir, "results")
sim_dir <- file.path(project_dir, "CountdownParadox_Manuscript_Simulations", "results")
fig_dir <- file.path(results_dir, "manuscript_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -- Color palette: color = method family, shape = coefficient ----------------
# Standard (biased) vs TV methods (valid); TV-AABC beta/gamma share green family
col_standard   <- "#D55E00"  # Vermillion: Standard/Countdown
col_tvc        <- "#0072B2"  # Blue: TV-BC
col_aabc_beta  <- "#009E73"  # Bluish green: TV-AABC beta
col_aabc_gamma <- "#66CC99"  # Light green: TV-AABC gamma

# Shapes
shp_standard   <- 16  # circle
shp_tvc        <- 17  # triangle
shp_aabc_beta  <- 15  # square
shp_aabc_gamma <- 18  # diamond

# Linetypes
lty_standard   <- "solid"
lty_tvc        <- "solid"
lty_aabc_beta  <- "solid"
lty_aabc_gamma <- "dashed"

# -- Named vectors for scale_*_manual (used across all figures) ---------------
method_colors <- c("Standard"       = col_standard,
                   "TV-BC"           = col_tvc,
                   "TV-AABC \u03b2"  = col_aabc_beta,
                   "TV-AABC \u03b3"  = col_aabc_gamma)
method_shapes <- c("Standard"       = shp_standard,
                   "TV-BC"           = shp_tvc,
                   "TV-AABC \u03b2"  = shp_aabc_beta,
                   "TV-AABC \u03b3"  = shp_aabc_gamma)
method_linetypes <- c("Standard"       = lty_standard,
                      "TV-BC"           = lty_tvc,
                      "TV-AABC \u03b2"  = lty_aabc_beta,
                      "TV-AABC \u03b3"  = lty_aabc_gamma)
method_levels <- c("Standard", "TV-BC", "TV-AABC \u03b2", "TV-AABC \u03b3")

# -- Theme (larger fonts for publication) -------------------
theme_manuscript <- theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    axis.title    = element_text(size = 13),
    axis.text     = element_text(size = 11),
    legend.position = "bottom",
    legend.title  = element_blank(),
    legend.text   = element_text(size = 11),
    plot.margin   = margin(10, 15, 10, 10)
  )

# -- Helper: save as PDF and PNG ---------------------------------------------
save_figure <- function(p, filename, width = 7, height = 5, dpi = 300) {
  pdf_path <- file.path(fig_dir, paste0(filename, ".pdf"))
  png_path <- file.path(fig_dir, paste0(filename, ".png"))
  ggsave(pdf_path, p, width = width, height = height, device = "pdf")
  ggsave(png_path, p, width = width, height = height, dpi = dpi, device = "png")
  cat(sprintf("Saved: %s (.pdf, .png)\n", filename))
}


# ==============================================================================
# FIGURE 1: TVC Schematic (single panel; degeneracy panel removed)
# ==============================================================================

cat("\n--- Figure 1: TVC Schematic ---\n")

# --- Left panel: TVC schematic with ~8 example subjects ---
schematic <- data.frame(
  subject  = 1:8,
  category = c("Already positive", "Already positive", "Already positive",
                "Transition", "Transition", "Transition",
                "Never positive", "Never positive"),
  entry    = c(55, 58, 52, 50, 55, 60, 48, 53),
  exit     = c(78, 85, 72, 82, 75, 88, 70, 80),
  eaoa     = c(48, 50, 45, 62, 65, 72, NA, NA),
  event    = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
)

# Build segment data for geom_segment
seg_list <- list()
for (i in seq_len(nrow(schematic))) {
  s <- schematic[i, ]
  if (s$category == "Already positive") {
    seg_list[[length(seg_list) + 1]] <- data.frame(
      subject = s$subject, xstart = s$entry, xend = s$exit,
      state = "After biomarker clock event")
  } else if (s$category == "Transition") {
    seg_list[[length(seg_list) + 1]] <- data.frame(
      subject = s$subject, xstart = s$entry, xend = s$eaoa,
      state = "No biomarker clock event yet")
    seg_list[[length(seg_list) + 1]] <- data.frame(
      subject = s$subject, xstart = s$eaoa, xend = s$exit,
      state = "After biomarker clock event")
  } else {
    seg_list[[length(seg_list) + 1]] <- data.frame(
      subject = s$subject, xstart = s$entry, xend = s$exit,
      state = "No biomarker clock event yet")
  }
}
seg_df <- do.call(rbind, seg_list)

event_pts <- schematic %>% filter(event) %>% select(subject, exit)

# ---- Panel (b): Time-varying analysis data structure (age scale, all subjects) ----
fig1_b <- ggplot(seg_df) +
  geom_segment(aes(x = xstart, xend = xend, y = subject, yend = subject,
                   color = state),
               linewidth = 3, lineend = "butt") +
  geom_point(data = event_pts, aes(x = exit, y = subject),
             shape = 4, size = 3, stroke = 1.5, color = "black") +
  scale_color_manual(
    values = c("No biomarker clock event yet" = "grey70",
               "After biomarker clock event" = col_standard),
    name = NULL
  ) +
  scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0.3, 9)) +
  geom_hline(yintercept = c(3.5, 6.5), linetype = "dotted", color = "grey80") +
  annotate("text", x = 49, y = 2, label = "Already\npositive",
           size = 2.8, fontface = "italic", color = "grey40", hjust = 0) +
  annotate("text", x = 49, y = 5, label = "Transition\nduring FU",
           size = 2.8, fontface = "italic", color = "grey40", hjust = 0) +
  annotate("text", x = 49, y = 7.5, label = "Never\npositive",
           size = 2.8, fontface = "italic", color = "grey40", hjust = 0) +
  annotate("text", x = 88, y = 8.8,
           label = expression(bold("\u00d7") * " = clinical onset"),
           size = 3, color = "black", hjust = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "Age (years)", y = NULL, title = "Time-varying analysis") +
  theme_manuscript +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.key.width = unit(1, "cm"),
        plot.margin = margin(10, 5, 10, 45))

# ---- Panel (a): Standard countdown analysis (same subjects) ----
# The standard analysis regresses the remaining time (age at onset - AABC) on the
# AABC, using biomarker-positive subjects only. The pre-positivity person-time is
# discarded (immortal time) and never-positive subjects are excluded entirely.
pos_sub <- schematic[!is.na(schematic$eaoa), ]   # S1-S6: positive
neg_sub <- schematic[is.na(schematic$eaoa), ]    # S7-S8: never positive

seg_a <- rbind(
  data.frame(subject = pos_sub$subject, xstart = pos_sub$entry, xend = pos_sub$eaoa,
             state = "Before positivity (discarded)"),
  data.frame(subject = pos_sub$subject, xstart = pos_sub$eaoa, xend = pos_sub$exit,
             state = "Remaining time (analyzed)")
)
seg_excl <- data.frame(subject = neg_sub$subject, xstart = neg_sub$entry,
                       xend = neg_sub$exit)
event_a <- pos_sub[pos_sub$event, c("subject", "exit")]

fig1_a <- ggplot() +
  # excluded never-positive subjects (faint, faded)
  geom_segment(data = seg_excl,
               aes(x = xstart, xend = xend, y = subject, yend = subject),
               color = "grey88", linewidth = 3, lineend = "butt", alpha = 0.6) +
  # positive subjects: discarded pre-positivity + analyzed remaining time
  geom_segment(data = seg_a,
               aes(x = xstart, xend = xend, y = subject, yend = subject, color = state),
               linewidth = 3, lineend = "butt") +
  # AABC (age at biomarker clock event) markers
  geom_point(data = pos_sub, aes(x = eaoa, y = subject),
             shape = 18, size = 2.8, color = "#7A3000") +
  geom_point(data = event_a, aes(x = exit, y = subject),
             shape = 4, size = 3, stroke = 1.5, color = "black") +
  scale_color_manual(
    values = c("Before positivity (discarded)" = "grey80",
               "Remaining time (analyzed)" = col_standard),
    name = NULL
  ) +
  scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0.3, 9)) +
  annotate("text", x = 49, y = 7.5, label = "Excluded\n(never positive)",
           size = 2.8, fontface = "italic", color = "grey55", hjust = 0) +
  annotate("text", x = 88, y = 8.8,
           label = expression(bold("\u00d7") * " = clinical onset"),
           size = 3, color = "black", hjust = 1) +
  annotate("point", x = 70, y = 8.2, shape = 18, size = 2.8, color = "#7A3000") +
  annotate("text", x = 71, y = 8.2, label = "= AABC",
           size = 3, color = "#7A3000", hjust = 0) +
  coord_cartesian(clip = "off") +
  labs(x = "Age (years)", y = NULL,
       title = "Standard countdown analysis") +
  theme_manuscript +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.key.width = unit(1, "cm"),
        plot.margin = margin(10, 5, 10, 45))

# ---- Figure 1: single panel — time-varying analysis data structure ----
# Only panel (b) is shown. fig1_a (the standard-analysis panel) is retained
# above but not composed: under the standard countdown analysis the gray
# pre-positivity segments are discarded and never-positive participants (whose
# follow-up is entirely gray) are excluded, as stated in the figure caption.
fig1 <- fig1_b + theme(legend.position = "bottom")

save_figure(fig1, "figure1_tvc_schematic", width = 7.5, height = 5.5)


# ==============================================================================
# FIGURE 2: Study 1 Simulation — Three panels
#   (a) Standard Type I error (0–100% scale)
#   (b) TV methods Type I error (zoomed ~0–15%)
#   (c) TV methods across sample sizes
# ==============================================================================

cat("\n--- Figure 2: Study 1 Simulation ---\n")

s1 <- read.csv(file.path(sim_dir, "study1", "summary_results.csv"),
               stringsAsFactors = FALSE)

# Exclude S8 (minimal overlap — TV methods not applicable). Then renumber to
# close the gap: original S9 → new S8, original S10 → new S9 (and n-variants).
s1 <- s1 %>% filter(scenario != "S8")
s1 <- s1 %>% mutate(
  scenario = case_when(
    scenario == "S9"        ~ "S8",
    scenario == "S10"       ~ "S9",
    scenario == "S10_n150"  ~ "S9_n150",
    scenario == "S10_n200"  ~ "S9_n200",
    scenario == "S10_n1000" ~ "S9_n1000",
    TRUE ~ scenario
  ),
  scenario_name = case_when(
    grepl("^S9: ",   scenario_name) ~ sub("^S9: ",   "S8: ", scenario_name),
    grepl("^S10: ",  scenario_name) ~ sub("^S10: ",  "S9: ", scenario_name),
    grepl("^S10_n", scenario_name) ~ sub("^S10_n", "S9_n", scenario_name),
    TRUE ~ scenario_name
  )
)

# Build long-form data with 4 method series
s1_long <- bind_rows(
  s1 %>%
    filter(method == "naive", !is.na(type1_error)) %>%
    mutate(t1e = type1_error, method_label = "Standard"),
  s1 %>%
    filter(method == "tvc_z_only", !is.na(type1_error)) %>%
    mutate(t1e = type1_error, method_label = "TV-BC"),
  s1 %>%
    filter(method == "tvc_interaction", !is.na(type1_error)) %>%
    mutate(t1e = type1_error, method_label = "TV-AABC \u03b2"),
  s1 %>%
    filter(method == "tvc_interaction", !is.na(type1_error_gamma)) %>%
    mutate(t1e = type1_error_gamma, method_label = "TV-AABC \u03b3")
) %>%
  mutate(scenario_short = gsub("_n", ", n=", sub(":.*", "", scenario_name)),
         method_label = factor(method_label, levels = method_levels))

# --- Distribution shape scenarios (S1-S5, S8, S9 in renumbered IDs) ---
distrib_scens <- c("S1", "S2", "S3", "S4", "S5", "S8", "S9")
distrib_labels <- c("S1" = "Unif\u2013Unif", "S2" = "Norm\u2013Norm",
                    "S3" = "Norm\u2013Weib", "S4" = "Gam\u2013Weib",
                    "S5" = "Beta\u2013Exp",  "S8" = "High Var(Z)",
                    "S9" = "BIOCARD\ncalibrated")

s1_distrib <- s1_long %>% filter(scenario %in% distrib_scens)
s1_distrib$scen_label <- factor(distrib_labels[s1_distrib$scenario],
                                levels = distrib_labels)

# --- Subset for TV methods (used in panels b and c) ---
# --- Panel (a): Standard only — bars, 0–100% scale ---
fig3a <- ggplot(s1_distrib %>% filter(method_label == "Standard"),
                aes(x = scen_label, y = t1e, fill = method_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  geom_col(width = 0.6) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     labels = percent) +
  labs(x = NULL, y = "Type I error rate",
       title = "Standard countdown analysis") +
  theme_manuscript +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# --- Panel (b): TV methods only — bars, 0–10% scale ---
fig3b <- ggplot(s1_distrib %>% filter(method_label != "Standard"),
                aes(x = scen_label, y = t1e, fill = method_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02),
                     labels = percent) +
  labs(x = NULL, y = "Type I error rate",
       title = "Time-varying methods") +
  theme_manuscript +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# --- Panel (c): TV methods across sample sizes — bars, grouped ---
# After renumbering: BIOCARD-calibrated variants are S9_* (was S10_*)
n_scens <- c("S1_n150", "S1_n200", "S1", "S1_n1000",
             "S9_n150", "S9_n200", "S9", "S9_n1000")
s1_nsize <- s1_long %>%
  filter(scenario %in% n_scens, method_label != "Standard")

n_map <- c("S1_n150" = "150", "S1_n200" = "200", "S1" = "500", "S1_n1000" = "1000",
           "S9_n150" = "150", "S9_n200" = "200", "S9" = "500", "S9_n1000" = "1000")
s1_nsize$n_label <- n_map[s1_nsize$scenario]
s1_nsize$base_scen <- ifelse(s1_nsize$scenario %in%
                               c("S1", "S1_n150", "S1_n200", "S1_n1000"),
                             "Unif\u2013Unif", "BIOCARD calibrated")

# Create a grouped x-axis: base_scen + n as facets
s1_nsize$n_label <- factor(s1_nsize$n_label, levels = c("150", "200", "500", "1000"))

fig3c <- ggplot(s1_nsize, aes(x = n_label, y = t1e, fill = method_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ base_scen, nrow = 1) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02),
                     labels = percent) +
  labs(x = "Sample size (n)", y = "Type I error rate",
       title = "Time-varying methods across sample sizes") +
  theme_manuscript +
  theme(strip.text = element_text(face = "bold", size = 11)) +
  guides(fill = "none")

# --- Study 1 panels (a/b/c) feed the merged manuscript Figure 2 (below) ---
# NOTE: We intentionally do NOT save a standalone 3-panel Study 1 figure. The
# manuscript Figure 2 is the merged 5-panel figure (figure2_combined, further
# below), which reuses panels fig3a/fig3b/fig3c. A separate
# "figure2_study1_natmed" file is not a manuscript figure, so it is not emitted.


# ==============================================================================
# FIGURE 3: Study 2 Simulation — Panel A (Type I error) + Panel B (Direction)
# ==============================================================================

cat("\n--- Figure 3: Study 2 Simulation ---\n")

s2 <- read.csv(file.path(sim_dir, "study2", "summary_results.csv"),
               stringsAsFactors = FALSE)

# -- Build 4 method series for Study 2 (same as Study 1) --
s2_long <- bind_rows(
  s2 %>% filter(method == "naive") %>%
    mutate(rej = reject_rate, method_label = "Standard"),
  s2 %>% filter(method == "tvc_z_only") %>%
    mutate(rej = reject_rate, method_label = "TV-BC"),
  s2 %>% filter(method == "tvc_interaction") %>%
    mutate(rej = reject_rate, method_label = "TV-AABC \u03b2"),
  s2 %>% filter(method == "tvc_interaction", !is.na(reject_rate_gamma)) %>%
    mutate(rej = reject_rate_gamma, method_label = "TV-AABC \u03b3")
) %>%
  mutate(method_label = factor(method_label,
    levels = c("Standard", "TV-BC", "TV-AABC \u03b2", "TV-AABC \u03b3")))

# --- Panel A: Type I error by spec sweep (4 sub-panels using null scenarios) ---

# Null scenario mapping to spec parameters
# n sweep (PED: base=500, S3h=200, S3i=1000, S3j=2000; BIO: base=500, S3h=150, S3i=1000)
null_n <- data.frame(
  scenario = c("S3a-ped", "S3h-ped-null", "S3i-ped-null", "S3j-ped-null",
               "S3a-bio", "S3h-bio-null", "S3i-bio-null"),
  spec_val = c(500, 200, 1000, 2000, 500, 150, 1000),
  spec_name = "Sample size",
  stringsAsFactors = FALSE
)

# Heterogeneity sweep (sigma_0: base=0.5, S3k=0.25, S3l=1.0)
null_het <- data.frame(
  scenario = c("S3a-ped", "S3k-ped-null", "S3l-ped-null",
               "S3a-bio", "S3k-bio-null", "S3l-bio-null"),
  spec_val = c(0.5, 0.25, 1.0, 0.5, 0.25, 1.0),
  spec_name = "Trajectory\nheterogeneity",
  stringsAsFactors = FALSE
)

# Disease rate sweep (lambda_0: PED base=0.03, S3m=0.02, S3n=0.05;
#                               BIO base=0.02, S3m=0.015, S3n=0.03)
null_lam <- data.frame(
  scenario = c("S3a-ped", "S3m-ped-null", "S3n-ped-null",
               "S3a-bio", "S3m-bio-null", "S3n-bio-null"),
  spec_val = c(0.03, 0.02, 0.05, 0.02, 0.015, 0.03),
  spec_name = "Baseline\nhazard",
  stringsAsFactors = FALSE
)

# Measurement precision sweep
# PED: base visit_interval=2/sigma_eps=0.25, S3o=freq(1), S3p=sparse(4), S3q=high_err(0.5)
# BIO: base visit_interval=4/sigma_eps=0.25, S3q=high_err(0.5)
# Use a categorical x-axis for measurement precision since it has mixed parameters
null_meas <- data.frame(
  scenario = c("S3o-ped-null", "S3a-ped", "S3p-ped-null", "S3q-ped-null",
               "S3a-bio", "S3q-bio-null"),
  spec_val = c(1, 2, 3, 4, 2, 4),  # ordinal: 1=frequent, 2=base, 3=sparse, 4=high error
  spec_name = "Measurement\nprecision",
  stringsAsFactors = FALSE
)
meas_labels <- c("1" = "Frequent\nvisits", "2" = "Base", "3" = "Sparse\nvisits",
                 "4" = "High\nmeas. error")

# Combine all spec sweeps
null_specs <- bind_rows(null_n, null_het, null_lam, null_meas)

# Join with s2_long to get rejection rates on nulls
panel_a_data <- s2_long %>%
  filter(beta_true == 0) %>%
  inner_join(null_specs, by = "scenario", relationship = "many-to-many") %>%
  mutate(spec_name = factor(spec_name,
    levels = c("Sample size",
               "Trajectory\nheterogeneity",
               "Baseline\nhazard",
               "Measurement\nprecision")))

# Add cohort to distinguish PED vs BIO at colliding spec_vals; each scenario
# gets its own x-position via a (P)/(B) suffix.
panel_a_data <- panel_a_data %>%
  mutate(cohort = factor(ifelse(grepl("-ped", scenario), "PED", "BIO"),
                         levels = c("PED", "BIO")))

# Create discrete x-axis labels per spec sweep, with cohort suffix on a new line
panel_a_data <- panel_a_data %>%
  mutate(spec_label_base = case_when(
    spec_name == "Sample size"               ~ as.character(spec_val),
    spec_name == "Trajectory\nheterogeneity" ~ as.character(spec_val),
    spec_name == "Baseline\nhazard"          ~ as.character(spec_val),
    spec_name == "Measurement\nprecision"    ~ c("1" = "Frequent", "2" = "Base",
                                                  "3" = "Sparse", "4" = "High err.")[as.character(spec_val)]
  ),
  spec_label = paste0(spec_label_base, " (", substr(as.character(cohort), 1, 1), ")"))

# Order spec_label within each facet by spec_val then cohort (PED before BIO)
panel_a_data <- panel_a_data %>%
  arrange(spec_name, spec_val, cohort) %>%
  mutate(spec_label = factor(spec_label, levels = unique(spec_label)))

# Sub-panels: facet by spec_name, discrete x-axis with dodged bars
fig4a <- ggplot(panel_a_data, aes(x = spec_label, y = rej,
                                  fill = method_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ spec_name, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     labels = percent) +
  labs(x = NULL, y = "Type I error rate",
       title = "Type I error across design specifications (null scenarios)") +
  theme_manuscript +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# --- Panel B: Direction accuracy by effect size (both families) ---
sweep_scens <- c("S3e-ped", "S3d-ped", "S3g-ped", "S3f-ped", "S3b-ped", "S3c-ped",
                 "S3e-bio", "S3d-bio", "S3g-bio", "S3f-bio", "S3b-bio", "S3c-bio")

# Direction accuracy for naive and tvc_z_only
s2_dir <- s2 %>%
  filter(scenario %in% sweep_scens,
         method %in% c("naive", "tvc_z_only"),
         beta_true != 0) %>%
  mutate(
    method_label = ifelse(method == "naive", "Standard", "TV-BC"),
    pct_correct = case_when(
      method == "naive" & beta_true > 0      ~ pct_hr_below_1 * 100,
      method == "naive" & beta_true < 0      ~ pct_hr_above_1 * 100,
      method == "tvc_z_only" & beta_true > 0 ~ pct_hr_above_1 * 100,
      method == "tvc_z_only" & beta_true < 0 ~ pct_hr_below_1 * 100),
    method_label = factor(method_label,
      levels = c("Standard", "TV-BC", "TV-AABC \u03b2", "TV-AABC \u03b3"))
  )

# Add TV-AABC beta direction accuracy (Z_tv coeff — same direction logic as TV-BC)
s2_dir_aabc_beta <- s2 %>%
  filter(scenario %in% sweep_scens,
         method == "tvc_interaction",
         beta_true != 0) %>%
  mutate(
    method_label = factor("TV-AABC \u03b2",
      levels = c("Standard", "TV-BC", "TV-AABC \u03b2", "TV-AABC \u03b3")),
    pct_correct = case_when(
      beta_true > 0 ~ pct_hr_above_1 * 100,
      beta_true < 0 ~ pct_hr_below_1 * 100)
  )

# Add TV-AABC gamma direction accuracy from per-replicate data
# gamma direction is INVERTED relative to beta: when beta_true > 0, earlier
# positivity (lower Z) → higher biomarker → higher hazard → gamma < 0 (HR_gamma < 1)
s2_rds_path <- file.path(sim_dir, "study2", "all_results.rds")
s2_rep <- readRDS(s2_rds_path)
s2_gamma_dir <- s2_rep %>%
  filter(scenario %in% sweep_scens,
         method == "tvc_interaction",
         beta_true != 0,
         !is.na(hr_gamma)) %>%
  group_by(scenario, beta_true) %>%
  summarise(
    pct_hr_above_1_gamma = mean(hr_gamma > 1, na.rm = TRUE),
    pct_hr_below_1_gamma = mean(hr_gamma < 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method_label = factor("TV-AABC \u03b3",
      levels = c("Standard", "TV-BC", "TV-AABC \u03b2", "TV-AABC \u03b3")),
    pct_correct = case_when(
      beta_true > 0 ~ pct_hr_below_1_gamma * 100,   # expect HR_gamma < 1
      beta_true < 0 ~ pct_hr_above_1_gamma * 100)    # expect HR_gamma > 1
  )

s2_dir_all <- bind_rows(s2_dir, s2_dir_aabc_beta, s2_gamma_dir)

# Convert beta_true to factor for bar grouping
s2_dir_all$beta_label <- factor(sprintf("%.1f", s2_dir_all$beta_true),
                                levels = sprintf("%.1f", sort(unique(s2_dir_all$beta_true))))

fig4b <- ggplot(s2_dir_all, aes(x = beta_label, y = pct_correct,
                                fill = method_label)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey60") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(-2, 105), breaks = seq(0, 100, 25),
                     labels = paste0(seq(0, 100, 25), "%")) +
  labs(x = expression("True " * beta[true] * " (biomarker effect)"),
       y = "Correct direction (%)",
       title = "Direction accuracy by effect size") +
  theme_manuscript +
  guides(fill = "none")

# --- Study 2 panels (d/e) feed the merged manuscript Figure 2 (below) ---
# NOTE: We intentionally do NOT save a standalone 2-panel Study 2 figure. Its
# panels fig4a/fig4b are merged into the manuscript Figure 2 (figure2_combined,
# below). Manuscript Figure 3 is figure_hr_aabc_panel, so a separate
# "figure3_study2_natmed" file is not a manuscript figure and is no longer
# emitted.


# ==============================================================================
# Combined Figure 2 (Study 1 + Study 2 simulation results merged into one figure)
#
# All 5 panels stacked vertically as a single column:
#   (a) Study 1 distribution shapes (Standard countdown)
#   (b) Study 1 distribution shapes (TV methods, magnified)
#   (c) Study 1 sample size (TV methods)
#   (d) Study 2 Type I error across design specifications
#   (e) Study 2 direction accuracy across effect sizes
#
# Design decisions:
#  - Vertical stack (single column) because each panel needs full width for
#    its x-axis labels (7-19 categories).
#  - Single shared legend at the bottom. This requires every panel that shows
#    a fill legend to declare the SAME 4-level scale (Standard, TV-BC,
#    TV-AABC β, TV-AABC γ). We rebuild fig3a and fig3b with `limits` set to
#    all 4 levels so patchwork::plot_layout(guides = "collect") can merge.
#  - Panels fig3c and fig4b already suppress their own fill legend.
# ==============================================================================
.method_levels <- c("Standard", "TV-BC", "TV-AABC β", "TV-AABC γ")

# Build the combined figure with a single shared legend AND compact panels
# so the entire 5-panel figure fits on one Nature Medicine page (~7.2 × 9 in
# usable). Compact strategy:
#  - Suppress fill legends on every panel except fig4a (canonical 4-level)
#  - Apply a compact theme via patchwork's `&` operator: smaller axis text
#    (7pt), smaller titles (8pt), tighter margins, smaller legend
#  - Reduce panel titles to one short line; rely on caption for full text
.compact_theme <- theme(
  plot.title = element_text(size = 8, face = "bold",
                            margin = margin(b = 2)),
  axis.title.x = element_text(size = 7, margin = margin(t = 1)),
  axis.title.y = element_text(size = 7, margin = margin(r = 1)),
  axis.text.x = element_text(size = 6.5),
  axis.text.y = element_text(size = 6.5),
  strip.text = element_text(size = 7, face = "bold",
                            margin = margin(t = 1, b = 1)),
  legend.text = element_text(size = 7),
  legend.key.size = unit(0.30, "cm"),
  legend.margin = margin(t = -2, b = 0),
  legend.box.margin = margin(t = -4, b = 0),
  plot.margin = margin(t = 2, r = 4, b = 0, l = 4),
  panel.spacing = unit(0.15, "lines")
)

fig3a_c <- fig3a + guides(fill = "none") + .compact_theme +
  theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1))
fig3b_c <- fig3b + guides(fill = "none") + .compact_theme +
  theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1))
fig3c_c <- fig3c + .compact_theme
fig4a_c <- fig4a +
  scale_fill_manual(values = method_colors,
                    limits = .method_levels,
                    drop = FALSE,
                    name = NULL) +
  .compact_theme +
  theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1))
fig4b_c <- fig4b + .compact_theme

combined_fig2 <- (fig3a_c / fig3b_c / fig3c_c / fig4a_c / fig4b_c) +
  plot_layout(ncol = 1, heights = c(1, 1, 1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(face = "bold", size = 9))) &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7))

save_figure(combined_fig2, "figure2_combined", width = 7.2, height = 9.2)


# ==============================================================================
# SUPPLEMENTARY FIGURE 1: Type I Error (all 15 Study 1 scenarios)
# ==============================================================================

cat("\n--- Supp Figure 1: Type I Error (4 method series, 15 scenarios) ---\n")

sf1_data <- s1_long

# Order scenarios
sf1_order <- sf1_data %>%
  distinct(scenario, scenario_short) %>%
  arrange(scenario) %>%
  pull(scenario_short)
sf1_data$scenario_short <- factor(sf1_data$scenario_short, levels = sf1_order)

# Panel (a): Standard only — 0–100% scale
sf1a <- ggplot(sf1_data %>% filter(method_label == "Standard"),
               aes(x = scenario_short, y = t1e, fill = method_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  geom_col(width = 0.6) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     labels = percent) +
  labs(x = NULL, y = "Type I error rate",
       title = "Standard countdown analysis") +
  theme_manuscript +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# Panel (b): TV methods only — 0–10% scale
sf1b <- ggplot(sf1_data %>% filter(method_label != "Standard"),
               aes(x = scenario_short, y = t1e, fill = method_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02),
                     labels = percent) +
  labs(x = NULL, y = "Type I error rate",
       title = "Time-varying methods") +
  theme_manuscript +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

supp_fig1 <- sf1a + sf1b +
  plot_layout(ncol = 2, guides = "collect", widths = c(1, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(face = "bold", size = 14))) &
  theme(legend.position = "bottom")

save_figure(supp_fig1, "supp_figure1_type1_error", width = 16, height = 6)


# ==============================================================================
# SUPPLEMENTARY FIGURE 2: HR Distribution Densities Under the Null
# ==============================================================================

cat("\n--- Supp Figure 2: HR Boxplots Under Null ---\n")

# Load per-replicate data for selected null scenarios
s1_all <- readRDS(file.path(sim_dir, "study1", "all_results.rds"))
s1_all <- s1_all %>% filter(scenario != "S8")
# Renumber to close the gap (S9 \u2192 S8, S10 \u2192 S9) for consistency with Fig 2 and Supp Table 1
s1_all <- s1_all %>% mutate(
  scenario = case_when(
    scenario == "S9"  ~ "S8",
    scenario == "S10" ~ "S9",
    TRUE ~ scenario
  )
)

# Select representative scenarios (note: original S10 is now S9 after relabel)
selected_scen <- c("S1", "S2", "S9")
sf2_data <- bind_rows(
  s1_all %>%
    filter(scenario %in% selected_scen, method == "naive",
           !is.na(hr), hr > 0, hr < 10) %>%
    mutate(hr_plot = hr, method_label = "Standard"),
  s1_all %>%
    filter(scenario %in% selected_scen, method == "tvc_z_only",
           !is.na(hr), hr > 0, hr < 10) %>%
    mutate(hr_plot = hr, method_label = "TV-BC"),
  s1_all %>%
    filter(scenario %in% selected_scen, method == "tvc_interaction",
           !is.na(hr), hr > 0, hr < 10) %>%
    mutate(hr_plot = hr, method_label = "TV-AABC \u03b2"),
  s1_all %>%
    filter(scenario %in% selected_scen, method == "tvc_interaction",
           !is.na(hr_gamma), hr_gamma > 0, hr_gamma < 10) %>%
    mutate(hr_plot = hr_gamma, method_label = "TV-AABC \u03b3")
) %>%
  mutate(
    method_label = factor(method_label, levels = method_levels),
    scenario_label = case_when(
      scenario == "S1" ~ "S1: Uniform\u2013Uniform",
      scenario == "S2" ~ "S2: Normal\u2013Normal",
      scenario == "S9" ~ "S9: BIOCARD-calibrated")
  )

supp_fig2 <- ggplot(sf2_data, aes(x = method_label, y = hr_plot,
                                   fill = method_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, width = 0.6) +
  facet_wrap(~ scenario_label, ncol = 3) +
  scale_fill_manual(values = method_colors) +
  coord_cartesian(ylim = c(0.8, 1.3)) +
  scale_y_continuous(breaks = seq(0.8, 1.3, 0.1)) +
  labs(x = NULL, y = "Hazard Ratio",
       title = "HR distributions under the null (1,000 replicates)") +
  theme_manuscript +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        legend.position = "none",
        plot.margin = margin(5, 10, 5, 5))

save_figure(supp_fig2, "supp_figure2_hr_boxplots", width = 12, height = 5)


# ==============================================================================
# Summary
# ==============================================================================

cat("\n=== All NatMed figures saved to:", fig_dir, "===\n")
cat("Files:\n")
for (f in sort(list.files(fig_dir, pattern = "figure"))) {
  cat("  ", f, "\n")
}
