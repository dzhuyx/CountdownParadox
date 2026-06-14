###############################################################################
## adni_event_types.R
##
## Among the ADNI analysis-cohort progressors (d == 1), tabulate whether the
## first impaired diagnosis was a direct CN -> dementia transition (no documented
## MCI) vs. the standard CN -> MCI path. This persists the manuscript statement
## (Results P029 / Methods P049):
##   ADNI amyloid PET:      7/106 (6.6%) direct-to-dementia
##   ADNI plasma p-tau217:  7/162 (4.3%) direct-to-dementia
##
## Computed on the analysis cohorts used for all other ADNI analyses (575 PET,
## 801 plasma). The denominator for the plasma direct-to-dementia count is the
## 162 progressors in the 801-subject plasma cohort.
##
## Input:  ADNI_2026_data/DXSUM_11Feb2026.csv
##         data/ADNI_pet_analysis_cohort.rda    (dat_complete,      n = 575)
##         data/ADNI_plasma_analysis_cohort.rda (dat_plasma_cohort, n = 801)
## Output: results/adni_event_types.csv
###############################################################################

rm(list = ls())
suppressMessages(library(dplyr))

# ---- Paths ---- #
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set CP_PROJECT_ROOT to the analysis root (the folder that holds data/ and results/).")
data_dir  <- file.path(project_root, "data")
data_2026 <- file.path(project_root, "ADNI_2026_data")
res_dir   <- file.path(project_root, "results")
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# ---- Full diagnosis trajectory per RID (to find each subject's first impaired Dx) ---- #
dat_dx <- read.csv(file.path(data_2026, "DXSUM_11Feb2026.csv"), na.strings = -4)
dat_dx$EXAMDATE <- as.Date(dat_dx$EXAMDATE, "%Y-%m-%d")
dat_dx$dx <- dat_dx$DIAGNOSIS   # 1 = CN, 2 = MCI, 3 = dementia

# Classify a subject's event path by the FIRST diagnosis in {2 (MCI), 3 (dementia)}
classify_event_path <- function(rid) {
  x <- dat_dx[dat_dx$RID == rid, ]
  x <- x[complete.cases(x[, c("dx", "EXAMDATE")]), ]
  x <- x[order(x$EXAMDATE), ]
  pos <- which(x$dx %in% c(2, 3))
  if (length(pos) == 0) return(NA_character_)
  if (x$dx[pos[1]] == 3) return("AD_direct")  # dementia without documented MCI
  if (x$dx[pos[1]] == 2) return("MCI_first")  # standard CN -> MCI path
  NA_character_
}

# ---- Summarize one analysis cohort ---- #
summarize_cohort <- function(rda_file, obj_name, cohort_label) {
  e <- new.env(); load(rda_file, envir = e)
  d <- get(obj_name, envir = e)
  ev_rids <- d$RID[d$d == 1]
  paths   <- sapply(ev_rids, classify_event_path)
  n_ad  <- sum(paths == "AD_direct", na.rm = TRUE)
  n_mci <- sum(paths == "MCI_first", na.rm = TRUE)
  data.frame(
    cohort        = cohort_label,
    n_events      = length(ev_rids),
    n_ad_direct   = n_ad,
    n_mci_first   = n_mci,
    pct_ad_direct = round(100 * n_ad / length(ev_rids), 1),
    stringsAsFactors = FALSE
  )
}

cat("================================================================\n")
cat("ADNI EVENT-TYPE SUMMARY (progressor first-impaired-diagnosis path)\n")
cat("================================================================\n\n")

res <- rbind(
  summarize_cohort(file.path(data_dir, "ADNI_pet_analysis_cohort.rda"),
                   "dat_complete",      "ADNI_Amyloid_PET"),
  summarize_cohort(file.path(data_dir, "ADNI_plasma_analysis_cohort.rda"),
                   "dat_plasma_cohort", "ADNI_Plasma_pTau217")
)

write.csv(res, file.path(res_dir, "adni_event_types.csv"), row.names = FALSE)
cat("Saved: results/adni_event_types.csv\n")
print(res, row.names = FALSE)
