###############################################################################
## R1_S6_beyond_curve.R — robustness check (A&D R1), pre-specified
## 2026-10-08 after the spaghetti/SILA plots showed AABCs projected far back for
## participants whose last biomarker value lies above the top of the fitted SILA
## trajectory (extrapolation beyond the curve). Rule (fixed before refitting):
## exclude from the TV analysis cohort the participants whose last value (SILA
## scale) exceeds max(tsila$val); refit TV-BC and TV-AABC with
## run_all_models_capture(). Counts and hazard ratios only.
## Repository copy of CountdownParadox_Analysis/R1_S6_beyond_curve_10082026.R; only the paths
## are changed. Needs the outputs of run_all.R (data/*.rda, results/full_coefficients_all_models.csv).
## Yuxin Zhu (drafted with Claude Code)
###############################################################################
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set CP_PROJECT_ROOT to the analysis root (the folder that holds data/ and results/).")
source(file.path(project_root, "conditional_survival_core.R"))
e1 <- new.env(); load(file.path(data_dir, "BIOCARD_CSF_SILA_intermediate.rda"), envir = e1)
e2 <- new.env(); load(file.path(data_dir, "BIOCARD_plasma_SILA_intermediate.rda"), envir = e2)
e3 <- new.env(); load(file.path(data_dir, "ADNI_SILA_intermediate_2026.rda"), envir = e3)
e4 <- new.env(); load(file.path(data_dir, "ADNI_plasma_SILA_intermediate.rda"), envir = e4)
spec <- list(CSF_AB42_AB40 = list(e1$resfit_last_ab, e1$res_sila_ab, "SUBJECT_ID"), CSF_pTau181 = list(e1$resfit_last_ptau, e1$res_sila_ptau, "SUBJECT_ID"),
             Plasma_pTau181 = list(e2$resfit_last_plasma, e2$res_sila_plasma, "SUBJECT_ID"), Amyloid_PET_FBP = list(e3$resfit_last, e3$res_sila, "RID"),
             Plasma_pTau217_Fuji = list(e4$resfit_last_fuj, e4$res_sila_fuj, "RID"))
hr <- function(fit, term) { s <- summary(fit); ci <- s$conf.int[term, ]; sprintf("%.2f (%.2f-%.2f), p=%.3g", ci[1], ci[3], ci[4], s$coefficients[term, "Pr(>|z|)"]) }
out <- list()
for (bm in names(fits)) { f <- fits[[bm]]; s <- spec[[bm]]; top <- max(s[[2]]$tsila$val)
  beyond <- s[[1]]$subid[s[[1]]$val > top]
  d <- f$dat; ex <- d[[s[[3]]]] %in% beyond
  r <- run_all_models_capture(d[!ex, ], f$cols$eaoa, f$cols$entry, f$cols$exit, f$cols$event, paste0(bm, "_S6"), f$cohort, f$cols$estpos)
  out[[bm]] <- data.frame(biomarker = bm, n_beyond_curve_in_cohort = sum(ex), events_excluded = sum(d[[f$cols$event]][ex]),
     min_AABC_main = min(f$aabc_pool), min_AABC_S6 = min(r$aabc_pool),
     main_TVBC = hr(f$fit_p2, "Z_tv"), S6_TVBC = hr(r$fit_p2, "Z_tv"),
     main_gamma = hr(f$fit_p3_bh, "AzZtv"), S6_gamma = hr(r$fit_p3_bh, "AzZtv")) }
o <- do.call(rbind, out); print(o, row.names = FALSE)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); od <- file.path(results_dir, paste0("R1_S6_", ts)); dir.create(od)
write.csv(o, file.path(od, "S6_beyond_curve.csv"), row.names = FALSE)
write.csv(data.frame(script = "R1_S6_beyond_curve.R", timestamp = ts, base_seed = NA, job_id = NA, r_version = R.version.string), file.path(od, "provenance.csv"), row.names = FALSE)
