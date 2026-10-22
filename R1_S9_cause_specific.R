###############################################################################
## R1_S9_cause_specific.R (A&D R1, Reviewer 2 Major 4a).
## Cause-specific analyses of AD-related MCI onset, as approved in
## CountdownParadox_Manuscript_A&D_R1/AD_R1_cause_specific_analysis_plan.md
## Section 9 answered by Daisy, 2026-10-17: Q1-Q3 (a), as planned. Q4 final (Daisy, 2026-10-20, after Corinne's
## reply): events with no recorded etiology are a third cause, unknown etiology (label "missing" below): they stay
## in the risk set and are censored at onset in the models of AD-related onset; no model is fitted for that cause,
## and the earlier coding (counted as AD-related) is not reported.
## Changed from R1_S9_cause_specific_10172026.R: only the AD-related event (AD only; was AD or missing) and the
## removal of the Q4 sensitivity fit (and this header).
##
##  (A) MCI events cross-classified by attributed cause (AD / other / missing)
##      and biomarker status at onset (the Table 2 groups: positive at onset =
##      positive at baseline or became positive before onset; became positive
##      only after onset; not positive at the last biomarker measurement).
##  (B) TV-BC and TV-AABC refitted with AD-related onset as the event; onset
##      from another cause censored at that onset (cause-specific hazard).
##      Everything else exactly as in the main analyses (same model runner).
## Cause is taken from the most recent diagnosis.
##   BIOCARD (consensus): AD primary = PROBADIF or POSSADIF; AD primary or
##     contributing = PROBAD or POSSAD (coding of BIOCARD_data_extraction.R);
##     a diagnosis with all four fields blank has a missing cause.
##   ADNI (site physician): MCI -> DXMDUE, dementia -> DXDDUE (1 = AD, 2 = other);
##     ADNI-1 dementia (no DXDDUE) -> DXAD = 1 means AD, else DXOTHDEM = 1 other.
##   Missing cause counts as AD-related (plan Q4) and is reported; where any event has a
##   missing cause, the models are refitted with it counted as other-cause (Q4 sensitivity).
## Gates: G0 (core), G1 (cause and status groups partition the events; every
##   AD-related event is an event; every event's most recent diagnosis is
##   impaired), G2 (the runner with the original event column reproduces the
##   main estimates exactly). Prints and writes aggregates only.
## Recycled: classify() and the refit pattern of R1_revision_analyses_10082026.R;
##   DXSUM handling of R1_S7/R1_S8; BIOCARD diagnosis reading of
##   BIOCARD_data_extraction.R.
## Repository copy of CountdownParadox_Analysis/R1_S9_cause_specific_10202026.R; only the paths
## are changed. Needs the outputs of run_all.R (data/*.rda, results/full_coefficients_all_models.csv),
## CP_BIOCARD_DIR (raw BIOCARD files) and ADNI_2026_data/DXSUM_11Feb2026.csv.
## Yuxin Zhu (analysis drafted with Claude Code), 2026-10-16; Q4 sensitivity 2026-10-17
###############################################################################
project_root <- Sys.getenv("CP_PROJECT_ROOT")
if (project_root == "") stop("CP_PROJECT_ROOT is not set. Run via run_all.R, or set CP_PROJECT_ROOT to the analysis root (the folder that holds data/ and results/).")
source(file.path(project_root, "conditional_survival_core.R"))
suppressPackageStartupMessages(library(readxl))

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(results_dir, paste0("R1_S9_", ts)); dir.create(out_dir, recursive = TRUE)
write.csv(data.frame(script = "R1_S9_cause_specific.R", timestamp = ts, base_seed = NA, job_id = NA,
                     r_version = R.version.string, survival_version = as.character(packageVersion("survival"))),
          file.path(out_dir, "provenance.csv"), row.names = FALSE)

id_col <- c(CSF_AB42_AB40 = "SUBJECT_ID", CSF_pTau181 = "SUBJECT_ID", Plasma_pTau181 = "SUBJECT_ID",
            Amyloid_PET_FBP = "RID", Plasma_pTau217_Fuji = "RID")
## ---- status at onset (recycled classify) ----
classify <- function(f) {
      d <- f$dat; cl <- f$cols
      A <- d[[cl$eaoa]]; en <- d[[cl$entry]]; ex <- d[[cl$exit]]; ev <- d[[cl$event]]
      ep <- if (!is.null(cl$estpos)) d[[cl$estpos]] else rep(NA, nrow(d))
      st <- ifelse(!is.na(ep) & !as.logical(ep), "a_not_positive_at_last_biomarker",
            ifelse(is.na(A), "e_estpos1_missing_AABC",
            ifelse(A <= en, "b_positive_at_baseline",
            ifelse(A < ex, "c_became_positive_during_followup",
            ifelse(ev == 1, "d1_positive_only_after_MCI_onset", "d2_positive_only_after_last_clinical_visit")))))
      data.frame(st = st, ev = ev, id = d[[id_col[[f$biomarker]]]])
}
grp <- function(st) ifelse(st %in% c("b_positive_at_baseline", "c_became_positive_during_followup"), "positive at onset",
                    ifelse(st == "d1_positive_only_after_MCI_onset", "became positive only after onset",
                    ifelse(st == "a_not_positive_at_last_biomarker", "not positive at last biomarker measurement", "other status")))

## ---- BIOCARD cause at the most recent diagnosis ----
biocard_dir <- Sys.getenv("CP_BIOCARD_DIR")
if (biocard_dir == "") stop("CP_BIOCARD_DIR is not set. Set it to the directory containing the raw BIOCARD data files.")
dx_b <- read_excel(path = file.path(biocard_dir, "BIOCARD_DiagnosisData_2024.09.08.xlsx"))
need_b <- c("SUBJECT_ID", "DIAGDATE", "PROBADIF", "POSSADIF", "PROBAD", "POSSAD")
stopifnot(all(need_b %in% names(dx_b)))
dx_b$DIAGDATE <- do.call("c", lapply(dx_b$DIAGDATE, function(x) as.Date(strsplit(as.character(x), " UTC")[[1]], "%Y-%m-%d")))
dx_b <- dx_b[!is.na(dx_b$DIAGDATE), ]
dx_b <- dx_b[order(dx_b$SUBJECT_ID, dx_b$DIAGDATE), ]
last_b <- dx_b[!duplicated(dx_b$SUBJECT_ID, fromLast = TRUE), ]
blank_b <- is.na(last_b$PROBADIF) & is.na(last_b$POSSADIF) & is.na(last_b$PROBAD) & is.na(last_b$POSSAD)
cause_b <- data.frame(id = last_b$SUBJECT_ID,
      AD_primary = ifelse(blank_b, "missing", ifelse(last_b$PROBADIF %in% 1 | last_b$POSSADIF %in% 1, "AD", "other")),
      AD_contrib = ifelse(blank_b, "missing", ifelse(last_b$PROBAD %in% 1 | last_b$POSSAD %in% 1, "AD", "other")))

## ---- ADNI cause at the most recent diagnosis ----
dx <- read.csv(file.path(project_root, "ADNI_2026_data", "DXSUM_11Feb2026.csv"), na.strings = -4)
need_a <- c("RID", "EXAMDATE", "DIAGNOSIS", "DXMDUE", "DXDDUE", "DXAD", "DXOTHDEM", "PHASE")
stopifnot(all(need_a %in% names(dx)))
dx$EXAMDATE <- as.Date(dx$EXAMDATE, "%Y-%m-%d"); dx <- dx[!is.na(dx$DIAGNOSIS) & !is.na(dx$EXAMDATE), ]
dx <- dx[order(dx$RID, dx$EXAMDATE), ]
last_a <- dx[!duplicated(dx$RID, fromLast = TRUE), ]
code <- ifelse(last_a$DIAGNOSIS == 2, last_a$DXMDUE, ifelse(last_a$DIAGNOSIS == 3, last_a$DXDDUE, NA))
adni1_dem <- last_a$DIAGNOSIS == 3 & is.na(last_a$DXDDUE)
code[adni1_dem] <- ifelse(last_a$DXAD[adni1_dem] %in% 1, 1, ifelse(last_a$DXOTHDEM[adni1_dem] %in% 1, 2, NA))
n_both_adni1 <- sum(adni1_dem & last_a$DXAD %in% 1 & last_a$DXOTHDEM %in% 1)
cause_a <- data.frame(id = last_a$RID, last_dx = last_a$DIAGNOSIS,
      suspected = ifelse(is.na(code), "missing", ifelse(code == 1, "AD", ifelse(code == 2, "other", "missing"))))

## ---- per analysis: (A) counts, (B) cause-specific refits ----
defs <- list(CSF_AB42_AB40 = c("AD_primary", "AD_contrib"), CSF_pTau181 = c("AD_primary", "AD_contrib"),
             Plasma_pTau181 = c("AD_primary", "AD_contrib"), Amyloid_PET_FBP = "suspected", Plasma_pTau217_Fuji = "suspected")
hr_row <- function(fit, term) { s <- summary(fit); ci <- s$conf.int[term, ]
      c(HR = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]), p = unname(s$coefficients[term, "Pr(>|z|)"])) }
cntA <- list(); hrB <- list()
for (bm in names(fits)) {
      f <- fits[[bm]]; cl <- f$cols; d <- f$dat; idc <- id_col[[bm]]
      x <- classify(f)
      cz <- if (f$cohort == "BIOCARD") cause_b else cause_a
      x <- merge(x, cz, by = "id", all.x = TRUE)
      ev <- x[x$ev == 1, ]
      ## G1: events' most recent diagnosis is impaired (ADNI) and causes are present for every event
      if (f$cohort == "ADNI") stopifnot(all(ev$last_dx %in% c(2, 3)))
      ## G2: the runner with the original event column reproduces the main fits
      r0 <- run_all_models_capture(d, cl$eaoa, cl$entry, cl$exit, cl$event, paste0(bm, "_G2"), f$cohort, cl$estpos)
      stopifnot(isTRUE(all.equal(coef(r0$fit_p2), coef(f$fit_p2), tolerance = 1e-10)),
                isTRUE(all.equal(coef(r0$fit_p3_bh), coef(f$fit_p3_bh), tolerance = 1e-10)))
      for (df in defs[[bm]]) {
            cause <- ev[[df]]; stopifnot(!any(is.na(cause)))
            tab <- as.data.frame(table(cause = factor(cause, levels = c("AD", "other", "missing")),
                                       status = factor(grp(ev$st), levels = c("positive at onset", "became positive only after onset",
                                                                              "not positive at last biomarker measurement", "other status"))))
            stopifnot(sum(tab$Freq) == sum(x$ev))                                   # G1: partition
            cntA[[paste(bm, df)]] <- data.frame(biomarker = bm, definition = df, tab)
            ## cause-specific event: AD-related onset only; other causes and unknown etiology ("missing") are censored at onset (plan Q4, 2026-10-20)
            ad_ids <- ev$id[cause == "AD"]
            d$event_cs <- ifelse(d[[cl$event]] == 1 & d[[idc]] %in% ad_ids, 1, 0)
            stopifnot(all(d$event_cs <= d[[cl$event]]), sum(d$event_cs) == length(ad_ids))   # G1
            r <- run_all_models_capture(d, cl$eaoa, cl$entry, cl$exit, "event_cs", paste0(bm, "_CS_", df), f$cohort, cl$estpos)
            hrB[[paste(bm, df)]] <- data.frame(biomarker = bm, definition = df, events_AD_related = sum(d$event_cs),
                  events_all = sum(d[[cl$event]]), events_missing_cause = sum(cause == "missing"),
                  term = c("TV-BC HR-beta", "TV-AABC HR-beta", "TV-AABC HR-gamma"),
                  rbind(hr_row(r$fit_p2, "Z_tv"), hr_row(r$fit_p3_bh, "Z_tv"), hr_row(r$fit_p3_bh, "AzZtv")))
      }
      hrB[[paste(bm, "all_cause")]] <- data.frame(biomarker = bm, definition = "all_cause_main", events_AD_related = NA,
            events_all = sum(d[[cl$event]]), events_missing_cause = NA,
            term = c("TV-BC HR-beta", "TV-AABC HR-beta", "TV-AABC HR-gamma"),
            rbind(hr_row(f$fit_p2, "Z_tv"), hr_row(f$fit_p3_bh, "Z_tv"), hr_row(f$fit_p3_bh, "AzZtv")))
}
A <- do.call(rbind, cntA); B <- do.call(rbind, hrB)
write.csv(A, file.path(out_dir, "S9A_events_by_cause_and_status.csv"), row.names = FALSE)
write.csv(B, file.path(out_dir, "S9B_cause_specific_hr.csv"), row.names = FALSE)
cat("\n[S9] G0 and G2 passed; G1 partitions passed.\n")
cat("[S9] ADNI-1 dementia diagnoses with both DXAD and DXOTHDEM = 1 (classified AD):", n_both_adni1, "\n")
cat("[S9] BIOCARD most recent diagnoses with no etiology recorded (all participants, not only events):", sum(blank_b), "\n")
cat("\n[S9] (A) events by cause and status at onset\n"); print(A[A$Freq > 0, ], row.names = FALSE)
cat("\n[S9] (B) hazard ratios\n"); print(B, row.names = FALSE, digits = 3)
cat("\nsaved", out_dir, "\n")
