###############################################################################
## tv_worked_example.R — Worked example of the time-varying (TV) analyses on
## simulated data
##
## Shows, without access to cohort data, how to structure the data and fit the
## two time-varying models of the manuscript, and why the standard countdown
## model misleads. Each simulated person has an age at study entry (the
## analytic baseline), an age at biomarker positivity (AABC) that may fall
## before entry, during follow-up or never, an age at MCI onset or censoring,
## and three covariates (sex, education, APOE e4).
##
## Truth in this simulation:
##   - biomarker positivity raises the hazard of onset by exp(0.8) = 2.23;
##   - the AABC does not modify that effect (HR-gamma = 1).
##
## Steps:
##   (1) build start-stop (counting-process) records on the age scale, by hand;
##   (2) build the same records with survival::tmerge() and check they agree;
##   (3) fit TV-BC and TV-AABC with survival::coxph();
##   (4) fit the standard countdown model for comparison.
## See Supplementary Note 6 of the manuscript (data structure, and the same
## models in SAS and Stata).
##
## Runtime: under a minute. Requires: R (>= 4.1) and the survival package.
## Usage: source("tv_worked_example.R") or Rscript tv_worked_example.R
##
## Author: Yuxin Zhu (drafted with Claude Code)
## Date: October 2026 (A&D revision, Reviewer 2 Major 1)
###############################################################################
library(survival)
set.seed(20261008)
n <- 1000
entry <- runif(n, 55, 75)                                # age at study entry (analytic baseline)
aabc  <- ifelse(runif(n) < 0.6, rnorm(n, 68, 8), Inf)    # 40% never become positive
exit_admin <- entry + runif(n, 3, 15)                    # end of follow-up if no event
sex <- rbinom(n, 1, 0.6); education <- rnorm(n); apoe4 <- rbinom(n, 1, 0.3)
lp_cov <- 0.5 * apoe4 - 0.1 * education                  # covariate effects (sex has none)

## onset hazard on the age scale: Gompertz-like baseline, x exp(0.8) after positivity
h0 <- function(a) 0.002 * exp(0.08 * (a - 55))
draw_onset <- function(e, z, x, lp) {                    # inversion on a fine age grid
  ages <- seq(e, x, by = 0.05); haz <- h0(ages) * exp(0.8 * (ages >= z) + lp)
  H <- cumsum(haz * 0.05); u <- -log(runif(1)); k <- which(H >= u)[1]
  if (is.na(k)) c(x, 0) else c(ages[k], 1)
}
sim <- t(mapply(draw_onset, entry, aabc, exit_admin, lp_cov))
d <- data.frame(id = 1:n, entry = entry, aabc = aabc, exit = sim[, 1], event = sim[, 2],
                sex = sex, education = education, apoe4 = apoe4)

## ---- (1) start-stop records on the age scale, by hand ----
## Positive at entry: one record with Z_tv = 1. Positive during follow-up: a
## not-yet-positive record up to the AABC (no event) and a positive record after
## it. Never positive during follow-up: one record with Z_tv = 0 (the AABC is not used).
mk <- function(r) {
  if (r$aabc <= r$entry)      data.frame(id = r$id, start = r$entry, stop = r$exit, event = r$event, Z_tv = 1, A = r$aabc)
  else if (r$aabc < r$exit)   rbind(data.frame(id = r$id, start = r$entry, stop = r$aabc, event = 0, Z_tv = 0, A = r$aabc),
                                    data.frame(id = r$id, start = r$aabc, stop = r$exit, event = r$event, Z_tv = 1, A = r$aabc))
  else                        data.frame(id = r$id, start = r$entry, stop = r$exit, event = r$event, Z_tv = 0, A = NA)
}
tv <- do.call(rbind, lapply(split(d, d$id), mk))
tv <- tv[tv$stop > tv$start, ]
tv <- merge(tv, d[, c("id", "sex", "education", "apoe4")], by = "id")
## standardize AABC among participants who contribute positive person-time;
## it enters only through its product with Z_tv, so it is 0 in not-yet-positive records
pos_A <- unique(tv[tv$Z_tv == 1, c("id", "A")])$A
tv$A_z <- ifelse(tv$Z_tv == 1, (tv$A - mean(pos_A)) / sd(pos_A), 0)
print(head(tv[order(tv$id), c("id", "start", "stop", "event", "Z_tv", "A_z")], 6))

## ---- (2) the same records with survival::tmerge() ----
tm <- tmerge(d[, c("id", "sex", "education", "apoe4")], d, id = id,
             mci = event(exit, event), tstart = entry, tstop = exit)
tm <- tmerge(tm, d[is.finite(d$aabc), ], id = id, Z_tv = tdc(aabc))
tm$Z_tv[is.na(tm$Z_tv)] <- 0
tm$A_z <- ifelse(tm$Z_tv == 1, (d$aabc[match(tm$id, d$id)] - mean(pos_A)) / sd(pos_A), 0)

## ---- (3) TV-BC and TV-AABC ----
fit_bc   <- coxph(Surv(start, stop, event) ~ Z_tv + sex + education + apoe4, data = tv)
fit_aabc <- coxph(Surv(start, stop, event) ~ Z_tv + Z_tv:A_z + sex + education + apoe4, data = tv)
fit_bc_tm   <- coxph(Surv(tstart, tstop, mci) ~ Z_tv + sex + education + apoe4, data = tm)
fit_aabc_tm <- coxph(Surv(tstart, tstop, mci) ~ Z_tv + Z_tv:A_z + sex + education + apoe4, data = tm)
stopifnot(isTRUE(all.equal(coef(fit_bc), coef(fit_bc_tm), tolerance = 1e-8)),
          isTRUE(all.equal(unname(coef(fit_aabc)), unname(coef(fit_aabc_tm)), tolerance = 1e-8)))
cat("\nRecords built by hand and by tmerge() give identical fits.\n")
print(summary(fit_bc)$conf.int); print(summary(fit_aabc)$conf.int)

## ---- (4) standard countdown model, for comparison ----
## restricted to participants positive before onset or censoring;
## outcome = remaining time from AABC to onset, predictor = standardized AABC
sd_ <- d[is.finite(d$aabc) & d$aabc < d$exit, ]
sd_$R <- sd_$exit - sd_$aabc; sd_$Z <- as.numeric(scale(sd_$aabc))
fit_std <- coxph(Surv(R, event) ~ Z + sex + education + apoe4, data = sd_)
print(summary(fit_std)$conf.int)

hr <- function(f, term) { ci <- summary(f)$conf.int[term, c(1, 3, 4)]; sprintf("%.2f (%.2f-%.2f)", ci[1], ci[2], ci[3]) }
cat("\nSummary (true values in brackets):\n",
    " TV-BC,   HR for positivity:            ", hr(fit_bc, "Z_tv"), " [2.23]\n",
    " TV-AABC, HR-beta for positivity:       ", hr(fit_aabc, "Z_tv"), " [2.23]\n",
    " TV-AABC, HR-gamma per SD of AABC:      ", hr(fit_aabc, "Z_tv:A_z"), " [1]\n",
    " Standard countdown, HR per SD of AABC: ", hr(fit_std, "Z"), " [1; any departure is the artifact]\n", sep = "")
