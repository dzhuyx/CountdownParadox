###############################################################################
## shared_utils.R — Shared helper functions for BIOCARD analysis scripts
##
## Sourced by: BIOCARD_data_extraction.R, BIOCARD_csf_sila.R,
##             BIOCARD_plasma_sila.R
##
## Functions:
##   FindJump() — detect consecutive equal values in a vector
##   GetSurv()  — extract survival information from a single BIOCARD subject
##
## Author: Yuxin Zhu
## Date: April 2026 (refactored from data_extraction.R)
###############################################################################

# FindJump: detect consecutive equal values (used by GetSurv)
FindJump <- function(x) {
      x <- c(x)
      temp <- NULL
      if (length(x) > 1) {
            for (i in 1:(length(x) - 1)) {
                  temp <- c(temp, x[i] == x[i + 1])
            }
      } else {
            temp <- FALSE
      }
      temp
}

# GetSurv: get survival information from single individual
GetSurv <- function(x) {
      result <- NULL
      x <- x[order(x$DIAGDATE, decreasing = F), ]
      x <- x[!is.na(x$DIAGDATE) & !is.na(x$DIAGNOSIS), ]
      if (all(!is.na(x$DIAGDATE) & !is.na(x$DIAGNOSIS))) {
            if (x[x$DIAGDATE == max(x$DIAGDATE), "DIAGNOSIS"] %in% c("NORMAL", "IMPAIRED NOT MCI")) {
                  d <- 0
                  baseline.age <- as.numeric(min(x$DIAGDATE - x$DOB[1])) / 365.25
                  onset.age <- as.numeric(max(x$DIAGDATE - x$DOB[1])) / 365.25
                  censor.age <- onset.age
                  result <- data.frame(SUBJECT_ID = x$SUBJECT_ID[1],
                                       LETTERCODE = x$LETTERCODE[1],
                                       baseline.age = baseline.age,
                                       onset.age = onset.age,
                                       d = d,
                                       diag = unlist(x[x$DIAGDATE == max(x$DIAGDATE), "DIAGNOSIS"]),
                                       DOB = x$DOB[1],
                                       censor.age = censor.age)
            } else if (x[x$DIAGDATE == max(x$DIAGDATE), "DIAGNOSIS"] %in% c("MCI", "DEMENTIA")) {
                  d <- 1
                  temp <- as.numeric(x$DIAGNOSIS %in% c("MCI", "DEMENTIA"))
                  onset.position <- nrow(x) - which(!FindJump(rev(temp)))[1] + 1
                  if (is.na(onset.position)) {
                        onset.position <- 1
                  }
                  baseline.age <- as.numeric(min(x$DIAGDATE - x$DOB[1])) / 365.25
                  onset.age <- x$DECAGE[onset.position]
                  censor.age <- as.numeric(max(x$DIAGDATE - x$DOB[1])) / 365.25
                  if (is.na(onset.age)) {
                        baseline.age <- as.numeric(min(x$DIAGDATE - x$DOB[1])) / 365.25
                        onset.age <- as.numeric(x$DIAGDATE[onset.position] - x$DOB[1])
                        onset.age <- onset.age / 365.25
                  }
                  result <- data.frame(SUBJECT_ID = x$SUBJECT_ID[1],
                                       LETTERCODE = x$LETTERCODE[1],
                                       baseline.age = baseline.age,
                                       onset.age = onset.age,
                                       d = d,
                                       diag = unlist(x[x$DIAGDATE == max(x$DIAGDATE), "DIAGNOSIS"]),
                                       DOB = x$DOB[1],
                                       censor.age = censor.age)
            } else {
                  print(x$SUBJECT_ID[1])
                  result <- NULL
            }
      }
      return(result)
}
