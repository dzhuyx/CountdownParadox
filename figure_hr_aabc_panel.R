###############################################################################
## 5-panel stacked: HR vs AABC for all biomarker-cohort combinations
##
## Manuscript order, shared x-axis covering all mean AABCs, independent y-axes
## per panel. Each panel has full annotations (mean AABC line + label,
## HR-β point + label, subtitle). Shared legend at top.
##
## Input:  results/fits_p2_p3_all.rds  (from extract_full_coefficients.R)
## Output: figure_hr_aabc_panel.pdf / .png
###############################################################################

rm(list = ls())
library(survival)
library(ggplot2)
library(patchwork)

project_root <- Sys.getenv("CP_PROJECT_ROOT", "/Users/daisyzhu/Documents/Research Projects/CountdownParadox_BiomarkerPositivity/CountdownParadox_Analysis")
results_dir <- file.path(project_root, "results")
out_dir     <- file.path(results_dir, "manuscript_figures")

all_fits <- readRDS(file.path(results_dir, "fits_p2_p3_all.rds"))

# Shared x-axis: covers all 5 mean AABCs (56.8–75.3) with tight margin.
X_MIN <- 55
X_MAX <- 80

# Manuscript order + short titles
panels <- list(
      list(key = "BIOCARD_CSF_AB42_AB40",    title = "BIOCARD CSF Aβ42/40"),
      list(key = "BIOCARD_CSF_pTau181",      title = "BIOCARD CSF p-tau181"),
      list(key = "BIOCARD_Plasma_pTau181",   title = "BIOCARD Plasma p-tau181"),
      list(key = "ADNI_Amyloid_PET_FBP",     title = "ADNI Amyloid PET"),
      list(key = "ADNI_Plasma_pTau217_Fuji", title = "ADNI Plasma p-tau217")
)

method_colors    <- c("TV-AABC" = "#1f77b4", "TV-BC" = "grey35")
method_linetypes <- c("TV-AABC" = "solid",   "TV-BC" = "dashed")

# ---- Panel builder: full annotations (matches prototype) ---- #
make_panel <- function(pn, show_x_axis = FALSE) {
      bc <- all_fits[[pn$key]]
      fit_p2 <- bc$fit_p2
      fit_p3 <- bc$fit_p3
      A_mean <- bc$A_mean
      A_sd   <- bc$A_sd

      # Curve range: intersect observed with shared x-axis
      a_lo <- max(X_MIN, bc$aabc_range[1])
      a_hi <- min(X_MAX, bc$aabc_range[2])

      # TV-BC
      p2_sum        <- summary(fit_p2)
      tvbc_logHR    <- p2_sum$coefficients["Z_tv", "coef"]
      tvbc_logHR_se <- p2_sum$coefficients["Z_tv", "se(coef)"]
      tvbc_HR    <- exp(tvbc_logHR)
      tvbc_HR_lo <- exp(tvbc_logHR - 1.96 * tvbc_logHR_se)
      tvbc_HR_hi <- exp(tvbc_logHR + 1.96 * tvbc_logHR_se)

      # TV-AABC
      p3_coef <- coef(fit_p3)
      p3_vcov <- vcov(fit_p3)
      beta  <- p3_coef["Z_tv"]
      gamma <- p3_coef["Z_tv:A_z"]
      v_bb  <- p3_vcov["Z_tv", "Z_tv"]
      v_gg  <- p3_vcov["Z_tv:A_z", "Z_tv:A_z"]
      v_bg  <- p3_vcov["Z_tv", "Z_tv:A_z"]

      A_grid   <- seq(a_lo, a_hi, length.out = 200)
      z_grid   <- (A_grid - A_mean) / A_sd
      logHR_A  <- beta + gamma * z_grid
      se_A     <- sqrt(v_bb + 2 * z_grid * v_bg + z_grid^2 * v_gg)

      df_curve <- data.frame(
            A  = A_grid,
            HR = exp(logHR_A),
            lo = exp(logHR_A - 1.96 * se_A),
            hi = exp(logHR_A + 1.96 * se_A)
      )

      beta_HR <- exp(beta)

      # Panel-specific y range
      y_min <- min(c(df_curve$lo, tvbc_HR_lo))
      y_max <- max(c(df_curve$hi, tvbc_HR_hi))

      df_lines <- rbind(
            data.frame(method = "TV-AABC", A = df_curve$A, HR = df_curve$HR),
            data.frame(method = "TV-BC",   A = df_curve$A, HR = tvbc_HR)
      )
      df_lines$method <- factor(df_lines$method, levels = c("TV-AABC", "TV-BC"))

      p <- ggplot() +
            # TV-BC horizontal CI band
            annotate("rect", xmin = a_lo, xmax = a_hi,
                      ymin = tvbc_HR_lo, ymax = tvbc_HR_hi,
                      fill = "grey70", alpha = 0.12) +
            # Mean AABC vertical line
            geom_vline(xintercept = A_mean, linetype = "dotted",
                        color = "grey45", linewidth = 0.5) +
            # TV-AABC CI ribbon
            geom_ribbon(data = df_curve,
                         aes(x = A, ymin = lo, ymax = hi),
                         fill = "#1f77b4", alpha = 0.25) +
            # Both lines (mapped -> legend)
            geom_line(data = df_lines,
                       aes(x = A, y = HR, color = method, linetype = method),
                       linewidth = 0.9) +
            # HR-β marker at (mean AABC, exp(β))
            annotate("point", x = A_mean, y = beta_HR,
                      color = "#1f77b4", fill = "white",
                      shape = 21, size = 2.8, stroke = 1.0) +
            # HR-β label
            annotate("text", x = A_mean + 0.7, y = beta_HR * 1.25,
                      label = sprintf("HR-β = %.2f", beta_HR),
                      hjust = 0, size = 4.0, color = "#1f77b4",
                      fontface = "bold") +
            # Mean AABC label at bottom
            annotate("text", x = A_mean, y = y_min,
                      label = sprintf("Mean AABC = %.1f", A_mean),
                      color = "grey40", size = 2.7, hjust = -0.05, vjust = -0.3) +
            scale_color_manual(values = method_colors, name = NULL) +
            scale_linetype_manual(values = method_linetypes, name = NULL) +
            scale_x_continuous(limits = c(X_MIN, X_MAX), breaks = seq(55, 80, 5)) +
            scale_y_continuous(trans = "log",
                                breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)) +
            labs(title = pn$title, x = NULL, y = NULL) +
            theme_bw() +
            theme(
                  plot.title       = element_text(face = "bold", size = 10.5),
                  panel.grid.minor = element_blank(),
                  legend.position  = "none",
                  plot.margin      = margin(4, 10, 2, 10)
            )

      if (!show_x_axis) {
            p <- p + theme(axis.text.x  = element_blank(),
                            axis.ticks.x = element_blank())
      }
      p
}

# ---- Build 5 panels ---- #
plots <- lapply(seq_along(panels), function(i) {
      make_panel(panels[[i]],
                  show_x_axis = (i == length(panels)))
})

# ---- Compose with shared legend at top and shared x/y axis labels ---- #
combined <- wrap_plots(plots, ncol = 1) +
      plot_layout(guides = "collect") +
      plot_annotation(
            caption = "Age at biomarker clock (AABC, years)",
            theme = theme(
                  plot.caption = element_text(hjust = 0.5, size = 11,
                                               margin = margin(t = 6))
            )
      ) &
      theme(legend.position  = "bottom",
             legend.key.width = unit(1.1, "cm"),
             legend.text      = element_text(size = 10))

# Shared y-axis label on the left
combined <- wrap_elements(panel = combined) +
      labs(tag = "Hazard ratio (vs not-yet-positive at same age)") +
      theme(
            plot.tag          = element_text(size = 11, angle = 90),
            plot.tag.position = "left"
      )

out_pdf <- file.path(out_dir, "figure_hr_aabc_panel.pdf")
out_png <- file.path(out_dir, "figure_hr_aabc_panel.png")
ggsave(out_pdf, combined, width = 9, height = 13, device = cairo_pdf)
ggsave(out_png, combined, width = 9, height = 13, dpi = 300)
cat(sprintf("\nSaved: %s\n", out_pdf))
cat(sprintf("Saved: %s\n", out_png))
