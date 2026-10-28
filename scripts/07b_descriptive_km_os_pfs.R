# ── Preamble ─────────────────────────────────────────────────
PROJECT_DIR <- normalizePath(file.path(dirname(
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    rstudioapi::getActiveDocumentContext()$path
  else {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- grep("^--file=", args, value = TRUE)
    if (length(f)) sub("^--file=", "", f[1]) else stop("Cannot detect script path.")
  }
), ".."), mustWork = FALSE)

sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

# ── Load checkpoint ────────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))

need <- c("survival", "survminer", "dplyr", "ggplot2", "gridExtra")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(survival); library(survminer); library(dplyr); library(ggplot2)

figure_dir <- file.path(paths$figures, "06_analysis5_survival")
table_dir  <- file.path(paths$tables,  "06_analysis5_survival")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)

## ---- 1. Assemble OS/PFS analysis sets ------------------------------------ ##

os_df <- clinical %>%
  filter(!is.na(os_months), !is.na(os_event)) %>%
  transmute(patient_id, time = os_months, event = os_event)

pfs_df <- clinical %>%
  filter(!is.na(pfs_months_main), !is.na(pfs_event_main)) %>%
  transmute(patient_id, time = pfs_months_main, event = pfs_event_main)

cat("CHECKPOINT B1: OS analysis set  ->", nrow(os_df),  "patients,",
    sum(os_df$event == 1),  "events\n")
cat("CHECKPOINT B1: PFS analysis set ->", nrow(pfs_df), "patients,",
    sum(pfs_df$event == 1), "events\n")

## ---- 2. Plain descriptive KM (no stratification) -------------------------- ##

fit_os  <- survival::survfit(Surv(time, event) ~ 1, data = os_df)
fit_pfs <- survival::survfit(Surv(time, event) ~ 1, data = pfs_df)

make_descriptive_km <- function(fit, df, label, colour) {
  survminer::ggsurvplot(
    fit, data = df, risk.table = TRUE, conf.int = TRUE,
    palette = colour, legend = "none",
    ggtheme = ggplot2::theme_bw(base_size = 11),
    title = paste0("Whole-cohort ", label),
    xlab = "Months", ylab = paste0("Probability of ", label),
    risk.table.y.text = FALSE, risk.table.height = 0.25,
    surv.median.line = "hv")
}

km_os  <- make_descriptive_km(fit_os,  os_df,  "overall survival",           "#2166AC")
km_pfs <- make_descriptive_km(fit_pfs, pfs_df, "progression-free survival",  "#B2182B")

# Individual saves (curve + risk table), consistent with 07_analysis5_survival.R
save_km <- function(km, stub) {
  cairo_pdf(paste0(stub, ".pdf"), width = 6.5, height = 5.5)
  print(km$plot); dev.off()
  png(paste0(stub, ".png"), width = 6.5, height = 5.5, units = "in", res = 300, bg = "white")
  print(km$plot); dev.off()
  cairo_pdf(paste0(stub, "_risk.pdf"), width = 6.5, height = 2.2)
  print(km$table); dev.off()
  png(paste0(stub, "_risk.png"), width = 6.5, height = 2.2, units = "in", res = 300, bg = "white")
  print(km$table); dev.off()
}
save_km(km_os,  file.path(figure_dir, "descriptive_km_os"))
save_km(km_pfs, file.path(figure_dir, "descriptive_km_pfs"))

# Combined A/B panel (matches the two-panel literature layout)
combined <- survminer::arrange_ggsurvplots(
  list(km_os, km_pfs), print = FALSE, ncol = 2, nrow = 1,
  labels = c("A", "B"))

cairo_pdf(file.path(figure_dir, "descriptive_km_os_pfs_panel.pdf"), width = 12, height = 6)
print(combined); dev.off()
png(file.path(figure_dir, "descriptive_km_os_pfs_panel.png"),
    width = 12, height = 6, units = "in", res = 300, bg = "white")
print(combined); dev.off()

cat("CHECKPOINT B2: descriptive KM panel saved to", figure_dir, "\n")

## ---- 3. Median survival summary table ------------------------------------- ##

summarise_fit <- function(fit, endpoint) {
  s <- summary(fit)$table
  tibble::tibble(
    endpoint    = endpoint,
    n           = s[["records"]],
    n_events    = s[["events"]],
    median_time = unname(s[["median"]]),
    lower_95    = unname(s[["0.95LCL"]]),
    upper_95    = unname(s[["0.95UCL"]]))
}
descriptive_summary <- dplyr::bind_rows(
  summarise_fit(fit_os,  "OS"),
  summarise_fit(fit_pfs, "PFS"))
write.csv(descriptive_summary,
          file.path(table_dir, "descriptive_km_summary.csv"), row.names = FALSE)
print(descriptive_summary)

## ---- 4. Multivariable clinical Cox PH model (no proteins) ----------------- ##
# Candidate covariates; script drops any that are constant, all-NA, or cause
# coxph to fail (e.g. a level with zero events) rather than halting.

candidate_covars <- c("age_at_diagnosis", "sex", "histology", "stage_collapsed", "new_group")

build_multivariable_cox <- function(base_df, clinical_df, time_col, event_col, covars) {
  df <- base_df %>%
    dplyr::left_join(
      clinical_df %>% dplyr::select(patient_id, dplyr::all_of(covars)),
      by = "patient_id")

  usable <- covars[sapply(covars, function(v) {
    x <- df[[v]]
    !all(is.na(x)) && dplyr::n_distinct(x, na.rm = TRUE) > 1
  })]
  dropped <- setdiff(covars, usable)
  if (length(dropped)) cat("  Dropping covariate(s) with no variation/all-NA:",
                            paste(dropped, collapse = ", "), "\n")

  df <- df[stats::complete.cases(df[, usable, drop = FALSE]), ]
  f <- as.formula(paste0("Surv(time, event) ~ ", paste(usable, collapse = " + ")))
  fit <- tryCatch(survival::coxph(f, data = df), error = function(e) {
    cat("  coxph failed:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(fit)) return(list(fit = NULL, tbl = tibble::tibble()))

  s <- summary(fit)
  tbl <- tibble::tibble(
    term       = rownames(s$coefficients),
    HR         = s$coefficients[, "exp(coef)"],
    lower_95   = s$conf.int[, "lower .95"],
    upper_95   = s$conf.int[, "upper .95"],
    p_value    = s$coefficients[, "Pr(>|z|)"],
    n          = s$n,
    n_events   = s$nevent)
  list(fit = fit, tbl = tbl, df = df)
}

cat("Multivariable Cox — OS:\n")
mv_os  <- build_multivariable_cox(os_df,  clinical, "time", "event", candidate_covars)
cat("Multivariable Cox — PFS:\n")
mv_pfs <- build_multivariable_cox(pfs_df, clinical, "time", "event", candidate_covars)

if (nrow(mv_os$tbl))
  write.csv(mv_os$tbl,  file.path(table_dir, "multivariable_cox_OS.csv"),  row.names = FALSE)
if (nrow(mv_pfs$tbl))
  write.csv(mv_pfs$tbl, file.path(table_dir, "multivariable_cox_PFS.csv"), row.names = FALSE)

print(mv_os$tbl)
print(mv_pfs$tbl)

## ---- 5. Forest plots for the multivariable models -------------------------- ##

make_clinical_forest <- function(tbl, title, path_stub) {
  if (nrow(tbl) == 0) return(invisible(NULL))
  tbl <- tbl %>% mutate(term = factor(term, levels = rev(term)))
  p <- ggplot(tbl, aes(x = HR, y = term)) +
    geom_point(size = 2.4) +
    geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_x_log10() +
    labs(x = "Hazard ratio (log scale)", y = NULL, title = title) +
    theme_bw(base_size = 11)
  cairo_pdf(paste0(path_stub, ".pdf"), width = 7, height = 0.4 * nrow(tbl) + 2)
  print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 7, height = 0.4 * nrow(tbl) + 2,
      units = "in", res = 300, bg = "white")
  print(p); dev.off()
}
make_clinical_forest(mv_os$tbl,  "Multivariable Cox PH — OS",  file.path(figure_dir, "multivariable_cox_forest_OS"))
make_clinical_forest(mv_pfs$tbl, "Multivariable Cox PH — PFS", file.path(figure_dir, "multivariable_cox_forest_PFS"))

cat("07b_descriptive_km_os_pfs.R complete. Tables in:", table_dir,
    "| Figures in:", figure_dir, "\n")
