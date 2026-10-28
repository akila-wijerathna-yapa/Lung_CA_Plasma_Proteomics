if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 03, 07 loaded")

while (dev.cur() > 1) dev.off()

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")
coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x)); op <- options(warn = -1)
  out <- as.numeric(raw_chr); options(op); out
}
get_expr <- function(id) { v <- expr_log2_rep[, id]; names(v) <- rownames(expr_log2_rep); v }
akila_theme <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, colour = "grey30"))
}

gzma_id <- feature_meta %>% dplyr::filter(Gene_Symbol == "GZMA") %>% dplyr::pull(nELISA_ID) %>% unique()
gzma_expr <- get_expr(gzma_id)

# ── CV percentile for GZMA ───────────────────────────────────────
message("\n=== 08w | GZMA technical CV percentile ===")
cv_tbl <- readr::read_csv(file.path(paths$tables, "01_qc", "protein_cv_distribution.csv"), show_col_types = FALSE)
gzma_cv <- cv_tbl$cv[cv_tbl$Gene_Symbol == "GZMA"][1]
gzma_pct <- mean(cv_tbl$cv < gzma_cv, na.rm = TRUE)
message("GZMA CV = ", round(gzma_cv, 3), " | percentile vs panel = ", round(gzma_pct * 100, 1), "%")

# ── Check 1: Surgery-exclusion sensitivity check (univariate OS) ──
message("\n=== 08w | Surgery-exclusion sensitivity check (univariate OS Cox) ===")

run_os_cox_univariate <- function(sample_df, expr_vec, label) {
  d <- sample_df %>%
    dplyr::mutate(expr = expr_vec[sample_id], time = as.numeric(os_months), event = coerce_event_indicator(os_event)) %>%
    dplyr::filter(!is.na(expr), !is.na(time), !is.na(event))
  if (nrow(d) < 20 || sum(d$event) < 8) return(tibble::tibble(subset = label, n = nrow(d), n_events = NA, hr = NA, p_value = NA))
  fit <- survival::coxph(survival::Surv(time, event) ~ expr, data = d)
  s <- summary(fit)
  tibble::tibble(
    subset = label, n = nrow(d), n_events = sum(d$event),
    hr = s$coefficients[1, "exp(coef)"], hr_low = s$conf.int[1, "lower .95"], hr_high = s$conf.int[1, "upper .95"],
    p_value = s$coefficients[1, "Pr(>|z|)"]
  )
}

sx_sensitivity <- dplyr::bind_rows(
  run_os_cox_univariate(survival_df, gzma_expr, "All patients"),
  run_os_cox_univariate(survival_df %>% dplyr::filter(new_group != "Sx"), gzma_expr, "Surgery excluded")
)
message("GZMA OS Cox: all patients vs. surgery-excluded:"); print(sx_sensitivity)

# ── Check 2: cytotoxic/immune-activation program correlation + OS ──
message("\n=== 08w | GZMA program coherence (GZMK, GNLY, CD8A, FASLG, IL2RA) ===")

program_genes <- feature_meta %>%
  dplyr::filter(Gene_Symbol %in% c("GZMK", "GNLY", "CD8A", "FASLG", "IL2RA")) %>%
  dplyr::distinct(Gene_Symbol, nELISA_ID)
message("Cytotoxic/activation program members found on panel: ", paste(program_genes$Gene_Symbol, collapse = ", "))

run_corr <- function(x, y, label) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 10) return(tibble::tibble(pair = label, n = sum(keep), rho = NA_real_, p_value = NA_real_))
  ct <- cor.test(x[keep], y[keep], method = "spearman")
  tibble::tibble(pair = label, n = sum(keep), rho = unname(ct$estimate), p_value = ct$p.value)
}

program_corrs <- purrr::pmap_dfr(program_genes, function(Gene_Symbol, nELISA_ID) {
  run_corr(gzma_expr, get_expr(nELISA_ID), paste0("GZMA vs ", Gene_Symbol))
})
program_os <- dplyr::bind_rows(
  run_os_cox_univariate(survival_df, gzma_expr, "GZMA"),
  purrr::pmap_dfr(program_genes, function(Gene_Symbol, nELISA_ID) run_os_cox_univariate(survival_df, get_expr(nELISA_ID), Gene_Symbol))
)
message("GZMA correlations with cytotoxic program members:"); print(program_corrs)
message("Independent OS associations across the program:"); print(program_os)

# ── Save combined summary ───────────────────────────────────────
safe_write_csv(
  dplyr::bind_rows(
    tibble::tibble(check = "cv_percentile", gzma_cv = gzma_cv, gzma_cv_percentile = gzma_pct),
    sx_sensitivity %>% dplyr::mutate(check = "sx_exclusion_sensitivity"),
    program_corrs %>% dplyr::mutate(check = "program_correlation"),
    program_os %>% dplyr::mutate(check = "program_os_association")
  ),
  file.path(OUT_TBL, "gzma_coherence_and_sx_sensitivity.csv")
)

# ── Figure 1: surgery-exclusion sensitivity forest plot ─────────
message("\n=== 08w | Building sensitivity forest plot ===")

p_sens <- ggplot(sx_sensitivity, aes(x = subset, y = hr)) +
  geom_pointrange(aes(ymin = hr_low, ymax = hr_high), size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "GZMA OS Hazard Ratio: All Patients vs. Surgery Excluded",
       subtitle = paste0("All: HR=", signif(sx_sensitivity$hr[1], 3), ", p=", signif(sx_sensitivity$p_value[1], 3),
                          " | Sx excluded: HR=", signif(sx_sensitivity$hr[2], 3), ", p=", signif(sx_sensitivity$p_value[2], 3)),
       x = NULL, y = "Hazard Ratio (95% CI)") +
  akila_theme()

stub_sens <- file.path(OUT_FIG, "gzma_sx_sensitivity_forest")
cairo_pdf(paste0(stub_sens, ".pdf"), width = 6.5, height = 5); print(p_sens); dev.off()
png(paste0(stub_sens, ".png"), width = 6.5, height = 5, units = "in", res = 300, bg = "white"); print(p_sens); dev.off()
message("Saved: ", stub_sens, ".pdf / .png")

# ── Figure 2: cytotoxic-program correlation scatterplots ────────
message("\n=== 08w | Building coherence correlation figure ===")

corr_plots <- purrr::pmap(program_genes, function(Gene_Symbol, nELISA_ID) {
  row <- program_corrs[program_corrs$pair == paste0("GZMA vs ", Gene_Symbol), ]
  d <- tibble::tibble(x = gzma_expr, y = get_expr(nELISA_ID)) %>% dplyr::filter(!is.na(x), !is.na(y))
  ggplot(d, aes(x, y)) +
    geom_point(alpha = 0.5, size = 1.4, colour = "#2f6fb0") +
    geom_smooth(method = "lm", se = TRUE, colour = "#c9772b", linewidth = 0.6) +
    labs(title = paste0("GZMA vs ", Gene_Symbol, " (\u03c1=", signif(row$rho, 3), ", p=", signif(row$p_value, 3), ")"),
         x = "GZMA (log2)", y = paste0(Gene_Symbol, " (log2)")) +
    akila_theme()
})
p_corr <- patchwork::wrap_plots(corr_plots, ncol = 3) +
  patchwork::plot_annotation(title = "GZMA vs. Cytotoxic/Immune-Activation Program Members")

stub_corr <- file.path(OUT_FIG, "gzma_program_coherence_scatterplots")
cairo_pdf(paste0(stub_corr, ".pdf"), width = 12, height = 5); print(p_corr); dev.off()
png(paste0(stub_corr, ".png"), width = 12, height = 5, units = "in", res = 300, bg = "white"); print(p_corr); dev.off()
