sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

ALTML_TBL_DIR <- file.path(paths$tables, "07_analysis6_machine_learning", "altml")
CV_TBL_DIR    <- file.path(paths$tables, "07_analysis6_machine_learning", "cv")
OUT_TBL       <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG       <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# All 9 binomial comparisons (5 response + 4 recurrence), with a clean
# display label and explicit ordering (response first, then recurrence,
# matching the pipeline's own run order).
COMPARISON_ORDER <- tibble::tribble(
  ~suffix,                  ~display_label,
  "01_CR_vs_PD",            "CR vs PD",
  "02_CRPR_vs_PD",          "CR/PR vs PD",
  "03_CRPR_vs_SDPD",        "CR/PR vs SD/PD",
  "04_CRPRSD_vs_PD",        "CR/PR/SD vs PD",
  "05_PR_vs_PD",            "PR vs PD",
  "06_CTx_recurrence",      "CTx recurrence",
  "07_Sx_recurrence",       "Sx recurrence",
  "08_IO_recurrence",       "IO recurrence",
  "09_CTxIO_recurrence",    "CTx+IO recurrence"
)

# ── Load in-sample results (altml) ─────────────────────────────
message("\n=== 08m | Loading in-sample (altml) results ===")

insample_df <- purrr::map_dfr(COMPARISON_ORDER$suffix, function(suf) {
  f <- file.path(ALTML_TBL_DIR, paste0("altml_", suf, "_all_methods.csv"))
  if (!file.exists(f)) {
    warning("Missing altml file for suffix: ", suf)
    return(NULL)
  }
  row1 <- readr::read_csv(f, show_col_types = FALSE, n_max = 1)
  tibble::tibble(
    suffix = suf,
    method = c("Elastic Net", "sPLS-DA", "Boruta+RF"),
    auc_insample = c(row1$en_model_auc, row1$splsda_train_accuracy, row1$rf_model_auc)
  )
})

message("In-sample rows loaded: ", nrow(insample_df))

# ── Load CV results ─────────────────────────────────────────────
message("\n=== 08m | Loading cross-validated (CV) results ===")

cv_csv_path <- file.path(CV_TBL_DIR, "cv_binomial_all_comparisons.csv")
if (!file.exists(cv_csv_path)) {
  stop("cv_binomial_all_comparisons.csv not found — run 08e_cv_validation.R first.")
}
cv_all <- readr::read_csv(cv_csv_path, show_col_types = FALSE) %>%
  dplyr::mutate(suffix = sub("^cv_", "", comparison))

message("CV rows loaded: ", nrow(cv_all))

# ── Join and reshape to long (phase) format ────────────────────
message("\n=== 08m | Joining in-sample + CV ===")

joined <- insample_df %>%
  dplyr::inner_join(
    cv_all %>% dplyr::select(suffix, method, auc_cv = cv_auc_mean, cv_auc_sd),
    by = c("suffix", "method")
  ) %>%
  dplyr::inner_join(COMPARISON_ORDER, by = "suffix") %>%
  dplyr::mutate(
    display_label = factor(display_label, levels = rev(COMPARISON_ORDER$display_label)),
    method         = factor(method, levels = c("Elastic Net", "sPLS-DA", "Boruta+RF")),
    overfit_gap    = auc_insample - auc_cv
  )

missing_pairs <- nrow(COMPARISON_ORDER) * 3 - nrow(joined)
if (missing_pairs != 0) {
  warning(missing_pairs, " expected comparison x method pairs did not join — check suffix naming consistency.")
}

safe_write_csv(
  joined %>% dplyr::select(display_label, method, auc_insample, auc_cv, cv_auc_sd, overfit_gap),
  file.path(OUT_TBL, "insample_vs_cv_binomial_summary.csv")
)

message("In-sample vs. CV summary (largest overfitting gaps first):")
print(joined %>%
        dplyr::arrange(dplyr::desc(overfit_gap)) %>%
        dplyr::select(display_label, method, auc_insample, auc_cv, overfit_gap) %>%
        head(10))

# ── Dumbbell plot ────────────────────────────────────────────────
message("\n=== 08m | Building in-sample vs CV dumbbell plot ===")

long_df <- dplyr::bind_rows(
  joined %>% dplyr::transmute(display_label, method, phase = "In-sample", auc = auc_insample),
  joined %>% dplyr::transmute(display_label, method, phase = "Cross-validated", auc = auc_cv)
) %>%
  dplyr::mutate(phase = factor(phase, levels = c("In-sample", "Cross-validated")))

p_dumbbell <- ggplot() +
  geom_segment(
    data = joined,
    aes(x = display_label, xend = display_label, y = auc_insample, yend = auc_cv),
    colour = "grey60", linewidth = 0.6
  ) +
  geom_point(
    data = long_df,
    aes(x = display_label, y = auc, colour = phase), size = 2.6
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  coord_flip(ylim = c(0, 1)) +
  facet_wrap(~ method, ncol = 3) +
  scale_colour_manual(values = c("In-sample" = "#c9772b", "Cross-validated" = "#2f6fb0"), name = NULL) +
  labs(
    title    = "Figure 6(a): In-Sample vs. Cross-Validated AUC, All Binomial Comparisons",
    subtitle = "Segment length = overfitting magnitude per method x comparison; dashed line = chance (AUC 0.5)",
    x = NULL, y = "AUC"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(size = 10.5, colour = "grey30"),
    strip.text        = element_text(face = "bold"),
    legend.position   = "top"
  )

stub_dumbbell <- file.path(OUT_FIG, "insample_vs_cv_binomial_dumbbell")

cairo_pdf(paste0(stub_dumbbell, ".pdf"), width = 11, height = 6.5)
print(p_dumbbell)
dev.off()

png(paste0(stub_dumbbell, ".png"), width = 11, height = 6.5, units = "in", res = 300, bg = "white")
print(p_dumbbell)
dev.off()

message("Dumbbell plot saved: ", stub_dumbbell, ".pdf / .png")
message("\n08m_insample_vs_cv_binomial_dumbbell.R complete")
