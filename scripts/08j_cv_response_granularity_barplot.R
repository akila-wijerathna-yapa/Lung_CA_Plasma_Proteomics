sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

while (dev.cur() > 1) dev.off()   # clear any stale devices

# ── Local constants ───────────────────────────────────────────
CV_TBL_DIR <- file.path(paths$tables,  "07_analysis6_machine_learning", "cv")
CV_FIG_DIR <- file.path(paths$figures, "07_analysis6_machine_learning", "cv")
dir.create(CV_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

CHANCE_LINE <- 0.5

# Response comparisons only, ordered by phenotype granularity:
# most extreme / cleanest contrast first, most graded / clinically
# realistic contrast last. This ordering is a judgment call about
# "distance between compared phenotypes", not a re-sort by AUC —
# the point of the figure is to show performance *given* this a
# priori ordering, not to cherry-pick an order that decays neatly.
RESPONSE_ORDER <- tibble::tribble(
  ~comparison,               ~contrast_label,                 ~granularity_rank,
  "cv_01_CR_vs_PD",          "CR vs PD\n(extreme)",                        1,
  "cv_05_PR_vs_PD",          "PR vs PD",                                   2,
  "cv_02_CRPR_vs_PD",        "CR/PR vs PD",                                3,
  "cv_04_CRPRSD_vs_PD",      "CR/PR/SD vs PD",                             4,
  "cv_03_CRPR_vs_SDPD",      "CR/PR vs SD/PD\n(clinical definition)",      5
)

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

akila_volcano_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 15),
      plot.subtitle     = element_text(size = 11, colour = "grey30"),
      axis.title        = element_text(face = "bold"),
      legend.position   = "top",
      legend.title      = element_blank()
    )
}

# ── Load CV binomial results (from CSV, not RData checkpoint) ───
# Read directly from the CSV 08e already wrote rather than chaining
# through checkpoint_08e.RData — avoids the checkpoint-dependency
# fragility that caused the earlier "analysis1_df not found" failure.
message("\n=== 08j | Loading CV binomial results ===")

cv_csv_path <- file.path(CV_TBL_DIR, "cv_binomial_all_comparisons.csv")
if (!file.exists(cv_csv_path)) {
  stop("cv_binomial_all_comparisons.csv not found at: ", cv_csv_path,
       "\nRun 08e_cv_validation.R (with the reverse=TRUE and checkpoint_03 fixes) first.")
}

cv_binomial_all <- readr::read_csv(cv_csv_path, show_col_types = FALSE)

message("Loaded ", nrow(cv_binomial_all), " rows across ",
        dplyr::n_distinct(cv_binomial_all$comparison), " comparisons")

# ── Filter to response comparisons and apply granularity order ──
message("\n=== 08j | Filtering to response comparisons ===")

plot_df <- cv_binomial_all %>%
  dplyr::inner_join(RESPONSE_ORDER, by = "comparison")

missing_comparisons <- setdiff(RESPONSE_ORDER$comparison, plot_df$comparison)
if (length(missing_comparisons) > 0) {
  warning("Expected response comparisons not found in cv_binomial_all_comparisons.csv: ",
          paste(missing_comparisons, collapse = ", "),
          " — check comparison names match those written by 08e_cv_validation.R")
}

plot_df <- plot_df %>%
  dplyr::mutate(
    contrast_label = factor(contrast_label,
                            levels = RESPONSE_ORDER$contrast_label[order(RESPONSE_ORDER$granularity_rank)]),
    method = factor(method, levels = c("Elastic Net", "sPLS-DA", "Boruta+RF"))
  ) %>%
  dplyr::arrange(granularity_rank, method)

safe_write_csv(plot_df, file.path(CV_TBL_DIR, "cv_response_granularity_ordered.csv"))

message("Response-only CV summary (ordered by phenotype granularity):")
print(plot_df %>% dplyr::select(contrast_label, method, cv_auc_mean, cv_auc_sd))

# ── Bar plot ───────────────────────────────────────────────────
message("\n=== 08j | Building CV AUC-by-granularity bar plot ===")

p_gran <- ggplot(plot_df, aes(x = contrast_label, y = cv_auc_mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68,
           colour = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = pmax(0, cv_auc_mean - cv_auc_sd), ymax = pmin(1, cv_auc_mean + cv_auc_sd)),
    position = position_dodge(width = 0.75), width = 0.18, linewidth = 0.4
  ) +
  geom_hline(yintercept = CHANCE_LINE, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  annotate("text", x = 0.6, y = CHANCE_LINE, label = "chance", vjust = -0.6,
           hjust = 0, size = 3.3, colour = "grey40", fontface = "italic") +
  scale_fill_manual(values = c("Elastic Net" = "#2f6fb0", "sPLS-DA" = "#2f9e5b",
                                "Boruta+RF" = "#c9772b")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "Response Classification: CV AUC by Phenotype Granularity",
    subtitle = "5-fold x 5-repeat cross-validated AUC, ordered from extreme to graded/clinical contrast",
    x = NULL,
    y = "Cross-validated AUC"
  ) +
  akila_volcano_theme() +
  theme(axis.text.x = element_text(size = 10))

stub_gran <- file.path(CV_FIG_DIR, "cv_response_granularity_barplot")

cairo_pdf(paste0(stub_gran, ".pdf"), width = 9, height = 5.8)
print(p_gran)
dev.off()

png(paste0(stub_gran, ".png"), width = 9, height = 5.8, units = "in", res = 300, bg = "white")
print(p_gran)
dev.off()

message("CV granularity bar plot saved: ", stub_gran, ".pdf / .png")
message("\n08j_cv_response_granularity_barplot.R complete")
