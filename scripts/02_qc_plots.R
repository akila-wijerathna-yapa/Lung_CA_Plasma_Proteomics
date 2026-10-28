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

# ── Load checkpoint ───────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")

# Close any stale graphics devices left from previous runs
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
MIN_PCA_FEATURES  <- 50    # minimum complete features required for PCA
MIN_BATCH_FEATURES <- 100  # minimum complete features for batch-correction plot

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

scale_rows_na <- function(mat) {
  # z-score each row; rows with zero variance become all-zero
  scaled <- t(apply(mat, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    mu  <- mean(v, na.rm = TRUE)
    sig <- sd(v,   na.rm = TRUE)
    if (!is.finite(sig) || sig == 0) return(rep(0, length(v)))
    (v - mu) / sig
  }))
  rownames(scaled) <- rownames(mat)
  colnames(scaled) <- colnames(mat)
  scaled
}

# Colorblind-friendly palette for the 7 treatment groups
GROUP_COLOURS <- c(
  "Sx"        = "#E69F00",
  "CTx"       = "#D55E00",
  "CTx + IO"  = "#009E73",
  "CTx + RTx" = "#0072B2",
  "IO"        = "#56B4E9",
  "TKI"       = "#CC79A7",
  "Other"     = "#999999",
  "Unmatched" = "#DDDDDD"
)

# ── Annotate pre_all with clinical columns ────────────────────
# pre_all from checkpoint contains all representative PRE samples,
# including those not matched to clinical. Add clinical columns here
# for QC stratification.
pre_all_qc <- pre_all %>%
  mutate(
    matched_clinical = ifelse(id_number %in% clinical$patient_id, "Matched", "Unmatched")
  ) %>%
  left_join(
    clinical %>% select(patient_id, new_group, response_4cat, response_binary),
    by = c("id_number" = "patient_id")
  ) %>%
  mutate(
    new_group_plot = dplyr::coalesce(new_group, "Unmatched")
  )

# ── Sample QC summary table ───────────────────────────────────
message("\n=== 02 | Sample QC summary ===")
sample_qc_summary <- pre_all_qc %>%
  mutate(
    sample_missing_rate    = rowMeans(is.na(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE])),
    median_log2_intensity  = apply(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE], 1, median, na.rm = TRUE),
    total_log2_intensity   = rowSums(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE], na.rm = TRUE)
  )
safe_write_csv(sample_qc_summary, file.path(paths$tables, "01_qc", "sample_qc_summary_pre.csv"))

# ── Protein CV table ──────────────────────────────────────────
protein_cv <- tibble(
  nELISA_ID    = colnames(expr_pgml_rep),
  cv           = apply(expr_pgml_rep[pre_all_qc$sample_id, , drop = FALSE], 2, function(x) {
    mu  <- mean(x, na.rm = TRUE)
    sig <- sd(x,   na.rm = TRUE)
    if (!is.finite(mu) || !is.finite(sig) || mu <= 0) return(NA_real_)
    sig / mu
  }),
  missing_rate = colMeans(is.na(expr_pgml_rep[pre_all_qc$sample_id, , drop = FALSE]))
) %>%
  left_join(feature_meta %>% select(nELISA_ID, Gene_Symbol, Protein_Label), by = "nELISA_ID")
safe_write_csv(protein_cv, file.path(paths$tables, "01_qc", "protein_cv_distribution.csv"))

# ── Plot 1: PCA ───────────────────────────────────────────────
message("\n=== 02 | Plot 1: PCA ===")
complete_pca_features <- names(which(
  colSums(is.na(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE])) == 0
))

if (length(complete_pca_features) >= MIN_PCA_FEATURES) {
  pca_obj <- prcomp(
    expr_log2_rep[pre_all_qc$sample_id, complete_pca_features, drop = FALSE],
    center = TRUE, scale. = TRUE
  )
  pca_df <- as_tibble(pca_obj$x[, 1:2], rownames = "sample_id") %>%
    left_join(pre_all_qc %>% select(sample_id, new_group_plot, plate_barcode), by = "sample_id")
  safe_write_csv(pca_df, file.path(paths$tables, "01_qc", "pca_pre_scores.csv"))

  pct1 <- round(100 * summary(pca_obj)$importance[2, 1], 1)
  pct2 <- round(100 * summary(pca_obj)$importance[2, 2], 1)

  p1 <- ggplot(pca_df, aes(PC1, PC2, colour = new_group_plot, shape = plate_barcode)) +
    geom_point(size = 3, alpha = 0.88) +
    scale_colour_manual(values = GROUP_COLOURS, na.value = "#DDDDDD") +
    labs(
      title    = paste0("QC PCA — baseline PRE after QC/dedup (n=", nrow(pca_df),
                        ", proteins=", length(complete_pca_features), ")"),
      x        = paste0("PC1 (", pct1, "%)"),
      y        = paste0("PC2 (", pct2, "%)"),
      colour   = "Treatment group",
      shape    = "Plate"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  stub <- file.path(paths$figures, "01_qc", "01_pca_new_group_by_plate")
  cairo_pdf(paste0(stub, ".pdf"), width = 8.2, height = 6)
  print(p1)
  dev.off()
  png(paste0(stub, ".png"), width = 8.2, height = 6, units = "in", res = 300, bg = "white")
  print(p1)
  dev.off()
  message("Plot 1 saved: PCA")
} else {
  message("Plot 1 skipped: fewer than ", MIN_PCA_FEATURES, " complete features for PCA")
}

# ── Plot 2: Sample correlation heatmap ────────────────────────
message("\n=== 02 | Plot 2: Sample correlation heatmap ===")
op   <- options(warn = -1)
corr_mat <- cor(
  t(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE]),
  use = "pairwise.complete.obs", method = "spearman"
)
options(op)

ann_corr <- pre_all_qc %>%
  transmute(sample_id, matched_clinical, plate_barcode) %>%
  tibble::column_to_rownames("sample_id")

safe_write_csv(
  as.data.frame(corr_mat) %>% tibble::rownames_to_column("sample_id"),
  file.path(paths$tables, "01_qc", "sample_correlation_matrix_pre.csv")
)

stub2 <- file.path(paths$figures, "01_qc", "02_sample_correlation_heatmap")
cairo_pdf(paste0(stub2, ".pdf"), width = 12, height = 11)
ph2 <- pheatmap::pheatmap(
  corr_mat,
  annotation_col    = ann_corr,
  annotation_row    = ann_corr,
  color             = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  clustering_method = "complete",
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  main              = paste0("Baseline PRE sample correlation (n=", nrow(corr_mat), ")"),
  silent            = TRUE
)
grid::grid.newpage(); grid::grid.draw(ph2$gtable)
dev.off()
png(paste0(stub2, ".png"), width = 12, height = 11, units = "in", res = 300, bg = "white")
ph2 <- pheatmap::pheatmap(
  corr_mat,
  annotation_col    = ann_corr,
  annotation_row    = ann_corr,
  color             = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  clustering_method = "complete",
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  main              = paste0("Baseline PRE sample correlation (n=", nrow(corr_mat), ")"),
  silent            = TRUE
)
grid::grid.newpage(); grid::grid.draw(ph2$gtable)
dev.off()
message("Plot 2 saved: sample correlation heatmap")

# ── Plot 3: Protein CV distribution ──────────────────────────
message("\n=== 02 | Plot 3: Protein CV distribution ===")
p3 <- ggplot(protein_cv %>% filter(is.finite(cv)), aes(cv)) +
  geom_histogram(bins = 50, fill = "#56B4E9", colour = "white") +
  labs(
    title = "Protein coefficient of variation — baseline PRE samples",
    x     = "CV (SD / mean in pg/mL)",
    y     = "Number of proteins"
  ) +
  theme_bw(base_size = 11)

stub3 <- file.path(paths$figures, "01_qc", "03_protein_cv_distribution")
cairo_pdf(paste0(stub3, ".pdf"), width = 7.5, height = 5.2)
print(p3)
dev.off()
png(paste0(stub3, ".png"), width = 7.5, height = 5.2, units = "in", res = 300, bg = "white")
print(p3)
dev.off()
message("Plot 3 saved: protein CV histogram")

# ── Plot 4: Missing-value heatmap ─────────────────────────────
message("\n=== 02 | Plot 4: Missing-value heatmap ===")
miss_sample_order   <- pre_all_qc$sample_id[
  order(rowMeans(is.na(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE])), decreasing = TRUE)
]
miss_protein_order  <- colnames(expr_log2_rep)[
  order(colMeans(is.na(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE])), decreasing = TRUE)
]
missing_mat <- t(is.na(expr_log2_rep[miss_sample_order, miss_protein_order, drop = FALSE]) * 1L)

stub4 <- file.path(paths$figures, "01_qc", "04_missing_value_heatmap")
cairo_pdf(paste0(stub4, ".pdf"), width = 13, height = 8.5)
ph4 <- pheatmap::pheatmap(
  missing_mat,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  color         = c("#f7fbff", "#D55E00"),
  legend_breaks = c(0, 1),
  legend_labels = c("Observed", "Missing"),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main          = "Missing-value map (samples \u00d7 proteins ordered by missing rate)",
  silent        = TRUE
)
grid::grid.newpage(); grid::grid.draw(ph4$gtable)
dev.off()
png(paste0(stub4, ".png"), width = 13, height = 8.5, units = "in", res = 300, bg = "white")
ph4 <- pheatmap::pheatmap(
  missing_mat,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  color         = c("#f7fbff", "#D55E00"),
  legend_breaks = c(0, 1),
  legend_labels = c("Observed", "Missing"),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main          = "Missing-value map (samples \u00d7 proteins ordered by missing rate)",
  silent        = TRUE
)
grid::grid.newpage(); grid::grid.draw(ph4$gtable)
dev.off()
message("Plot 4 saved: missing-value heatmap")

# ── Plot 5: Per-sample boxplots ───────────────────────────────
message("\n=== 02 | Plot 5: Per-sample boxplots ===")

sample_id_order <- sample_qc_summary$sample_id[order(sample_qc_summary$total_log2_intensity)]

# Matrix: proteins × samples, columns ordered by total intensity
mat5 <- t(expr_log2_rep[sample_id_order, , drop = FALSE])

plot_boxplot5 <- function() {
  boxplot(mat5,
          outline  = FALSE,
          col      = "#009E7366",
          border   = "#004D40",
          xaxt     = "n",
          main     = "Per-sample log2 concentration distributions (ordered by total intensity)",
          xlab     = "Samples (ordered by total intensity)",
          ylab     = "log2 concentration (pg/mL)",
          cex.main = 0.85,
          las      = 1)
}

stub5 <- file.path(paths$figures, "01_qc", "05_total_intensity_boxplots")

cairo_pdf(paste0(stub5, ".pdf"), width = 12, height = 6.2)
plot_boxplot5()
dev.off()

png(paste0(stub5, ".png"), width = 12, height = 6.2, units = "in", res = 300, bg = "white")
plot_boxplot5()
dev.off()
message("Plot 5 saved: per-sample boxplots")

# ── Plot 6: Batch-correction density (QC visual only) ─────────
message("\n=== 02 | Plot 6: Batch-correction density ===")
batch_features <- names(which(
  colSums(is.na(expr_log2_rep[pre_all_qc$sample_id, , drop = FALSE])) == 0
))

if (length(batch_features) >= MIN_BATCH_FEATURES) {
  before_df <- as_tibble(
    expr_log2_rep[pre_all_qc$sample_id, batch_features, drop = FALSE],
    rownames = "sample_id"
  ) %>%
    tidyr::pivot_longer(-sample_id, names_to = "nELISA_ID", values_to = "log2_value") %>%
    left_join(pre_all_qc %>% select(sample_id, plate_barcode), by = "sample_id") %>%
    mutate(stage = "Before correction")

  corrected <- limma::removeBatchEffect(
    t(expr_log2_rep[pre_all_qc$sample_id, batch_features, drop = FALSE]),
    batch = pre_all_qc$plate_barcode
  )

  after_df <- as_tibble(t(corrected), rownames = "sample_id") %>%
    tidyr::pivot_longer(-sample_id, names_to = "nELISA_ID", values_to = "log2_value") %>%
    left_join(pre_all_qc %>% select(sample_id, plate_barcode), by = "sample_id") %>%
    mutate(stage = "After correction")

  density_df <- bind_rows(before_df, after_df) %>%
    mutate(stage = factor(stage, levels = c("Before correction", "After correction")))

  p6 <- ggplot(
    density_df,
    aes(log2_value, colour = plate_barcode, group = interaction(plate_barcode, sample_id))
  ) +
    geom_density(alpha = 0.15, linewidth = 0.3) +
    facet_wrap(~stage, ncol = 1) +
    labs(
      title  = "Density distributions before and after plate correction (QC visual only)",
      x      = "log2 concentration (pg/mL)",
      y      = "Density",
      colour = "Plate"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  stub6 <- file.path(paths$figures, "01_qc", "06_density_before_after_batch_correction")
  cairo_pdf(paste0(stub6, ".pdf"), width = 8.2, height = 7.2)
  print(p6)
  dev.off()
  png(paste0(stub6, ".png"), width = 8.2, height = 7.2, units = "in", res = 300, bg = "white")
  print(p6)
  dev.off()
  message("Plot 6 saved: batch-correction density")
} else {
  message("Plot 6 skipped: fewer than ", MIN_BATCH_FEATURES, " complete features")
}

message("\n02_qc_plots.R complete")
