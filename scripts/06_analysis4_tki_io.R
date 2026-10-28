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

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
TOP_BAR      <- 25     # proteins to show in ranked bar plot
TOP_HEATMAP  <- 50     # proteins for heatmap

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  left_join(df, feature_meta, by = "nELISA_ID")
}

scale_rows_na <- function(mat) {
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

save_ggplot <- function(p, path_stub, width, height) {
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height)
  print(p)
  dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white")
  print(p)
  dev.off()
}

# ── Analysis 4 data frame ─────────────────────────────────────
message("\n=== 06 | Analysis 4 TKI vs IO setup ===")
analysis4_df <- analysis1_df %>% filter(new_group %in% c("TKI", "IO"))
n_io  <- sum(analysis4_df$new_group == "IO")
n_tki <- sum(analysis4_df$new_group == "TKI")
message("IO: ", n_io, " | TKI: ", n_tki)

if (nrow(analysis4_df) == 0 || length(unique(analysis4_df$new_group)) < 2) {
  message("Analysis 4: insufficient samples — skipped")
  analysis4_desc <- tibble()
} else {

  mat4 <- expr_log2_rep[analysis4_df$sample_id, , drop = FALSE]

  # ── Descriptive logFC table ─────────────────────────────────
  analysis4_desc <- tibble(
    nELISA_ID = colnames(mat4),
    mean_IO   = colMeans(mat4[analysis4_df$new_group == "IO",  , drop = FALSE], na.rm = TRUE),
    mean_TKI  = colMeans(mat4[analysis4_df$new_group == "TKI", , drop = FALSE], na.rm = TRUE),
    logFC     = colMeans(mat4[analysis4_df$new_group == "TKI", , drop = FALSE], na.rm = TRUE) -
                colMeans(mat4[analysis4_df$new_group == "IO",  , drop = FALSE], na.rm = TRUE),
    n_IO      = colSums(!is.na(mat4[analysis4_df$new_group == "IO",  , drop = FALSE])),
    n_TKI     = colSums(!is.na(mat4[analysis4_df$new_group == "TKI", , drop = FALSE]))
  ) %>%
    add_feature_annotation(feature_meta) %>%
    mutate(analysis_type = "descriptive_logFC_only") %>%
    arrange(desc(abs(logFC)))

  safe_write_csv(analysis4_desc,
    file.path(paths$tables, "05_analysis4_tki_vs_io",
              "analysis4_tki_vs_io_logfc_only.csv"))
  message("Descriptive logFC table saved: ", nrow(analysis4_desc), " proteins")

  # ── Ranked bar: top absolute logFC ───────────────────────────
  bar_df <- analysis4_desc %>%
    slice_head(n = min(TOP_BAR, nrow(.))) %>%
    mutate(
      direction   = ifelse(logFC >= 0, "Higher TKI", "Higher IO"),
      label_order = factor(Protein_Label, levels = rev(Protein_Label))
    )
  p_bar <- ggplot(bar_df, aes(x = logFC, y = label_order, fill = direction)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("Higher TKI" = "#CC79A7", "Higher IO" = "#56B4E9")) +
    labs(
      title = paste0("Analysis 4: TKI vs IO — top ", nrow(bar_df),
                     " absolute logFC (n IO=", n_io, ", n TKI=", n_tki, ")"),
      subtitle = "Descriptive only — group sizes too small for formal testing",
      x = "Mean log2FC (TKI \u2212 IO)", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top",
          plot.subtitle = element_text(color = "grey50", face = "italic"))
  save_ggplot(p_bar,
    file.path(paths$figures, "05_analysis4_tki_vs_io",
              "analysis4_tki_vs_io_logfc_ranked_bar"),
    width = 8.5, height = 7)
  message("Ranked bar saved")

  # ── Heatmap: top absolute logFC proteins ─────────────────────
  heat_ids4 <- analysis4_desc %>%
    slice_head(n = min(TOP_HEATMAP, nrow(.))) %>%
    pull(nELISA_ID)
  heat_ann4 <- analysis4_df %>%
    select(sample_id, new_group, plate_barcode) %>%
    tibble::column_to_rownames("sample_id")
  heat_mat4   <- t(mat4[, heat_ids4, drop = FALSE])
  heat_scaled4 <- scale_rows_na(heat_mat4)

  stub_heat4 <- file.path(paths$figures, "05_analysis4_tki_vs_io",
                           "analysis4_tki_vs_io_heatmap")
  cairo_pdf(paste0(stub_heat4, ".pdf"), width = 10, height = 10.5)
  ph4 <- pheatmap::pheatmap(
    heat_scaled4,
    annotation_col = heat_ann4,
    show_colnames  = FALSE,
    color          = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    na_col         = "grey90",
    main           = paste0("Analysis 4: Top ", length(heat_ids4),
                            " absolute logFC proteins (TKI vs IO)"),
    silent         = TRUE
  )
  grid::grid.newpage(); grid::grid.draw(ph4$gtable)
  dev.off()
  png(paste0(stub_heat4, ".png"), width = 10, height = 10.5,
      units = "in", res = 300, bg = "white")
  ph4 <- pheatmap::pheatmap(
    heat_scaled4,
    annotation_col = heat_ann4,
    show_colnames  = FALSE,
    color          = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    na_col         = "grey90",
    main           = paste0("Analysis 4: Top ", length(heat_ids4),
                            " absolute logFC proteins (TKI vs IO)"),
    silent         = TRUE
  )
  grid::grid.newpage(); grid::grid.draw(ph4$gtable)
  dev.off()
  message("Heatmap saved")
}

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 06 | Saving checkpoint_06.RData ===")
save(
  analysis4_df,
  analysis4_desc,
  file = file.path(paths$rdata, "checkpoint_06.RData")
)
message("06_analysis4_tki_io.R complete")
