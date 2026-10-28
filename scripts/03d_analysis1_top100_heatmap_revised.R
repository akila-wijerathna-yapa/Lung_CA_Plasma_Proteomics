sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")

while (dev.cur() > 1) dev.off()

# ── Local constants (match 03_analysis1_seven_group.R) ────────
TOP_N_HEATMAP <- 100

# ── Local function (same as 03_analysis1_seven_group.R) ───────
scale_rows_na <- function(mat) {
  scaled <- t(apply(mat, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    mu <- mean(v, na.rm = TRUE); sdv <- sd(v, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) return(rep(0, length(v)))
    (v - mu) / sdv
  }))
  dimnames(scaled) <- dimnames(mat)
  scaled
}

# ── Rebuild top-100 selection (identical to original script) ──
message("\n=== 03d | Rebuilding top-100 omnibus selection ===")

if (!exists("analysis1_omnibus") || nrow(analysis1_omnibus) == 0) {
  stop("analysis1_omnibus not found — confirm checkpoint_03.RData is up to date.")
}

top100_ids <- analysis1_omnibus %>%
  dplyr::arrange(adj_p_value, dplyr::desc(f_stat)) %>%
  dplyr::slice_head(n = min(TOP_N_HEATMAP, nrow(.))) %>%
  dplyr::pull(nELISA_ID)

message("Top-", length(top100_ids), " proteins selected")

if (length(top100_ids) == 0) stop("No proteins selected for heatmap.")

# ── Build gene-symbol row labels ────────────────────────────────
message("\n=== 03d | Mapping nELISA_ID -> Gene_Symbol for row labels ===")

id_to_label <- feature_meta %>%
  dplyr::distinct(nELISA_ID, .keep_all = TRUE) %>%
  dplyr::mutate(
    row_label = dplyr::case_when(
      !is.na(Gene_Symbol)   & Gene_Symbol   != "" ~ Gene_Symbol,
      !is.na(Protein_Label) & Protein_Label != "" ~ Protein_Label,
      TRUE ~ nELISA_ID
    )
  ) %>%
  dplyr::select(nELISA_ID, row_label) %>%
  tibble::deframe()

row_labels <- unname(id_to_label[top100_ids])
row_labels[is.na(row_labels)] <- top100_ids[is.na(row_labels)]   # fallback if any ID is unmapped
row_labels <- make.unique(row_labels)                             # guard against duplicate symbols

n_unmapped <- sum(row_labels == top100_ids & !top100_ids %in% id_to_label)
if (n_unmapped > 0) message("Note: ", n_unmapped, " protein(s) had no Gene_Symbol/Protein_Label — nELISA_ID shown instead.")

# ── Column annotation, plate_barcode removed ───────────────────
message("\n=== 03d | Building column annotation (plate_barcode removed) ===")

heat_ann <- analysis1_df %>%
  dplyr::select(sample_id, new_group, response_binary) %>%
  tibble::column_to_rownames("sample_id")

# ── Build scaled matrix with gene-symbol row names ──────────────
heat_mat           <- t(expr_log2_rep[analysis1_df$sample_id, top100_ids, drop = FALSE])
rownames(heat_mat) <- row_labels
heat_scaled        <- scale_rows_na(heat_mat)

# ── Save heatmap (same style/palette as original) ───────────────
message("\n=== 03d | Saving revised heatmap ===")

stub_heat <- file.path(paths$figures, "02_analysis1_seven_group", "analysis1_top100_heatmap")

cairo_pdf(paste0(stub_heat, ".pdf"), width = 12, height = 14)
ph_a1 <- pheatmap::pheatmap(
  heat_scaled,
  annotation_col = heat_ann,
  show_colnames  = FALSE,
  show_rownames  = TRUE,
  fontsize_row   = 6,
  color          = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  na_col         = "grey90",
  main           = paste0("Analysis 1: Top ", length(top100_ids), " omnibus proteins"),
  silent         = TRUE
)
grid::grid.newpage(); grid::grid.draw(ph_a1$gtable)
dev.off()

png(paste0(stub_heat, ".png"), width = 12, height = 14, units = "in", res = 300, bg = "white")
ph_a1 <- pheatmap::pheatmap(
  heat_scaled,
  annotation_col = heat_ann,
  show_colnames  = FALSE,
  show_rownames  = TRUE,
  fontsize_row   = 6,
  color          = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  na_col         = "grey90",
  main           = paste0("Analysis 1: Top ", length(top100_ids), " omnibus proteins"),
  silent         = TRUE
)
grid::grid.newpage(); grid::grid.draw(ph_a1$gtable)
dev.off()

message("Revised heatmap saved: ", stub_heat, ".pdf / .png")
message("\n03d_analysis1_top100_heatmap_revised.R complete")
