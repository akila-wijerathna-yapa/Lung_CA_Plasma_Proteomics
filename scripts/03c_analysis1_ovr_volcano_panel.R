sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
suppressPackageStartupMessages({ library(ggrepel); library(patchwork) })

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")

while (dev.cur() > 1) dev.off()

# ── Local constants (match 03_analysis1_seven_group.R) ────────
FDR_CUTOFF   <- 0.1
LFC_THR      <- 0.5
LABELS_PER_PANEL <- 10   # fewer than the single-plot version (18) — panels are smaller

if (!exists("analysis1_ovr_results") || length(analysis1_ovr_results) == 0) {
  stop("analysis1_ovr_results not found — confirm checkpoint_03.RData is up to date.")
}

# ── Local function: single mini-volcano panel ─────────────────
# Same status/sig_level logic as make_volcano() in 03_analysis1_seven_group.R,
# simplified for compact multi-panel display (fewer labels, smaller theme).
make_mini_volcano <- function(df, group_label,
                               lfc_thr = LFC_THR, padj_thr = FDR_CUTOFF,
                               label_col = "Protein_Label",
                               n_labels = LABELS_PER_PANEL) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  plot_df <- df %>%
    dplyr::mutate(
      neg_log10 = -log10(pmax(adj_p_value, 1e-300)),
      status = dplyr::case_when(
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC >=  lfc_thr ~ paste0("Higher in ", group_label),
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC <= -lfc_thr ~ "Higher in rest",
        TRUE ~ "Not significant"
      )
    )

  n_sig <- sum(plot_df$status != "Not significant")

  ann <- plot_df %>%
    dplyr::filter(status != "Not significant") %>%
    dplyr::arrange(adj_p_value, dplyr::desc(abs(logFC))) %>%
    dplyr::slice_head(n = n_labels)

  label_vals <- if (label_col %in% names(ann)) ann[[label_col]] else ann$nELISA_ID

  ggplot(plot_df, aes(x = logFC, y = neg_log10)) +
    geom_point(aes(colour = status), size = 1.1, alpha = 0.75) +
    { if (nrow(ann) > 0) ggrepel::geom_text_repel(
        data = ann %>% dplyr::mutate(.lbl = label_vals),
        aes(label = .lbl), size = 2.6, max.overlaps = 30,
        segment.size = 0.2, min.segment.length = 0
      ) } +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
    scale_colour_manual(values = stats::setNames(
      c("grey75", "#d55e00", "#0072b2"),
      c("Not significant", paste0("Higher in ", group_label), "Higher in rest")
    ), guide = "none") +
    labs(
      title = paste0(group_label, " vs rest  (n sig. = ", n_sig, ")"),
      x = paste0("log2FC (", group_label, " / rest)"),
      y = expression(-log[10](adj~p))
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 10),
      axis.title       = element_text(size = 8.5)
    )
}

# ── Build one mini-volcano per treatment group ─────────────────
message("\n=== 03c | Building 7-panel volcano grid ===")

panels <- lapply(names(analysis1_ovr_results), function(grp) {
  make_mini_volcano(analysis1_ovr_results[[grp]], group_label = grp)
})
panels <- Filter(Negate(is.null), panels)

message("Panels built: ", length(panels), " of ", length(analysis1_ovr_results), " groups")

combined <- patchwork::wrap_plots(panels, ncol = 3) +
  patchwork::plot_annotation(
    title    = "Supplementary Figure: One-vs-Rest Differential Expression by Treatment Group",
    subtitle = paste0("Each panel: 577 proteins tested, FDR<", FDR_CUTOFF,
                       " dashed line, top ", LABELS_PER_PANEL, " labeled per panel"),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, colour = "grey30"))
  )

stub_panel <- file.path(paths$figures, "02_analysis1_seven_group", "analysis1_ovr_volcano_panel")

cairo_pdf(paste0(stub_panel, ".pdf"), width = 13, height = 9.5)
print(combined)
dev.off()

png(paste0(stub_panel, ".png"), width = 13, height = 9.5, units = "in", res = 300, bg = "white")
print(combined)
dev.off()

message("7-panel volcano grid saved: ", stub_panel, ".pdf / .png")
message("\n03c_analysis1_ovr_volcano_panel.R complete")
