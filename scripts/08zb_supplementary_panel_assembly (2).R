if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
source(file.path(PROJECT_DIR, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(cowplot)
  library(magick)
  library(ggplot2)
})

ML_FIG   <- file.path(paths$figures, "07_analysis6_machine_learning")
CV_FIG   <- file.path(ML_FIG, "cv")
OUT_DIR  <- file.path(paths$figures, "08_supplementary_panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Helper: safe PNG -> ggdraw panel ────────────────────────────────────────
safe_panel <- function(path, label = NULL) {
  if (!file.exists(path)) {
    warning("Missing source figure: ", path)
    return(ggdraw() + draw_label(paste("MISSING:\n", basename(path)),
                                  size = 10, color = "red"))
  }
  p <- ggdraw() + draw_image(image_read(path))
  if (!is.null(label)) p <- p + draw_plot_label(label, size = 14, fontface = "bold")
  p
}

save_panel <- function(plot_obj, stub, width, height) {
  cairo_pdf(paste0(stub, ".pdf"), width = width, height = height)
  print(plot_obj); dev.off()
  png(paste0(stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white")
  print(plot_obj); dev.off()
  message("Saved: ", stub, ".pdf/.png")
}

# ── Sup Fig S1 — OS (GDF15 cross-method coherence) ─────────────────────────
# Panel B now sources 08zd's output (GDF15/GZMA, Cooper-integrated consensus)
# instead of 08f's stale consensus_OS_facet_grid.png (GDF11, pre-Cooper list).
p1 <- safe_panel(file.path(ML_FIG, "gdf15_coherence_check.png"), "A")
p2 <- safe_panel(file.path(ML_FIG, "consensus_km", "gdf15_gzma_OS_facet_grid.png"), "B")

s1 <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 0.7))
save_panel(s1, file.path(OUT_DIR, "Sup_Fig_S1_OS"), width = 10, height = 11)

# ── Sup Fig S2 — Recurrence (CEACAM5 robustness) ───────────────────────────
# Panel B: recurrence_binary_sensitivity endpoint, TKI-sensitivity (08zc).
# Panel C (new): stage-stratified/stage-adjusted sensitivity check (08zf) —
# CEACAM5 recurrence OR attenuates from 1.79 (unadjusted) to 1.66 (adjusted
# for stage) but remains significant (p=0.017), addressing the stage
# confound flagged by 08w's stage-progression trend test (CEACAM5 FDR=0.036).
p3 <- safe_panel(file.path(ML_FIG, "ceacam5_coherence_check.png"), "A")
p4 <- safe_panel(file.path(ML_FIG, "ceacam5_recurrence_tki_sensitivity_boxplot.png"), "B")
p5_s2 <- safe_panel(file.path(ML_FIG, "ceacam5_recurrence_by_stage_boxplot.png"), "C")

s2_top <- plot_grid(p3, p4, ncol = 2, rel_widths = c(1, 1))
s2 <- plot_grid(s2_top, p5_s2, ncol = 1, rel_heights = c(1, 0.9))
save_panel(s2, file.path(OUT_DIR, "Sup_Fig_S2_Recurrence"), width = 14, height = 11)

# ── Sup Fig S3 — Response (SHH phenotype-extreme dependency) ──────────────
p5 <- safe_panel(file.path(ML_FIG, "shh_pathway_coherence_check.png"), "A")
p6 <- safe_panel(file.path(CV_FIG, "cv_response_granularity_barplot.png"), "B")

s3 <- plot_grid(p5, p6, ncol = 1, rel_heights = c(1.3, 1))
save_panel(s3, file.path(OUT_DIR, "Sup_Fig_S3_Response"), width = 11, height = 10.5)

message("\n=== 08za complete: 3 supplementary panels written to ", OUT_DIR, " ===")
