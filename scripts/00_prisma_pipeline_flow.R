library(ggplot2)
library(grid)

# ---- EDIT FOR YOUR MACHINE --------------------------------------------------
PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA"
OUTPUT_PDF  <- file.path(PROJECT_DIR, "figures", "00_overview", "00_prisma_pipeline_flow.pdf")
# -----------------------------------------------------------------------------

dir.create(dirname(OUTPUT_PDF), recursive = TRUE, showWarnings = FALSE)

PAGE_W <- 100
PAGE_H <- 350
TITLE_H <- 6
PAD <- 2.4
LINE_H <- 4.3

FILL_PALETTE <- c(
  source    = "#dbe9f6",
  qc        = "#d7eedb",
  analysis  = "#fde7c7",
  biomarker = "#e2dcf4",
  consensus = "#f8d9de",
  report    = "#ececec"
)

make_box <- function(x, y_top, w, h, fill, title, detail = NULL) {
  rect <- data.frame(xmin = x, xmax = x + w, ymin = y_top - h, ymax = y_top, fill_group = fill)
  ttl  <- data.frame(x = x + PAD, y = y_top - TITLE_H / 2, label = title)
  det  <- NULL
  if (!is.null(detail)) {
    det <- data.frame(
      x = x + PAD, y = y_top - TITLE_H - PAD - LINE_H / 2,
      label = paste0(paste0("- ", detail), collapse = "\n")
    )
  }
  list(rect = rect, title = ttl, detail = det)
}

boxes <- list(
  # ---- trunk: identification / screening -----------------------------------
  A     = make_box(12, 330, 40, 20, "source", "Clinical cohort identified",
                    "n = 173 patients (9 raw 1L_Treatment_Class categories)"),
  A_n   = make_box(56, 330, 44, 20, "source", "Source",
                    "Clinical_Akila_FINAL_cleaned_20260610.xlsx"),
  B     = make_box(12, 298, 40, 20, "qc", "New_Group collapse: 9 -> 7",
                    "n = 173 (unchanged)"),
  B_n   = make_box(56, 298, 44, 20, "qc", "Merged",
                    c("Classes 7-9 (n=6): Chemo+IO+RT, RT only,",
                      "RT+IO sequential -> 'Other'")),
  C     = make_box(12, 266, 40, 20, "qc", "Proteomics QC applied",
                    "n = 173 retained"),
  C_n   = make_box(56, 266, 44, 20, "qc", "Excluded / adjusted",
                    c("Excluded: sample_comments = 'qc_fail'",
                      "No imputation; LLOD/ULOD/bead-count -> NA")),
  D     = make_box(12, 234, 40, 24, "qc", "Baseline (pre) matched cohort: n = 173",
                    c("973 proteins across 3 nELISA assay plates",
                      "Group n: CTx 44, Sx 40, IO 28,",
                      "CTx+IO 26, CTx+RTx 15, TKI 13, Other 6")),
  D_n   = make_box(56, 234, 44, 24, "qc", "Model note",
                    c("Plate included as covariate in all",
                      "models (not batch-corrected)")),

  # ---- left branch: cross-sectional contrasts -------------------------------
  L1    = make_box(1, 192, 30, 30, "analysis", "A1 - Seven-group omnibus",
                    c("limma F-test, ~0+group+plate", "Min 3 samples/group/feature",
                      "FDR<0.10: 0 significant hits")),
  L2    = make_box(1, 154, 30, 28, "analysis", "A3/A4 - Response & TKI vs IO",
                    c("A3: PR/CR vs SD/PD, n=173 (prelim.)",
                      "A4: TKI n=13 vs IO n=28 (logFC only)")),

  # ---- middle branch: paired analysis ---------------------------------------
  M1    = make_box(35, 192, 30, 32, "analysis", "A2 - Paired pre/post (n=29)",
                    c("limma: ~pair_id + timepoint", "Surgery n=14 (>=5 pairs/feature)",
                      "IO-containing n=9: logFC only")),

  # ---- right branch: survival + biomarker discovery -------------------------
  R1    = make_box(69, 192, 30, 28, "analysis", "A5 - Survival (n=172)",
                    c("Cox: min n=18, min events=6", "BH-FDR; icenReg for interval-cens.")),
  R2    = make_box(69, 158, 30, 32, "biomarker", "A6 - ML discovery cohorts",
                    c("Comp. A: IO/CTx+IO resp., n=57", "Comp. B: Chemo+other resp., n=88",
                      "Surgery recurrence, n=40")),
  R3    = make_box(69, 120, 30, 32, "biomarker", "LASSO + STABL + alt. ML",
                    c("Stability threshold >= 0.60", "STABL: FDP-controlled, P<0.05",
                      "CV-corrected AUC (was in-sample)")),
  R4    = make_box(69, 82, 30, 32, "biomarker", "Consensus scoring (0-10)",
                    c("FDR<0.05=2, Cox OS=2, PFS=2", "LASSO=2, Stability=1, AUC=1",
                      "Priority: score >= 6/10")),

  # ---- convergence -----------------------------------------------------------
  CONV  = make_box(15, 34, 70, 32, "consensus", "09 - Pathway/network + Reporting (10-11)",
                    c("fgsea GO/KEGG/Reactome (FDR<0.05) + STRING PPI",
                      "Consensus: SHH (5/5), CEACAM5 (4/5) -> validation list",
                      "Poster (3 panels) + self-contained HTML/PDF report"))
)

rect_df   <- do.call(rbind, lapply(boxes, `[[`, "rect"))
title_df  <- do.call(rbind, lapply(boxes, `[[`, "title"))
detail_df <- do.call(rbind, Filter(Negate(is.null), lapply(boxes, `[[`, "detail")))

arrows <- rbind(
  data.frame(x = 32, xend = 32, y = 310, yend = 298),           # A -> B
  data.frame(x = 32, xend = 32, y = 278, yend = 266),           # B -> C
  data.frame(x = 32, xend = 32, y = 246, yend = 234),           # C -> D
  data.frame(x = 32, xend = 16, y = 210, yend = 192),           # D -> L1
  data.frame(x = 32, xend = 50, y = 210, yend = 192),           # D -> M1
  data.frame(x = 32, xend = 84, y = 210, yend = 192),           # D -> R1
  data.frame(x = 16, xend = 16, y = 162, yend = 154),           # L1 -> L2
  data.frame(x = 84, xend = 84, y = 164, yend = 158),           # R1 -> R2
  data.frame(x = 84, xend = 84, y = 126, yend = 120),           # R2 -> R3
  data.frame(x = 84, xend = 84, y = 88,  yend = 82),            # R3 -> R4
  data.frame(x = 16, xend = 16, y = 126, yend = 34),            # L2 -> convergence
  data.frame(x = 50, xend = 50, y = 160, yend = 34),            # M1 -> convergence
  data.frame(x = 84, xend = 84, y = 50,  yend = 34)             # R4 -> convergence
)

connectors <- rbind(  # trunk box -> its side note
  data.frame(x = 52, xend = 56, y = 320, yend = 320),
  data.frame(x = 52, xend = 56, y = 288, yend = 288),
  data.frame(x = 52, xend = 56, y = 256, yend = 256),
  data.frame(x = 52, xend = 56, y = 222, yend = 222)
)

leg_labels <- c("Data source", "QC / cohort", "Analysis", "Biomarker discovery", "Consensus / reporting")
leg_fill   <- c("source", "qc", "analysis", "biomarker", "consensus")
leg_y      <- 345 - seq_along(leg_labels) * 4.2
leg_rect   <- data.frame(xmin = 68, xmax = 71, ymin = leg_y - 1.3, ymax = leg_y + 1.3, f = leg_fill)
leg_text   <- data.frame(x = 73, y = leg_y, label = leg_labels)

p <- ggplot() +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_group),
            color = "grey30", linewidth = 0.4) +
  geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.11, "inches")), color = "grey30", linewidth = 0.5) +
  geom_segment(data = connectors, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.08, "inches")), color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_text(data = title_df, aes(x = x, y = y, label = label),
            hjust = 0, vjust = 0.5, fontface = "bold", size = 3.5, family = "sans") +
  geom_text(data = detail_df, aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 2.75, lineheight = 1.1, family = "sans") +
  geom_rect(data = leg_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = f),
            color = "grey40", linewidth = 0.3) +
  geom_text(data = leg_text, aes(x = x, y = y, label = label), hjust = 0, size = 2.6, color = "grey30") +
  scale_fill_manual(values = FILL_PALETTE, guide = "none") +
  annotate("text", x = 12, y = 345, label = "nELISA NSCLC study & analysis flow (PRISMA-style)",
           hjust = 0, size = 5, fontface = "bold") +
  coord_cartesian(xlim = c(0, PAGE_W), ylim = c(0, PAGE_H), clip = "off") +
  theme_void()

cairo_pdf(OUTPUT_PDF, width = 12, height = 17)
print(p)
dev.off()

message("PRISMA-style flow diagram written to: ", OUTPUT_PDF)
