library(ggplot2)
library(grid)

# ---- EDIT THESE TWO PATHS FOR YOUR MACHINE ---------------------------------
PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA"
OUTPUT_PDF  <- file.path(PROJECT_DIR, "figures", "00_overview", "00_pipeline_overview_with_cutoffs_and_outcomes.pdf")
# -----------------------------------------------------------------------------

dir.create(dirname(OUTPUT_PDF), recursive = TRUE, showWarnings = FALSE)

PAGE_UNIT_HEIGHT <- 240   # fixed coordinate budget per page (constant across pages)
PAGE_UNIT_WIDTH  <- 100
BOX_WIDTH        <- 86
X_LEFT           <- 7
TITLE_H          <- 7
PAD              <- 3
LINE_H           <- 4.6
GAP              <- 7
TOP_MARGIN       <- 16

FILL_PALETTE <- c(
  source    = "#dbe9f6",
  qc        = "#d7eedb",
  analysis  = "#fde7c7",
  biomarker = "#e2dcf4",
  consensus = "#f8d9de",
  report    = "#ececec"
)

# ---- generic box-stack builder ---------------------------------------------
build_flow_page <- function(boxes, page_title, page_subtitle = NULL) {
  n <- length(boxes)
  heights <- vapply(boxes, function(b) TITLE_H + 2 * PAD + length(b$detail) * LINE_H, numeric(1))

  y_top <- numeric(n); y_bottom <- numeric(n)
  cursor <- PAGE_UNIT_HEIGHT - TOP_MARGIN
  for (i in seq_len(n)) {
    y_top[i] <- cursor
    y_bottom[i] <- cursor - heights[i]
    cursor <- y_bottom[i] - GAP
  }

  rect_df <- data.frame(
    xmin = X_LEFT, xmax = X_LEFT + BOX_WIDTH,
    ymin = y_bottom, ymax = y_top,
    fill_group = vapply(boxes, function(b) b$fill, character(1))
  )
  title_df <- data.frame(
    x = X_LEFT + PAD, y = y_top - TITLE_H / 2,
    label = vapply(boxes, function(b) b$title, character(1))
  )
  detail_df <- data.frame(
    x = X_LEFT + PAD,
    y = y_top - TITLE_H - PAD - LINE_H / 2,
    label = vapply(boxes, function(b) paste0(paste0("- ", b$detail), collapse = "\n"), character(1))
  )

  p <- ggplot() +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_group),
      color = "grey30", linewidth = 0.4
    ) +
    geom_text(
      data = title_df, aes(x = x, y = y, label = label),
      hjust = 0, vjust = 0.5, fontface = "bold", size = 4.1, family = "sans"
    ) +
    geom_text(
      data = detail_df, aes(x = x, y = y, label = label),
      hjust = 0, vjust = 1, size = 3.1, lineheight = 1.15, family = "sans"
    ) +
    scale_fill_manual(values = FILL_PALETTE, guide = "none") +
    annotate("text", x = X_LEFT, y = PAGE_UNIT_HEIGHT - 6, label = page_title,
             hjust = 0, size = 5.3, fontface = "bold") +
    coord_cartesian(xlim = c(0, PAGE_UNIT_WIDTH), ylim = c(0, PAGE_UNIT_HEIGHT), clip = "off") +
    theme_void()

  if (n > 1) {
    arrow_df <- data.frame(
      x = X_LEFT + BOX_WIDTH / 2, xend = X_LEFT + BOX_WIDTH / 2,
      y = y_bottom[-n], yend = y_top[-1]
    )
    p <- p + geom_segment(
      data = arrow_df, aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(0.12, "inches")), color = "grey30", linewidth = 0.5
    )
  }

  if (!is.null(page_subtitle)) {
    p <- p + annotate("text", x = X_LEFT, y = PAGE_UNIT_HEIGHT - 11, label = page_subtitle,
                       hjust = 0, size = 3.3, color = "grey30")
  }

  # compact legend, same on every page
  leg_labels <- c("Data source", "QC / cohort", "Analysis", "Biomarker discovery", "Consensus", "Reporting")
  leg_fill   <- c("source", "qc", "analysis", "biomarker", "consensus", "report")
  leg_x0     <- X_LEFT + BOX_WIDTH - 3
  leg_y      <- PAGE_UNIT_HEIGHT - 6 - seq_along(leg_labels) * 4.4
  leg_rect <- data.frame(xmin = leg_x0 - 30, xmax = leg_x0 - 27, ymin = leg_y - 1.4, ymax = leg_y + 1.4, f = leg_fill)
  leg_text <- data.frame(x = leg_x0 - 25, y = leg_y, label = leg_labels)
  p <- p +
    geom_rect(data = leg_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = f), color = "grey40", linewidth = 0.3) +
    geom_text(data = leg_text, aes(x = x, y = y, label = label), hjust = 0, size = 2.7, color = "grey30")

  p
}

# ---- PAGE 1: data sources, QC, cohort definition ----------------------------
page1_boxes <- list(
  list(fill = "source", title = "Data sources", detail = c(
    "Clinical_Akila_FINAL_cleaned_20260610.xlsx (design n = 173; REAL clinical n = 172)",
    "Protein_pg_ml_Akila (Nomic Omni 1000 nELISA, ~973 proteins, 3 assay plates)",
    "sample_IDs_before_after.xlsx (pre/post sample mapping)"
  )),
  list(fill = "qc", title = "QC filtering (00-02)", detail = c(
    "Remove samples flagged sample_comments == 'qc_fail'",
    "No imputation for missing values, LLOD, or ULOD",
    "'Low bead count' values recoded to NA",
    "Batch correction (limma::removeBatchEffect) used for QC visualization ONLY",
    "Inferential models use raw log2 matrix with plate as covariate, not batch-corrected"
  )),
  list(fill = "qc", title = "Cohort & group definition", detail = c(
    "Design target: 173 baseline / 29 post. REAL: 154 baseline-matched / 28 post-matched",
    "1L_Treatment_Class (9 raw classes) collapsed to 7 New_Group levels",
    "Classes 7-9 (Chemo+IO+RT, RT only, RT+IO sequential) merged into 'Other'",
    "REAL group n (matched, n=154): CTx 36, Sx 36, IO 26, CTx+IO 22, CTx+RTx 15, TKI 13, Other 6",
    "REAL: 27 patients have both pre & post samples (13 of these are surgery pairs)"
  ))
)
p1 <- build_flow_page(page1_boxes, "nELISA NSCLC pipeline - 1/4: data, QC, cohort (design vs. REAL)")

# ---- PAGE 2: parallel contrasts (Analyses 1-5) with cutoffs -----------------
page2_boxes <- list(
  list(fill = "analysis", title = "A1 - Seven-group omnibus", detail = c(
    "limma lmFit/eBayes, design ~ 0 + group + plate_barcode",
    "Feature kept only if >= 3 non-missing samples per group",
    "Omnibus F-test (topTableF), ranked by adj.P.Val",
    "One-vs-rest logFC contrasts computed per group",
    "REAL: 972 features tested; top hit CSF1R, adj.p = 4.1e-288 (strong group signal, as expected)"
  )),
  list(fill = "analysis", title = "A2 - Paired pre vs post", detail = c(
    "limma paired design ~ pair_id + timepoint",
    "Feature kept only if >= 5 pairs have non-missing values each side",
    "REAL: surgery n = 13 pairs (formally tested); top hit adj.p = 0.195 - no significant hits",
    "REAL: remaining 14 pairs (logFC only) - CTx+IO 4, CTx 3, CTx+RTx 3, IO 2, TKI 1, Other 1"
  )),
  list(fill = "analysis", title = "A3 / A4 - Response & TKI vs IO", detail = c(
    "A3: PR/CR vs SD/PD, all patients - preliminary, explicitly underpowered",
    "REAL primary: PR 57 vs PD 39 (n=96), top hit adj.p = 0.835 - no significant hits",
    "REAL secondary: PR/CR 78 vs SD/PD 76 (n=154), top hit adj.p = 1.00 - no significant hits",
    "A4: TKI (n=13) vs IO (REAL n=26; design doc listed 28) - descriptive logFC only"
  )),
  list(fill = "analysis", title = "A5 - Survival (OS / PFS)", detail = c(
    "Univariate Cox per protein: coxph(Surv(time, event) ~ scale(x))",
    "Minimum n = 18, minimum events = 6 per model (delta-OS threshold)",
    "BH-FDR correction across all tested proteins",
    "icenReg sensitivity analysis for interval-censored PFS cases",
    "REAL OS: n=154, events=87; top hit GDF15, adj.p = 1.8e-5 (FDR<0.10 CONFIRMED)",
    "REAL PFS (main): n=140, events=112; top hits GDF2/CEACAM5 tied, adj.p = 0.094 (borderline)"
  )),
  list(fill = "analysis", title = "FDR reality check", detail = c(
    "Primary significance cutoff: adj.P.Val < 0.10",
    "REAL: 0 hits at FDR<0.10 in A2 paired (min adj.p=0.195) and A3 response (min adj.p=0.835/1.00)",
    "REAL: A1 omnibus and A5-OS are NOT underpowered - CSF1R and GDF15 both pass FDR<0.10 easily",
    "Power limit applies to response/recurrence signal only, not omnibus or OS survival"
  ))
)
p2 <- build_flow_page(page2_boxes, "nELISA NSCLC pipeline - 2/4: parallel contrasts, cutoffs & REAL outcomes")

# ---- PAGE 3: biomarker discovery methods and thresholds ---------------------
page3_boxes <- list(
  list(fill = "biomarker", title = "08 - LASSO + stability selection", detail = c(
    "cv.glmnet(family='binomial', alpha=1) on complete-case features only",
    "Selected: non-zero coefficient at lambda.min",
    "Stability: 100 resamples (80% subsample each), stability_fraction = selection frequency",
    "Stable-feature flag threshold: stability_fraction >= 0.60",
    "REAL comparison A (PR vs PD): 0 selected at lambda.min despite stability up to 0.99",
    "REAL comparison B (PR/CR vs SD/PD): 11 selected; top = SHH (1.00), MMP7 (0.93), GDF11 (0.91)"
  )),
  list(fill = "biomarker", title = "08b - STABL discovery", detail = c(
    "Adaptive LASSO within the STABL framework (Hedou et al., Nat Biotechnol 2024)",
    "Artificial feature injection drives false discovery proportion (FDP) control",
    "Data-driven stability threshold optimized per comparison",
    "Retained if: exceeds stability threshold AND model P < 0.05",
    "REAL recurrence tau* (low): IL27 Sx=0.31, CTRC IO=0.29, KNG1 CTx+IO=0.30, AZU1 Sx pre/post=0.41",
    "REAL response tau* (high): SHH CRPR-vs-SDPD=0.99; CEACAM5 (3 comparisons) 0.96-0.99",
    "REAL Cox hits: GDF15 (OS, tau*=0.99), CEACAM5 (PFS, tau*=0.99)"
  )),
  list(fill = "biomarker", title = "08d - Alternative ML methods", detail = c(
    "Elastic Net: glmnet(alpha = 0.5)",
    "sPLS-DA: mixOmics",
    "Boruta / Random Forest variable importance",
    "REAL: all 4 methods TRUE for SHH (2 comparisons) and CEACAM5 (2 comparisons)"
  )),
  list(fill = "biomarker", title = "08e - Cross-validated performance", detail = c(
    "Repeated stratified k-fold CV AUC (replaces in-sample AUC)",
    "Boruta/RF originally reported in-sample AUC = 1.000 (overfit); now CV-corrected",
    "Reported metric: mean CV AUC +/- SD across repeats",
    "REAL best single-protein AUC: AKR1B10=0.69 (A), AR=0.66 (B) - neither clears 0.70 cutoff"
  )),
  list(fill = "biomarker", title = "Consensus scoring rubric (max 10 pts)", detail = c(
    "limma FDR < 0.05 in primary comparison: 2 pts",
    "Cox OS HR significant: 2 pts | Cox PFS HR significant: 2 pts",
    "Selected by LASSO: 2 pts | Stability >= 60%: 1 pt | Single-protein AUC >= 0.70: 1 pt",
    "Priority validation candidates: total score >= 6 / 10",
    "REAL: only 1 protein scored >= 6 overall - CEACAM5 (PR vs PD), score = 6",
    "Driven by OS+PFS Cox significance; nominal limma p=0.04 but adj.p=0.995"
  ))
)
p3 <- build_flow_page(page3_boxes, "nELISA NSCLC pipeline - 3/4: biomarker discovery cutoffs & REAL outcomes")

# ---- PAGE 4: consensus, pathway/network, reporting --------------------------
page4_boxes <- list(
  list(fill = "consensus", title = "Cross-method consensus (08f-08i)", detail = c(
    "Cross-method overlap tracks 4 ML methods (EN, sPLS-DA, Boruta, STABL); LASSO tracked separately",
    "REAL: SHH and CEACAM5 both TRUE on all 4 methods in 2 comparisons each",
    "Data note: n_methods_selected reads 3 despite 4 TRUE flags - check 08d/08f counting logic",
    "REAL: only CEACAM5 (PR vs PD) clears the full 10-pt consensus rubric (score >= 6)",
    "KM stratification uses median split (High/Low); boxplots keep continuous distribution"
  )),
  list(fill = "consensus", title = "Survival signal flags", detail = c(
    "REAL: GDF15 - OS top hit, adj.p = 1.8e-5 - FDR<0.10 CONFIRMED",
    "REAL: GDF2 / CEACAM5 tied - PFS top hits, adj.p = 0.094 (borderline FDR<0.10)",
    "OS signal is FDR-robust; PFS signal remains fragile and flagged for independent validation"
  )),
  list(fill = "consensus", title = "09 - Pathway & network analysis", detail = c(
    "fgsea: GO BP/MF/CC, KEGG, Reactome on ranked gene lists (BH-FDR < 0.05)",
    "STRING PPI network per contrast, built when enough genes are available",
    "REAL: signal is sparse - most collections show 0 significant pathways per contrast",
    "REAL: GO_MF is the exception (up to 16 significant pathways, e.g. CTx vs rest)"
  )),
  list(fill = "report", title = "10-11 - Poster & reporting outputs", detail = c(
    "Poster: 3 panels (IO/CTx+IO responders, Chemo+other responders, Surgery recurrence)",
    "REAL: Poster1 PR/CR21 vs SD/PD31; Poster2 PR/CR24 vs SD/PD29; Poster3 no-recur20 vs recur16",
    "Self-contained HTML report (inlined CSS/JS) + PDF report",
    "All figures saved as PDF (vector, for submission) and PNG (300 dpi, for slides)"
  ))
)
p4 <- build_flow_page(page4_boxes, "nELISA NSCLC pipeline - 4/4: consensus, pathway, reporting - REAL outcomes")

# ---- write all pages to one PDF (direct device calls, no wrapper function) -
cairo_pdf(OUTPUT_PDF, width = 12, height = 16.7, onefile = TRUE)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

message("Pipeline overview PDF written to: ", OUTPUT_PDF)
