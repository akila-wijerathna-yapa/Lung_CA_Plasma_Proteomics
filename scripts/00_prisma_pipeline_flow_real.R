library(ggplot2)
library(grid)

# ---- EDIT FOR YOUR MACHINE --------------------------------------------------
PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA"
OUTPUT_PDF  <- file.path(PROJECT_DIR, "figures", "00_overview", "00_prisma_pipeline_flow_real.pdf")
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
  # ---- trunk: identification / screening (real counts from cohort_summary.csv)
  A     = make_box(12, 330, 40, 20, "source", "Clinical cohort identified",
                    "n = 172 patients"),
  A_n   = make_box(56, 330, 44, 20, "source", "Source",
                    "Clinical_Akila_FINAL_cleaned_20260610.xlsx"),
  B     = make_box(12, 298, 40, 20, "qc", "Proteomics QC filtering",
                    c("Raw rows = 223 -> 213 after qc_fail filter",
                      "(184 PRE + 29 POST retained)")),
  B_n   = make_box(56, 298, 44, 20, "qc", "Excluded",
                    "10 rows removed: sample_comments = 'qc_fail'"),
  C     = make_box(12, 266, 40, 20, "qc", "Duplicate resolution",
                    c("PRE: 184 rows -> 175 unique samples",
                      "POST: 29 rows -> 29 unique (no dups)")),
  C_n   = make_box(56, 266, 44, 20, "qc", "Method",
                    "See duplicate_sample_choice_table.csv"),
  D     = make_box(12, 234, 40, 24, "qc", "Clinically matched cohort",
                    c("PRE matched = 154 patients",
                      "POST matched = 28 patients")),
  D_n   = make_box(56, 234, 44, 24, "qc", "Pairing",
                    c("27 patients with both pre & post",
                      "13 of these are surgery pairs")),

  # ---- left branch: cross-sectional contrasts (real n's) --------------------
  L1    = make_box(1, 192, 30, 30, "analysis", "A1 - Seven-group omnibus (n=154)",
                    c("CTx 36, Sx 36, IO 26, CTx+IO 22",
                      "CTx+RTx 15, TKI 13, Other 6")),
  L2    = make_box(1, 154, 30, 32, "analysis", "A3/A4 - Response & TKI vs IO",
                    c("A3: PR 57 vs PD 39 (n=96)",
                      "A3b: PR/CR 78 vs SD/PD 76 (n=154)",
                      "A4: TKI n=13 vs IO n=26")),

  # ---- middle branch: paired analysis (real pairs) --------------------------
  M1    = make_box(35, 192, 30, 30, "analysis", "A2 - Paired pre/post (27 pairs)",
                    c("Surgery 13, CTx+IO 4, CTx 3",
                      "CTx+RTx 3, IO 2, TKI 1, Other 1")),

  # ---- right branch: survival + biomarker discovery (real outcomes) ---------
  R1    = make_box(69, 192, 30, 32, "analysis", "A5 - Survival (real n/events)",
                    c("OS: n=154, events=87, cens=67",
                      "PFS: n=140, events=112, cens=28",
                      "EFS: n=154, events=131, cens=23")),
  R2    = make_box(69, 156, 30, 28, "biomarker", "A6 - ML comparisons run (12)",
                    c("5 response, 4 recurrence",
                      "1 pre/post, 2 Cox (OS, PFS)")),
  R3    = make_box(69, 124, 30, 36, "biomarker", "STABL selected (real tau*)",
                    c("tau* range 0.29-0.99 across 12 runs",
                      "1-2 features retained per comparison",
                      "Sel./comparison: LASSO 0-11, STABL 1-2",
                      "EN 0-39, sPLS-DA 27-30, Boruta 0-9")),
  R4    = make_box(69, 84, 30, 40, "biomarker", "Consensus - real outcome",
                    c("Only 1 priority candidate: CEACAM5",
                      "(PR vs PD), score 6/10",
                      "SHH & CEACAM5: 4/4 ML methods in 2 comps",
                      "Full biomarker lists: see companion CSV export")),

  # ---- convergence: real pathway + poster outcomes ---------------------------
  CONV  = make_box(15, 34, 70, 32, "consensus", "09 - Pathway/network + poster (real)",
                    c("Pathway signal sparse (mostly 0 sig; GO_MF up to 16)",
                      "Poster1 IO/CTx+IO: PR/CR 21 vs SD/PD 31",
                      "Poster2 Chemo+other: PR/CR 24 vs SD/PD 29",
                      "Poster3 Surgery: No recur 20 vs Recur 16"))
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
  data.frame(x = 84, xend = 84, y = 160, yend = 156),           # R1 -> R2
  data.frame(x = 84, xend = 84, y = 128, yend = 124),           # R2 -> R3
  data.frame(x = 84, xend = 84, y = 88,  yend = 84),            # R3 -> R4
  data.frame(x = 16, xend = 16, y = 122, yend = 34),            # L2 -> convergence
  data.frame(x = 50, xend = 50, y = 162, yend = 34),            # M1 -> convergence
  data.frame(x = 84, xend = 84, y = 44,  yend = 34)             # R4 -> convergence
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
  annotate("text", x = 12, y = 345, label = "nELISA NSCLC study & analysis flow - REAL OUTCOMES",
           hjust = 0, size = 5, fontface = "bold") +
  coord_cartesian(xlim = c(0, PAGE_W), ylim = c(0, PAGE_H), clip = "off") +
  theme_void()

cairo_pdf(OUTPUT_PDF, width = 12, height = 17)
print(p)
dev.off()

message("Real-outcome PRISMA flow diagram written to: ", OUTPUT_PDF)
