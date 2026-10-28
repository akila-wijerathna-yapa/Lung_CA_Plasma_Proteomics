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

suppressPackageStartupMessages({
  if (!requireNamespace("multcompView", quietly = TRUE)) install.packages("multcompView")
  library(multcompView)
})

# ── Load checkpoints ───────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01, 03 loaded")
while (dev.cur() > 1) dev.off()

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning", "stage_progression")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning", "stage_progression")
dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ── Local constants ────────────────────────────────────────────
ANOVA_FDR  <- 0.05
PANEL_NCOL <- 5

# ── Local functions (self-contained) ───────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

save_ggplot <- function(p, path_stub, width, height) {
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white"); print(p); dev.off()
}

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  dplyr::left_join(df, feature_meta, by = "nELISA_ID")
}

# Maps raw Stage_Group strings to an ordinal rank 1-4 — same convention
# as 08w_stage_progression_biomarkers.R, reproduced here (self-contained).
stage_to_rank <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    startsWith(x, "IV")  ~ 4L,
    startsWith(x, "III") ~ 3L,
    startsWith(x, "II")  ~ 2L,
    startsWith(x, "I")   ~ 1L,
    TRUE ~ NA_integer_
  )
}

# One-way ANOVA + Tukey HSD + compact-letter-display for one protein.
# Uses an UNORDERED factor for the model fit (avoids polynomial-contrast
# complications with TukeyHSD on an ordered factor); level order "1..4"
# is preserved for display purposes only.
run_anova_letters <- function(expr_vec, stage_f) {
  ok <- !is.na(expr_vec)
  stage_plain <- droplevels(factor(as.character(stage_f[ok]), levels = c("1", "2", "3", "4")))
  df <- data.frame(expr = expr_vec[ok], stage = stage_plain)
  empty_letters <- stats::setNames(rep(NA_character_, 4), c("1", "2", "3", "4"))
  if (nlevels(df$stage) < 3 || nrow(df) < 20 || any(table(df$stage) < 3)) {
    return(list(anova_p = NA_real_, letters = empty_letters))
  }
  fit <- stats::aov(expr ~ stage, data = df)
  anova_p <- summary(fit)[[1]][["Pr(>F)"]][1]
  tuk <- stats::TukeyHSD(fit)
  pvals <- tuk$stage[, "p adj"]
  cld <- tryCatch(multcompView::multcompLetters(pvals)$Letters, error = function(e) NULL)
  full_letters <- empty_letters
  if (!is.null(cld)) full_letters[names(cld)] <- cld
  list(anova_p = anova_p, letters = full_letters)
}

# ============================================================
# Rebuild stage sample frame (must match 08w exactly)
# ============================================================
message("\n=== 08z | Rebuilding stage variables (matches 08w) ===")

stage_df <- analysis1_df %>%
  dplyr::mutate(stage_rank = stage_to_rank(stage_group)) %>%
  dplyr::filter(!is.na(stage_rank), !is.na(stage_collapsed))

stage_rank_f <- factor(stage_df$stage_rank, levels = 1:4, ordered = TRUE)
expr_mat <- expr_log2_rep[stage_df$sample_id, , drop = FALSE]
message("Per-stage n: "); print(table(stage_rank_f))

# ── Load the Primary (FDR<=0.05) trend protein list from 08w ─────
trend_path <- file.path(OUT_TBL, "stage_ordinal_trend_primary_FDR05.csv")
if (!file.exists(trend_path)) {
  stop("Cannot find ", trend_path, " -- run 08w_stage_progression_biomarkers.R first.")
}
trend_primary <- readr::read_csv(trend_path, show_col_types = FALSE)
message("Loaded ", nrow(trend_primary), " Primary (FDR<=0.05) trend proteins from 08w.")

# ============================================================
# One-way ANOVA + Tukey HSD + compact letters, per protein
# ============================================================
message("\n=== 08z | Running one-way ANOVA + Tukey HSD + letters ===")

anova_results <- purrr::map_dfr(seq_len(nrow(trend_primary)), function(i) {
  pid <- trend_primary$nELISA_ID[i]
  res <- run_anova_letters(expr_mat[, pid], stage_rank_f)
  tibble::tibble(
    nELISA_ID  = pid,
    direction  = trend_primary$direction[i],
    jt_p_value = trend_primary$jt_p_value[i],
    anova_p    = res$anova_p,
    letter_I   = res$letters["1"],
    letter_II  = res$letters["2"],
    letter_III = res$letters["3"],
    letter_IV  = res$letters["4"]
  )
})

anova_results <- anova_results %>%
  dplyr::mutate(anova_adj_p = p.adjust(anova_p, method = "BH")) %>%
  dplyr::arrange(anova_p)
anova_results <- add_feature_annotation(anova_results, feature_meta)
safe_write_csv(anova_results, file.path(OUT_TBL, "stage_anova_letters_results.csv"))
message("ANOVA significant (adj p<", ANOVA_FDR, "): ",
        sum(anova_results$anova_adj_p < ANOVA_FDR, na.rm = TRUE), " / ", nrow(anova_results))

# ============================================================
# Faceted boxplots with significance letters, split by direction
# (mirrors the "Higher in III-IV" / "Higher in I-II" split from 08w)
# ============================================================
message("\n=== 08z | Building faceted ANOVA + letters panels ===")

make_anova_facet_panel <- function(results_df, direction_label, ncol = PANEL_NCOL) {
  sub <- results_df %>% dplyr::filter(direction == direction_label, !is.na(anova_p))
  if (nrow(sub) == 0) return(NULL)

  box_long <- purrr::map_dfr(seq_len(nrow(sub)), function(i) {
    pid <- sub$nELISA_ID[i]
    lbl <- coalesce(sub$Gene_Symbol[i], pid)
    tibble::tibble(
      facet_label = sprintf("%s (ANOVA p=%.2g)", lbl, sub$anova_p[i]),
      stage = stage_rank_f,
      expr  = expr_mat[, pid]
    )
  }) %>%
    dplyr::mutate(facet_label = factor(facet_label, levels = unique(facet_label)))

  # Per-protein, per-stage letter positioned just above that stage's max
  letter_long <- purrr::map_dfr(seq_len(nrow(sub)), function(i) {
    pid <- sub$nELISA_ID[i]
    lbl <- coalesce(sub$Gene_Symbol[i], pid)
    facet_lbl <- sprintf("%s (ANOVA p=%.2g)", lbl, sub$anova_p[i])
    v <- expr_mat[, pid]
    stage_max <- tapply(v, stage_rank_f, max, na.rm = TRUE)
    protein_range <- diff(range(v, na.rm = TRUE))
    letters_vec <- c(sub$letter_I[i], sub$letter_II[i], sub$letter_III[i], sub$letter_IV[i])
    tibble::tibble(
      facet_label = facet_lbl,
      stage  = factor(1:4, levels = 1:4, ordered = TRUE),
      y_pos  = as.numeric(stage_max) + 0.08 * protein_range,
      letter = letters_vec
    )
  }) %>%
    dplyr::mutate(facet_label = factor(facet_label, levels = levels(box_long$facet_label)))

  ggplot2::ggplot(box_long, ggplot2::aes(stage, expr, fill = stage)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.3, size = 0.5) +
    ggplot2::geom_text(data = letter_long, ggplot2::aes(x = stage, y = y_pos, label = letter),
                       inherit.aes = FALSE, fontface = "bold", size = 3.2, na.rm = TRUE) +
    ggplot2::scale_x_discrete(labels = c("I", "II", "III", "IV")) +
    ggplot2::scale_fill_manual(values = c("#fbb4ae", "#fdcdac", "#b3cde3", "#80cdc1"), guide = "none") +
    ggplot2::facet_wrap(~facet_label, scales = "free_y", ncol = ncol) +
    ggplot2::labs(title = paste0(direction_label, "  (one-way ANOVA + Tukey HSD letters, n=", nrow(sub), ")"),
                  subtitle = "Groups sharing a letter are not significantly different (Tukey HSD, alpha=0.05)",
                  x = "Stage", y = "log2 concentration") +
    ggplot2::theme_minimal(base_size = 9)
}

n_up   <- sum(anova_results$direction == "Higher in III-IV" & !is.na(anova_results$anova_p))
n_down <- sum(anova_results$direction == "Higher in I-II"   & !is.na(anova_results$anova_p))
rows_up   <- ceiling(n_up   / PANEL_NCOL)
rows_down <- ceiling(n_down / PANEL_NCOL)

p_up <- make_anova_facet_panel(anova_results, "Higher in III-IV")
if (!is.null(p_up)) {
  save_ggplot(p_up, file.path(OUT_FIG, "stage_anova_higher_III_IV_panel"),
              width = 3 * PANEL_NCOL, height = 3 * max(rows_up, 1))
} else {
  message("No proteins in 'Higher in III-IV' direction -- figure skipped.")
}

p_down <- make_anova_facet_panel(anova_results, "Higher in I-II")
if (!is.null(p_down)) {
  save_ggplot(p_down, file.path(OUT_FIG, "stage_anova_higher_I_II_panel"),
              width = 3 * PANEL_NCOL, height = 3 * max(rows_down, 1))
} else {
  message("No proteins in 'Higher in I-II' direction -- figure skipped.")
}

message("\n08z_stage_progression_anova_letters.R complete")
