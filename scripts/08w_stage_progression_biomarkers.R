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
  if (!requireNamespace("clinfun", quietly = TRUE)) install.packages("clinfun")
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(clinfun)
  library(glmnet)
  library(ggrepel)
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
RNG_SEED     <- 42
EN_ALPHA     <- 0.5
MIN_PER_GRP  <- 10   # min samples per stage class for univariate tests
FDR_CUTOFF   <- 0.1
LFC_THR      <- 0.5  # |log2FC| threshold, matches 03_analysis1_seven_group.R convention

# ── Local functions (self-contained) ───────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

save_ggplot <- function(p, path_stub, width, height) {
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white"); print(p); dev.off()
}

# Project-standard volcano theme/plot (matches 03_analysis1_seven_group.R /
# 05_analysis3_response.R exactly): bordered theme_bw, dual Primary/Secondary
# FDR shapes, group-direction colour coding, bottom horizontal legend.
akila_volcano_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(color = "black", size = 12),
      axis.title        = element_text(color = "black", size = 14, face = "bold"),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
      legend.title      = element_text(face = "bold", size = 10),
      legend.text       = element_text(size = 9),
      plot.margin       = ggplot2::margin(68, 56, 46, 56)
    )
}

make_volcano <- function(df, title, path_stub,
                         lfc_thr = LFC_THR, padj_thr = FDR_CUTOFF,
                         label_col = "Protein_Label",
                         group_a_label = "Group A", group_b_label = "Group B") {
  if (nrow(df) == 0) return(invisible(NULL))
  x_label       <- paste0("Fold change log2(", group_b_label, " / ", group_a_label, ")")
  secondary_lbl <- paste0("Secondary (FDR<=", formatC(padj_thr, format = "f", digits = 2), ")")
  left_label    <- paste("Higher in", group_a_label)
  right_label   <- paste("Higher in", group_b_label)

  plot_df <- df %>%
    dplyr::mutate(
      neg_log10 = -log10(pmax(adj_p_value, 1e-300)),
      status = dplyr::case_when(
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC >=  lfc_thr ~ right_label,
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC <= -lfc_thr ~ left_label,
        TRUE ~ "Not significant"
      ),
      sig_level = dplyr::case_when(
        status == "Not significant" ~ "Not significant",
        !is.na(adj_p_value) & adj_p_value <= 0.05 ~ "Primary (FDR<=0.05)",
        TRUE ~ secondary_lbl
      )
    )

  ann <- plot_df %>%
    dplyr::filter(sig_level == "Primary (FDR<=0.05)") %>%
    dplyr::arrange(adj_p_value, dplyr::desc(abs(logFC)))
  if (nrow(ann) == 0) {
    ann <- plot_df %>%
      dplyr::filter(status != "Not significant") %>%
      dplyr::arrange(adj_p_value, dplyr::desc(abs(logFC))) %>%
      dplyr::slice_head(n = 18)
  }

  x_rng         <- range(plot_df$logFC, na.rm = TRUE)
  y_thr         <- -log10(padj_thr)
  y_max         <- max(plot_df$neg_log10, y_thr, na.rm = TRUE)
  x_pad         <- max(0.35, diff(x_rng) * 0.22)
  y_pad         <- max(0.18, y_max * 0.16)
  x_panel_left  <- x_rng[1] - 0.45 * x_pad
  x_panel_right <- x_rng[2] + 0.45 * x_pad
  y_panel_top   <- y_max  + 0.12 * y_pad
  y_label       <- y_panel_top + 0.9 * y_pad

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = logFC, y = neg_log10)) +
    ggplot2::geom_point(ggplot2::aes(color = status, shape = sig_level), size = 2.4, alpha = 0.82, stroke = 0.7) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, 0, lfc_thr),
               linetype = c("dashed", "solid", "dashed"),
               linewidth = c(0.6, 0.4, 0.6), color = "black") +
    ggplot2::geom_hline(yintercept = y_thr, linetype = "dashed", linewidth = 0.6, color = "black") +
    ggplot2::annotate("text", x = x_panel_left - 0.18 * x_pad, y = y_label,
             label = paste0("Higher in\n", group_a_label),
             hjust = 0, vjust = 0, size = 4.8, fontface = "bold", color = "#2166AC") +
    ggplot2::annotate("text", x = x_panel_right + 0.18 * x_pad, y = y_label,
             label = paste0("Higher in\n", group_b_label),
             hjust = 1, vjust = 0, size = 4.8, fontface = "bold", color = "#B2182B") +
    ggplot2::annotate("text", x = x_rng[1] + 0.02 * diff(x_rng), y = y_thr,
             label = paste0(" P Adj. = ", formatC(padj_thr, format = "f", digits = 2)),
             hjust = 0, vjust = -0.6, size = 3.8, fontface = "italic") +
    ggplot2::scale_color_manual(values = c(
      stats::setNames("#2166AC", left_label), stats::setNames("#B2182B", right_label),
      "Not significant" = "grey85"), drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("Not significant" = 1, "Primary (FDR<=0.05)" = 19,
                                  stats::setNames(17, secondary_lbl)),
                       breaks = c("Not significant", "Primary (FDR<=0.05)", secondary_lbl),
                       drop = FALSE) +
    ggrepel::geom_text_repel(data = ann, ggplot2::aes(label = .data[[label_col]]),
                             size = 3.1, fontface = "bold", max.overlaps = Inf,
                             box.padding = 0.5, point.padding = 0.2,
                             segment.color = "grey60", segment.size = 0.2,
                             show.legend = FALSE) +
    ggplot2::labs(title = title, x = x_label, y = "P value (-log10)",
         color = NULL, shape = "Significance") +
    ggplot2::guides(color = "none", shape = ggplot2::guide_legend(override.aes = list(color = "black"))) +
    ggplot2::coord_cartesian(xlim = c(x_panel_left, x_panel_right),
                    ylim = c(0, y_panel_top), clip = "off") +
    akila_volcano_theme()
  cairo_pdf(paste0(path_stub, ".pdf"), width = 8.6, height = 8.2); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 8.6, height = 8.2, units = "in", res = 300, bg = "white"); print(p); dev.off()
}

# MA plot (avg expression vs logFC), same colour/shape/status logic and
# theme as make_volcano() above, for the same primary/secondary FDR tiers.
make_ma_plot <- function(df, title, path_stub,
                         lfc_thr = LFC_THR, padj_thr = FDR_CUTOFF,
                         label_col = "Protein_Label",
                         group_a_label = "Group A", group_b_label = "Group B") {
  if (nrow(df) == 0) return(invisible(NULL))
  secondary_lbl <- paste0("Secondary (FDR<=", formatC(padj_thr, format = "f", digits = 2), ")")
  left_label    <- paste("Higher in", group_a_label)
  right_label   <- paste("Higher in", group_b_label)

  plot_df <- df %>%
    dplyr::mutate(
      status = dplyr::case_when(
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC >=  lfc_thr ~ right_label,
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC <= -lfc_thr ~ left_label,
        TRUE ~ "Not significant"
      ),
      sig_level = dplyr::case_when(
        status == "Not significant" ~ "Not significant",
        !is.na(adj_p_value) & adj_p_value <= 0.05 ~ "Primary (FDR<=0.05)",
        TRUE ~ secondary_lbl
      )
    )

  ann <- plot_df %>%
    dplyr::filter(sig_level == "Primary (FDR<=0.05)") %>%
    dplyr::arrange(adj_p_value, dplyr::desc(abs(logFC)))
  if (nrow(ann) == 0) {
    ann <- plot_df %>%
      dplyr::filter(status != "Not significant") %>%
      dplyr::arrange(adj_p_value, dplyr::desc(abs(logFC))) %>%
      dplyr::slice_head(n = 18)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = avg_expr, y = logFC)) +
    ggplot2::geom_point(ggplot2::aes(color = status, shape = sig_level), size = 2.4, alpha = 0.82, stroke = 0.7) +
    ggplot2::geom_hline(yintercept = c(-lfc_thr, 0, lfc_thr),
               linetype = c("dashed", "solid", "dashed"),
               linewidth = c(0.6, 0.4, 0.6), color = "black") +
    ggplot2::scale_color_manual(values = c(
      stats::setNames("#2166AC", left_label), stats::setNames("#B2182B", right_label),
      "Not significant" = "grey85"), drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("Not significant" = 1, "Primary (FDR<=0.05)" = 19,
                                  stats::setNames(17, secondary_lbl)),
                       breaks = c("Not significant", "Primary (FDR<=0.05)", secondary_lbl),
                       drop = FALSE) +
    { if (nrow(ann) > 0) ggrepel::geom_text_repel(
        data = ann, ggplot2::aes(label = .data[[label_col]]),
        size = 3.1, fontface = "bold", max.overlaps = Inf,
        box.padding = 0.5, point.padding = 0.2,
        segment.color = "grey60", segment.size = 0.2, show.legend = FALSE) } +
    ggplot2::labs(title = title, x = "Average log2 expression (AveExpr)",
         y = paste0("log2 fold change (", group_b_label, " / ", group_a_label, ")"),
         color = NULL, shape = "Significance") +
    ggplot2::guides(color = "none", shape = ggplot2::guide_legend(override.aes = list(color = "black"))) +
    akila_volcano_theme()
  cairo_pdf(paste0(path_stub, ".pdf"), width = 8.6, height = 7.4); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 8.6, height = 7.4, units = "in", res = 300, bg = "white"); print(p); dev.off()
}

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  dplyr::left_join(df, feature_meta, by = "nELISA_ID")
}

# Maps raw Stage_Group strings (which may carry sub-stage letters, e.g.
# "IIIA", "IVB") to an ordinal rank 1-4. Checked longest-prefix-first
# ("IV" before "III" before "II" before "I") to avoid "IV" matching "I".
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

run_binary_limma <- function(mat_samples_by_features, sample_df, group_col,
                             group_levels, adjust_cols = character(), min_per_group = 3) {
  keep <- !is.na(sample_df[[group_col]]) & sample_df[[group_col]] %in% group_levels
  df   <- sample_df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  grp  <- factor(df[[group_col]], levels = group_levels)
  if (length(unique(grp)) < 2 || any(table(grp) < min_per_group)) return(NULL)
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  feature_keep <- colSums(!is.na(mat_sub[grp == group_levels[1], , drop = FALSE])) >= min_per_group &
    colSums(!is.na(mat_sub[grp == group_levels[2], , drop = FALSE])) >= min_per_group
  mat_sub <- mat_sub[, feature_keep, drop = FALSE]
  if (ncol(mat_sub) == 0) return(NULL)
  design_df <- data.frame(primary_group = relevel(grp, ref = group_levels[1]))
  for (col in adjust_cols) {
    vals <- df[[col]]
    if (is.numeric(vals)) {
      if (length(unique(vals[!is.na(vals)])) > 1) design_df[[col]] <- vals
    } else {
      fac <- factor(vals)
      if (length(unique(fac[!is.na(fac)])) > 1) design_df[[col]] <- fac
    }
  }
  design <- model.matrix(~ ., data = design_df)
  if (qr(design)$rank < ncol(design)) return(NULL)  # guard: rank-deficient (confounded) design
  coef_name <- grep("^primary_group", colnames(design), value = TRUE)[1]
  fit <- limma::lmFit(t(mat_sub), design)
  fit <- limma::eBayes(fit)
  limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("nELISA_ID") %>%
    tibble::as_tibble() %>%
    dplyr::transmute(nELISA_ID, logFC, t_stat = t, p_value = P.Value,
                      adj_p_value = adj.P.Val, avg_expr = AveExpr)
}

# Forced-covariate design: clinical covariates unpenalized, proteins
# penalized. Stage itself is NEVER included here (it is the outcome).
build_clinical_covariates_nostage <- function(sample_ids, clinical_df = analysis1_df) {
  cd  <- clinical_df[match(sample_ids, clinical_df$sample_id), ]
  age <- as.numeric(difftime(cd$start_date, cd$dob, units = "days")) / 365.25
  mat <- cbind(
    age        = age,
    sex_female = as.integer(cd$sex == "FEMALE"),
    hist_adeno = as.integer(cd$histology == "Adenocarcinoma"),
    tx_chemo   = as.integer(grepl("Chemo",   cd$treatment_class, ignore.case = TRUE)),
    tx_io      = as.integer(grepl("IO",      cd$treatment_class, ignore.case = TRUE)),
    tx_rt      = as.integer(grepl("RT",      cd$treatment_class, ignore.case = TRUE)),
    tx_surgery = as.integer(grepl("Surgery", cd$treatment_class, ignore.case = TRUE)),
    tx_tki     = as.integer(grepl("TKI",     cd$treatment_class, ignore.case = TRUE))
  )
  for (j in seq_len(ncol(mat))) {
    mat[is.na(mat[, j]), j] <- stats::median(mat[, j], na.rm = TRUE)
  }
  col_sd <- apply(mat, 2, sd, na.rm = TRUE)
  mat <- mat[, col_sd > 0, drop = FALSE]
  rownames(mat) <- sample_ids
  mat
}


# Drop proteins with >20% missingness, median-impute any remainder.
# Matches the missingness-handling convention used in
# 09_holdout_validation_forced_covariates.R before glmnet fitting.
impute_for_glmnet <- function(X, max_missing_frac = 0.20) {
  miss_frac <- colMeans(is.na(X))
  X <- X[, miss_frac <= max_missing_frac, drop = FALSE]
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[, j])) X[is.na(X[, j]), j] <- stats::median(X[, j], na.rm = TRUE)
  }
  X
}

run_elastic_net_binomial <- function(X, y, X_clin) {
  X <- impute_for_glmnet(X)
  X_std <- scale(X); X_std[is.nan(X_std)] <- 0
  X_design <- cbind(X_clin[rownames(X_std), , drop = FALSE], X_std)
  pf <- c(rep(0, ncol(X_clin)), rep(1, ncol(X_std)))
  set.seed(RNG_SEED)
  cv_fit <- glmnet::cv.glmnet(X_design, y, family = "binomial", alpha = EN_ALPHA,
                              penalty.factor = pf, nfolds = min(5, min(table(y))))
  coef_all <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
  names(coef_all) <- colnames(X_design)
  coef_v <- coef_all[colnames(X)]
  prob <- as.numeric(predict(cv_fit, newx = X_design, s = "lambda.min", type = "response"))
  op <- options(warn = -1)
  auc_val <- tryCatch(as.numeric(pROC::auc(pROC::roc(y, prob, quiet = TRUE))), error = function(e) NA_real_)
  options(op)
  tibble::tibble(nELISA_ID = colnames(X), en_coef = coef_v,
                 en_selected = coef_v != 0, en_model_auc = auc_val)
}

run_elastic_net_ordinal <- function(X, y_rank, X_clin) {
  X <- impute_for_glmnet(X)
  X_std <- scale(X); X_std[is.nan(X_std)] <- 0
  X_design <- cbind(X_clin[rownames(X_std), , drop = FALSE], X_std)
  pf <- c(rep(0, ncol(X_clin)), rep(1, ncol(X_std)))
  set.seed(RNG_SEED)
  cv_fit <- glmnet::cv.glmnet(X_design, y_rank, family = "gaussian", alpha = EN_ALPHA,
                              penalty.factor = pf, nfolds = 5)
  coef_all <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
  names(coef_all) <- colnames(X_design)
  coef_v <- coef_all[colnames(X)]
  pred <- as.numeric(predict(cv_fit, newx = X_design, s = "lambda.min"))
  rho  <- suppressWarnings(cor(pred, y_rank, method = "spearman"))
  tibble::tibble(nELISA_ID = colnames(X), en_coef = coef_v,
                 en_selected = coef_v != 0, en_model_spearman_rho = rho)
}

# ============================================================
# Build sample-level stage frame from analysis1_df (144-patient cohort)
# ============================================================
message("\n=== 08w | Building stage variables ===")

stage_df <- analysis1_df %>%
  dplyr::mutate(stage_rank = stage_to_rank(stage_group)) %>%
  dplyr::filter(!is.na(stage_rank), !is.na(stage_collapsed))

message("Stage_collapsed counts:"); print(table(stage_df$stage_collapsed, useNA = "ifany"))
message("Stage_rank (I-IV) counts:"); print(table(stage_df$stage_rank, useNA = "ifany"))
safe_write_csv(stage_df %>% dplyr::count(stage_collapsed, stage_rank),
               file.path(OUT_TBL, "stage_variable_counts.csv"))

expr_mat <- expr_log2_rep[stage_df$sample_id, , drop = FALSE]

# ============================================================
# PART A - Binary: I-II vs III-IV
# ============================================================
message("\n=== 08w | PART A - Binary I-II vs III-IV ===")

# A1a. Univariate, primary model (plate_barcode only)
binA_primary <- run_binary_limma(expr_mat, stage_df, "stage_collapsed",
                                  c("I_II", "III_IV"),
                                  adjust_cols = c("plate_barcode"),
                                  min_per_group = MIN_PER_GRP)
binA_primary <- add_feature_annotation(binA_primary, feature_meta)
safe_write_csv(binA_primary, file.path(OUT_TBL, "stage_binary_limma_primary.csv"))

# A1b. Sensitivity model, additionally adjusted for treatment group
#      (guarded: returns NULL if design is rank-deficient, i.e. stage
#      and new_group are perfectly confounded in this cohort)
binA_sens <- run_binary_limma(expr_mat, stage_df, "stage_collapsed",
                               c("I_II", "III_IV"),
                               adjust_cols = c("plate_barcode", "new_group"),
                               min_per_group = MIN_PER_GRP)
if (is.null(binA_sens)) {
  message("A1b sensitivity model skipped: stage_collapsed and new_group are rank-deficient ",
          "(likely near-perfectly confounded, e.g. surgery patients concentrated in one stage class).")
} else {
  binA_sens <- add_feature_annotation(binA_sens, feature_meta)
  safe_write_csv(binA_sens, file.path(OUT_TBL, "stage_binary_limma_sensitivity_newgroup.csv"))
}

# A2. Multivariate elastic net, forced clinical covariates (no stage, no new_group)
X_clin_A <- build_clinical_covariates_nostage(stage_df$sample_id)
y_A <- as.integer(factor(stage_df$stage_collapsed, levels = c("I_II", "III_IV"))) - 1L
en_A <- run_elastic_net_binomial(expr_mat[stage_df$sample_id, , drop = FALSE], y_A, X_clin_A)
en_A <- add_feature_annotation(en_A, feature_meta)
safe_write_csv(en_A, file.path(OUT_TBL, "stage_binary_elasticnet.csv"))
message("Binary EN model AUC (in-sample): ", round(en_A$en_model_auc[1], 3),
        " | proteins selected: ", sum(en_A$en_selected))

# A3. Volcano + MA plots for primary univariate model (project-standard style)
if (!is.null(binA_primary) && nrow(binA_primary) > 0) {
  make_volcano(binA_primary, "Stage I-II vs III-IV (plate-adjusted limma)",
               file.path(OUT_FIG, "stage_binary_volcano"),
               group_a_label = "I-II", group_b_label = "III-IV")
  make_ma_plot(binA_primary, "Stage I-II vs III-IV — MA plot (plate-adjusted limma)",
               file.path(OUT_FIG, "stage_binary_ma_plot"),
               group_a_label = "I-II", group_b_label = "III-IV")
}

# ============================================================
# PART B - Ordinal: I -> II -> III -> IV
# ============================================================
message("\n=== 08w | PART B - Ordinal trend I->IV ===")

stage_rank_f <- factor(stage_df$stage_rank, levels = 1:4, ordered = TRUE)
message("Per-stage n: "); print(table(stage_rank_f))

# B1. Per-protein Jonckheere-Terpstra trend test + Spearman rho
run_trend_tests <- function(mat, rank_f) {
  purrr::map_dfr(colnames(mat), function(pid) {
    v <- mat[, pid]
    ok <- !is.na(v)
    if (sum(ok) < 20 || length(unique(rank_f[ok])) < 3) {
      return(tibble::tibble(nELISA_ID = pid, jt_p_value = NA_real_, spearman_rho = NA_real_, spearman_p = NA_real_))
    }
    jt_p <- tryCatch(clinfun::jonckheere.test(v[ok], as.integer(rank_f[ok]))$p.value,
                      error = function(e) NA_real_)
    sp   <- tryCatch(cor.test(v[ok], as.integer(rank_f[ok]), method = "spearman", exact = FALSE),
                      error = function(e) NULL)
    tibble::tibble(nELISA_ID = pid, jt_p_value = jt_p,
                   spearman_rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
                   spearman_p   = if (is.null(sp)) NA_real_ else sp$p.value)
  })
}

ordB_trend <- run_trend_tests(expr_mat, stage_rank_f) %>%
  dplyr::mutate(jt_adj_p = p.adjust(jt_p_value, method = "BH"),
                spearman_adj_p = p.adjust(spearman_p, method = "BH")) %>%
  dplyr::arrange(jt_p_value)
ordB_trend <- add_feature_annotation(ordB_trend, feature_meta)
safe_write_csv(ordB_trend, file.path(OUT_TBL, "stage_ordinal_trend_tests.csv"))
message("Proteins with JT-trend FDR < ", FDR_CUTOFF, ": ",
        sum(ordB_trend$jt_adj_p < FDR_CUTOFF, na.rm = TRUE))

# B2. Multivariate elastic net, ordinal rank as continuous outcome proxy
X_clin_B <- build_clinical_covariates_nostage(stage_df$sample_id)
en_B <- run_elastic_net_ordinal(expr_mat[stage_df$sample_id, , drop = FALSE],
                                 stage_df$stage_rank, X_clin_B)
en_B <- add_feature_annotation(en_B, feature_meta)
safe_write_csv(en_B, file.path(OUT_TBL, "stage_ordinal_elasticnet.csv"))
message("Ordinal EN model Spearman rho (in-sample): ", round(en_B$en_model_spearman_rho[1], 3),
        " | proteins selected: ", sum(en_B$en_selected))

# B3. Faceted panels of ordinal-trend proteins at Primary (FDR<=0.05),
#     split by trend direction: increasing with stage ("Higher in III-IV")
#     vs decreasing with stage ("Higher in I-II") — direction taken from
#     the sign of the Spearman rho vs ordinal stage rank, same
#     "Higher in <group>" / Primary FDR<=0.05 convention as the binary
#     volcano/MA plots above. Saved as two separate figures, each with a
#     mean-expression trend line overlaid per facet.
PRIMARY_FDR <- 0.05

trend_primary <- ordB_trend %>%
  dplyr::filter(!is.na(jt_adj_p), jt_adj_p <= PRIMARY_FDR) %>%
  dplyr::mutate(direction = ifelse(spearman_rho > 0, "Higher in III-IV", "Higher in I-II")) %>%
  dplyr::arrange(direction, jt_p_value)

message("Ordinal trend proteins at Primary (FDR<=", PRIMARY_FDR, "): ", nrow(trend_primary),
        " | Higher in III-IV: ", sum(trend_primary$direction == "Higher in III-IV"),
        " | Higher in I-II: ", sum(trend_primary$direction == "Higher in I-II"))
safe_write_csv(trend_primary, file.path(OUT_TBL, "stage_ordinal_trend_primary_FDR05.csv"))

make_trend_facet_panel <- function(trend_df, direction_label, ncol = 5) {
  sub <- trend_df %>% dplyr::filter(direction == direction_label)
  if (nrow(sub) == 0) return(NULL)
  box_long <- purrr::map_dfr(seq_len(nrow(sub)), function(i) {
    pid <- sub$nELISA_ID[i]
    lbl <- coalesce(sub$Gene_Symbol[i], pid)
    tibble::tibble(
      facet_label = sprintf("%s (JT p=%.2g)", lbl, sub$jt_p_value[i]),
      stage = stage_rank_f,
      expr  = expr_mat[, pid]
    )
  }) %>%
    dplyr::mutate(facet_label = factor(facet_label, levels = unique(facet_label)),
                  stage_num = as.numeric(stage))

  ggplot2::ggplot(box_long, ggplot2::aes(stage, expr, fill = stage)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.6) +
    ggplot2::stat_summary(ggplot2::aes(x = stage_num, y = expr, group = facet_label), fun = mean,
                          geom = "line", linewidth = 0.9, colour = "black", inherit.aes = FALSE,
                          data = box_long) +
    ggplot2::stat_summary(ggplot2::aes(x = stage_num, y = expr, group = facet_label), fun = mean,
                          geom = "point", shape = 18, size = 2.2, colour = "black",
                          inherit.aes = FALSE, data = box_long) +
    ggplot2::scale_x_discrete(labels = c("I", "II", "III", "IV")) +
    ggplot2::scale_fill_manual(values = c("#fbb4ae", "#fdcdac", "#b3cde3", "#80cdc1"), guide = "none") +
    ggplot2::facet_wrap(~facet_label, scales = "free_y", ncol = ncol) +
    ggplot2::labs(title = paste0(direction_label, "  (Primary FDR<=", PRIMARY_FDR, ", n=", nrow(sub), ")"),
                  subtitle = "Black diamond/line = mean expression per stage",
                  x = "Stage", y = "log2 concentration") +
    ggplot2::theme_minimal(base_size = 9)
}

PANEL_NCOL <- 5
p_up   <- make_trend_facet_panel(trend_primary, "Higher in III-IV", ncol = PANEL_NCOL)
p_down <- make_trend_facet_panel(trend_primary, "Higher in I-II",   ncol = PANEL_NCOL)
n_up   <- sum(trend_primary$direction == "Higher in III-IV")
n_down <- sum(trend_primary$direction == "Higher in I-II")
rows_up   <- ceiling(n_up   / PANEL_NCOL)
rows_down <- ceiling(n_down / PANEL_NCOL)

if (!is.null(p_up)) {
  save_ggplot(p_up, file.path(OUT_FIG, "stage_ordinal_higher_III_IV_panel"),
              width = 3 * PANEL_NCOL, height = 3 * max(rows_up, 1))
} else {
  message("No proteins pass Primary (FDR<=", PRIMARY_FDR, ") in the 'Higher in III-IV' direction — figure skipped.")
}
if (!is.null(p_down)) {
  save_ggplot(p_down, file.path(OUT_FIG, "stage_ordinal_higher_I_II_panel"),
              width = 3 * PANEL_NCOL, height = 3 * max(rows_down, 1))
} else {
  message("No proteins pass Primary (FDR<=", PRIMARY_FDR, ") in the 'Higher in I-II' direction — figure skipped.")
}

message("\n08w_stage_progression_biomarkers.R complete")
