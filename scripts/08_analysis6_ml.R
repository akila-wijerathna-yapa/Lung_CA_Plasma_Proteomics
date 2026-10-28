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
load(file.path(paths$rdata, "checkpoint_05.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 05, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
RNG_SEED           <- 42
N_RESAMPLES        <- 100   # stability selection bootstrap iterations
STABILITY_THR      <- 0.6   # fraction threshold for "stable" selection
AUC_THR            <- 0.70  # per-protein AUC threshold for scoring
LIMMA_P_THR        <- 0.05  # nominal p threshold for limma score
COX_FDR_THR        <- 0.10  # FDR threshold for Cox score
SCORE_PRIORITY_THR <- 6L    # minimum composite score for priority candidates

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  left_join(df, feature_meta, by = "nELISA_ID")
}

run_lasso_with_stability <- function(mat_samples_by_features, sample_df, outcome_col,
                                     positive_label, comparison_label,
                                     table_dir, figure_dir, n_resamples = N_RESAMPLES) {
  df <- sample_df %>% filter(!is.na(.data[[outcome_col]]))
  y  <- ifelse(df[[outcome_col]] == positive_label, 1, 0)
  if (length(unique(y)) < 2 || min(table(y)) < 5) return(tibble())

  mat_sub          <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  complete_features <- colnames(mat_sub)[colSums(is.na(mat_sub)) == 0]
  mat_complete     <- mat_sub[, complete_features, drop = FALSE]
  if (ncol(mat_complete) < 20) return(tibble())

  x <- scale(mat_complete)
  x[!is.finite(x)] <- 0   # zero-variance columns produce NaN after scale(); replace with 0
  op <- options(warn = -1)
  cv_fit <- glmnet::cv.glmnet(x = x, y = y, family = "binomial", alpha = 1,
                               nfolds = min(5, length(y)))
  options(op)
  coef_mat <- coef(cv_fit, s = "lambda.min")
  selected <- setdiff(rownames(coef_mat)[as.numeric(coef_mat[, 1]) != 0], "(Intercept)")

  # Stability selection
  stability <- setNames(rep(0, ncol(x)), colnames(x))
  set.seed(RNG_SEED)
  for (i in seq_len(n_resamples)) {
    idx <- sort(sample(seq_len(nrow(x)), size = max(ceiling(0.8 * nrow(x)), 10)))
    x_idx <- x[idx, , drop = FALSE]
    x_idx[!is.finite(x_idx)] <- 0
    fit_i <- tryCatch(
      glmnet::glmnet(x_idx, y[idx], family = "binomial", alpha = 1),
      error = function(e) NULL)
    if (is.null(fit_i)) next
    cf_i       <- coef(fit_i, s = fit_i$lambda[min(20, length(fit_i$lambda))])
    selected_i <- intersect(rownames(cf_i)[as.numeric(cf_i[, 1]) != 0], names(stability))
    stability[selected_i] <- stability[selected_i] + 1
  }
  stability <- stability / n_resamples

  prob    <- as.numeric(predict(cv_fit, newx = x, s = "lambda.min", type = "response"))
  roc_obj <- pROC::roc(y, prob, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  message("  ", comparison_label, " LASSO — features: ", ncol(x),
          " | selected: ", length(selected),
          " | model AUC: ", round(auc_val, 3))

  # LASSO CV plot (base graphics — direct device)
  cairo_pdf(file.path(figure_dir, paste0(comparison_label, "_lasso_cv.pdf")), width = 8, height = 6)
  plot(cv_fit); title(main = paste(comparison_label, "LASSO cross-validation"))
  dev.off()
  png(file.path(figure_dir, paste0(comparison_label, "_lasso_cv.png")),
      width = 8, height = 6, units = "in", res = 300, bg = "white")
  plot(cv_fit); title(main = paste(comparison_label, "LASSO cross-validation"))
  dev.off()

  # ROC curve
  roc_df <- tibble(specificity = roc_obj$specificities, sensitivity = roc_obj$sensitivities)
  p_roc  <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(color = "#2166AC", linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    labs(title    = paste(comparison_label, "LASSO ROC"),
         subtitle = paste("AUC =", round(auc_val, 3)),
         x = "False positive rate", y = "True positive rate") +
    theme_minimal(base_size = 11)
  cairo_pdf(file.path(figure_dir, paste0(comparison_label, "_lasso_roc.pdf")),
            width = 6.5, height = 5)
  print(p_roc); dev.off()
  png(file.path(figure_dir, paste0(comparison_label, "_lasso_roc.png")),
      width = 6.5, height = 5, units = "in", res = 300, bg = "white")
  print(p_roc); dev.off()

  # Stability bar for top stable features
  stable_df <- tibble(nELISA_ID = names(stability), stability_fraction = stability) %>%
    filter(stability_fraction > 0) %>%
    left_join(feature_meta %>% select(nELISA_ID, Protein_Label), by = "nELISA_ID") %>%
    arrange(desc(stability_fraction)) %>%
    slice_head(n = min(30, nrow(.))) %>%
    mutate(selected    = nELISA_ID %in% selected,
           label_order = factor(Protein_Label, levels = rev(Protein_Label)))
  if (nrow(stable_df) > 0) {
    p_stab <- ggplot(stable_df, aes(stability_fraction, label_order, fill = selected)) +
      geom_col(width = 0.75) +
      geom_vline(xintercept = STABILITY_THR, linetype = "dashed", color = "#B2182B") +
      scale_fill_manual(values = c("TRUE" = "#1565C0", "FALSE" = "#90A4AE"),
                        labels = c("TRUE" = "LASSO selected", "FALSE" = "Not selected")) +
      labs(title    = paste(comparison_label, "— stability selection"),
           subtitle = paste0("Red dashed line = stability threshold (", STABILITY_THR, ")"),
           x = "Selection frequency (100 resamples)", y = NULL, fill = NULL) +
      theme_minimal(base_size = 11) + theme(legend.position = "top")
    cairo_pdf(file.path(figure_dir, paste0(comparison_label, "_stability_bar.pdf")),
              width = 9, height = 7.5)
    print(p_stab); dev.off()
    png(file.path(figure_dir, paste0(comparison_label, "_stability_bar.png")),
        width = 9, height = 7.5, units = "in", res = 300, bg = "white")
    print(p_stab); dev.off()
  }

  tibble(
    nELISA_ID          = colnames(x),
    lasso_selected     = colnames(x) %in% selected,
    stability_fraction = unname(stability[colnames(x)]),
    lasso_model_auc    = auc_val
  )
}

single_protein_auc <- function(mat_samples_by_features, sample_df, outcome_col, positive_label) {
  df <- sample_df %>% filter(!is.na(.data[[outcome_col]]))
  y  <- ifelse(df[[outcome_col]] == positive_label, 1, 0)
  if (length(unique(y)) < 2) return(tibble())
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  aucs <- lapply(colnames(mat_sub), function(feature) {
    x    <- mat_sub[, feature]
    keep <- !is.na(x)
    if (sum(keep) < 15 || length(unique(y[keep])) < 2) return(NULL)
    roc_obj <- tryCatch(pROC::roc(y[keep], x[keep], quiet = TRUE, direction = "<"),
                        error = function(e) NULL)
    if (is.null(roc_obj)) return(NULL)
    tibble(nELISA_ID = feature, single_protein_auc = as.numeric(pROC::auc(roc_obj)))
  })
  bind_rows(aucs)
}

# ── ML specifications ─────────────────────────────────────────
# Use exploratory (plate_barcode only) results as the limma reference
# since the adjusted analysis (module 05) had no power.
message("\n=== 08 | ML setup ===")
ml_specs <- list(
  analysis6a_pr_vs_pd = list(
    df              = analysis3_df,
    outcome         = "response_prpd",
    positive        = "PR",
    comparison_name = "Analysis 6A: PR vs PD extreme contrast",
    limma_ref       = analysis3_primary_expl
  ),
  analysis6b_prcr_vs_sdpd = list(
    df              = analysis3_df,
    outcome         = "response_binary2",
    positive        = "PR/CR",
    comparison_name = "Analysis 6B: PR/CR vs SD/PD",
    limma_ref       = analysis3_secondary_expl
  )
)

# Cox lookup from module 07
cox_lookup <- full_join(
  os_res  %>% select(nELISA_ID, os_cox_adj_p  = adj_p_value),
  pfs_res %>% select(nELISA_ID, pfs_cox_adj_p = adj_p_value),
  by = "nELISA_ID"
)

tbl7  <- file.path(paths$tables,  "07_analysis6_machine_learning")
fig7  <- file.path(paths$figures, "07_analysis6_machine_learning")
ml_outputs     <- list()
priority_tables <- list()

# ── LASSO + AUC + scoring loop ────────────────────────────────
for (nm in names(ml_specs)) {
  spec  <- ml_specs[[nm]]
  ml_df <- spec$df %>% filter(!is.na(.data[[spec$outcome]]))
  n_pos <- sum(ml_df[[spec$outcome]] == spec$positive,  na.rm = TRUE)
  n_neg <- sum(ml_df[[spec$outcome]] != spec$positive,  na.rm = TRUE)
  message("\n--- ", spec$comparison_name, " ---")
  message("  ", spec$positive, ": n=", n_pos, " | other: n=", n_neg)
  if (nrow(ml_df) < 20 || n_pos < 5 || n_neg < 5) {
    message("  Skipped: insufficient class sizes (need >= 5 per class)")
    next
  }

  lasso_tbl <- run_lasso_with_stability(
    expr_log2_rep, ml_df, spec$outcome, spec$positive, nm, tbl7, fig7)
  if (nrow(lasso_tbl) == 0) { message("  LASSO returned 0 rows — skipped"); next }

  message("  Running single-protein AUC screen...")
  auc_tbl <- single_protein_auc(expr_log2_rep, ml_df, spec$outcome, spec$positive)
  if (nrow(auc_tbl) == 0) { message("  AUC table empty — skipped"); next }

  safe_write_csv(lasso_tbl %>%
    add_feature_annotation(feature_meta) %>%
    arrange(desc(stability_fraction), desc(lasso_selected)),
    file.path(tbl7, paste0(nm, "_lasso_stability.csv")))
  safe_write_csv(auc_tbl %>%
    add_feature_annotation(feature_meta) %>%
    arrange(desc(single_protein_auc)),
    file.path(tbl7, paste0(nm, "_single_protein_auc.csv")))

  # Composite biomarker score
  limma_ref <- spec$limma_ref
  if (is.null(limma_ref) || nrow(limma_ref) == 0) limma_ref <- tibble(nELISA_ID = character(), adj_p_value = numeric(), p_value = numeric())
  score_tbl <- limma_ref %>%
    transmute(nELISA_ID,
              Protein_Label    = dplyr::coalesce(Protein_Label, nELISA_ID),
              limma_adj_p      = adj_p_value,
              limma_p          = p_value) %>%
    left_join(cox_lookup,  by = "nELISA_ID") %>%
    left_join(lasso_tbl,   by = "nELISA_ID") %>%
    left_join(auc_tbl,     by = "nELISA_ID") %>%
    left_join(feature_meta %>% select(nELISA_ID, Gene_Symbol, Protein_Name, Uniprot_ID),
              by = "nELISA_ID") %>%
    mutate(
      score_limma     = ifelse(!is.na(limma_p) & limma_p < LIMMA_P_THR,     2L, 0L),
      score_os        = ifelse(!is.na(os_cox_adj_p)  & os_cox_adj_p  < COX_FDR_THR, 2L, 0L),
      score_pfs       = ifelse(!is.na(pfs_cox_adj_p) & pfs_cox_adj_p < COX_FDR_THR, 2L, 0L),
      score_lasso     = ifelse(isTRUE(lasso_selected), 2L, 0L),
      score_stability = ifelse(!is.na(stability_fraction) & stability_fraction >= STABILITY_THR, 1L, 0L),
      score_auc       = ifelse(!is.na(single_protein_auc) & single_protein_auc >= AUC_THR,       1L, 0L),
      biomarker_score = score_limma + score_os + score_pfs + score_lasso + score_stability + score_auc,
      priority_candidate = biomarker_score >= SCORE_PRIORITY_THR
    ) %>%
    arrange(desc(biomarker_score), limma_p)
  safe_write_csv(score_tbl, file.path(tbl7, paste0(nm, "_biomarker_scoring.csv")))
  ml_outputs[[nm]] <- score_tbl
  priority_tables[[nm]] <- score_tbl %>% filter(priority_candidate)
  message("  Biomarker scoring done — priority candidates: ",
          sum(score_tbl$priority_candidate, na.rm = TRUE),
          " | top score: ", max(score_tbl$biomarker_score, na.rm = TRUE))

  # Top AUC bar plot
  top_auc <- score_tbl %>%
    filter(!is.na(single_protein_auc)) %>%
    arrange(desc(single_protein_auc), desc(biomarker_score)) %>%
    slice_head(n = min(15, nrow(.))) %>%
    mutate(label = factor(Protein_Label, levels = rev(Protein_Label)))
  if (nrow(top_auc) > 0) {
    p_auc <- ggplot(top_auc, aes(single_protein_auc, label, fill = biomarker_score)) +
      geom_col(width = 0.75) +
      scale_fill_gradient(low = "#BBDEFB", high = "#0D47A1") +
      geom_vline(xintercept = AUC_THR, linetype = "dashed", color = "#B2182B") +
      labs(title    = paste(spec$comparison_name, "— top single-protein AUC"),
           subtitle = paste0("Red dashed = AUC threshold (", AUC_THR, ")"),
           x = "AUC", y = NULL, fill = "Composite\nscore") +
      theme_minimal(base_size = 11)
    cairo_pdf(file.path(fig7, paste0(nm, "_top_auc.pdf")), width = 8, height = 6.5)
    print(p_auc); dev.off()
    png(file.path(fig7, paste0(nm, "_top_auc.png")), width = 8, height = 6.5,
        units = "in", res = 300, bg = "white")
    print(p_auc); dev.off()
  }
}

# ── Priority candidates summary ───────────────────────────────
message("\n=== 08 | Priority candidates summary ===")
priority_all <- bind_rows(priority_tables, .id = "comparison")
safe_write_csv(priority_all, file.path(tbl7, "analysis6_priority_candidates_all.csv"))
message("Priority candidates across all comparisons: ", nrow(priority_all))

top_scoring_all <- bind_rows(lapply(names(ml_specs), function(nm) {
  score_path <- file.path(tbl7, paste0(nm, "_biomarker_scoring.csv"))
  if (!file.exists(score_path)) return(NULL)
  readr::read_csv(score_path, show_col_types = FALSE) %>%
    mutate(comparison = nm) %>%
    arrange(desc(biomarker_score), desc(single_protein_auc), desc(stability_fraction), limma_p) %>%
    slice_head(n = 20)
}))
safe_write_csv(top_scoring_all, file.path(tbl7, "analysis6_top20_scoring_candidates_all.csv"))

if (nrow(top_scoring_all) > 0) {
  top_summary <- top_scoring_all %>%
    filter(!is.na(Protein_Label)) %>%
    arrange(desc(biomarker_score), desc(single_protein_auc)) %>%
    slice_head(n = min(20, nrow(.))) %>%
    mutate(label = factor(Protein_Label, levels = rev(unique(Protein_Label))))
  p_top <- ggplot(top_summary,
                  aes(biomarker_score, label,
                      color = comparison, size = dplyr::coalesce(single_protein_auc, 0.5))) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c(
      "analysis6a_pr_vs_pd"       = "#E69F00",
      "analysis6b_prcr_vs_sdpd"   = "#56B4E9"),
      labels = c("analysis6a_pr_vs_pd" = "PR vs PD",
                 "analysis6b_prcr_vs_sdpd" = "PR/CR vs SD/PD")) +
    scale_size_continuous(range = c(3, 8)) +
    labs(title    = "Analysis 6: Top scoring biomarker candidates",
         subtitle = paste0("Score = limma(p<0.05) + OS Cox(FDR<0.1) + PFS Cox(FDR<0.1) +",
                           " LASSO selected + stability(\u22650.6) + AUC(\u22650.70)"),
         x = "Composite biomarker score", y = NULL, color = "Comparison", size = "Single-protein AUC") +
    theme_minimal(base_size = 11) + theme(legend.position = "right")
  cairo_pdf(file.path(fig7, "analysis6_top_candidates_dotplot.pdf"), width = 10, height = 7)
  print(p_top); dev.off()
  png(file.path(fig7, "analysis6_top_candidates_dotplot.png"), width = 10, height = 7,
      units = "in", res = 300, bg = "white")
  print(p_top); dev.off()
}

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 08 | Saving checkpoint_08.RData ===")
save(
  ml_outputs,
  priority_all,
  top_scoring_all,
  file = file.path(paths$rdata, "checkpoint_08.RData")
)
message("08_analysis6_ml.R complete")
