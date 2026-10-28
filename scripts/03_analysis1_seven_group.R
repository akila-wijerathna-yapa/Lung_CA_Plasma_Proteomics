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

# ── Load checkpoint ───────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")
while (dev.cur() > 1) dev.off()   # clear any stale devices

# ── Local constants ───────────────────────────────────────────
FDR_CUTOFF    <- 0.1    # volcano significance threshold
LFC_THR       <- 0.5    # |log2FC| threshold for volcano
MIN_GROUP_SIZE <- 3     # minimum samples per group for limma
TOP_N_HEATMAP  <- 100   # top omnibus proteins in heatmap

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  ifelse(nzchar(out), out, "item")
}

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

run_omnibus_limma <- function(mat_samples_by_features, sample_df, group_col,
                              adjust_cols = character(), min_group_size = 3) {
  keep <- !is.na(sample_df[[group_col]])
  df   <- sample_df[keep, , drop = FALSE]
  grp  <- factor(df[[group_col]])
  grp_counts <- table(grp)
  if (length(grp_counts) < 3 || any(grp_counts < min_group_size)) return(NULL)
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  feature_keep <- sapply(seq_len(ncol(mat_sub)), function(j) {
    all(tapply(!is.na(mat_sub[, j]), grp, sum) >= min_group_size)
  })
  mat_sub <- mat_sub[, feature_keep, drop = FALSE]
  if (ncol(mat_sub) == 0) return(NULL)
  design_df <- data.frame(omnibus_group = grp)
  for (col in adjust_cols) {
    vals <- df[[col]]
    if (is.numeric(vals)) {
      if (length(unique(vals[!is.na(vals)])) > 1) design_df[[col]] <- vals
    } else {
      fac <- factor(vals)
      if (length(unique(fac[!is.na(fac)])) > 1) design_df[[col]] <- fac
    }
  }
  design <- model.matrix(~ 0 + omnibus_group + ., data = design_df)
  fit <- limma::lmFit(t(mat_sub), design)
  fit <- limma::eBayes(fit)
  limma::topTableF(fit, number = Inf, sort.by = "F") %>%
    tibble::rownames_to_column("nELISA_ID") %>%
    as_tibble() %>%
    transmute(nELISA_ID, f_stat = F, p_value = P.Value, adj_p_value = adj.P.Val)
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
  design    <- model.matrix(~ ., data = design_df)
  coef_name <- grep("^primary_group", colnames(design), value = TRUE)[1]
  fit <- limma::lmFit(t(mat_sub), design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("nELISA_ID") %>%
    as_tibble() %>%
    transmute(nELISA_ID, logFC, t_stat = t, p_value = P.Value,
              adj_p_value = adj.P.Val, avg_expr = AveExpr)
  means    <- lapply(group_levels, function(gl) colMeans(mat_sub[grp == gl, , drop = FALSE], na.rm = TRUE))
  n_nonmiss <- lapply(group_levels, function(gl) colSums(!is.na(mat_sub[grp == gl, , drop = FALSE])))
  tt %>% mutate(
    group_a      = group_levels[1],
    group_b      = group_levels[2],
    mean_group_a = unname(means[[1]][match(nELISA_ID, names(means[[1]]))]),
    mean_group_b = unname(means[[2]][match(nELISA_ID, names(means[[2]]))]),
    n_group_a    = unname(n_nonmiss[[1]][match(nELISA_ID, names(n_nonmiss[[1]]))]),
    n_group_b    = unname(n_nonmiss[[2]][match(nELISA_ID, names(n_nonmiss[[2]]))])
  )
}

akila_volcano_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text          = element_text(color = "black", size = 12),
      axis.title         = element_text(color = "black", size = 14, face = "bold"),
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.background  = element_rect(color = "black", fill = "white", linewidth = 0.5),
      legend.title       = element_text(face = "bold", size = 10),
      legend.text        = element_text(size = 9),
      plot.margin        = margin(68, 56, 46, 56)
    )
}

make_volcano <- function(df, title, path_stub,
                         lfc_thr = LFC_THR, padj_thr = FDR_CUTOFF,
                         label_col = "Protein_Label",
                         group_a_label = NULL, group_b_label = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
  group_a_label <- group_a_label %||% unique(na.omit(df$group_a))[1] %||% "Group A"
  group_b_label <- group_b_label %||% unique(na.omit(df$group_b))[1] %||% "Group B"
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  group_a_label <- group_a_label %||% "Group A"
  group_b_label <- group_b_label %||% "Group B"
  x_label       <- paste0("Fold change log2(", group_b_label, " / ", group_a_label, ")")
  secondary_lbl <- paste0("Secondary (FDR<=", formatC(padj_thr, format = "f", digits = 2), ")")
  left_label    <- paste("Higher in", group_a_label)
  right_label   <- paste("Higher in", group_b_label)

  plot_df <- df %>%
    mutate(
      neg_log10 = -log10(pmax(adj_p_value, 1e-300)),
      status = case_when(
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC >=  lfc_thr ~ right_label,
        !is.na(adj_p_value) & adj_p_value < padj_thr & logFC <= -lfc_thr ~ left_label,
        TRUE ~ "Not significant"
      ),
      sig_level = case_when(
        status == "Not significant" ~ "Not significant",
        !is.na(adj_p_value) & adj_p_value <= 0.05 ~ "Primary (FDR<=0.05)",
        TRUE ~ secondary_lbl
      )
    )

  ann <- plot_df %>%
    filter(status != "Not significant") %>%
    arrange(adj_p_value, desc(abs(logFC))) %>%
    slice_head(n = 18)
  if (nrow(ann) == 0) {
    ann <- plot_df %>% arrange(adj_p_value, desc(abs(logFC))) %>% slice_head(n = 18)
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

  p <- ggplot(plot_df, aes(x = logFC, y = neg_log10)) +
    geom_point(aes(color = status, shape = sig_level), size = 2.4, alpha = 0.82, stroke = 0.7) +
    geom_vline(xintercept = c(-lfc_thr, 0, lfc_thr),
               linetype = c("dashed", "solid", "dashed"),
               linewidth = c(0.6, 0.4, 0.6), color = "black") +
    geom_hline(yintercept = y_thr, linetype = "dashed", linewidth = 0.6, color = "black") +
    annotate("text", x = x_panel_left - 0.18 * x_pad, y = y_label,
             label = paste0("Higher in\n", group_a_label),
             hjust = 0, vjust = 0, size = 4.8, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_panel_right + 0.18 * x_pad, y = y_label,
             label = paste0("Higher in\n", group_b_label),
             hjust = 1, vjust = 0, size = 4.8, fontface = "bold", color = "#B2182B") +
    annotate("text", x = x_rng[1] + 0.02 * diff(x_rng), y = y_thr,
             label = paste0(" P Adj. = ", formatC(padj_thr, format = "f", digits = 2)),
             hjust = 0, vjust = -0.6, size = 3.8, fontface = "italic", color = "black") +
    scale_color_manual(values = c(
      setNames("#2166AC", left_label),
      setNames("#B2182B", right_label),
      "Not significant" = "grey85"
    ), drop = FALSE) +
    scale_shape_manual(values = c(
      "Not significant" = 1,
      "Primary (FDR<=0.05)" = 19,
      setNames(17, secondary_lbl)
    ), breaks = c("Not significant", "Primary (FDR<=0.05)", secondary_lbl), drop = FALSE) +
    ggrepel::geom_text_repel(
      data = ann, aes(label = .data[[label_col]]),
      size = 3.1, fontface = "bold", max.overlaps = 30,
      box.padding = 0.5, point.padding = 0.2,
      segment.color = "grey60", segment.size = 0.2, show.legend = FALSE
    ) +
    labs(title = title, x = x_label, y = "P value (-log10)",
         color = NULL, shape = "Significance") +
    guides(color = "none", shape = guide_legend(override.aes = list(color = "black"))) +
    coord_cartesian(xlim = c(x_panel_left, x_panel_right),
                    ylim = c(0, y_panel_top), clip = "off") +
    akila_volcano_theme()

  cairo_pdf(paste0(path_stub, ".pdf"), width = 8.6, height = 8.2)
  print(p)
  dev.off()
  png(paste0(path_stub, ".png"), width = 8.6, height = 8.2, units = "in", res = 300, bg = "white")
  print(p)
  dev.off()
}

# ── Analysis 1 data frame ─────────────────────────────────────
message("\n=== 03 | Analysis 1 setup ===")
analysis1_df <- pre_matched   # baseline PRE samples matched to clinical

analysis1_counts <- analysis1_df %>% count(new_group, name = "n")
safe_write_csv(analysis1_counts,
               file.path(paths$tables, "02_analysis1_seven_group", "analysis1_group_counts.csv"))
message("Group counts:")
print(analysis1_counts)

# ── Omnibus F-test ────────────────────────────────────────────
message("\n=== 03 | Omnibus F-test ===")
analysis1_omnibus <- run_omnibus_limma(
  expr_log2_rep, analysis1_df, "new_group",
  adjust_cols    = c("plate_barcode"),
  min_group_size = MIN_GROUP_SIZE
)
analysis1_omnibus <- add_feature_annotation(analysis1_omnibus, feature_meta)
safe_write_csv(analysis1_omnibus,
               file.path(paths$tables, "02_analysis1_seven_group", "analysis1_omnibus_f_test.csv"))
message("Omnibus F-test: ", nrow(analysis1_omnibus), " proteins tested")
message("FDR<0.1: ", sum(analysis1_omnibus$adj_p_value < 0.1, na.rm = TRUE))

# ── Top 100 heatmap ───────────────────────────────────────────
message("\n=== 03 | Top 100 heatmap ===")
top100_ids <- analysis1_omnibus %>%
  arrange(adj_p_value, desc(f_stat)) %>%
  slice_head(n = min(TOP_N_HEATMAP, nrow(.))) %>%
  pull(nELISA_ID)

if (length(top100_ids) > 0) {
  heat_ann <- analysis1_df %>%
    select(sample_id, new_group, plate_barcode, response_binary) %>%
    tibble::column_to_rownames("sample_id")
  heat_mat    <- t(expr_log2_rep[analysis1_df$sample_id, top100_ids, drop = FALSE])
  heat_scaled <- scale_rows_na(heat_mat)

  stub_heat <- file.path(paths$figures, "02_analysis1_seven_group", "analysis1_top100_heatmap")

  cairo_pdf(paste0(stub_heat, ".pdf"), width = 12, height = 14)
  ph_a1 <- pheatmap::pheatmap(
    heat_scaled,
    annotation_col = heat_ann,
    show_colnames  = FALSE,
    show_rownames  = TRUE,
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
    color          = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    na_col         = "grey90",
    main           = paste0("Analysis 1: Top ", length(top100_ids), " omnibus proteins"),
    silent         = TRUE
  )
  grid::grid.newpage(); grid::grid.draw(ph_a1$gtable)
  dev.off()
  message("Heatmap saved")
}

# ── One-vs-rest contrasts ─────────────────────────────────────
message("\n=== 03 | One-vs-rest contrasts ===")
analysis1_ovr_results <- list()

for (grp in sort(unique(analysis1_df$new_group))) {
  df_tmp <- analysis1_df %>%
    mutate(group_vs_rest = ifelse(new_group == grp, grp, "Rest"))

  res <- run_binary_limma(
    expr_log2_rep, df_tmp, "group_vs_rest",
    group_levels   = c("Rest", grp),
    adjust_cols    = c("plate_barcode"),
    min_per_group  = MIN_GROUP_SIZE
  )
  if (is.null(res) || nrow(res) == 0) {
    message("  Skipped (insufficient n): ", grp)
    next
  }
  res   <- add_feature_annotation(res, feature_meta)
  label <- paste0("analysis1_", sanitize_filename(tolower(grp)), "_vs_rest")
  analysis1_ovr_results[[label]] <- res

  safe_write_csv(res,
    file.path(paths$tables, "02_analysis1_seven_group", paste0(label, ".csv")))

  make_volcano(
    res,
    title     = paste("Analysis 1:", grp, "vs rest"),
    path_stub = file.path(paths$figures, "02_analysis1_seven_group", label)
  )

  n_sig <- sum(res$adj_p_value < FDR_CUTOFF, na.rm = TRUE)
  message("  ", grp, " vs rest — proteins: ", nrow(res), " | FDR<0.1: ", n_sig)
}

message("One-vs-rest complete: ", length(analysis1_ovr_results), " contrasts saved")

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 03 | Saving checkpoint_03.RData ===")
save(
  analysis1_df,
  analysis1_omnibus,
  analysis1_ovr_results,
  file = file.path(paths$rdata, "checkpoint_03.RData")
)
message("03_analysis1_seven_group.R complete")
