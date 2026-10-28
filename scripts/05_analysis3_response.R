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
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
FDR_CUTOFF        <- 0.1
LFC_THR           <- 0.5
MIN_PER_GROUP_PRI <- 8    # primary: PR vs PD (smaller classes)
MIN_PER_GROUP_SEC <- 15   # secondary: PR/CR vs SD/PD (larger classes)

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  left_join(df, feature_meta, by = "nELISA_ID")
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
  means     <- lapply(group_levels, function(gl) colMeans(mat_sub[grp == gl, , drop = FALSE], na.rm = TRUE))
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
      plot.margin       = margin(68, 56, 46, 56)
    )
}

make_volcano <- function(df, title, path_stub,
                         lfc_thr = LFC_THR, padj_thr = FDR_CUTOFF,
                         label_col = "Protein_Label",
                         group_a_label = NULL, group_b_label = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
  `%||%`        <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  group_a_label <- group_a_label %||% unique(na.omit(df$group_a))[1] %||% "Group A"
  group_b_label <- group_b_label %||% unique(na.omit(df$group_b))[1] %||% "Group B"
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
  if (nrow(ann) == 0)
    ann <- plot_df %>% arrange(adj_p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  x_rng         <- range(plot_df$logFC, na.rm = TRUE)
  y_thr         <- -log10(padj_thr)
  y_max         <- max(plot_df$neg_log10, y_thr, na.rm = TRUE)
  x_pad         <- max(0.35, diff(x_rng) * 0.22)
  y_pad         <- max(0.18, y_max * 0.16)
  x_panel_left  <- x_rng[1] - 0.45 * x_pad
  x_panel_right <- x_rng[2] + 0.45 * x_pad
  y_panel_top   <- y_max + 0.12 * y_pad
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
             hjust = 0, vjust = -0.6, size = 3.8, fontface = "italic") +
    scale_color_manual(values = c(
      setNames("#2166AC", left_label), setNames("#B2182B", right_label),
      "Not significant" = "grey85"), drop = FALSE) +
    scale_shape_manual(values = c("Not significant" = 1, "Primary (FDR<=0.05)" = 19,
                                  setNames(17, secondary_lbl)),
                       breaks = c("Not significant", "Primary (FDR<=0.05)", secondary_lbl),
                       drop = FALSE) +
    ggrepel::geom_text_repel(data = ann, aes(label = .data[[label_col]]),
                             size = 3.1, fontface = "bold", max.overlaps = 30,
                             box.padding = 0.5, point.padding = 0.2,
                             segment.color = "grey60", segment.size = 0.2,
                             show.legend = FALSE) +
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

# ── Build analysis 3 data frame ───────────────────────────────
message("\n=== 05 | Analysis 3 setup ===")
analysis3_df <- analysis1_df %>%
  mutate(
    response_prpd   = ifelse(response_4cat %in% c("PR", "PD"), response_4cat, NA_character_),
    response_binary2 = response_binary
  )

n_pr  <- sum(analysis3_df$response_prpd == "PR",   na.rm = TRUE)
n_pd  <- sum(analysis3_df$response_prpd == "PD",   na.rm = TRUE)
n_prcr <- sum(analysis3_df$response_binary2 == "PR/CR", na.rm = TRUE)
n_sdpd <- sum(analysis3_df$response_binary2 == "SD/PD", na.rm = TRUE)
message("Primary   — PR: ", n_pr,   " | PD: ", n_pd)
message("Secondary — PR/CR: ", n_prcr, " | SD/PD: ", n_sdpd)

response_counts <- tibble(
  analysis     = c("primary", "primary", "secondary", "secondary"),
  group        = c("PR", "PD", "PR/CR", "SD/PD"),
  n            = c(n_pr, n_pd, n_prcr, n_sdpd)
)
safe_write_csv(response_counts,
               file.path(paths$tables, "04_analysis3_pr_vs_pd", "analysis3_response_counts.csv"))

# ── Primary: PR vs PD ─────────────────────────────────────────
message("\n=== 05 | Primary: PR vs PD ===")
analysis3_primary <- run_binary_limma(
  expr_log2_rep, analysis3_df, "response_prpd",
  group_levels  = c("PD", "PR"),
  adjust_cols   = c("plate_barcode", "new_group"),
  min_per_group = MIN_PER_GROUP_PRI
)

if (!is.null(analysis3_primary) && nrow(analysis3_primary) > 0) {
  analysis3_primary <- add_feature_annotation(analysis3_primary, feature_meta)
  safe_write_csv(analysis3_primary,
    file.path(paths$tables, "04_analysis3_pr_vs_pd", "analysis3_primary_pr_vs_pd.csv"))
  make_volcano(
    analysis3_primary,
    title     = "Analysis 3 primary: PR vs PD",
    path_stub = file.path(paths$figures, "04_analysis3_pr_vs_pd",
                          "analysis3_primary_pr_vs_pd")
  )
  message("Primary PR vs PD: ", nrow(analysis3_primary),
          " proteins | FDR<0.1: ", sum(analysis3_primary$adj_p_value < FDR_CUTOFF, na.rm = TRUE))
} else {
  message("Primary PR vs PD: insufficient data (need >= ", MIN_PER_GROUP_PRI, " per group)")
  analysis3_primary <- tibble()
}

# ── Secondary: PR/CR vs SD/PD ─────────────────────────────────
message("\n=== 05 | Secondary: PR/CR vs SD/PD ===")
analysis3_secondary <- run_binary_limma(
  expr_log2_rep, analysis3_df, "response_binary2",
  group_levels  = c("SD/PD", "PR/CR"),
  adjust_cols   = c("plate_barcode", "new_group"),
  min_per_group = MIN_PER_GROUP_SEC
)

if (!is.null(analysis3_secondary) && nrow(analysis3_secondary) > 0) {
  analysis3_secondary <- add_feature_annotation(analysis3_secondary, feature_meta)
  safe_write_csv(analysis3_secondary,
    file.path(paths$tables, "04_analysis3_pr_vs_pd", "analysis3_secondary_prcr_vs_sdpd.csv"))
  make_volcano(
    analysis3_secondary,
    title     = "Analysis 3 secondary: PR/CR vs SD/PD",
    path_stub = file.path(paths$figures, "04_analysis3_pr_vs_pd",
                          "analysis3_secondary_prcr_vs_sdpd")
  )
  message("Secondary PR/CR vs SD/PD: ", nrow(analysis3_secondary),
          " proteins | FDR<0.1: ", sum(analysis3_secondary$adj_p_value < FDR_CUTOFF, na.rm = TRUE))
} else {
  message("Secondary PR/CR vs SD/PD: insufficient data (need >= ", MIN_PER_GROUP_SEC, " per group)")
  analysis3_secondary <- tibble()
}

# ── Overlap: proteins significant in both ─────────────────────
message("\n=== 05 | Overlap primary + secondary ===")
if (nrow(analysis3_primary) > 0 && nrow(analysis3_secondary) > 0) {
  overlap_tbl <- full_join(
    analysis3_primary %>%
      filter(adj_p_value < FDR_CUTOFF) %>%
      transmute(nELISA_ID, Protein_Label, Protein_Name, Gene_Symbol,
                primary_logFC = logFC, primary_p = p_value, primary_fdr = adj_p_value),
    analysis3_secondary %>%
      filter(adj_p_value < FDR_CUTOFF) %>%
      transmute(nELISA_ID, Protein_Label, Protein_Name, Gene_Symbol,
                secondary_logFC = logFC, secondary_p = p_value, secondary_fdr = adj_p_value),
    by = c("nELISA_ID", "Protein_Label", "Protein_Name", "Gene_Symbol")
  ) %>%
    mutate(
      shared_hit     = !is.na(primary_fdr) & !is.na(secondary_fdr),
      overlap_basis  = "Strict overlap: FDR<0.10 in both response analyses",
      same_direction = sign(primary_logFC) == sign(secondary_logFC)
    ) %>%
    arrange(desc(shared_hit), dplyr::coalesce(primary_fdr, secondary_fdr))

  # Fall back to nominal p<0.05 if no strict FDR hits
  if (nrow(overlap_tbl) == 0 || !any(overlap_tbl$shared_hit, na.rm = TRUE)) {
    overlap_tbl <- inner_join(
      analysis3_primary %>%
        filter(p_value < 0.05) %>%
        transmute(nELISA_ID, Protein_Label, Protein_Name, Gene_Symbol,
                  primary_logFC = logFC, primary_p = p_value, primary_fdr = adj_p_value),
      analysis3_secondary %>%
        filter(p_value < 0.05) %>%
        transmute(nELISA_ID, Protein_Label, Protein_Name, Gene_Symbol,
                  secondary_logFC = logFC, secondary_p = p_value, secondary_fdr = adj_p_value),
      by = c("nELISA_ID", "Protein_Label", "Protein_Name", "Gene_Symbol")
    ) %>%
      mutate(
        shared_hit     = FALSE,
        overlap_basis  = "Exploratory overlap: nominal p<0.05 in both response analyses",
        same_direction = sign(primary_logFC) == sign(secondary_logFC)
      ) %>%
      arrange(primary_p, secondary_p)
  }

  safe_write_csv(overlap_tbl,
    file.path(paths$tables, "04_analysis3_pr_vs_pd",
              "analysis3_overlap_primary_secondary.csv"))
  message("Overlap table: ", nrow(overlap_tbl), " proteins | strict shared hits: ",
          sum(overlap_tbl$shared_hit, na.rm = TRUE))
}

# ── Exploratory pass: plate_barcode only (no group adjustment) ─
# Rationale: new_group adjustment consumes 6 df and removes signal
# correlated with treatment. This pass trades confounding control
# for statistical power. Results are labelled exploratory.
message("\n=== 05 | Exploratory: plate_barcode adjustment only ===")

analysis3_primary_expl <- run_binary_limma(
  expr_log2_rep, analysis3_df, "response_prpd",
  group_levels  = c("PD", "PR"),
  adjust_cols   = c("plate_barcode"),
  min_per_group = MIN_PER_GROUP_PRI
)
if (!is.null(analysis3_primary_expl) && nrow(analysis3_primary_expl) > 0) {
  analysis3_primary_expl <- add_feature_annotation(analysis3_primary_expl, feature_meta)
  safe_write_csv(analysis3_primary_expl,
    file.path(paths$tables, "04_analysis3_pr_vs_pd",
              "analysis3_primary_pr_vs_pd_EXPLORATORY_plate_only.csv"))
  make_volcano(
    analysis3_primary_expl,
    title     = "Analysis 3 primary EXPLORATORY (plate adj. only): PR vs PD",
    path_stub = file.path(paths$figures, "04_analysis3_pr_vs_pd",
                          "analysis3_primary_pr_vs_pd_exploratory")
  )
  message("Exploratory primary — FDR<0.1: ",
          sum(analysis3_primary_expl$adj_p_value < FDR_CUTOFF, na.rm = TRUE),
          " | nominal p<0.05: ",
          sum(analysis3_primary_expl$p_value < 0.05, na.rm = TRUE))
}

analysis3_secondary_expl <- run_binary_limma(
  expr_log2_rep, analysis3_df, "response_binary2",
  group_levels  = c("SD/PD", "PR/CR"),
  adjust_cols   = c("plate_barcode"),
  min_per_group = MIN_PER_GROUP_SEC
)
if (!is.null(analysis3_secondary_expl) && nrow(analysis3_secondary_expl) > 0) {
  analysis3_secondary_expl <- add_feature_annotation(analysis3_secondary_expl, feature_meta)
  safe_write_csv(analysis3_secondary_expl,
    file.path(paths$tables, "04_analysis3_pr_vs_pd",
              "analysis3_secondary_prcr_vs_sdpd_EXPLORATORY_plate_only.csv"))
  make_volcano(
    analysis3_secondary_expl,
    title     = "Analysis 3 secondary EXPLORATORY (plate adj. only): PR/CR vs SD/PD",
    path_stub = file.path(paths$figures, "04_analysis3_pr_vs_pd",
                          "analysis3_secondary_prcr_vs_sdpd_exploratory")
  )
  message("Exploratory secondary — FDR<0.1: ",
          sum(analysis3_secondary_expl$adj_p_value < FDR_CUTOFF, na.rm = TRUE),
          " | nominal p<0.05: ",
          sum(analysis3_secondary_expl$p_value < 0.05, na.rm = TRUE))
}

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 05 | Saving checkpoint_05.RData ===")
save(
  analysis3_df,
  analysis3_primary,
  analysis3_secondary,
  analysis3_primary_expl,
  analysis3_secondary_expl,
  file = file.path(paths$rdata, "checkpoint_05.RData")
)
message("05_analysis3_response.R complete")
