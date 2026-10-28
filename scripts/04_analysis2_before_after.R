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
FDR_CUTOFF   <- 0.1
LFC_THR      <- 0.5
MIN_PAIRS    <- 6     # minimum pairs for surgery paired limma
TOP_SPAG     <- 6     # top proteins for spaghetti plot
TOP_PANEL    <- 4     # top proteins for paired panel plot
TOP_WILCOX   <- 6     # top proteins for individual Wilcoxon plots

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  ifelse(nzchar(out), out, "item")
}

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  left_join(df, feature_meta, by = "nELISA_ID")
}

run_paired_limma <- function(mat_samples_by_features, sample_df,
                             pair_id_col = "patient_id",
                             time_col    = "pre_post",
                             time_levels = c("pre", "post"),
                             min_pairs   = 5) {
  pair_counts <- sample_df %>%
    count(.data[[pair_id_col]], .data[[time_col]]) %>%
    tidyr::pivot_wider(names_from = .data[[time_col]], values_from = n, values_fill = 0)
  keep_ids <- pair_counts %>%
    filter(.data[[time_levels[1]]] > 0, .data[[time_levels[2]]] > 0) %>%
    pull(.data[[pair_id_col]])
  df <- sample_df %>%
    filter(.data[[pair_id_col]] %in% keep_ids) %>%
    mutate(timepoint = factor(.data[[time_col]], levels = time_levels))
  if (length(unique(df[[pair_id_col]])) < min_pairs) return(NULL)
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  feature_keep <- colSums(!is.na(mat_sub[df$timepoint == time_levels[1], , drop = FALSE])) >= min_pairs &
    colSums(!is.na(mat_sub[df$timepoint == time_levels[2], , drop = FALSE])) >= min_pairs
  mat_sub <- mat_sub[, feature_keep, drop = FALSE]
  if (ncol(mat_sub) == 0) return(NULL)
  design_df <- data.frame(
    pair_id   = factor(df[[pair_id_col]]),
    timepoint = factor(df$timepoint, levels = time_levels)
  )
  design    <- model.matrix(~ pair_id + timepoint, data = design_df)
  coef_name <- grep("^timepoint", colnames(design), value = TRUE)[1]
  fit <- limma::lmFit(t(mat_sub), design)
  fit <- limma::eBayes(fit)
  limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("nELISA_ID") %>%
    as_tibble() %>%
    transmute(nELISA_ID, logFC, t_stat = t, p_value = P.Value,
              adj_p_value = adj.P.Val, avg_expr = AveExpr) %>%
    mutate(group_a = time_levels[1], group_b = time_levels[2])
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
  group_a_label  <- group_a_label %||% unique(na.omit(df$group_a))[1] %||% "Group A"
  group_b_label  <- group_b_label %||% unique(na.omit(df$group_b))[1] %||% "Group B"
  `%||%`         <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  group_a_label  <- group_a_label %||% "Group A"
  group_b_label  <- group_b_label %||% "Group B"
  x_label        <- paste0("Fold change log2(", group_b_label, " / ", group_a_label, ")")
  secondary_lbl  <- paste0("Secondary (FDR<=", formatC(padj_thr, format = "f", digits = 2), ")")
  left_label     <- paste("Higher in", group_a_label)
  right_label    <- paste("Higher in", group_b_label)
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

save_ggplot <- function(p, path_stub, width, height) {
  # Direct device save for any ggplot object — avoids top-level echo issue
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height)
  print(p)
  dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white")
  print(p)
  dev.off()
}

# ── Build paired sample sets ──────────────────────────────────
message("\n=== 04 | Build paired sample sets ===")
paired_pre  <- pre_matched  %>% filter(patient_id %in% paired_ids)
paired_post <- post_matched %>% filter(patient_id %in% paired_ids)
paired_summary <- paired_clinical %>% count(new_group, name = "n_pairs")
safe_write_csv(paired_summary,
               file.path(paths$tables, "03_analysis2_before_after", "paired_group_counts.csv"))
message("Paired counts:"); print(paired_summary)

# ── Surgery paired limma ──────────────────────────────────────
message("\n=== 04 | Surgery paired limma ===")
surgery_pair_ids <- paired_clinical %>% filter(new_group == "Sx") %>% pull(patient_id)

surgery_long <- bind_rows(
  paired_pre  %>% filter(patient_id %in% surgery_pair_ids),
  paired_post %>% filter(patient_id %in% surgery_pair_ids)
) %>% mutate(pre_post = factor(pre_post, levels = c("pre", "post")))

analysis2_surgery <- run_paired_limma(
  expr_log2_rep, surgery_long,
  pair_id_col = "patient_id",
  time_col    = "pre_post",
  time_levels = c("pre", "post"),
  min_pairs   = MIN_PAIRS
)

if (!is.null(analysis2_surgery) && nrow(analysis2_surgery) > 0) {
  analysis2_surgery <- add_feature_annotation(analysis2_surgery, feature_meta)
  safe_write_csv(analysis2_surgery,
    file.path(paths$tables, "03_analysis2_before_after", "analysis2_surgery_paired_limma.csv"))
  message("Surgery paired: ", nrow(analysis2_surgery), " proteins | FDR<0.1: ",
          sum(analysis2_surgery$adj_p_value < FDR_CUTOFF, na.rm = TRUE))

  make_volcano(
    analysis2_surgery,
    title     = paste0("Analysis 2: Surgery paired PRE vs POST (n=", length(surgery_pair_ids), ")"),
    path_stub = file.path(paths$figures, "03_analysis2_before_after",
                          "analysis2_surgery_paired_volcano"),
    group_a_label = "PRE", group_b_label = "POST"
  )

  # Biology cross-reference for known surgery proteins
  surgery_biology_reference <- tibble::tribble(
    ~Gene_Symbol, ~cross_reference_theme,                         ~rationale_short,
    "THBS4",  "ECM remodelling / wound healing",        "Thrombospondin-family matrix protein linked to tissue repair.",
    "POSTN",  "ECM remodelling / wound healing",        "Periostin: classic stromal wound-healing ECM protein.",
    "CHI3L1", "Inflammation / tissue remodelling",      "CHI3L1: inflammatory tissue injury and repair responses.",
    "ANGPTL2","Inflammation / tissue remodelling",      "ANGPTL2: inflammatory stromal activation and remodelling.",
    "CCL8",   "Chemokine / immune-cell trafficking",    "CCL8: chemokine for leukocyte recruitment.",
    "CXCL3",  "Chemokine / neutrophil trafficking",     "CXCL3: neutrophil-directed chemokine, postoperative signalling.",
    "CXCL5",  "Chemokine / neutrophil trafficking",     "CXCL5: neutrophil chemokine, acute inflammatory recruitment.",
    "CXCL6",  "Chemokine / neutrophil trafficking",     "CXCL6: granulocyte trafficking and inflammatory activation.",
    "GP6",    "Platelet activation / coagulation",      "GP6: platelet activation receptor, perioperative biology.",
    "HSPB1",  "Cellular stress response",               "HSPB1: stress-response chaperone, acute tissue stress.",
    "HPSE",   "Matrix remodelling / injury response",   "Heparanase: ECM turnover and injury-associated remodelling.",
    "CD69",   "Early immune activation",                "CD69: early activation marker, postoperative immune response."
  )
  surgery_biology_tbl <- analysis2_surgery %>%
    mutate(direction_post_vs_pre = ifelse(logFC > 0, "Higher POST", "Lower POST")) %>%
    inner_join(surgery_biology_reference, by = "Gene_Symbol") %>%
    arrange(adj_p_value, p_value) %>%
    transmute(Protein_Name, Gene_Symbol, direction_post_vs_pre,
              p_value, adj_p_value, cross_reference_theme, rationale_short)
  safe_write_csv(surgery_biology_tbl,
    file.path(paths$tables, "03_analysis2_before_after",
              "analysis2_surgery_biological_cross_reference.csv"))
} else {
  message("Surgery paired limma: insufficient pairs (need >= ", MIN_PAIRS, ")")
  analysis2_surgery <- tibble()
}

# ── Surgery spaghetti and paired plots ────────────────────────
message("\n=== 04 | Surgery spaghetti + paired Wilcoxon plots ===")
if (!is.null(analysis2_surgery) && nrow(analysis2_surgery) > 0 &&
    length(surgery_pair_ids) >= 4) {

  top_surgery <- analysis2_surgery %>%
    arrange(adj_p_value) %>%
    slice_head(n = min(TOP_SPAG, nrow(.))) %>%
    pull(nELISA_ID)

  spaghetti_base <- bind_rows(
    paired_pre  %>% filter(patient_id %in% surgery_pair_ids) %>% mutate(timepoint = "PRE"),
    paired_post %>% filter(patient_id %in% surgery_pair_ids) %>% mutate(timepoint = "POST")
  ) %>% select(sample_id, patient_id, timepoint)

  spaghetti_df <- spaghetti_base %>%
    left_join(
      as.data.frame(expr_log2_rep[spaghetti_base$sample_id, top_surgery, drop = FALSE]) %>%
        tibble::rownames_to_column("sample_id") %>%
        pivot_longer(-sample_id, names_to = "nELISA_ID", values_to = "log2_value"),
      by = "sample_id"
    ) %>%
    left_join(
      as.data.frame(expr_pgml_rep[spaghetti_base$sample_id, top_surgery, drop = FALSE]) %>%
        tibble::rownames_to_column("sample_id") %>%
        pivot_longer(-sample_id, names_to = "nELISA_ID", values_to = "pgml_value"),
      by = c("sample_id", "nELISA_ID")
    ) %>%
    left_join(feature_meta %>% select(nELISA_ID, Protein_Label), by = "nELISA_ID")

  # Spaghetti: top 6 proteins
  p_spag <- ggplot(spaghetti_df %>% filter(!is.na(log2_value)),
                   aes(timepoint, log2_value, group = patient_id)) +
    geom_line(alpha = 0.5, color = "#6a51a3") +
    geom_point(size = 1.7, color = "#54278f") +
    facet_wrap(~Protein_Label, scales = "free_y") +
    labs(title = "Top surgery paired proteins", x = NULL, y = "log2 concentration") +
    theme_minimal(base_size = 11)
  save_ggplot(p_spag,
    file.path(paths$figures, "03_analysis2_before_after", "analysis2_surgery_top_spaghetti"),
    width = 11, height = 7.5)

  # Paired panel plots (top 4 proteins, Wilcoxon p-values)
  panel_targets <- analysis2_surgery %>%
    arrange(adj_p_value) %>%
    slice_head(n = min(TOP_PANEL, nrow(.))) %>%
    pull(nELISA_ID)
  panel_df <- spaghetti_df %>% filter(nELISA_ID %in% panel_targets)
  panel_pvals <- bind_rows(lapply(split(panel_df, panel_df$nELISA_ID), function(d) {
    wide  <- d %>% select(patient_id, timepoint, log2_value) %>%
      tidyr::pivot_wider(names_from = timepoint, values_from = log2_value)
    p_val <- tryCatch(stats::wilcox.test(wide$PRE, wide$POST, paired = TRUE)$p.value,
                      error = function(e) NA_real_)
    tibble(nELISA_ID = unique(d$nELISA_ID), Protein_Label = unique(d$Protein_Label),
           paired_wilcoxon_p = p_val)
  }))
  panel_df <- panel_df %>%
    left_join(panel_pvals, by = c("nELISA_ID", "Protein_Label")) %>%
    mutate(
      facet_label = paste0(Protein_Label, "\nPaired Wilcoxon p=",
                           scales::number(paired_wilcoxon_p, accuracy = 0.0001)),
      timepoint   = factor(timepoint, levels = c("PRE", "POST"))
    )
  p_panel <- ggplot(panel_df %>% filter(!is.na(log2_value)),
                    aes(timepoint, log2_value, group = patient_id)) +
    geom_boxplot(width = 0.55, fill = "#deebf7", color = "#08519c", outlier.shape = NA) +
    geom_line(alpha = 0.45, color = "grey45") +
    geom_point(size = 1.6, color = "#08519c") +
    facet_wrap(~facet_label, scales = "free_y") +
    labs(title = "Top surgery paired proteins: paired panels", x = NULL, y = "log2 concentration") +
    theme_bw(base_size = 11)
  save_ggplot(p_panel,
    file.path(paths$figures, "03_analysis2_before_after", "analysis2_surgery_top_paired_panels"),
    width = 10.5, height = 7.2)

  # Top 6 individual Wilcoxon plots (pg/mL)
  wilcox_targets <- analysis2_surgery %>%
    arrange(adj_p_value, p_value, desc(abs(logFC))) %>%
    slice_head(n = min(TOP_WILCOX, nrow(.))) %>%
    pull(nELISA_ID)
  wilcox_df <- spaghetti_df %>%
    filter(nELISA_ID %in% wilcox_targets) %>%
    left_join(
      feature_meta %>% transmute(nELISA_ID,
        feature_display = dplyr::coalesce(Protein_Name, Protein_Label, nELISA_ID)),
      by = "nELISA_ID"
    )
  wilcox_pvals <- bind_rows(lapply(split(wilcox_df, wilcox_df$nELISA_ID), function(d) {
    wide  <- d %>% select(patient_id, timepoint, pgml_value) %>%
      tidyr::pivot_wider(names_from = timepoint, values_from = pgml_value) %>%
      filter(is.finite(PRE), is.finite(POST))
    p_val <- tryCatch(
      stats::wilcox.test(wide$PRE, wide$POST, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_)
    tibble(nELISA_ID = unique(d$nELISA_ID), paired_wilcoxon_p = p_val)
  }))
  wilcox_df <- wilcox_df %>%
    filter(is.finite(pgml_value)) %>%
    left_join(wilcox_pvals, by = "nELISA_ID") %>%
    mutate(
      timepoint     = factor(timepoint, levels = c("PRE", "POST")),
      feature_display = dplyr::coalesce(feature_display, Protein_Label, nELISA_ID),
      facet_label   = paste0(feature_display, "\nPaired Wilcoxon p=",
                             scales::number(paired_wilcoxon_p, accuracy = 0.0001))
    )
  safe_write_csv(
    wilcox_df %>% distinct(nELISA_ID, feature_display, paired_wilcoxon_p) %>% arrange(paired_wilcoxon_p),
    file.path(paths$tables, "03_analysis2_before_after",
              "analysis2_surgery_top6_paired_wilcox_summary.csv"))

  # Multi-panel overview (all 6 together)
  p_wilcox_multi <- ggplot(wilcox_df, aes(timepoint, pgml_value)) +
    geom_boxplot(aes(fill = timepoint), alpha = 0.35, outlier.shape = NA,
                 width = 0.55, color = "grey30") +
    geom_line(aes(group = patient_id), alpha = 0.45, color = "grey65", linewidth = 0.35) +
    geom_point(aes(color = timepoint), size = 1.9) +
    scale_fill_manual(values  = c(PRE = "#f4b6b2", POST = "#8dd3e7")) +
    scale_color_manual(values = c(PRE = "#d95f5f", POST = "#1f9fb5")) +
    facet_wrap(~facet_label, ncol = 3, scales = "free_y") +
    labs(title = "Pre/Post Treatment Protein Changes (Paired Wilcoxon)",
         x = NULL, y = "pg/mL") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          plot.title  = element_text(size = 11, face = "bold"),
          strip.text  = element_text(face = "bold", size = 9))
  save_ggplot(p_wilcox_multi,
    file.path(paths$figures, "03_analysis2_before_after",
              "analysis2_surgery_top6_paired_wilcox"),
    width = 12.2, height = 8.6)

  # Individual protein plots
  for (feature_i in wilcox_targets) {
    df_i <- wilcox_df %>% filter(nELISA_ID == feature_i)
    if (!nrow(df_i)) next
    fname_i <- unique(df_i$feature_display)[1]
    p_i <- ggplot(df_i, aes(timepoint, pgml_value)) +
      geom_boxplot(aes(fill = timepoint), alpha = 0.35, outlier.shape = NA,
                   width = 0.55, color = "grey30") +
      geom_line(aes(group = patient_id), alpha = 0.5, color = "grey65", linewidth = 0.4) +
      geom_point(aes(color = timepoint), size = 2.1) +
      ggrepel::geom_text_repel(aes(label = sample_id), size = 2.4,
                               max.overlaps = Inf, segment.alpha = 0.2) +
      scale_fill_manual(values  = c(PRE = "#f4b6b2", POST = "#8dd3e7")) +
      scale_color_manual(values = c(PRE = "#d95f5f", POST = "#1f9fb5")) +
      labs(title    = paste0("Surgery paired PRE vs POST — ", fname_i),
           subtitle = paste0("Paired Wilcoxon p = ",
                             scales::number(unique(df_i$paired_wilcoxon_p)[1], accuracy = 0.0001)),
           x = NULL, y = "pg/mL") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none",
            plot.title = element_text(size = 11, face = "bold"))
    save_ggplot(p_i,
      file.path(paths$figures, "03_analysis2_before_after",
                paste0("analysis2_surgery_paired_wilcox_", sanitize_filename(fname_i))),
      width = 7.4, height = 5.6)
  }
  message("Surgery plots saved")
}

# ── IO-containing pairs: descriptive delta logFC ──────────────
message("\n=== 04 | IO-containing delta logFC (descriptive) ===")
delta_io_tbl <- tibble()
io_containing_ids <- paired_clinical %>%
  filter(treatment_class %in% c("Immunotherapy", "Chemo + IO",
                                 "Chemo + IO + RT", "RT + IO (sequential)", "IO + RT")) %>%
  pull(patient_id)

if (length(io_containing_ids) >= 4) {
  pre_io    <- paired_pre  %>% filter(patient_id %in% io_containing_ids) %>% arrange(patient_id)
  post_io   <- paired_post %>% filter(patient_id %in% io_containing_ids) %>% arrange(patient_id)
  shared_io <- intersect(pre_io$patient_id, post_io$patient_id)
  pre_io    <- pre_io  %>% filter(patient_id %in% shared_io) %>% arrange(patient_id)
  post_io   <- post_io %>% filter(patient_id %in% shared_io) %>% arrange(patient_id)

  delta_io     <- expr_log2_rep[post_io$sample_id, , drop = FALSE] -
                  expr_log2_rep[pre_io$sample_id,  , drop = FALSE]
  delta_io_tbl <- tibble(
    nELISA_ID    = colnames(delta_io),
    mean_logFC   = colMeans(delta_io, na.rm = TRUE),
    median_logFC = apply(delta_io, 2, median, na.rm = TRUE),
    n_pairs      = colSums(!is.na(delta_io))
  ) %>%
    add_feature_annotation(feature_meta) %>%
    arrange(desc(abs(mean_logFC)))
  safe_write_csv(delta_io_tbl,
    file.path(paths$tables, "03_analysis2_before_after",
              "analysis2_io_containing_logfc_only.csv"))

  io_comp_summary <- pre_io %>%
    count(treatment_class, new_group, response_binary, name = "n_pairs") %>%
    arrange(desc(n_pairs), treatment_class, response_binary)
  io_comp_detail <- pre_io %>%
    transmute(patient_id, treatment_class, new_group, response_4cat,
              response_binary, pfs_censor_class, recurrence_basis,
              pre_sample_id = sample_id, pre_plate_barcode = plate_barcode) %>%
    left_join(
      post_io %>% transmute(patient_id, post_sample_id = sample_id,
                             post_plate_barcode = plate_barcode),
      by = "patient_id"
    ) %>%
    arrange(treatment_class, patient_id)
  safe_write_csv(io_comp_summary,
    file.path(paths$tables, "03_analysis2_before_after",
              "analysis2_io_containing_subclass_composition_summary.csv"))
  safe_write_csv(io_comp_detail,
    file.path(paths$tables, "03_analysis2_before_after",
              "analysis2_io_containing_subclass_composition_detail.csv"))

  # Ranked bar: top absolute delta logFC
  top_io <- delta_io_tbl %>%
    arrange(desc(abs(mean_logFC))) %>%
    slice_head(n = min(25, nrow(.))) %>%
    mutate(direction  = ifelse(mean_logFC >= 0, "Higher POST", "Higher PRE"),
           label_order = factor(Protein_Label, levels = rev(Protein_Label)))
  p_io_bar <- ggplot(top_io, aes(x = mean_logFC, y = label_order, fill = direction)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("Higher POST" = "#B2182B", "Higher PRE" = "#2166AC")) +
    labs(title = paste0("Analysis 2: IO-containing paired logFC only (n=", length(shared_io), ")"),
         x = "Mean log2 fold change (POST - PRE)", y = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  save_ggplot(p_io_bar,
    file.path(paths$figures, "03_analysis2_before_after", "analysis2_io_containing_logfc"),
    width = 8.5, height = 7)

  message("IO delta logFC: ", nrow(delta_io_tbl), " proteins | pairs n=", length(shared_io))
} else {
  message("IO-containing delta: fewer than 4 pairs — skipped")
}

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 04 | Saving checkpoint_04.RData ===")
save(
  paired_pre, paired_post,
  surgery_pair_ids,
  analysis2_surgery,
  delta_io_tbl,
  file = file.path(paths$rdata, "checkpoint_04.RData")
)
message("04_analysis2_before_after.R complete")
