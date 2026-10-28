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

# ── Load checkpoint ────────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")
while (dev.cur() > 1) dev.off()

# ── Output directories (self-contained, not reliant on 00_setup) ─
tbl_dir <- file.path(paths$tables,  "10_exploratory_subgroups")
fig_dir <- file.path(paths$figures, "10_exploratory_subgroups")
tbl_pfs <- file.path(paths$tables,  "11_pfs_risk")
fig_pfs <- file.path(paths$figures, "11_pfs_risk")
invisible(lapply(c(tbl_dir, fig_dir, tbl_pfs, fig_pfs), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

# ── Local constants ────────────────────────────────────────────
PVAL_CUTOFF   <- 0.1   # nominal p-value threshold (NOT FDR-adjusted)
LFC_THR       <- 0.5
MIN_PER_GROUP <- 5     # minimum n per arm for binary limma contrasts
MIN_PAIRS     <- 5     # minimum pairs for surgery paired limma
FOREST_N      <- 20
KM_TOP_N      <- 4
MIN_N_COX     <- 30
MIN_EVENTS_COX<- 10

# ── Local utility functions ────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  ifelse(nzchar(out), out, "item")
}

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  left_join(df, feature_meta, by = "nELISA_ID") %>%
    mutate(Protein_Label_Name = dplyr::coalesce(Protein_Name, nELISA_ID))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

# ── Differential expression engines ────────────────────────────
run_binary_limma <- function(mat_samples_by_features, sample_df, group_col,
                             group_levels, adjust_cols = character(),
                             min_per_group = MIN_PER_GROUP) {
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

run_paired_limma <- function(mat_samples_by_features, sample_df,
                             pair_id_col = "patient_id",
                             time_col    = "pre_post",
                             time_levels = c("pre", "post"),
                             min_pairs   = MIN_PAIRS) {
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
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("nELISA_ID") %>%
    as_tibble() %>%
    transmute(nELISA_ID, logFC, t_stat = t, p_value = P.Value,
              adj_p_value = adj.P.Val, avg_expr = AveExpr) %>%
    mutate(group_a = time_levels[1], group_b = time_levels[2])
  n_pre  <- colSums(!is.na(mat_sub[df$timepoint == time_levels[1], , drop = FALSE]))
  n_post <- colSums(!is.na(mat_sub[df$timepoint == time_levels[2], , drop = FALSE]))
  tt %>% mutate(
    n_group_a = unname(n_pre[match(nELISA_ID, names(n_pre))]),
    n_group_b = unname(n_post[match(nELISA_ID, names(n_post))])
  )
}

# ── Plot theme + volcano + MA ─────────────────────────────────
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

save_ggplot <- function(p, path_stub, width, height) {
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height)
  print(p); dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

make_volcano <- function(df, title, path_stub,
                         lfc_thr = LFC_THR, p_thr = PVAL_CUTOFF,
                         label_col = "Protein_Label_Name",
                         group_a_label = NULL, group_b_label = NULL) {
  # NOTE: thresholding/coloring below uses NOMINAL p_value (uncorrected for
  # multiple testing). The true BH-adjusted value is still carried in the
  # data as adj_p_value for reference (see CSV output) but is NOT used here.
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  group_a_label <- group_a_label %||% (unique(na.omit(df$group_a))[1] %||% "Group A")
  group_b_label <- group_b_label %||% (unique(na.omit(df$group_b))[1] %||% "Group B")
  secondary_lbl <- paste0("Secondary (nominal p<=", formatC(p_thr, format = "f", digits = 2), ")")
  left_label    <- paste("Higher in", group_a_label)
  right_label   <- paste("Higher in", group_b_label)
  plot_df <- df %>%
    mutate(
      neg_log10 = -log10(pmax(p_value, 1e-300)),
      status = case_when(
        !is.na(p_value) & p_value < p_thr & logFC >=  lfc_thr ~ right_label,
        !is.na(p_value) & p_value < p_thr & logFC <= -lfc_thr ~ left_label,
        TRUE ~ "Not significant"
      ),
      sig_level = case_when(
        status == "Not significant" ~ "Not significant",
        !is.na(p_value) & p_value <= 0.05 ~ "Primary (nominal p<=0.05)",
        TRUE ~ secondary_lbl
      )
    )
  ann <- plot_df %>%
    filter(status != "Not significant") %>%
    arrange(p_value, desc(abs(logFC))) %>%
    slice_head(n = 18)
  if (nrow(ann) == 0)
    ann <- plot_df %>% arrange(p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  x_rng         <- range(plot_df$logFC, na.rm = TRUE)
  y_thr         <- -log10(p_thr)
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
             label = paste0(" Nominal P = ", formatC(p_thr, format = "f", digits = 2)),
             hjust = 0, vjust = -0.6, size = 3.8, fontface = "italic") +
    scale_color_manual(values = c(
      setNames("#2166AC", left_label), setNames("#B2182B", right_label),
      "Not significant" = "grey85"), drop = FALSE) +
    scale_shape_manual(values = c("Not significant" = 1, "Primary (nominal p<=0.05)" = 19,
                                 setNames(17, secondary_lbl)),
                       breaks = c("Not significant", "Primary (nominal p<=0.05)", secondary_lbl), drop = FALSE) +
    ggrepel::geom_text_repel(data = ann, aes(label = .data[[label_col]]),
                             size = 3.1, fontface = "bold", max.overlaps = 30,
                             box.padding = 0.5, point.padding = 0.2,
                             segment.color = "grey60", segment.size = 0.2, show.legend = FALSE) +
    labs(title = title, x = paste0("Fold change log2(", group_b_label, " / ", group_a_label, ")"),
         y = "-log10(nominal p value)", color = NULL, shape = "Significance") +
    guides(color = "none", shape = guide_legend(override.aes = list(color = "black"))) +
    coord_cartesian(xlim = c(x_panel_left, x_panel_right), ylim = c(0, y_panel_top), clip = "off") +
    akila_volcano_theme()
  save_ggplot(p, path_stub, width = 8.6, height = 8.2)
}

make_ma_plot <- function(df, title, path_stub,
                         lfc_thr = LFC_THR, p_thr = PVAL_CUTOFF,
                         label_col = "Protein_Label_Name",
                         group_a_label = NULL, group_b_label = NULL) {
  # NOTE: thresholding/coloring below uses NOMINAL p_value (uncorrected).
  # True BH-adjusted value is carried as adj_p_value in the data/CSV only.
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  group_a_label <- group_a_label %||% (unique(na.omit(df$group_a))[1] %||% "Group A")
  group_b_label <- group_b_label %||% (unique(na.omit(df$group_b))[1] %||% "Group B")
  left_label  <- paste("Higher in", group_a_label)
  right_label <- paste("Higher in", group_b_label)
  plot_df <- df %>%
    mutate(
      status = case_when(
        !is.na(p_value) & p_value < p_thr & logFC >=  lfc_thr ~ right_label,
        !is.na(p_value) & p_value < p_thr & logFC <= -lfc_thr ~ left_label,
        TRUE ~ "Not significant"
      )
    )
  ann <- plot_df %>%
    filter(status != "Not significant") %>%
    arrange(p_value, desc(abs(logFC))) %>%
    slice_head(n = 18)
  p <- ggplot(plot_df, aes(x = avg_expr, y = logFC)) +
    geom_point(aes(color = status), size = 2.2, alpha = 0.8, stroke = 0.6) +
    geom_hline(yintercept = c(-lfc_thr, 0, lfc_thr),
               linetype = c("dashed", "solid", "dashed"),
               linewidth = c(0.6, 0.4, 0.6), color = "black") +
    scale_color_manual(values = c(
      setNames("#2166AC", left_label), setNames("#B2182B", right_label),
      "Not significant" = "grey80"), drop = FALSE) +
    ggrepel::geom_text_repel(data = ann, aes(label = .data[[label_col]]),
                             size = 3.1, fontface = "bold", max.overlaps = 30,
                             box.padding = 0.5, point.padding = 0.2,
                             segment.color = "grey60", segment.size = 0.2, show.legend = FALSE) +
    labs(title = title, x = "Average log2 expression (A)",
         y = paste0("log2 fold change (M): ", group_b_label, " / ", group_a_label),
         color = NULL) +
    akila_volcano_theme()
  save_ggplot(p, path_stub, width = 8.6, height = 8.2)
}

run_exploratory_contrast <- function(result_df, feature_meta, tables_dir, figures_dir,
                                     file_stub, title_stub,
                                     group_a_label = NULL, group_b_label = NULL) {
  if (is.null(result_df) || nrow(result_df) == 0) {
    message(title_stub, ": insufficient data — skipped")
    return(tibble())
  }
  result_df <- add_feature_annotation(result_df, feature_meta)
  safe_write_csv(result_df, file.path(tables_dir, paste0(file_stub, ".csv")))
  make_volcano(result_df, paste0(title_stub, " — Volcano"),
              file.path(figures_dir, paste0(file_stub, "_volcano")),
              group_a_label = group_a_label, group_b_label = group_b_label)
  make_ma_plot(result_df, paste0(title_stub, " — MA plot"),
              file.path(figures_dir, paste0(file_stub, "_ma")),
              group_a_label = group_a_label, group_b_label = group_b_label)
  message(title_stub, ": ", nrow(result_df), " proteins | nominal p<", PVAL_CUTOFF, ": ",
          sum(result_df$p_value < PVAL_CUTOFF, na.rm = TRUE),
          " | BH-adjusted p<", PVAL_CUTOFF, ": ",
          sum(result_df$adj_p_value < PVAL_CUTOFF, na.rm = TRUE))
  result_df
}

# ============================================================
# (A) Surgery: Post vs Pre (paired)
# ============================================================
message("\n=== A | Surgery: Post vs Pre ===")
surgery_pair_ids <- paired_clinical %>% filter(new_group == "Sx") %>% pull(patient_id)
paired_pre_sx    <- pre_matched  %>% filter(patient_id %in% surgery_pair_ids)
paired_post_sx   <- post_matched %>% filter(patient_id %in% surgery_pair_ids)
surgery_sample_df <- bind_rows(
  paired_pre_sx  %>% mutate(pre_post = "pre"),
  paired_post_sx %>% mutate(pre_post = "post")
)
message("Surgery pairs available: ", length(intersect(paired_pre_sx$patient_id, paired_post_sx$patient_id)))

surgery_ppost_result <- run_paired_limma(
  expr_log2_rep, surgery_sample_df,
  pair_id_col = "patient_id", time_col = "pre_post",
  time_levels = c("pre", "post"), min_pairs = MIN_PAIRS
)
surgery_ppost_result <- run_exploratory_contrast(
  surgery_ppost_result, feature_meta, tbl_dir, fig_dir,
  file_stub = "A_surgery_post_vs_pre",
  title_stub = "Surgery: Post vs Pre",
  group_a_label = "Pre-surgery", group_b_label = "Post-surgery"
)

# ============================================================
# (B) Immunotherapy response: PR/CR vs SD/PD (new_group == "IO")
# ============================================================
message("\n=== B | Immunotherapy response ===")
# Assumption: restricted to single-agent immunotherapy (new_group == "IO");
# chemo-immunotherapy combinations are excluded to isolate the IO effect.
io_df <- pre_matched %>%
  filter(new_group == "IO") %>%
  mutate(response_binary2 = response_binary)
message("IO patients with response data: ",
        sum(!is.na(io_df$response_binary2)))

io_response_result <- run_binary_limma(
  expr_log2_rep, io_df, "response_binary2",
  group_levels  = c("SD/PD", "PR/CR"),
  adjust_cols   = c("plate_barcode"),
  min_per_group = MIN_PER_GROUP
)
io_response_result <- run_exploratory_contrast(
  io_response_result, feature_meta, tbl_dir, fig_dir,
  file_stub = "B_immunotherapy_response_prcr_vs_sdpd",
  title_stub = "Immunotherapy response: PR/CR vs SD/PD",
  group_a_label = "SD/PD", group_b_label = "PR/CR"
)

# ============================================================
# (C) Surgery recurrence: Recurrence vs No recurrence (baseline PRE)
# ============================================================
message("\n=== C | Surgery recurrence profile ===")
sx_df <- pre_matched %>%
  filter(new_group == "Sx") %>%
  mutate(recurrence_binary2 = recurrence_binary)
message("Surgery patients with recurrence status: ",
        sum(!is.na(sx_df$recurrence_binary2)))

sx_recurrence_result <- run_binary_limma(
  expr_log2_rep, sx_df, "recurrence_binary2",
  group_levels  = c("No recurrence", "Recurrence"),
  adjust_cols   = c("plate_barcode"),
  min_per_group = MIN_PER_GROUP
)
sx_recurrence_result <- run_exploratory_contrast(
  sx_recurrence_result, feature_meta, tbl_dir, fig_dir,
  file_stub = "C_surgery_recurrence_vs_no_recurrence",
  title_stub = "Surgery: Recurrence vs No recurrence",
  group_a_label = "No recurrence", group_b_label = "Recurrence"
)

# ============================================================
# (D) Chemotherapy response: PR/CR vs SD/PD (new_group == "CTx")
# ============================================================
message("\n=== D | Chemotherapy response ===")
ctx_df <- pre_matched %>%
  filter(new_group == "CTx") %>%
  mutate(response_binary2 = response_binary)
message("Chemotherapy patients with response data: ",
        sum(!is.na(ctx_df$response_binary2)))

ctx_response_result <- run_binary_limma(
  expr_log2_rep, ctx_df, "response_binary2",
  group_levels  = c("SD/PD", "PR/CR"),
  adjust_cols   = c("plate_barcode"),
  min_per_group = MIN_PER_GROUP
)
ctx_response_result <- run_exploratory_contrast(
  ctx_response_result, feature_meta, tbl_dir, fig_dir,
  file_stub = "D_chemotherapy_response_prcr_vs_sdpd",
  title_stub = "Chemotherapy response: PR/CR vs SD/PD",
  group_a_label = "SD/PD", group_b_label = "PR/CR"
)

# ============================================================
# (E) Progression Free Survival risk: forest + KM plots
# ============================================================
message("\n=== E | PFS univariate Cox risk screen ===")

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op  <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out[tolower(raw_chr) %in% c("dead", "deceased", "event", "yes", "true")] <- 1
  out[tolower(raw_chr) %in% c("alive", "censored", "no", "false")]        <- 0
  out[!is.finite(out) | !(out %in% c(0, 1))] <- NA_real_
  as.integer(out)
}

run_univariate_cox <- function(mat_samples_by_features, sample_df, time_col, event_col,
                               min_n = MIN_N_COX, min_events = MIN_EVENTS_COX) {
  op  <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  df <- sample_df %>%
    mutate(time_value  = time_val,
           event_value = coerce_event_indicator(.data[[event_col]])) %>%
    filter(!is.na(time_value), !is.na(event_value))
  if (nrow(df) < min_n || sum(df$event_value == 1, na.rm = TRUE) < min_events)
    return(tibble())
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  results <- lapply(colnames(mat_sub), function(feature) {
    x    <- mat_sub[, feature]
    keep <- !is.na(x) & is.finite(x) & !is.na(df$time_value) & !is.na(df$event_value)
    if (sum(keep) < min_n || sum(df$event_value[keep] == 1, na.rm = TRUE) < min_events)
      return(NULL)
    fit <- tryCatch(
      survival::coxph(survival::Surv(df$time_value[keep], df$event_value[keep]) ~ scale(x[keep])),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    tibble(
      nELISA_ID = feature,
      beta      = unname(s$coef[1, "coef"]),
      HR        = unname(s$coef[1, "exp(coef)"]),
      HR_low    = unname(s$conf.int[1, "lower .95"]),
      HR_high   = unname(s$conf.int[1, "upper .95"]),
      z         = unname(s$coef[1, "z"]),
      p_value   = unname(s$coef[1, "Pr(>|z|)"]),
      n         = sum(keep),
      n_events  = sum(df$event_value[keep] == 1, na.rm = TRUE)
    )
  })
  bind_rows(results) %>% mutate(adj_p_value = p.adjust(p_value, method = "BH"))
}

make_forest_plot <- function(df, title, path_stub, n = FOREST_N) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  plot_df <- df %>%
    filter(is.finite(HR), is.finite(HR_low), is.finite(HR_high)) %>%
    arrange(p_value) %>%
    slice_head(n = min(n, nrow(.))) %>%
    mutate(label     = factor(paste0(Protein_Label_Name, " (events=", n_events, ")"),
                              levels = rev(paste0(Protein_Label_Name, " (events=", n_events, ")"))),
           direction = ifelse(HR >= 1, "Higher risk", "Lower risk"))
  if (nrow(plot_df) == 0) return(invisible(NULL))
  p <- ggplot(plot_df, aes(x = HR, y = label, color = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbarh(aes(xmin = HR_low, xmax = HR_high), height = 0.2, linewidth = 0.5) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("Higher risk" = "#B2182B", "Lower risk" = "#2166AC")) +
    scale_x_log10() +
    labs(title = title, x = "Hazard ratio (log scale)", y = NULL, color = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  save_ggplot(p, path_stub, width = 8.2, height = 6.8)
}

run_km_curves <- function(mat_samples_by_features, sample_df, result_df,
                         time_col, event_col, endpoint_label, figure_dir, table_dir,
                         top_n = KM_TOP_N) {
  if (is.null(result_df) || nrow(result_df) == 0) return(tibble())
  op  <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  df <- sample_df %>%
    mutate(time_value  = time_val,
           event_value = coerce_event_indicator(.data[[event_col]])) %>%
    filter(!is.na(time_value), !is.na(event_value),
           sample_id %in% rownames(mat_samples_by_features))
  if (nrow(df) < 15 || sum(df$event_value == 1, na.rm = TRUE) < 6) return(tibble())
  candidates <- result_df %>%
    filter(nELISA_ID %in% colnames(mat_samples_by_features), is.finite(HR), is.finite(p_value)) %>%
    arrange(p_value, adj_p_value)
  if (nrow(candidates) == 0) return(tibble())
  km_rows <- list(); saved <- 0L
  for (i in seq_len(nrow(candidates))) {
    if (saved >= top_n) break
    feature_i <- candidates$nELISA_ID[i]
    value_i   <- mat_samples_by_features[df$sample_id, feature_i]
    keep      <- is.finite(value_i) & !is.na(df$time_value) & !is.na(df$event_value)
    if (sum(keep) < 15) next
    km_df <- df[keep, , drop = FALSE] %>%
      mutate(expr_value = value_i[keep],
             expr_group = factor(ifelse(expr_value > median(expr_value, na.rm = TRUE), "High", "Low"),
                                 levels = c("Low", "High")))
    if (any(table(km_df$expr_group) < 4) || length(unique(km_df$event_value)) < 2 ||
        sum(km_df$event_value == 1, na.rm = TRUE) < 5) next
    fit <- tryCatch(
      survival::survfit(survival::Surv(time_value, event_value) ~ expr_group, data = km_df),
      error = function(e) NULL)
    diff_fit <- tryCatch(
      survival::survdiff(survival::Surv(time_value, event_value) ~ expr_group, data = km_df),
      error = function(e) NULL)
    if (is.null(fit)) next
    logrank_p <- if (!is.null(diff_fit))
      stats::pchisq(diff_fit$chisq, df = max(length(diff_fit$n) - 1, 1), lower.tail = FALSE) else NA_real_
    title_stub <- paste0(endpoint_label, " KM — ", candidates$Protein_Label_Name[i])
    km_plot <- survminer::ggsurvplot(
      fit, data = km_df, pval = TRUE, risk.table = FALSE, conf.int = FALSE,
      palette = c("#2166AC", "#B2182B"), legend.title = "Expression",
      legend.labs = c("Low", "High"), ggtheme = ggplot2::theme_bw(base_size = 11),
      title = title_stub,
      subtitle = paste0("Median split | HR=", round(candidates$HR[i], 2),
                        " | Cox adj.P=", scales::number(candidates$adj_p_value[i], accuracy = 0.001)),
      xlab = "Months", ylab = "Progression-free survival probability")
    km_risk <- survminer::ggsurvplot(
      fit, data = km_df, pval = FALSE, risk.table = TRUE, conf.int = FALSE,
      palette = c("#2166AC", "#B2182B"), legend = "none",
      ggtheme = ggplot2::theme_bw(base_size = 11), risk.table.y.text = FALSE,
      xlab = "Months", ylab = "PFS probability")
    saved <- saved + 1L
    path_stub_km <- file.path(figure_dir,
      paste0("E_pfs_km_", sprintf("%02d_", saved), sanitize_filename(candidates$Protein_Label_Name[i])))
    cairo_pdf(paste0(path_stub_km, ".pdf"), width = 7.8, height = 5.5)
    print(km_plot$plot); dev.off()
    png(paste0(path_stub_km, ".png"), width = 7.8, height = 5.5, units = "in", res = 300, bg = "white")
    print(km_plot$plot); dev.off()
    cairo_pdf(paste0(path_stub_km, "_risk.pdf"), width = 7.8, height = 2.6)
    print(km_risk$table); dev.off()
    png(paste0(path_stub_km, "_risk.png"), width = 7.8, height = 2.6, units = "in", res = 300, bg = "white")
    print(km_risk$table); dev.off()
    km_rows[[saved]] <- tibble(
      endpoint = endpoint_label, rank = saved,
      nELISA_ID = feature_i, Protein_Label_Name = candidates$Protein_Label_Name[i],
      HR = candidates$HR[i], p_value = candidates$p_value[i],
      adj_p_value = candidates$adj_p_value[i], logrank_p_value = logrank_p,
      n = nrow(km_df), n_events = sum(km_df$event_value == 1, na.rm = TRUE),
      n_low = sum(km_df$expr_group == "Low"), n_high = sum(km_df$expr_group == "High"))
  }
  km_tbl <- bind_rows(km_rows)
  if (nrow(km_tbl))
    safe_write_csv(km_tbl, file.path(table_dir, "E_pfs_km_summary.csv"))
  km_tbl
}

pfs_df  <- pre_matched   # baseline PRE samples matched to clinical (whole cohort)
pfs_res <- run_univariate_cox(expr_log2_rep, pfs_df, "pfs_months_main", "pfs_event_main")
pfs_res <- add_feature_annotation(pfs_res, feature_meta)

if (nrow(pfs_res) > 0) {
  safe_write_csv(pfs_res, file.path(tbl_pfs, "E_pfs_univariate_cox.csv"))
  make_forest_plot(pfs_res, "PFS risk: univariate Cox top proteins",
                   file.path(fig_pfs, "E_pfs_forest"))
  pfs_km <- run_km_curves(expr_log2_rep, pfs_df, pfs_res,
                          "pfs_months_main", "pfs_event_main", "PFS",
                          fig_pfs, tbl_pfs)
  message("PFS Cox: ", nrow(pfs_res), " proteins | nominal p<", PVAL_CUTOFF, ": ",
          sum(pfs_res$p_value < PVAL_CUTOFF, na.rm = TRUE),
          " | BH-adjusted p<", PVAL_CUTOFF, ": ",
          sum(pfs_res$adj_p_value < PVAL_CUTOFF, na.rm = TRUE),
          " | KM plots saved: ", nrow(pfs_km))
} else {
  message("PFS Cox: insufficient data (need >= ", MIN_N_COX, " samples / ",
          MIN_EVENTS_COX, " events)")
}

# ── Save checkpoint ────────────────────────────────────────────
message("\n=== 13 | Saving checkpoint_13.RData ===")
save(
  surgery_ppost_result, io_response_result,
  sx_recurrence_result, ctx_response_result,
  pfs_res,
  file = file.path(paths$rdata, "checkpoint_13.RData")
)
message("13_exploratory_subgroup_analyses.R complete")
