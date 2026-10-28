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
message("Checkpoints 01, 03 loaded")

gene_sets <- NULL
ck09 <- file.path(paths$rdata, "checkpoint_09.RData")
if (file.exists(ck09)) {
  load(ck09)
  message("Checkpoint 09 loaded (gene_sets reused)")
}

while (dev.cur() > 1) dev.off()   # clear any stale devices

# ── Local constants ───────────────────────────────────────────
GROUP_LABEL     <- "CTx + IO"
OVR_KEY         <- "analysis1_ctx_io_vs_rest"   # correct key (see bug note above)
FDR_CUTOFF      <- 0.10
FGSEA_FDR       <- 0.25
MIN_GENE_SET    <- 5
MAX_GENE_SET    <- 500
TOP_PLOT        <- 10
EMAP_TOP        <- 10
RNG_SEED        <- 42
SUBGROUP_MIN_N       <- 10   # Cox / Wilcoxon floor within subgroup
SUBGROUP_MIN_EVENTS  <- 4
KM_MIN_N        <- 10
KM_MIN_EVENTS   <- 4
FOREST_N        <- 20   # top-N proteins shown in HR forest plots
HEADLINE_GENES  <- c("GDF15", "GZMA", "NOTUM", "SHH", "CEACAM5", "PLAUR")

TBL_DIR <- file.path(paths$tables,  "11_ctxio_subgroup_deep_dive")
FIG_DIR <- file.path(paths$figures, "11_ctxio_subgroup_deep_dive")
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

PW_TBL_BASE <- file.path(paths$tables,  ANALYSIS8_DIR)
PW_FIG_BASE <- file.path(paths$figures, ANALYSIS8_DIR)

# ── Local functions ────────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  out <- substr(out, 1, 60)
  ifelse(nzchar(out), out, "item")
}

add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  left_join(df, feature_meta, by = "nELISA_ID")
}

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op  <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out[!is.finite(out) | !(out %in% c(0, 1))] <- NA_real_
  as.integer(out)
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

# Simpler theme for boxplots (matches 08k_ceacam5_tki_sensitivity_check.R style) —
# deliberately distinct from akila_volcano_theme(): boxplots keep gridlines for
# readability and don't need the wide margin/bottom-legend machinery volcanoes need.
akila_boxplot_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 15),
      plot.subtitle     = element_text(size = 11, colour = "grey30"),
      axis.title        = element_text(face = "bold")
    )
}

get_gene_id <- function(gene_symbol, feature_meta) {
  ids <- feature_meta %>% filter(Gene_Symbol == gene_symbol) %>% pull(nELISA_ID) %>% unique()
  if (length(ids) != 1) {
    message("  WARN: expected exactly one nELISA_ID for ", gene_symbol,
            ", found ", length(ids), " — skipping")
    return(NA_character_)
  }
  ids
}

# Generic two-group Wilcoxon (mirrors 08k's run_ceacam5_wilcox, generalised)
run_group_wilcox <- function(sample_df, expr_vec, group_col, protein_label, min_per_group = 3) {
  df <- sample_df %>%
    mutate(expr = expr_vec[sample_id]) %>%
    filter(!is.na(.data[[group_col]]), !is.na(expr))
  lv <- levels(factor(df[[group_col]]))
  if (length(lv) != 2) return(tibble())
  n1 <- sum(df[[group_col]] == lv[1]); n2 <- sum(df[[group_col]] == lv[2])
  if (n1 < min_per_group || n2 < min_per_group) {
    return(tibble(protein = protein_label, group_col = group_col,
                  level_1 = lv[1], level_2 = lv[2], n1 = n1, n2 = n2,
                  median_1 = NA_real_, median_2 = NA_real_, wilcox_p = NA_real_,
                  status = "SKIPPED: n<3 in a cell"))
  }
  wt <- wilcox.test(stats::as.formula(paste("expr ~", group_col)), data = df)
  tibble(protein = protein_label, group_col = group_col,
         level_1 = lv[1], level_2 = lv[2], n1 = n1, n2 = n2,
         median_1 = median(df$expr[df[[group_col]] == lv[1]]),
         median_2 = median(df$expr[df[[group_col]] == lv[2]]),
         wilcox_p = wt$p.value, status = "OK")
}

# Binary limma (copied from 03_analysis1_seven_group.R) — used here for the
# within-subgroup recurrence and response contrasts (Parts 5-6 below).
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

# Volcano plot (copied from 03_analysis1_seven_group.R, unchanged)
make_volcano <- function(df, title, path_stub, lfc_thr = 0.5, padj_thr = FDR_CUTOFF,
                         label_col = "Protein_Label", group_a_label = NULL, group_b_label = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
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
  if (nrow(ann) == 0) ann <- plot_df %>% arrange(adj_p_value, desc(abs(logFC))) %>% slice_head(n = 18)

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
      setNames("#2166AC", left_label), setNames("#B2182B", right_label), "Not significant" = "grey85"
    ), drop = FALSE) +
    scale_shape_manual(values = c(
      "Not significant" = 1, "Primary (FDR<=0.05)" = 19, setNames(17, secondary_lbl)
    ), breaks = c("Not significant", "Primary (FDR<=0.05)", secondary_lbl), drop = FALSE) +
    ggrepel::geom_text_repel(
      data = ann, aes(label = .data[[label_col]]),
      size = 3.1, fontface = "bold", max.overlaps = 30,
      box.padding = 0.5, point.padding = 0.2,
      segment.color = "grey60", segment.size = 0.2, show.legend = FALSE
    ) +
    labs(title = title, x = x_label, y = "P value (-log10)", color = NULL, shape = "Significance") +
    guides(color = "none", shape = guide_legend(override.aes = list(color = "black"))) +
    coord_cartesian(xlim = c(x_panel_left, x_panel_right), ylim = c(0, y_panel_top), clip = "off") +
    akila_volcano_theme()

  cairo_pdf(paste0(path_stub, ".pdf"), width = 8.6, height = 8.2); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 8.6, height = 8.2, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# HR forest plot (copied from 07_analysis5_survival.R, unchanged)
make_forest_plot <- function(df, title, path_stub, n = FOREST_N) {
  if (nrow(df) == 0) return(invisible(NULL))
  plot_df <- df %>%
    filter(is.finite(HR), is.finite(HR_low), is.finite(HR_high)) %>%
    arrange(p_value) %>%
    slice_head(n = min(n, nrow(.))) %>%
    mutate(label     = factor(paste0(Protein_Label, " (events=", n_events, ")"),
                              levels = rev(paste0(Protein_Label, " (events=", n_events, ")"))),
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
  cairo_pdf(paste0(path_stub, ".pdf"), width = 8.2, height = 6.8); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 8.2, height = 6.8, units="in", res=300, bg="white")
  print(p); dev.off()
}

# Survival volcano: log2(HR) vs -log10(p) (copied from 07_analysis5_survival.R, unchanged)
make_survival_volcano <- function(df, title, path_stub, label_col = "Protein_Label",
                                  n = FOREST_N, padj_thr = FDR_CUTOFF) {
  if (nrow(df) == 0) return(invisible(NULL))
  secondary_label <- paste0("Secondary (FDR<=", formatC(padj_thr, format="f", digits=2), ")")
  plot_df <- df %>%
    filter(is.finite(HR), HR > 0, is.finite(p_value)) %>%
    mutate(
      log2_hr   = log2(HR),
      neg_log10 = -log10(pmax(p_value, 1e-300)),
      direction = case_when(
        p_value < 0.05 & log2_hr > 0 ~ "Higher risk",
        p_value < 0.05 & log2_hr < 0 ~ "Protective",
        TRUE ~ "Not significant"),
      sig_level = case_when(
        direction == "Not significant" ~ "Not significant",
        p_value <= 0.05               ~ "Primary (P<=0.05)",
        TRUE                          ~ secondary_label))
  if (nrow(plot_df) == 0) return(invisible(NULL))
  ann <- plot_df %>% filter(direction != "Not significant") %>%
    arrange(p_value, desc(abs(log2_hr))) %>% slice_head(n = min(n, nrow(.)))
  if (nrow(ann) == 0)
    ann <- plot_df %>% arrange(p_value, desc(abs(log2_hr))) %>% slice_head(n = min(n, nrow(.)))
  x_rng <- range(plot_df$log2_hr, na.rm=TRUE)
  y_thr <- -log10(0.05)
  y_max <- max(plot_df$neg_log10, y_thr, na.rm=TRUE)
  x_pad <- max(0.22, diff(x_rng)*0.24); y_pad <- max(0.16, y_max*0.16)
  x_pl  <- x_rng[1]-0.45*x_pad; x_pr <- x_rng[2]+0.45*x_pad
  y_pt  <- y_max+0.12*y_pad;     y_lbl <- y_pt+0.9*y_pad
  p <- ggplot(plot_df, aes(x = log2_hr, y = neg_log10)) +
    geom_point(aes(color = direction, shape = sig_level), size=2.4, alpha=0.82, stroke=0.7) +
    geom_vline(xintercept=0, linetype="solid",  linewidth=0.4, color="black") +
    geom_hline(yintercept=y_thr, linetype="dashed", linewidth=0.6, color="black") +
    annotate("text", x=x_pl-0.18*x_pad, y=y_lbl, label="Protective",
             hjust=0, vjust=0, size=4.8, fontface="bold", color="#2166AC") +
    annotate("text", x=x_pr+0.18*x_pad, y=y_lbl, label="Higher risk",
             hjust=1, vjust=0, size=4.8, fontface="bold", color="#B2182B") +
    annotate("text", x=x_rng[1]+0.02*diff(x_rng), y=y_thr,
             label=paste0(" P = 0.05"), hjust=0, vjust=-0.6, size=3.8, fontface="italic") +
    scale_color_manual(values=c("Higher risk"="#B2182B","Protective"="#2166AC","Not significant"="grey85"), drop=FALSE) +
    scale_shape_manual(values=c("Not significant"=1,"Primary (P<=0.05)"=19,setNames(17,secondary_label)),
                       breaks=c("Not significant","Primary (P<=0.05)",secondary_label), drop=FALSE) +
    ggrepel::geom_text_repel(data=ann, aes(label=.data[[label_col]]),
                             size=3.1, fontface="bold", max.overlaps=30,
                             box.padding=0.5, point.padding=0.2,
                             segment.color="grey60", segment.size=0.2, show.legend=FALSE) +
    labs(title=title, x="log2(Hazard Ratio)", y="-log10(p value)", color=NULL, shape="Significance") +
    guides(color="none", shape=guide_legend(override.aes=list(color="black"))) +
    coord_cartesian(xlim=c(x_pl,x_pr), ylim=c(0,y_pt), clip="off") +
    akila_volcano_theme()
  cairo_pdf(paste0(path_stub, ".pdf"), width=8.6, height=8.2); print(p); dev.off()
  png(paste0(path_stub, ".png"), width=8.6, height=8.2, units="in", res=300, bg="white")
  print(p); dev.off()
}

run_univariate_cox <- function(mat_samples_by_features, sample_df, time_col, event_col,
                               min_n = SUBGROUP_MIN_N, min_events = SUBGROUP_MIN_EVENTS) {
  op  <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  df <- sample_df %>%
    mutate(time_value = time_val, event_value = coerce_event_indicator(.data[[event_col]])) %>%
    filter(!is.na(time_value), !is.na(event_value))
  if (nrow(df) < min_n || sum(df$event_value == 1, na.rm = TRUE) < min_events) return(tibble())
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  results <- lapply(colnames(mat_sub), function(feature) {
    x    <- mat_sub[, feature]
    keep <- !is.na(x) & is.finite(x) & !is.na(df$time_value) & !is.na(df$event_value)
    if (sum(keep) < min_n || sum(df$event_value[keep] == 1, na.rm = TRUE) < min_events) return(NULL)
    fit <- tryCatch(
      survival::coxph(survival::Surv(df$time_value[keep], df$event_value[keep]) ~ scale(x[keep])),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    tibble(nELISA_ID = feature, beta = unname(s$coef[1, "coef"]), HR = unname(s$coef[1, "exp(coef)"]),
           HR_low = unname(s$conf.int[1, "lower .95"]), HR_high = unname(s$conf.int[1, "upper .95"]),
           z = unname(s$coef[1, "z"]), p_value = unname(s$coef[1, "Pr(>|z|)"]),
           n = sum(keep), n_events = sum(df$event_value[keep] == 1, na.rm = TRUE))
  })
  bind_rows(results) %>% mutate(adj_p_value = p.adjust(p_value, method = "BH"))
}

plot_km_median_split <- function(sample_df, expr_vec, time_col, event_col, protein_label,
                                 subtitle_extra, path_stub, min_n = KM_MIN_N, min_events = KM_MIN_EVENTS) {
  op <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  event_val <- coerce_event_indicator(sample_df[[event_col]])
  keep <- is.finite(expr_vec) & !is.na(time_val) & !is.na(event_val)
  if (sum(keep) < min_n) { message("    SKIP KM (n<", min_n, "): ", protein_label); return(invisible(NULL)) }
  df <- tibble::tibble(time_value = time_val[keep], event_value = event_val[keep], expr_value = expr_vec[keep]) %>%
    dplyr::mutate(expr_group = factor(ifelse(expr_value > median(expr_value, na.rm = TRUE), "High", "Low"),
                                      levels = c("Low", "High")))
  if (any(table(df$expr_group) < 3) || sum(df$event_value == 1) < min_events) {
    message("    SKIP KM (small stratum/events): ", protein_label); return(invisible(NULL))
  }
  fit <- survival::survfit(survival::Surv(time_value, event_value) ~ expr_group, data = df)
  diff_fit <- tryCatch(survival::survdiff(survival::Surv(time_value, event_value) ~ expr_group, data = df),
                       error = function(e) NULL)
  logrank_p <- if (!is.null(diff_fit))
    stats::pchisq(diff_fit$chisq, df = max(length(diff_fit$n) - 1, 1), lower.tail = FALSE) else NA_real_
  km_plot <- survminer::ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    palette = c("#2166AC", "#B2182B"), legend.title = "Expression", legend.labs = c("Low", "High"),
    title = paste0("KM \u2014 ", protein_label, " (CTx+IO subgroup, n=22)"),
    subtitle = paste0(subtitle_extra, " | logrank p=", scales::number(logrank_p, accuracy = 0.001)),
    xlab = "Months", ylab = "Survival probability", ggtheme = ggplot2::theme_bw(base_size = 11))
  tryCatch({
    cairo_pdf(paste0(path_stub, ".pdf"), width = 7.8, height = 6.5); print(km_plot); dev.off()
    png(paste0(path_stub, ".png"), width = 7.8, height = 6.5, units = "in", res = 300, bg = "white")
    print(km_plot); dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("    ERROR saving KM plot: ", conditionMessage(e), " -> close open PDF viewers and re-run.")
  })
  tibble::tibble(protein = protein_label, endpoint = NA_character_, n = nrow(df),
                 n_events = sum(df$event_value == 1), n_low = sum(df$expr_group == "Low"),
                 n_high = sum(df$expr_group == "High"), logrank_p = logrank_p)
}

# ── fgsea machinery (copied from 09_pathway_network.R, unchanged) ─
clean_pathway_name <- function(x) {
  x %>%
    stringr::str_replace_all("^HALLMARK_", "") %>% stringr::str_replace_all("^REACTOME_", "") %>%
    stringr::str_replace_all("^GOBP_|^GO:BP_", "") %>% stringr::str_replace_all("^GOMF_|^GO:MF_", "") %>%
    stringr::str_replace_all("^GOCC_|^GO:CC_", "") %>% stringr::str_replace_all("^KEGG_(MEDICUS|LEGACY)_", "") %>%
    stringr::str_replace_all("^KEGG_", "") %>% stringr::str_replace_all("_", " ") %>% stringr::str_squish()
}

get_msigdbr_table <- function(collection, subcollection = NULL) {
  fn_args <- names(formals(msigdbr::msigdbr))
  args <- list(species = "Homo sapiens")
  if ("collection" %in% fn_args) {
    args$collection <- collection; if (!is.null(subcollection)) args$subcollection <- subcollection
  } else {
    args$category <- collection; if (!is.null(subcollection)) args$subcategory <- subcollection
  }
  do.call(msigdbr::msigdbr, args)
}

build_gene_sets <- function() {
  kegg_medicus <- tryCatch(get_msigdbr_table("C2", "CP:KEGG_MEDICUS"), error = function(e) tibble())
  kegg_legacy  <- tryCatch(get_msigdbr_table("C2", "CP:KEGG_LEGACY"),  error = function(e) tibble())
  sets <- list(
    HALLMARK = get_msigdbr_table("H"), GO_BP = get_msigdbr_table("C5", "GO:BP"),
    GO_MF = get_msigdbr_table("C5", "GO:MF"), GO_CC = get_msigdbr_table("C5", "GO:CC"),
    KEGG = bind_rows(kegg_medicus %>% mutate(gs_name = paste0("KEGG_MEDICUS_", gs_name)),
                     kegg_legacy  %>% mutate(gs_name = paste0("KEGG_LEGACY_",  gs_name))),
    REACTOME = get_msigdbr_table("C2", "CP:REACTOME"))
  out <- lapply(sets, function(df) {
    if (is.null(df) || nrow(df) == 0) return(list())
    df %>% select(gs_name, gene_symbol) %>%
      filter(!is.na(gs_name), !is.na(gene_symbol), gene_symbol != "") %>% distinct() %>%
      split(.$gs_name) %>% lapply(function(tbl) unique(tbl$gene_symbol))
  })
  out[lengths(out) > 0]
}

prepare_rank_vector <- function(result_df, stat_col = "t_stat") {
  rank_df <- result_df %>%
    filter(!is.na(Gene_Symbol), !is.na(.data[[stat_col]]), is.finite(.data[[stat_col]])) %>%
    group_by(Gene_Symbol) %>%
    summarise(stat = .data[[stat_col]][which.max(abs(.data[[stat_col]]))], .groups = "drop") %>%
    filter(!is.na(stat), is.finite(stat))
  stats <- rank_df$stat; names(stats) <- rank_df$Gene_Symbol
  sort(stats, decreasing = TRUE)
}

run_fgsea_one_collection <- function(pathways_use, stats, label, collection_name, coll_tbl_dir, coll_fig_dir) {
  pathways_sub <- lapply(pathways_use, function(gs) intersect(gs, names(stats)))
  pathways_sub <- pathways_sub[lengths(pathways_sub) >= MIN_GENE_SET]
  if (!length(pathways_sub)) return(invisible(NULL))
  set.seed(RNG_SEED)
  fg <- tryCatch(fgsea::fgsea(pathways = pathways_sub, stats = stats, minSize = MIN_GENE_SET,
                              maxSize = MAX_GENE_SET, eps = 0),
                 error = function(e) { message("  fgsea error (", collection_name, "): ", conditionMessage(e)); NULL })
  if (is.null(fg) || nrow(fg) == 0) return(invisible(NULL))
  fg <- fg %>% as_tibble() %>%
    mutate(pathway_clean = clean_pathway_name(pathway),
           leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))) %>%
    arrange(padj, desc(abs(NES)))
  sig_fg <- fg %>% filter(!is.na(padj), padj < FGSEA_FDR)
  lbl_short <- sanitize_filename(label); col_short <- tolower(collection_name)
  safe_write_csv(fg, file.path(coll_tbl_dir, paste0(lbl_short, "_", col_short, "_fgsea.csv")))
  safe_write_csv(sig_fg, file.path(coll_tbl_dir, paste0(lbl_short, "_", col_short, "_fgsea_sig.csv")))
  safe_write_csv(
    fg %>% select(pathway, pathway_clean, NES, pval, padj, leadingEdge) %>%
      tidyr::separate_rows(leadingEdge, sep = ";") %>% rename(Gene_Symbol = leadingEdge) %>%
      filter(!is.na(Gene_Symbol), Gene_Symbol != ""),
    file.path(coll_tbl_dir, paste0(lbl_short, "_", col_short, "_leading_edge.csv")))

  top_plot <- fg %>% filter(!is.na(padj), is.finite(NES)) %>% arrange(padj, desc(abs(NES))) %>%
    slice_head(n = min(TOP_PLOT, nrow(.))) %>%
    mutate(rank = row_number(), abs_nes = abs(NES), neglog_fdr = -log10(pmax(padj, 1e-300)),
           core_gene_count = vapply(strsplit(dplyr::coalesce(leadingEdge, ""), ";", fixed = TRUE),
                                    function(x) sum(nzchar(x)), integer(1)),
           pathway_plot = paste0(rank, ". ", pathway_clean),
           pathway_plot = factor(stringr::str_wrap(pathway_plot, width = 48),
                                 levels = rev(stringr::str_wrap(pathway_plot, width = 48))))
  if (nrow(top_plot)) {
    subtitle_str <- if (nrow(sig_fg) >= 1) paste(gsub("_", " ", label), "| ranked by adjusted p") else
      paste(gsub("_", " ", label), "| top pathways by p; none reached FDR<0.25")
    x_pad <- max(top_plot$abs_nes, na.rm = TRUE) * 0.22
    p <- ggplot(top_plot, aes(abs_nes, pathway_plot)) +
      geom_segment(aes(x = 0, xend = abs_nes, yend = pathway_plot), color = "grey82", linewidth = 0.7) +
      geom_point(aes(size = core_gene_count, fill = neglog_fdr), shape = 21, color = "#53306d", stroke = 0.7) +
      scale_fill_gradientn(colors = c("#f3e5f5", "#d499eb", "#a334d0", "#6f1597"), name = expression(-log[10](FDR))) +
      scale_size_continuous(range = c(4, 13), name = "Core genes") +
      labs(title = paste(collection_name, "top pathways"), subtitle = stringr::str_wrap(subtitle_str, width = 68),
           x = "Normalized Enrichment Score (|NES|)", y = NULL) +
      coord_cartesian(xlim = c(0, max(top_plot$abs_nes, na.rm = TRUE) + x_pad), clip = "off") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right", axis.text.y = element_text(face = "bold"),
            panel.grid.major.y = element_blank(), plot.title = element_text(size = 12.5, face = "bold"),
            plot.subtitle = element_text(size = 10.2), plot.margin = ggplot2::margin(8, 10, 8, 8))
    tryCatch({
      cairo_pdf(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_fgsea.pdf")), width = 9, height = 6)
      print(p); dev.off()
      png(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_fgsea.png")), width = 9, height = 6,
          units = "in", res = 300, bg = "white"); print(p); dev.off()
    }, error = function(e) { if (dev.cur() > 1) dev.off(); message("    WARN lollipop: ", conditionMessage(e)) })
  }

  enrich_candidates <- bind_rows(
    fg %>% filter(NES > 0) %>% slice_min(order_by = padj, n = 1, with_ties = FALSE),
    fg %>% filter(NES < 0) %>% slice_min(order_by = padj, n = 1, with_ties = FALSE)
  ) %>% distinct(pathway, .keep_all = TRUE)
  for (i in seq_len(nrow(enrich_candidates))) {
    pw_i <- enrich_candidates$pathway[i]
    if (!pw_i %in% names(pathways_sub)) next
    p_enrich <- fgsea::plotEnrichment(pathways_sub[[pw_i]], stats) +
      labs(title = paste(label, "-", collection_name),
           subtitle = paste0(clean_pathway_name(pw_i), " | NES=", round(enrich_candidates$NES[i], 2),
                             ", FDR=", scales::number(enrich_candidates$padj[i], accuracy = 0.001)),
           x = "Ranked proteins", y = "Enrichment score") + theme_bw(base_size = 11)
    fname_enrich <- sanitize_filename(clean_pathway_name(pw_i))
    tryCatch({
      cairo_pdf(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_enrich_", fname_enrich, ".pdf")),
                width = 8.6, height = 5.4); print(p_enrich); dev.off()
      png(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_enrich_", fname_enrich, ".png")),
          width = 8.6, height = 5.4, units = "in", res = 300, bg = "white"); print(p_enrich); dev.off()
    }, error = function(e) { if (dev.cur() > 1) dev.off(); message("    WARN enrich plot: ", conditionMessage(e)) })
  }

  emap_nodes <- if (nrow(sig_fg) >= 3) sig_fg else fg %>% arrange(padj, desc(abs(NES)))
  emap_nodes <- emap_nodes %>% slice_head(n = min(EMAP_TOP, nrow(.))) %>%
    mutate(leading_edge_list = strsplit(leadingEdge, ";", fixed = TRUE), core_gene_count = lengths(leading_edge_list))
  if (nrow(emap_nodes) >= 3 && requireNamespace("igraph", quietly = TRUE)) {
    edge_list <- bind_rows(lapply(combn(seq_len(nrow(emap_nodes)), 2, simplify = FALSE), function(idx) {
      ga <- unique(emap_nodes$leading_edge_list[[idx[1]]]); gb <- unique(emap_nodes$leading_edge_list[[idx[2]]])
      shared <- length(intersect(ga, gb)); union_n <- length(unique(c(ga, gb)))
      jaccard <- ifelse(union_n > 0, shared / union_n, 0)
      if (!is.finite(jaccard) || jaccard < 0.15) return(NULL)
      tibble(from = emap_nodes$pathway_clean[idx[1]], to = emap_nodes$pathway_clean[idx[2]],
             weight = jaccard, shared_genes = shared)
    }))
    vertex_df <- emap_nodes %>% transmute(name = pathway_clean, NES, neglog_fdr = -log10(pmax(padj, 1e-300)), core_gene_count)
    if (is.null(edge_list) || nrow(edge_list) == 0)
      edge_list <- tibble(from = character(), to = character(), weight = numeric(), shared_genes = integer())
    g_emap <- igraph::graph_from_data_frame(edge_list, directed = FALSE, vertices = vertex_df)
    emap_title <- if (nrow(sig_fg) >= 3) paste(label, "-", collection_name, "EMAP (FDR<0.25 pathways)") else
      paste(label, "-", collection_name, "EMAP top pathways (none FDR<0.25)")
    tryCatch({
      cairo_pdf(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_emap.pdf")), width = 9.2, height = 7.2)
      set.seed(RNG_SEED)
      lay <- if (igraph::ecount(g_emap) > 0) igraph::layout_with_fr(g_emap) else igraph::layout_in_circle(g_emap)
      vcols <- ifelse(igraph::V(g_emap)$NES >= 0, "#B2182B", "#2166AC")
      vsizes <- scales::rescale(igraph::V(g_emap)$core_gene_count, to = c(14, 28))
      ewidth <- if (igraph::ecount(g_emap) > 0) scales::rescale(igraph::E(g_emap)$weight, to = c(1.2, 4.5)) else numeric()
      plot(g_emap, layout = lay, vertex.color = vcols, vertex.size = vsizes, vertex.label = igraph::V(g_emap)$name,
           vertex.label.cex = 0.72, vertex.label.family = "sans", vertex.frame.color = "grey30",
           edge.width = ewidth, edge.color = "grey70", main = emap_title)
      legend("topleft", legend = c("Positive NES", "Negative NES"), col = c("#B2182B", "#2166AC"), pch = 16, bty = "n", cex = 0.9)
      dev.off()
      png(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_emap.png")), width = 9.2, height = 7.2,
          units = "in", res = 300, bg = "white")
      set.seed(RNG_SEED)
      plot(g_emap, layout = lay, vertex.color = vcols, vertex.size = vsizes, vertex.label = igraph::V(g_emap)$name,
           vertex.label.cex = 0.72, vertex.label.family = "sans", vertex.frame.color = "grey30",
           edge.width = ewidth, edge.color = "grey70", main = emap_title)
      legend("topleft", legend = c("Positive NES", "Negative NES"), col = c("#B2182B", "#2166AC"), pch = 16, bty = "n", cex = 0.9)
      dev.off()
    }, error = function(e) { if (dev.cur() > 1) dev.off(); message("    WARN EMAP: ", conditionMessage(e)) })
  }
  message("    ", collection_name, ": ", nrow(fg), " pathways | FDR<0.25: ", nrow(sig_fg))
  invisible(fg)
}

run_fgsea_suite <- function(result_df, label, base_tbl_dir, base_fig_dir, stat_col = "t_stat") {
  if (!length(gene_sets)) { message("  No gene sets available — skipped"); return(invisible(NULL)) }
  stats <- prepare_rank_vector(result_df, stat_col = stat_col)
  if (length(stats) < 20) { message("  Fewer than 20 ranked genes — skipped"); return(invisible(NULL)) }
  message("  Ranked list: ", length(stats), " genes")
  for (collection_name in names(gene_sets)) {
    coll_tbl_dir <- file.path(base_tbl_dir, sanitize_filename(label), tolower(collection_name))
    coll_fig_dir <- file.path(base_fig_dir, sanitize_filename(label), tolower(collection_name))
    dir.create(coll_tbl_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(coll_fig_dir, recursive = TRUE, showWarnings = FALSE)
    run_fgsea_one_collection(gene_sets[[collection_name]], stats, label, collection_name, coll_tbl_dir, coll_fig_dir)
  }
}

run_string_ppi <- function(result_df, label, base_tbl_dir, base_fig_dir, max_genes = 30, nominal_fallback = TRUE) {
  lbl_dir <- file.path(base_tbl_dir, sanitize_filename(label))
  dir.create(lbl_dir, recursive = TRUE, showWarnings = FALSE)
  genes_fdr <- result_df %>% filter(!is.na(Gene_Symbol), adj_p_value < 0.1) %>%
    arrange(adj_p_value, p_value) %>% distinct(Gene_Symbol, .keep_all = TRUE) %>% slice_head(n = max_genes)
  genes_tbl <- genes_fdr; selection_rule <- "FDR<0.10"
  if (nominal_fallback && nrow(genes_tbl) < 5 && "p_value" %in% names(result_df)) {
    genes_nom <- result_df %>% filter(!is.na(Gene_Symbol), p_value < 0.05) %>% arrange(p_value) %>%
      distinct(Gene_Symbol, .keep_all = TRUE) %>% slice_head(n = max_genes)
    if (nrow(genes_nom) >= 5) { genes_tbl <- genes_nom; selection_rule <- "Nominal p<0.05" }
  }
  safe_write_csv(genes_tbl %>% transmute(Gene_Symbol, selection_rule = selection_rule, adj_p_value, p_value),
                file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_input_genes.csv")))
  ppi_status <- tibble(label = label, genes = nrow(genes_tbl), rule = selection_rule,
                       string_attempted = FALSE, status = "Input genes saved (STRING API call disabled here)")
  safe_write_csv(ppi_status, file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_status.csv")))
}

# ============================================================
# PART 0 — Subgroup snapshot
# ============================================================
message("\n=== 08y | PART 0: CTx+IO subgroup snapshot ===")

ctxio_df <- analysis1_df %>% filter(new_group == GROUP_LABEL)
message("CTx+IO baseline n = ", nrow(ctxio_df))

snapshot <- tibble(
  n_total               = nrow(ctxio_df),
  n_os_evaluable        = sum(!is.na(ctxio_df$os_months)  & !is.na(ctxio_df$os_event)),
  n_os_events           = sum(coerce_event_indicator(ctxio_df$os_event) == 1, na.rm = TRUE),
  n_pfs_evaluable       = sum(!is.na(ctxio_df$pfs_months_main) & !is.na(ctxio_df$pfs_event_main)),
  n_pfs_events          = sum(coerce_event_indicator(ctxio_df$pfs_event_main) == 1, na.rm = TRUE),
  n_response_evaluable  = sum(!is.na(ctxio_df$response_binary)),
  n_recurrence_evaluable = sum(!is.na(ctxio_df$recurrence_binary_sensitivity))
)
safe_write_csv(snapshot, file.path(TBL_DIR, "ctxio_subgroup_snapshot.csv"))
print(snapshot)

response_counts    <- ctxio_df %>% count(response_binary, name = "n")
response4_counts   <- ctxio_df %>% count(response_4cat, name = "n")
recurrence_counts  <- ctxio_df %>% count(recurrence_binary_sensitivity, name = "n")
safe_write_csv(response_counts,   file.path(TBL_DIR, "ctxio_response_binary_counts.csv"))
safe_write_csv(response4_counts,  file.path(TBL_DIR, "ctxio_response_4cat_counts.csv"))
safe_write_csv(recurrence_counts, file.path(TBL_DIR, "ctxio_recurrence_counts.csv"))
message("Response (binary):"); print(response_counts)
message("Recurrence:"); print(recurrence_counts)
message("NOTE: recurrence cell counts are too imbalanced/small for a defensible binary",
        " group comparison — PFS Cox/KM (Part 3) is the time-to-event proxy used instead.")

# Baseline characteristics vs rest (descriptive only)
char_cols <- c("sex", "stage_collapsed", "histology", "pdl1_tps")
have_cols <- intersect(char_cols, names(analysis1_df))
if (length(have_cols) > 0) {
  char_df <- analysis1_df %>%
    mutate(is_ctxio = ifelse(new_group == GROUP_LABEL, "CTx+IO", "Rest")) %>%
    select(is_ctxio, all_of(have_cols))
  safe_write_csv(char_df %>% count(is_ctxio, across(all_of(setdiff(have_cols, "pdl1_tps")))),
                file.path(TBL_DIR, "ctxio_baseline_characteristics_categorical.csv"))
  if ("pdl1_tps" %in% have_cols) {
    pdl1_summary <- char_df %>% group_by(is_ctxio) %>%
      summarise(median_pdl1 = median(pdl1_tps, na.rm = TRUE),
                q1 = quantile(pdl1_tps, 0.25, na.rm = TRUE),
                q3 = quantile(pdl1_tps, 0.75, na.rm = TRUE), n = sum(!is.na(pdl1_tps)), .groups = "drop")
    safe_write_csv(pdl1_summary, file.path(TBL_DIR, "ctxio_baseline_pdl1_summary.csv"))
  }
  message("Baseline characteristics (descriptive) saved.")
} else {
  message("NOTE: none of sex/stage_collapsed/histology/pdl1_tps found in analysis1_df —",
          " check column names against 01_data_loading.R build_clinical() before re-running this block.")
}

# ============================================================
# PART 1 — One-vs-rest DE recap + headline biomarker extraction
# (limma itself already ran correctly in 03_analysis1_seven_group.R;
#  we only recap + pull headline proteins here)
# ============================================================
message("\n=== 08y | PART 1: One-vs-rest DE recap ===")

ovr_res <- analysis1_ovr_results[[OVR_KEY]]
if (is.null(ovr_res) || nrow(ovr_res) == 0) {
  stop("Expected '", OVR_KEY, "' in analysis1_ovr_results (checkpoint_03) but found none — ",
       "re-run 03_analysis1_seven_group.R first.")
}
message("CTx+IO vs rest — proteins tested: ", nrow(ovr_res),
        " | FDR<0.10: ", sum(ovr_res$adj_p_value < FDR_CUTOFF, na.rm = TRUE))
safe_write_csv(ovr_res %>% arrange(p_value), file.path(TBL_DIR, "ctxio_vs_rest_de_recap.csv"))

headline_ids <- setNames(vapply(HEADLINE_GENES, get_gene_id, character(1), feature_meta = feature_meta),
                         HEADLINE_GENES)
headline_ids <- headline_ids[!is.na(headline_ids)]

headline_de <- ovr_res %>% filter(nELISA_ID %in% headline_ids) %>%
  arrange(p_value)
safe_write_csv(headline_de, file.path(TBL_DIR, "ctxio_headline_biomarker_de.csv"))
message("Headline biomarkers in CTx+IO vs rest DE:"); print(headline_de %>% select(Gene_Symbol, logFC, p_value, adj_p_value))

# Headline biomarker boxplots (CTx+IO vs Rest), one file per protein
for (gene in names(headline_ids)) {
  id <- headline_ids[[gene]]
  plot_df <- analysis1_df %>%
    mutate(expr = expr_log2_rep[sample_id, id], grp = ifelse(new_group == GROUP_LABEL, "CTx+IO", "Rest")) %>%
    filter(!is.na(expr))
  stat_row <- headline_de %>% filter(Gene_Symbol == gene)
  subtitle_txt <- if (nrow(stat_row)) paste0("logFC=", round(stat_row$logFC, 2), " | p=", signif(stat_row$p_value, 3),
                                             " | FDR=", signif(stat_row$adj_p_value, 3)) else ""
  p <- ggplot(plot_df, aes(x = grp, y = expr, fill = grp)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    scale_fill_manual(values = c("CTx+IO" = "#B2182B", "Rest" = "#a6bddb"), guide = "none") +
    labs(title = paste0(gene, " — CTx+IO (n=22) vs Rest"), subtitle = subtitle_txt,
         x = NULL, y = "log2 concentration") +
    akila_boxplot_theme()
  stub <- file.path(FIG_DIR, paste0("headline_boxplot_", gene))
  cairo_pdf(paste0(stub, ".pdf"), width = 6, height = 5.5); print(p); dev.off()
  png(paste0(stub, ".png"), width = 6, height = 5.5, units = "in", res = 300, bg = "white"); print(p); dev.off()
}
message("Headline biomarker boxplots saved (", length(headline_ids), " proteins).")

# ============================================================
# PART 2 — Pathway enrichment + PPI input (the missing analysis;
# bug-fix retroactive run)
# ============================================================
message("\n=== 08y | PART 2: fgsea + STRING PPI input (bug-fix retroactive run) ===")

if (is.null(gene_sets) || !length(gene_sets)) {
  message("checkpoint_09 gene_sets not found — rebuilding from msigdbr (one-time cost)")
  gene_sets <- tryCatch(build_gene_sets(), error = function(e) { message("Gene set build failed: ", conditionMessage(e)); list() })
}
message("Collections available: ", paste(names(gene_sets), collapse = ", "))

run_fgsea_suite(ovr_res, OVR_KEY, PW_TBL_BASE, PW_FIG_BASE, stat_col = "t_stat")
run_string_ppi(ovr_res, OVR_KEY, PW_TBL_BASE, PW_FIG_BASE)

message("fgsea + PPI input now saved under: ", file.path(PW_TBL_BASE, OVR_KEY),
        " (previously silently skipped due to the key-name bug documented at the top of this script).")

# ============================================================
# PART 3 — Subgroup survival (OS, PFS): Cox screen + headline KM
# ============================================================
message("\n=== 08y | PART 3: Subgroup survival (n=22, hypothesis-generating only) ===")

os_cox_sub  <- add_feature_annotation(run_univariate_cox(expr_log2_rep, ctxio_df, "os_months", "os_event"), feature_meta)
pfs_cox_sub <- add_feature_annotation(run_univariate_cox(expr_log2_rep, ctxio_df, "pfs_months_main", "pfs_event_main"), feature_meta)
efs_cox_sub <- add_feature_annotation(run_univariate_cox(expr_log2_rep, ctxio_df, "efs_months", "efs_event"), feature_meta)

safe_write_csv(os_cox_sub  %>% arrange(p_value), file.path(TBL_DIR, "ctxio_os_univariate_cox.csv"))
safe_write_csv(pfs_cox_sub %>% arrange(p_value), file.path(TBL_DIR, "ctxio_pfs_univariate_cox.csv"))
safe_write_csv(efs_cox_sub %>% arrange(p_value), file.path(TBL_DIR, "ctxio_efs_univariate_cox.csv"))
message("OS Cox (subgroup):  ", nrow(os_cox_sub),  " proteins tested (nominal p only — FDR not meaningful at n=22)")
message("PFS Cox (subgroup): ", nrow(pfs_cox_sub), " proteins tested (nominal p only — FDR not meaningful at n=22)")
message("EFS Cox (subgroup): ", nrow(efs_cox_sub), " proteins tested (nominal p only — FDR not meaningful at n=22)")

# HR forest plots + survival volcanoes — OS, PFS, EFS (top 20 by nominal p each)
message("\n=== 08y | PART 3b: OS/PFS/EFS forest + survival volcano plots ===")
survival_endpoints <- list(
  OS  = list(df = os_cox_sub,  label = "OS"),
  PFS = list(df = pfs_cox_sub, label = "PFS"),
  EFS = list(df = efs_cox_sub, label = "EFS")
)
for (ep in names(survival_endpoints)) {
  spec <- survival_endpoints[[ep]]
  if (is.null(spec$df) || nrow(spec$df) == 0) {
    message("  ", ep, ": no Cox results — forest/volcano skipped (check n/events floor)")
    next
  }
  make_forest_plot(
    spec$df,
    title     = paste0("CTx+IO (n=22): ", spec$label, " Cox — top proteins by nominal p"),
    path_stub = file.path(FIG_DIR, paste0("ctxio_", tolower(spec$label), "_forest"))
  )
  make_survival_volcano(
    spec$df,
    title     = paste0("CTx+IO (n=22): ", spec$label, " survival volcano"),
    path_stub = file.path(FIG_DIR, paste0("ctxio_", tolower(spec$label), "_volcano"))
  )
  message("  ", ep, ": forest + volcano saved (", nrow(spec$df), " proteins tested, ",
          "nominal p<0.05: ", sum(spec$df$p_value < 0.05, na.rm = TRUE), ")")
}

km_summary_list <- list()
for (gene in names(headline_ids)) {
  id <- headline_ids[[gene]]
  res_os <- plot_km_median_split(ctxio_df, expr_log2_rep[ctxio_df$sample_id, id], "os_months", "os_event",
                                 gene, "CTx+IO subgroup | OS", file.path(FIG_DIR, paste0("ctxio_km_OS_", gene)))
  if (!is.null(res_os)) km_summary_list[[paste0(gene, "_OS")]] <- res_os %>% mutate(endpoint = "OS")
  res_pfs <- plot_km_median_split(ctxio_df, expr_log2_rep[ctxio_df$sample_id, id], "pfs_months_main", "pfs_event_main",
                                  gene, "CTx+IO subgroup | PFS", file.path(FIG_DIR, paste0("ctxio_km_PFS_", gene)))
  if (!is.null(res_pfs)) km_summary_list[[paste0(gene, "_PFS")]] <- res_pfs %>% mutate(endpoint = "PFS")
}
km_summary <- bind_rows(km_summary_list)
safe_write_csv(km_summary, file.path(TBL_DIR, "ctxio_headline_km_summary.csv"))
message("Headline biomarker KM checks within CTx+IO subgroup:")
print(km_summary)

# ============================================================
# PART 4 — Response comparison (headline biomarkers), if cell
# counts allow (n>=3 per response category)
# ============================================================
message("\n=== 08y | PART 4: Response comparison within CTx+IO (if n allows) ===")

n_prcr <- sum(ctxio_df$response_binary == "PR/CR", na.rm = TRUE)
n_sdpd <- sum(ctxio_df$response_binary == "SD/PD", na.rm = TRUE)
message("Response cells — PR/CR: ", n_prcr, " | SD/PD: ", n_sdpd)

response_wilcox <- tibble()
if (n_prcr >= 3 && n_sdpd >= 3) {
  response_wilcox <- bind_rows(lapply(names(headline_ids), function(gene) {
    id <- headline_ids[[gene]]
    run_group_wilcox(ctxio_df, expr_log2_rep[, id], "response_binary", gene)
  }))
  message("Response Wilcoxon tests run (nominal p, n=22 total — treat as descriptive):")
  print(response_wilcox)
} else {
  message("SKIPPED: response cell(s) below n=3 — reporting counts only, no test run.")
  response_wilcox <- tibble(protein = names(headline_ids), status = "SKIPPED: response cell n<3",
                            n_prcr = n_prcr, n_sdpd = n_sdpd)
}
safe_write_csv(response_wilcox, file.path(TBL_DIR, "ctxio_response_headline_wilcoxon.csv"))

# ============================================================
# PART 5 — Recurrence vs non-recurrence, WITHIN CTx+IO (n=22):
# full protein-wide limma DE + volcano
# ============================================================
message("\n=== 08y | PART 5: Recurrence vs non-recurrence DE (within CTx+IO) ===")

rec_counts <- ctxio_df %>% count(recurrence_binary_sensitivity)
message("Recurrence cells: "); print(rec_counts)
message("NOTE: unlike CTx+RTx (2 vs 13), CTx+IO's recurrence split is reasonably balanced",
        " (see counts above), so the standard pipeline min_per_group=3 floor is used here —",
        " no relaxation needed. Still n=22 total, so treat as hypothesis-generating for a",
        " manuscript, not a confirmed biomarker without external validation.")

recurrence_de <- run_binary_limma(
  expr_log2_rep, ctxio_df, "recurrence_binary_sensitivity",
  group_levels  = c("No recurrence", "Recurrence"),
  min_per_group = 3
)

if (is.null(recurrence_de) || nrow(recurrence_de) == 0) {
  message("Recurrence DE returned no results (insufficient non-missing values per feature).")
  recurrence_de <- tibble()
} else {
  recurrence_de <- add_feature_annotation(recurrence_de, feature_meta)
  message("Recurrence DE (within CTx+IO): ", nrow(recurrence_de), " proteins tested | ",
          "nominal p<0.05: ", sum(recurrence_de$p_value < 0.05, na.rm = TRUE),
          " | FDR<0.10: ", sum(recurrence_de$adj_p_value < FDR_CUTOFF, na.rm = TRUE))
  safe_write_csv(recurrence_de %>% arrange(p_value),
                file.path(TBL_DIR, "ctxio_recurrence_de_within_subgroup.csv"))

  make_volcano(
    recurrence_de,
    title     = "CTx+IO (n=22): Recurrence vs No recurrence",
    path_stub = file.path(FIG_DIR, "ctxio_recurrence_vs_no_recurrence_volcano"),
    group_a_label = "No recurrence", group_b_label = "Recurrence"
  )
  message("Recurrence volcano saved.")

  run_fgsea_suite(recurrence_de, "ctxio_recurrence_vs_no_recurrence", PW_TBL_BASE, PW_FIG_BASE, stat_col = "t_stat")
  run_string_ppi(recurrence_de, "ctxio_recurrence_vs_no_recurrence", PW_TBL_BASE, PW_FIG_BASE)
  message("Recurrence fgsea + PPI input saved under: ",
          file.path(PW_TBL_BASE, "ctxio_recurrence_vs_no_recurrence"))
}

# ============================================================
# PART 6 — Responder vs non-responder, WITHIN CTx+IO
# (counts determined at runtime below): full protein-wide limma
# DE + volcano + fgsea
# ============================================================
message("\n=== 08y | PART 6: Response (PR/CR vs SD/PD) DE (within CTx+IO) ===")
message("Cell sizes (PR/CR=", n_prcr, ", SD/PD=", n_sdpd, ") checked against",
        " min_per_group=3 below. Still n=22 total, so treat as hypothesis-generating for a",
        " manuscript, not a confirmed biomarker without external validation.")

response_de <- run_binary_limma(
  expr_log2_rep, ctxio_df, "response_binary",
  group_levels  = c("SD/PD", "PR/CR"),
  min_per_group = 3
)

if (is.null(response_de) || nrow(response_de) == 0) {
  message("Response DE returned no results.")
  response_de <- tibble()
} else {
  response_de <- add_feature_annotation(response_de, feature_meta)
  message("Response DE (within CTx+IO): ", nrow(response_de), " proteins tested | ",
          "nominal p<0.05: ", sum(response_de$p_value < 0.05, na.rm = TRUE),
          " | FDR<0.10: ", sum(response_de$adj_p_value < FDR_CUTOFF, na.rm = TRUE))
  safe_write_csv(response_de %>% arrange(p_value),
                file.path(TBL_DIR, "ctxio_response_de_within_subgroup.csv"))

  make_volcano(
    response_de,
    title     = "CTx+IO (n=22): PR/CR vs SD/PD",
    path_stub = file.path(FIG_DIR, "ctxio_response_prcr_vs_sdpd_volcano"),
    group_a_label = "SD/PD", group_b_label = "PR/CR"
  )
  message("Response volcano saved.")

  run_fgsea_suite(response_de, "ctxio_response_prcr_vs_sdpd", PW_TBL_BASE, PW_FIG_BASE, stat_col = "t_stat")
  run_string_ppi(response_de, "ctxio_response_prcr_vs_sdpd", PW_TBL_BASE, PW_FIG_BASE)
  message("Response fgsea + PPI input saved under: ",
          file.path(PW_TBL_BASE, "ctxio_response_prcr_vs_sdpd"))
}

# ============================================================
# PART 7 — New biomarker candidates: cross-method ranking
# Combines: one-vs-rest DE (Part 1), OS/PFS Cox (Part 3),
# recurrence DE (Part 5), response DE (Part 6) into one table,
# scored by how many analyses flag each protein at nominal p<0.05.
# This is the "look beyond the six headline proteins" deliverable —
# ranked by n_methods_flagged, independent of whether a protein
# was already known from whole-cohort analyses.
# ============================================================
message("\n=== 08y | PART 7: Cross-method new-candidate ranking (CTx+IO-specific) ===")

flag_tbl <- function(df, p_col = "p_value", suffix) {
  pcol_name <- paste0("p_", suffix)
  fcol_name <- paste0("flag_", suffix)
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(nELISA_ID = character(), !!pcol_name := numeric(), !!fcol_name := logical()))
  }
  df %>% transmute(nELISA_ID, !!pcol_name := .data[[p_col]], !!fcol_name := .data[[p_col]] < 0.05)
}

candidate_tbl <- feature_meta %>% distinct(nELISA_ID, Gene_Symbol, Protein_Label, Uniprot_ID) %>%
  left_join(flag_tbl(ovr_res,          suffix = "ovr_vs_rest"), by = "nELISA_ID") %>%
  left_join(flag_tbl(os_cox_sub,       suffix = "os_cox"),      by = "nELISA_ID") %>%
  left_join(flag_tbl(pfs_cox_sub,      suffix = "pfs_cox"),     by = "nELISA_ID") %>%
  left_join(flag_tbl(efs_cox_sub,      suffix = "efs_cox"),     by = "nELISA_ID") %>%
  left_join(flag_tbl(recurrence_de,    suffix = "recurrence"),  by = "nELISA_ID") %>%
  left_join(flag_tbl(response_de,      suffix = "response"),    by = "nELISA_ID")

flag_cols <- grep("^flag_", names(candidate_tbl), value = TRUE)
candidate_tbl <- candidate_tbl %>%
  mutate(n_methods_flagged = rowSums(across(all_of(flag_cols), ~ coalesce(., FALSE))),
         is_headline_six   = Gene_Symbol %in% HEADLINE_GENES) %>%
  arrange(desc(n_methods_flagged), p_ovr_vs_rest)

safe_write_csv(candidate_tbl, file.path(TBL_DIR, "ctxio_new_candidate_biomarkers_ranked.csv"))

new_candidates_top <- candidate_tbl %>%
  filter(!is_headline_six, n_methods_flagged >= 2) %>%
  arrange(desc(n_methods_flagged), p_ovr_vs_rest)

