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
if (file.exists(ck09)) { load(ck09); message("Checkpoint 09 loaded (gene_sets reused)") }

while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
FDR_CUTOFF     <- 0.10
FGSEA_FDR      <- 0.25
MIN_GENE_SET   <- 5
MAX_GENE_SET   <- 500
TOP_PLOT       <- 10
EMAP_TOP       <- 10
RNG_SEED       <- 42
FOREST_N       <- 20
MIN_GROUP_SIZE <- 3     # DE contrasts
STAGE_COX_MIN_N      <- 30   # stratified (early-only / late-only) Cox floor
STAGE_COX_MIN_EVENTS <- 10
INTERACTION_MIN_N      <- 40  # whole-cohort interaction model floor
INTERACTION_MIN_EVENTS <- 15
KM_MIN_N       <- 10
KM_MIN_EVENTS  <- 4
HEADLINE_GENES <- c("GDF15", "GZMA", "NOTUM", "SHH", "CEACAM5", "PLAUR")

TBL_DIR <- file.path(paths$tables,  "13_stage_early_vs_late")
FIG_DIR <- file.path(paths$figures, "13_stage_early_vs_late")
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

get_gene_id <- function(gene_symbol, feature_meta) {
  ids <- feature_meta %>% filter(Gene_Symbol == gene_symbol) %>% pull(nELISA_ID) %>% unique()
  if (length(ids) != 1) {
    message("  WARN: expected exactly one nELISA_ID for ", gene_symbol, ", found ", length(ids), " — skipping")
    return(NA_character_)
  }
  ids
}

# Full volcano theme (matches 03_analysis1_seven_group.R / 07_analysis5_survival.R exactly)
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

# Simpler theme for boxplots (matches 08k_ceacam5_tki_sensitivity_check.R)
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

# ── Binary limma (copied from 03_analysis1_seven_group.R) ──────
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
  design <- tryCatch(model.matrix(~ ., data = design_df), error = function(e) NULL)
  if (is.null(design) || qr(design)$rank < ncol(design)) {
    message("    WARN: design matrix rank-deficient (likely near-perfect collinearity with adjust_cols) — skipped")
    return(NULL)
  }
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
    group_a      = group_levels[1], group_b = group_levels[2],
    mean_group_a = unname(means[[1]][match(nELISA_ID, names(means[[1]]))]),
    mean_group_b = unname(means[[2]][match(nELISA_ID, names(means[[2]]))]),
    n_group_a    = unname(n_nonmiss[[1]][match(nELISA_ID, names(n_nonmiss[[1]]))]),
    n_group_b    = unname(n_nonmiss[[2]][match(nELISA_ID, names(n_nonmiss[[2]]))])
  )
}

# ── Volcano (copied from 03_analysis1_seven_group.R) ────────────
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
        TRUE ~ "Not significant"),
      sig_level = case_when(
        status == "Not significant" ~ "Not significant",
        !is.na(adj_p_value) & adj_p_value <= 0.05 ~ "Primary (FDR<=0.05)",
        TRUE ~ secondary_lbl))
  ann <- plot_df %>% filter(status != "Not significant") %>%
    arrange(adj_p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  if (nrow(ann) == 0) ann <- plot_df %>% arrange(adj_p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  x_rng <- range(plot_df$logFC, na.rm = TRUE); y_thr <- -log10(padj_thr)
  y_max <- max(plot_df$neg_log10, y_thr, na.rm = TRUE)
  x_pad <- max(0.35, diff(x_rng) * 0.22); y_pad <- max(0.18, y_max * 0.16)
  x_panel_left  <- x_rng[1] - 0.45 * x_pad; x_panel_right <- x_rng[2] + 0.45 * x_pad
  y_panel_top   <- y_max  + 0.12 * y_pad;   y_label <- y_panel_top + 0.9 * y_pad
  p <- ggplot(plot_df, aes(x = logFC, y = neg_log10)) +
    geom_point(aes(color = status, shape = sig_level), size = 2.4, alpha = 0.82, stroke = 0.7) +
    geom_vline(xintercept = c(-lfc_thr, 0, lfc_thr), linetype = c("dashed", "solid", "dashed"),
               linewidth = c(0.6, 0.4, 0.6), color = "black") +
    geom_hline(yintercept = y_thr, linetype = "dashed", linewidth = 0.6, color = "black") +
    annotate("text", x = x_panel_left - 0.18 * x_pad, y = y_label, label = paste0("Higher in\n", group_a_label),
             hjust = 0, vjust = 0, size = 4.8, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_panel_right + 0.18 * x_pad, y = y_label, label = paste0("Higher in\n", group_b_label),
             hjust = 1, vjust = 0, size = 4.8, fontface = "bold", color = "#B2182B") +
    annotate("text", x = x_rng[1] + 0.02 * diff(x_rng), y = y_thr,
             label = paste0(" P Adj. = ", formatC(padj_thr, format = "f", digits = 2)),
             hjust = 0, vjust = -0.6, size = 3.8, fontface = "italic", color = "black") +
    scale_color_manual(values = c(setNames("#2166AC", left_label), setNames("#B2182B", right_label),
                                  "Not significant" = "grey85"), drop = FALSE) +
    scale_shape_manual(values = c("Not significant" = 1, "Primary (FDR<=0.05)" = 19, setNames(17, secondary_lbl)),
                       breaks = c("Not significant", "Primary (FDR<=0.05)", secondary_lbl), drop = FALSE) +
    ggrepel::geom_text_repel(data = ann, aes(label = .data[[label_col]]), size = 3.1, fontface = "bold",
                             max.overlaps = 30, box.padding = 0.5, point.padding = 0.2,
                             segment.color = "grey60", segment.size = 0.2, show.legend = FALSE) +
    labs(title = title, x = x_label, y = "P value (-log10)", color = NULL, shape = "Significance") +
    guides(color = "none", shape = guide_legend(override.aes = list(color = "black"))) +
    coord_cartesian(xlim = c(x_panel_left, x_panel_right), ylim = c(0, y_panel_top), clip = "off") +
    akila_volcano_theme()
  cairo_pdf(paste0(path_stub, ".pdf"), width = 8.6, height = 8.2); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 8.6, height = 8.2, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# ── HR forest plot (copied from 07_analysis5_survival.R) ────────
make_forest_plot <- function(df, title, path_stub, n = FOREST_N) {
  if (nrow(df) == 0) return(invisible(NULL))
  plot_df <- df %>% filter(is.finite(HR), is.finite(HR_low), is.finite(HR_high)) %>%
    arrange(p_value) %>% slice_head(n = min(n, nrow(.))) %>%
    mutate(label = factor(paste0(Protein_Label, " (events=", n_events, ")"),
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

# ── Survival volcano (copied from 07_analysis5_survival.R) ──────
make_survival_volcano <- function(df, title, path_stub, label_col = "Protein_Label",
                                  n = FOREST_N, padj_thr = FDR_CUTOFF) {
  if (nrow(df) == 0) return(invisible(NULL))
  secondary_label <- paste0("Secondary (FDR<=", formatC(padj_thr, format="f", digits=2), ")")
  plot_df <- df %>% filter(is.finite(HR), HR > 0, is.finite(p_value)) %>%
    mutate(log2_hr = log2(HR), neg_log10 = -log10(pmax(p_value, 1e-300)),
           direction = case_when(p_value < 0.05 & log2_hr > 0 ~ "Higher risk",
                                  p_value < 0.05 & log2_hr < 0 ~ "Protective", TRUE ~ "Not significant"),
           sig_level = case_when(direction == "Not significant" ~ "Not significant",
                                  p_value <= 0.05 ~ "Primary (P<=0.05)", TRUE ~ secondary_label))
  if (nrow(plot_df) == 0) return(invisible(NULL))
  ann <- plot_df %>% filter(direction != "Not significant") %>%
    arrange(p_value, desc(abs(log2_hr))) %>% slice_head(n = min(n, nrow(.)))
  if (nrow(ann) == 0) ann <- plot_df %>% arrange(p_value, desc(abs(log2_hr))) %>% slice_head(n = min(n, nrow(.)))
  x_rng <- range(plot_df$log2_hr, na.rm=TRUE); y_thr <- -log10(0.05)
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
                               min_n = STAGE_COX_MIN_N, min_events = STAGE_COX_MIN_EVENTS) {
  op  <- options(warn = -1); time_val <- as.numeric(sample_df[[time_col]]); options(op)
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

# ── Protein x stage interaction Cox screen (NEW — the "predictive"/
# effect-modification test: does a protein's hazard association
# differ significantly BETWEEN early and late stage?) ───────────
run_cox_interaction <- function(mat_samples_by_features, sample_df, time_col, event_col,
                                min_n = INTERACTION_MIN_N, min_events = INTERACTION_MIN_EVENTS) {
  op <- options(warn = -1); time_val <- as.numeric(sample_df[[time_col]]); options(op)
  df <- sample_df %>%
    mutate(time_value = time_val, event_value = coerce_event_indicator(.data[[event_col]]),
           stage_f = factor(stage_collapsed, levels = c("I_II", "III_IV"))) %>%
    filter(!is.na(time_value), !is.na(event_value), !is.na(stage_f))
  if (nrow(df) < min_n || sum(df$event_value == 1, na.rm = TRUE) < min_events) return(tibble())
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  results <- lapply(colnames(mat_sub), function(feature) {
    x <- mat_sub[, feature]
    keep <- !is.na(x) & is.finite(x)
    if (sum(keep) < min_n) return(NULL)
    dsub <- df[keep, , drop = FALSE]
    if (sum(dsub$event_value == 1) < min_events || length(unique(dsub$stage_f)) < 2) return(NULL)
    xs <- as.numeric(scale(x[keep]))
    fit <- tryCatch(
      survival::coxph(survival::Surv(dsub$time_value, dsub$event_value) ~ xs * dsub$stage_f),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s  <- summary(fit); cn <- rownames(s$coef)
    int_row <- grep(":", cn, value = TRUE)
    main_row <- grep("^xs$", cn, value = TRUE)
    if (!length(int_row) || !length(main_row)) return(NULL)
    ir <- int_row[1]; mr <- main_row[1]
    tibble(
      nELISA_ID        = feature,
      main_HR_early    = unname(exp(s$coef[mr, "coef"])),
      main_p_early     = unname(s$coef[mr, "Pr(>|z|)"]),
      HR               = unname(exp(s$coef[ir, "coef"])),      # interaction HR (named for volcano/forest reuse)
      HR_low           = unname(s$conf.int[ir, "lower .95"]),
      HR_high          = unname(s$conf.int[ir, "upper .95"]),
      p_value          = unname(s$coef[ir, "Pr(>|z|)"]),       # interaction p (named for volcano/forest reuse)
      n                = nrow(dsub), n_events = sum(dsub$event_value == 1)
    )
  })
  bind_rows(results) %>% mutate(adj_p_value = p.adjust(p_value, method = "BH"))
}

plot_km_median_split <- function(sample_df, expr_vec, time_col, event_col, protein_label,
                                 subtitle_extra, path_stub, min_n = KM_MIN_N, min_events = KM_MIN_EVENTS) {
  op <- options(warn = -1); time_val <- as.numeric(sample_df[[time_col]]); options(op)
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
    title = paste0("KM \u2014 ", protein_label),
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
  tibble::tibble(protein = protein_label, n = nrow(df), n_events = sum(df$event_value == 1),
                 n_low = sum(df$expr_group == "Low"), n_high = sum(df$expr_group == "High"), logrank_p = logrank_p)
}

# ── fgsea machinery (copied from 09_pathway_network.R, unchanged) ─
clean_pathway_name <- function(x) {
  x %>% stringr::str_replace_all("^HALLMARK_", "") %>% stringr::str_replace_all("^REACTOME_", "") %>%
    stringr::str_replace_all("^GOBP_|^GO:BP_", "") %>% stringr::str_replace_all("^GOMF_|^GO:MF_", "") %>%
    stringr::str_replace_all("^GOCC_|^GO:CC_", "") %>% stringr::str_replace_all("^KEGG_(MEDICUS|LEGACY)_", "") %>%
    stringr::str_replace_all("^KEGG_", "") %>% stringr::str_replace_all("_", " ") %>% stringr::str_squish()
}
get_msigdbr_table <- function(collection, subcollection = NULL) {
  fn_args <- names(formals(msigdbr::msigdbr)); args <- list(species = "Homo sapiens")
  if ("collection" %in% fn_args) { args$collection <- collection; if (!is.null(subcollection)) args$subcollection <- subcollection }
  else { args$category <- collection; if (!is.null(subcollection)) args$subcategory <- subcollection }
  do.call(msigdbr::msigdbr, args)
}
build_gene_sets <- function() {
  kegg_medicus <- tryCatch(get_msigdbr_table("C2", "CP:KEGG_MEDICUS"), error = function(e) tibble())
  kegg_legacy  <- tryCatch(get_msigdbr_table("C2", "CP:KEGG_LEGACY"),  error = function(e) tibble())
  sets <- list(HALLMARK = get_msigdbr_table("H"), GO_BP = get_msigdbr_table("C5", "GO:BP"),
               GO_MF = get_msigdbr_table("C5", "GO:MF"), GO_CC = get_msigdbr_table("C5", "GO:CC"),
               KEGG = bind_rows(kegg_medicus %>% mutate(gs_name = paste0("KEGG_MEDICUS_", gs_name)),
                                kegg_legacy  %>% mutate(gs_name = paste0("KEGG_LEGACY_",  gs_name))),
               REACTOME = get_msigdbr_table("C2", "CP:REACTOME"))
  out <- lapply(sets, function(df) {
    if (is.null(df) || nrow(df) == 0) return(list())
    df %>% select(gs_name, gene_symbol) %>% filter(!is.na(gs_name), !is.na(gene_symbol), gene_symbol != "") %>%
      distinct() %>% split(.$gs_name) %>% lapply(function(tbl) unique(tbl$gene_symbol))
  })
  out[lengths(out) > 0]
}
prepare_rank_vector <- function(result_df, stat_col = "t_stat") {
  rank_df <- result_df %>% filter(!is.na(Gene_Symbol), !is.na(.data[[stat_col]]), is.finite(.data[[stat_col]])) %>%
    group_by(Gene_Symbol) %>% summarise(stat = .data[[stat_col]][which.max(abs(.data[[stat_col]]))], .groups = "drop") %>%
    filter(!is.na(stat), is.finite(stat))
  stats <- rank_df$stat; names(stats) <- rank_df$Gene_Symbol; sort(stats, decreasing = TRUE)
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
# PART 0 — Snapshot: stage counts + the confounding cross-tabs
# ============================================================
message("\n=== 13 | PART 0: Stage snapshot ===")

stage_df <- analysis1_df %>% filter(!is.na(stage_collapsed))
message("Stage-evaluable baseline n = ", nrow(stage_df), " (of ", nrow(analysis1_df), " total)")

stage_counts <- stage_df %>% count(stage_collapsed, name = "n")
safe_write_csv(stage_counts, file.path(TBL_DIR, "stage_counts.csv"))
print(stage_counts)

stage_by_group <- stage_df %>% count(new_group, stage_collapsed) %>%
  tidyr::pivot_wider(names_from = stage_collapsed, values_from = n, values_fill = 0)
safe_write_csv(stage_by_group, file.path(TBL_DIR, "stage_by_treatment_group_crosstab.csv"))
message("\nStage x treatment-group crosstab (this is the confounding — note how concentrated",
        " Surgery is in I_II):")
print(stage_by_group)

stage_by_response <- stage_df %>% filter(!is.na(response_binary)) %>% count(stage_collapsed, response_binary) %>%
  tidyr::pivot_wider(names_from = response_binary, values_from = n, values_fill = 0)
safe_write_csv(stage_by_response, file.path(TBL_DIR, "stage_by_response_crosstab.csv"))
message("\nStage x response crosstab:"); print(stage_by_response)

surv_snapshot <- tibble(
  n_os_evaluable  = sum(!is.na(stage_df$os_months)  & !is.na(stage_df$os_event)),
  n_pfs_evaluable = sum(!is.na(stage_df$pfs_months_main) & !is.na(stage_df$pfs_event_main)),
  n_efs_evaluable = sum(!is.na(stage_df$efs_months) & !is.na(stage_df$efs_event))
)
safe_write_csv(surv_snapshot, file.path(TBL_DIR, "stage_survival_evaluable_snapshot.csv"))
print(surv_snapshot)

headline_ids <- setNames(vapply(HEADLINE_GENES, get_gene_id, character(1), feature_meta = feature_meta), HEADLINE_GENES)
headline_ids <- headline_ids[!is.na(headline_ids)]

# ============================================================
# PART 1 — Early vs late stage DE: crude, treatment-adjusted,
# and Surgery-excluded (see confounding warning in header)
# ============================================================
message("\n=== 13 | PART 1: Early (I_II) vs Late (III_IV) DE — three versions ===")

# 1a. Crude — reflects mostly WHO GETS SURGERY, flagged explicitly
de_crude <- run_binary_limma(expr_log2_rep, stage_df, "stage_collapsed",
                             group_levels = c("I_II", "III_IV"), min_per_group = MIN_GROUP_SIZE)
if (!is.null(de_crude)) {
  de_crude <- add_feature_annotation(de_crude, feature_meta)
  message("Crude stage DE: ", nrow(de_crude), " proteins | FDR<0.10: ", sum(de_crude$adj_p_value < FDR_CUTOFF, na.rm = TRUE),
          " — CAUTION: unadjusted for treatment group, expect this to largely reflect Surgery-vs-rest.")
  safe_write_csv(de_crude %>% arrange(p_value), file.path(TBL_DIR, "stage_de_crude.csv"))
  make_volcano(de_crude, title = "Stage I-II vs III-IV — CRUDE (confounded with treatment, see caution)",
              path_stub = file.path(FIG_DIR, "stage_volcano_crude"),
              group_a_label = "I-II", group_b_label = "III-IV")
}

# 1b. Adjusted for treatment group (new_group) — may be rank-deficient
# given severe collinearity; run_binary_limma checks this and skips cleanly.
de_adjusted <- run_binary_limma(expr_log2_rep, stage_df, "stage_collapsed",
                                group_levels = c("I_II", "III_IV"),
                                adjust_cols = c("new_group", "plate_barcode"), min_per_group = MIN_GROUP_SIZE)
if (!is.null(de_adjusted)) {
  de_adjusted <- add_feature_annotation(de_adjusted, feature_meta)
  message("Treatment-adjusted stage DE: ", nrow(de_adjusted), " proteins | FDR<0.10: ",
          sum(de_adjusted$adj_p_value < FDR_CUTOFF, na.rm = TRUE))
  safe_write_csv(de_adjusted %>% arrange(p_value), file.path(TBL_DIR, "stage_de_adjusted_for_treatment.csv"))
  make_volcano(de_adjusted, title = "Stage I-II vs III-IV — adjusted for treatment group",
              path_stub = file.path(FIG_DIR, "stage_volcano_adjusted_for_treatment"),
              group_a_label = "I-II", group_b_label = "III-IV")
} else {
  message("Treatment-adjusted stage DE not estimable (design rank-deficient — near-perfect collinearity",
          " between stage and treatment group). Use the Surgery-excluded version (1c) instead.")
}

# 1c. Surgery excluded — removes the single biggest confound directly
de_no_surgery <- run_binary_limma(expr_log2_rep, stage_df %>% filter(new_group != "Sx"), "stage_collapsed",
                                  group_levels = c("I_II", "III_IV"), min_per_group = MIN_GROUP_SIZE)
if (!is.null(de_no_surgery)) {
  de_no_surgery <- add_feature_annotation(de_no_surgery, feature_meta)
  message("Surgery-excluded stage DE: ", nrow(de_no_surgery), " proteins | FDR<0.10: ",
          sum(de_no_surgery$adj_p_value < FDR_CUTOFF, na.rm = TRUE),
          " — this is the most defensible version for stage-intrinsic biology.")
  safe_write_csv(de_no_surgery %>% arrange(p_value), file.path(TBL_DIR, "stage_de_surgery_excluded.csv"))
  make_volcano(de_no_surgery, title = "Stage I-II vs III-IV — Surgery patients excluded",
              path_stub = file.path(FIG_DIR, "stage_volcano_surgery_excluded"),
              group_a_label = "I-II", group_b_label = "III-IV")
}

# ============================================================
# PART 2 — fgsea + PPI input for crude and Surgery-excluded DE
# ============================================================
message("\n=== 13 | PART 2: fgsea + STRING PPI input ===")
if (is.null(gene_sets) || !length(gene_sets)) {
  message("checkpoint_09 gene_sets not found — rebuilding from msigdbr (one-time cost)")
  gene_sets <- tryCatch(build_gene_sets(), error = function(e) { message("Gene set build failed: ", conditionMessage(e)); list() })
}
if (!is.null(de_crude)) {
  run_fgsea_suite(de_crude, "stage_early_vs_late_crude", PW_TBL_BASE, PW_FIG_BASE, stat_col = "t_stat")
  run_string_ppi(de_crude, "stage_early_vs_late_crude", PW_TBL_BASE, PW_FIG_BASE)
}
if (!is.null(de_no_surgery)) {
  run_fgsea_suite(de_no_surgery, "stage_early_vs_late_surgery_excluded", PW_TBL_BASE, PW_FIG_BASE, stat_col = "t_stat")
  run_string_ppi(de_no_surgery, "stage_early_vs_late_surgery_excluded", PW_TBL_BASE, PW_FIG_BASE)
}
message("fgsea/PPI saved under: ", file.path(PW_TBL_BASE, "stage_early_vs_late_crude"), " and ",
        file.path(PW_TBL_BASE, "stage_early_vs_late_surgery_excluded"))

# ============================================================
# PART 3 — How stage behaves for OS/PFS/EFS: sanity KM, stage as
# a single Cox covariate, then full stratified (early-only /
# late-only) protein-wide Cox screens with forest + volcano
# ============================================================
message("\n=== 13 | PART 3: Stage vs survival ===")

# 3a. KM sanity check — stage itself vs OS/PFS/EFS (expected: worse for III_IV)
km_stage_summary <- list()
for (ep in list(list(time = "os_months", event = "os_event", tag = "OS"),
                list(time = "pfs_months_main", event = "pfs_event_main", tag = "PFS"),
                list(time = "efs_months", event = "efs_event", tag = "EFS"))) {
  df_ep <- stage_df %>% filter(!is.na(.data[[ep$time]]), !is.na(.data[[ep$event]]))
  time_val <- as.numeric(df_ep[[ep$time]]); event_val <- coerce_event_indicator(df_ep[[ep$event]])
  keep <- !is.na(time_val) & !is.na(event_val)
  df_r <- df_ep[keep, ] %>% mutate(time_value = time_val[keep], event_value = event_val[keep],
                                   stage_f = factor(stage_collapsed, levels = c("I_II", "III_IV")))
  fit_r <- survival::survfit(survival::Surv(time_value, event_value) ~ stage_f, data = df_r)
  diff_r <- tryCatch(survival::survdiff(survival::Surv(time_value, event_value) ~ stage_f, data = df_r), error = function(e) NULL)
  p_r <- if (!is.null(diff_r)) stats::pchisq(diff_r$chisq, df = max(length(diff_r$n) - 1, 1), lower.tail = FALSE) else NA_real_
  km_r <- survminer::ggsurvplot(fit_r, data = df_r, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    palette = c("#2166AC", "#B2182B"), legend.title = "Stage", legend.labs = c("I-II", "III-IV"),
    title = paste0("KM \u2014 Stage sanity check (", ep$tag, ")"),
    subtitle = paste0("I-II vs III-IV | logrank p=", scales::number(p_r, accuracy = 0.001)),
    xlab = "Months", ylab = "Survival probability", ggtheme = ggplot2::theme_bw(base_size = 11))
  stub_r <- file.path(FIG_DIR, paste0("stage_sanity_km_", ep$tag))
  cairo_pdf(paste0(stub_r, ".pdf"), width = 7.8, height = 6.5); print(km_r); dev.off()
  png(paste0(stub_r, ".png"), width = 7.8, height = 6.5, units = "in", res = 300, bg = "white"); print(km_r); dev.off()
  km_stage_summary[[ep$tag]] <- tibble(endpoint = ep$tag, n = nrow(df_r), n_events = sum(df_r$event_value == 1), logrank_p = p_r)

  # 3b. stage as a single Cox covariate (unadjusted HR for stage itself)
  fit_cox <- tryCatch(survival::coxph(survival::Surv(time_value, event_value) ~ stage_f, data = df_r), error = function(e) NULL)
  if (!is.null(fit_cox)) {
    sc <- summary(fit_cox)
    km_stage_summary[[paste0(ep$tag, "_cox")]] <- tibble(
      endpoint = paste0(ep$tag, "_stage_HR"), n = nrow(df_r), n_events = sum(df_r$event_value == 1),
      HR = unname(sc$coef[1, "exp(coef)"]), HR_low = unname(sc$conf.int[1, "lower .95"]),
      HR_high = unname(sc$conf.int[1, "upper .95"]), p_value = unname(sc$coef[1, "Pr(>|z|)"]), logrank_p = NA_real_)
  }
}
km_stage_tbl <- bind_rows(km_stage_summary)
safe_write_csv(km_stage_tbl, file.path(TBL_DIR, "stage_survival_sanity_and_cox_hr.csv"))
message("Stage sanity KM + stage-as-covariate Cox HR:"); print(km_stage_tbl)

# 3c. Stratified protein-wide Cox screens: EARLY-ONLY and LATE-ONLY,
# for OS, PFS, EFS — forest + volcano each. This is the "early
# prognostic markers" deliverable.
message("\n=== 13 | PART 3c: Stratified (early-only / late-only) Cox screens ===")
early_df <- stage_df %>% filter(stage_collapsed == "I_II")
late_df  <- stage_df %>% filter(stage_collapsed == "III_IV")
message("Early-stage (I_II) n = ", nrow(early_df), " | Late-stage (III_IV) n = ", nrow(late_df))

stratified_cox <- list()
for (strat in list(list(df = early_df, tag = "early"), list(df = late_df, tag = "late"))) {
  for (ep in list(list(time = "os_months", event = "os_event", label = "OS"),
                  list(time = "pfs_months_main", event = "pfs_event_main", label = "PFS"),
                  list(time = "efs_months", event = "efs_event", label = "EFS"))) {
    res <- add_feature_annotation(run_univariate_cox(expr_log2_rep, strat$df, ep$time, ep$event), feature_meta)
    key <- paste0(strat$tag, "_", ep$label)
    stratified_cox[[key]] <- res
    if (nrow(res) == 0) { message("  ", key, ": skipped (insufficient n/events)"); next }
    safe_write_csv(res %>% arrange(p_value), file.path(TBL_DIR, paste0("stage_", key, "_cox.csv")))
    make_forest_plot(res, title = paste0(strat$tag, "-stage (n=", nrow(strat$df), "): ", ep$label, " Cox — top proteins"),
                     path_stub = file.path(FIG_DIR, paste0("stage_", key, "_forest")))
    make_survival_volcano(res, title = paste0(strat$tag, "-stage (n=", nrow(strat$df), "): ", ep$label, " survival volcano"),
                          path_stub = file.path(FIG_DIR, paste0("stage_", key, "_volcano")))
    message("  ", key, ": ", nrow(res), " proteins tested | nominal p<0.05: ", sum(res$p_value < 0.05, na.rm = TRUE))
  }
}

# ============================================================
# PART 4 — Protein x stage interaction models (whole cohort):
# the formally correct "predictive"/effect-modification test
# ============================================================
message("\n=== 13 | PART 4: Protein x stage interaction Cox models ===")

interaction_results <- list()
for (ep in list(list(time = "os_months", event = "os_event", label = "OS"),
                list(time = "pfs_months_main", event = "pfs_event_main", label = "PFS"),
                list(time = "efs_months", event = "efs_event", label = "EFS"))) {
  res_int <- add_feature_annotation(run_cox_interaction(expr_log2_rep, stage_df, ep$time, ep$event), feature_meta)
  interaction_results[[ep$label]] <- res_int
  if (nrow(res_int) == 0) { message("  ", ep$label, " interaction: skipped (insufficient n/events)"); next }
  safe_write_csv(res_int %>% arrange(p_value), file.path(TBL_DIR, paste0("stage_interaction_", ep$label, "_cox.csv")))
  make_forest_plot(res_int, title = paste0("Protein x Stage interaction — ", ep$label,
                                           " (HR>1: protein's hazard association is stronger in III-IV)"),
                   path_stub = file.path(FIG_DIR, paste0("stage_interaction_", ep$label, "_forest")))
  make_survival_volcano(res_int, title = paste0("Protein x Stage interaction volcano — ", ep$label),
                        path_stub = file.path(FIG_DIR, paste0("stage_interaction_", ep$label, "_volcano")))
  message("  ", ep$label, " interaction: ", nrow(res_int), " proteins tested | nominal p<0.05: ",
          sum(res_int$p_value < 0.05, na.rm = TRUE),
          " — interaction models need power in BOTH main effect and interaction term;",
          " treat as exploratory even more so than the stratified screens above.")
}

# ============================================================
# PART 5 — Headline biomarkers by stage: DE extraction, boxplots,
# and stage-stratified KM
# ============================================================
message("\n=== 13 | PART 5: Headline biomarkers by stage ===")

headline_stage_de <- if (!is.null(de_crude)) de_crude %>% filter(nELISA_ID %in% headline_ids) %>% arrange(p_value) else tibble()
safe_write_csv(headline_stage_de, file.path(TBL_DIR, "stage_headline_biomarker_de_crude.csv"))
message("Headline biomarkers in crude stage DE:")
print(headline_stage_de %>% select(Gene_Symbol, logFC, p_value, adj_p_value))

for (gene in names(headline_ids)) {
  id <- headline_ids[[gene]]
  plot_df <- stage_df %>%
    mutate(expr = expr_log2_rep[sample_id, id],
           stage_lbl = ifelse(stage_collapsed == "I_II", "I-II", "III-IV")) %>%
    filter(!is.na(expr))
  stat_row <- headline_stage_de %>% filter(Gene_Symbol == gene)
  subtitle_txt <- if (nrow(stat_row)) paste0("Crude: logFC=", round(stat_row$logFC, 2), " | p=", signif(stat_row$p_value, 3),
                                             " | FDR=", signif(stat_row$adj_p_value, 3)) else ""
  p <- ggplot(plot_df, aes(x = stage_lbl, y = expr, fill = stage_lbl)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    scale_fill_manual(values = c("I-II" = "#2166AC", "III-IV" = "#B2182B"), guide = "none") +
    labs(title = paste0(gene, " — Stage I-II vs III-IV"), subtitle = subtitle_txt, x = NULL, y = "log2 concentration") +
    akila_boxplot_theme()
  stub <- file.path(FIG_DIR, paste0("headline_boxplot_", gene))
  cairo_pdf(paste0(stub, ".pdf"), width = 6, height = 5.5); print(p); dev.off()
  png(paste0(stub, ".png"), width = 6, height = 5.5, units = "in", res = 300, bg = "white"); print(p); dev.off()

  for (ep in list(list(time = "os_months", event = "os_event", label = "OS"),
                  list(time = "pfs_months_main", event = "pfs_event_main", label = "PFS"),
                  list(time = "efs_months", event = "efs_event", label = "EFS"))) {
    plot_km_median_split(early_df, expr_log2_rep[early_df$sample_id, id], ep$time, ep$event, gene,
                         paste0("Early stage (I-II) | ", ep$label),
                         file.path(FIG_DIR, paste0("headline_km_early_", ep$label, "_", gene)))
    plot_km_median_split(late_df, expr_log2_rep[late_df$sample_id, id], ep$time, ep$event, gene,
                         paste0("Late stage (III-IV) | ", ep$label),
                         file.path(FIG_DIR, paste0("headline_km_late_", ep$label, "_", gene)))
  }
}
message("Headline biomarker boxplots + stage-stratified KM saved (", length(headline_ids), " proteins).")

# ============================================================
# PART 6 — Candidate ranking: early-specific prognostic markers
# and predictive (interaction-significant) markers
# ============================================================
message("\n=== 13 | PART 6: Candidate ranking ===")

flag_tbl <- function(df, p_col = "p_value", suffix) {
  pcol_name <- paste0("p_", suffix)
  fcol_name <- paste0("flag_", suffix)
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(nELISA_ID = character(), !!pcol_name := numeric(), !!fcol_name := logical()))
  }
  df %>% transmute(nELISA_ID, !!pcol_name := .data[[p_col]], !!fcol_name := .data[[p_col]] < 0.05)
}

candidate_tbl <- feature_meta %>% distinct(nELISA_ID, Gene_Symbol, Protein_Label, Uniprot_ID) %>%
  left_join(flag_tbl(stratified_cox[["early_OS"]],  suffix = "early_os"),  by = "nELISA_ID") %>%
  left_join(flag_tbl(stratified_cox[["early_PFS"]], suffix = "early_pfs"), by = "nELISA_ID") %>%
  left_join(flag_tbl(stratified_cox[["early_EFS"]], suffix = "early_efs"), by = "nELISA_ID") %>%
  left_join(flag_tbl(stratified_cox[["late_OS"]],   suffix = "late_os"),   by = "nELISA_ID") %>%
  left_join(flag_tbl(stratified_cox[["late_PFS"]],  suffix = "late_pfs"),  by = "nELISA_ID") %>%
  left_join(flag_tbl(stratified_cox[["late_EFS"]],  suffix = "late_efs"),  by = "nELISA_ID") %>%
  left_join(flag_tbl(interaction_results[["OS"]],   suffix = "interact_os"),  by = "nELISA_ID") %>%
  left_join(flag_tbl(interaction_results[["PFS"]],  suffix = "interact_pfs"), by = "nELISA_ID") %>%
  left_join(flag_tbl(interaction_results[["EFS"]],  suffix = "interact_efs"), by = "nELISA_ID")

flag_cols_early <- grep("^flag_early_", names(candidate_tbl), value = TRUE)
flag_cols_late  <- grep("^flag_late_",  names(candidate_tbl), value = TRUE)
flag_cols_int   <- grep("^flag_interact_", names(candidate_tbl), value = TRUE)

candidate_tbl <- candidate_tbl %>%
  mutate(
    n_early_flags = rowSums(across(all_of(flag_cols_early), ~ coalesce(., FALSE))),
    n_late_flags  = rowSums(across(all_of(flag_cols_late),  ~ coalesce(., FALSE))),
    n_interaction_flags = rowSums(across(all_of(flag_cols_int), ~ coalesce(., FALSE))),
    is_headline_six = Gene_Symbol %in% HEADLINE_GENES
  )
safe_write_csv(candidate_tbl %>% arrange(desc(n_early_flags), desc(n_interaction_flags)),
              file.path(TBL_DIR, "stage_candidate_biomarkers_full.csv"))

# Early-specific prognostic: significant in >=1 early-stage endpoint, NOT significant
# in the SAME endpoint(s) in late stage (i.e. the signal doesn't just carry over)
early_specific <- candidate_tbl %>%
  filter(n_early_flags >= 1,
         (flag_early_os  & !coalesce(flag_late_os,  FALSE)) |
         (flag_early_pfs & !coalesce(flag_late_pfs, FALSE)) |
         (flag_early_efs & !coalesce(flag_late_efs, FALSE))) %>%
  arrange(desc(n_early_flags), p_early_os)
safe_write_csv(early_specific, file.path(TBL_DIR, "stage_early_specific_prognostic_candidates.csv"))
message("Early-specific prognostic candidates (significant early, not the same endpoint in late): ", nrow(early_specific))
print(early_specific %>% select(Gene_Symbol, n_early_flags, p_early_os, p_early_pfs, p_early_efs,
                                p_late_os, p_late_pfs, p_late_efs) %>% slice_head(n = min(20, nrow(.))))

# Predictive: significant protein x stage interaction on any endpoint
predictive_candidates <- candidate_tbl %>%
  filter(n_interaction_flags >= 1) %>%
  arrange(desc(n_interaction_flags), p_interact_os)
safe_write_csv(predictive_candidates, file.path(TBL_DIR, "stage_predictive_interaction_candidates.csv"))
message("\nPredictive (stage-interacting) candidates: ", nrow(predictive_candidates))
print(predictive_candidates %>% select(Gene_Symbol, n_interaction_flags, p_interact_os, p_interact_pfs, p_interact_efs) %>%
        slice_head(n = min(20, nrow(.))))
