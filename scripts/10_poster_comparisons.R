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
FDR_CUTOFF    <- 0.1
LFC_THR       <- 0.5
MIN_PER_GROUP <- 4    # poster subgroups are small

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
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(color = "black", size = 14, face = "bold"),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9), plot.margin = ggplot2::margin(68, 56, 46, 56))
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
  secondary_lbl <- paste0("Secondary (FDR<=", formatC(padj_thr, format="f", digits=2), ")")
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
  if (nrow(ann) == 0)
    ann <- plot_df %>% arrange(adj_p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  x_rng <- range(plot_df$logFC, na.rm=TRUE)
  y_thr <- -log10(padj_thr)
  y_max <- max(plot_df$neg_log10, y_thr, na.rm=TRUE)
  x_pad <- max(0.35, diff(x_rng)*0.22); y_pad <- max(0.18, y_max*0.16)
  x_pl  <- x_rng[1]-0.45*x_pad; x_pr <- x_rng[2]+0.45*x_pad
  y_pt  <- y_max+0.12*y_pad;     y_lbl <- y_pt+0.9*y_pad
  p <- ggplot(plot_df, aes(x=logFC, y=neg_log10)) +
    geom_point(aes(color=status, shape=sig_level), size=2.4, alpha=0.82, stroke=0.7) +
    geom_vline(xintercept=c(-lfc_thr,0,lfc_thr), linetype=c("dashed","solid","dashed"),
               linewidth=c(0.6,0.4,0.6), color="black") +
    geom_hline(yintercept=y_thr, linetype="dashed", linewidth=0.6, color="black") +
    annotate("text", x=x_pl-0.18*x_pad, y=y_lbl, label=paste0("Higher in\n",group_a_label),
             hjust=0, vjust=0, size=4.8, fontface="bold", color="#2166AC") +
    annotate("text", x=x_pr+0.18*x_pad, y=y_lbl, label=paste0("Higher in\n",group_b_label),
             hjust=1, vjust=0, size=4.8, fontface="bold", color="#B2182B") +
    annotate("text", x=x_rng[1]+0.02*diff(x_rng), y=y_thr,
             label=paste0(" P Adj. = ",formatC(padj_thr,format="f",digits=2)),
             hjust=0, vjust=-0.6, size=3.8, fontface="italic") +
    scale_color_manual(values=c(setNames("#2166AC",left_label),
                                setNames("#B2182B",right_label),
                                "Not significant"="grey85"), drop=FALSE) +
    scale_shape_manual(values=c("Not significant"=1,"Primary (FDR<=0.05)"=19,
                                setNames(17,secondary_lbl)),
                       breaks=c("Not significant","Primary (FDR<=0.05)",secondary_lbl), drop=FALSE) +
    ggrepel::geom_text_repel(data=ann, aes(label=.data[[label_col]]),
                             size=3.1, fontface="bold", max.overlaps=30,
                             box.padding=0.5, point.padding=0.2,
                             segment.color="grey60", segment.size=0.2, show.legend=FALSE) +
    labs(title=title, x=x_label, y="P value (-log10)", color=NULL, shape="Significance") +
    guides(color="none", shape=guide_legend(override.aes=list(color="black"))) +
    coord_cartesian(xlim=c(x_pl,x_pr), ylim=c(0,y_pt), clip="off") +
    akila_volcano_theme()
  tryCatch({
    cairo_pdf(paste0(path_stub,".pdf"), width=8.6, height=8.2); print(p); dev.off()
    png(paste0(path_stub,".png"), width=8.6, height=8.2, units="in", res=300, bg="white")
    print(p); dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("    WARN: volcano save failed — ", conditionMessage(e),
            "\n    -> Close any open PDF viewers and re-run.")
  })
}

# ── Poster specifications ─────────────────────────────────────
message("\n=== 10 | Poster comparisons setup ===")
IO_CLASSES   <- c("Immunotherapy","Chemo + IO","Chemo + IO + RT",
                  "RT + IO (sequential)","IO + RT")
CHEMO_CLASSES <- c("Chemotherapy","Chemo + RT","RT only")

poster_specs <- list(
  poster_01_io_and_ctxio_response = list(
    df           = analysis1_df %>%
      filter(treatment_class %in% IO_CLASSES) %>%
      mutate(poster_outcome = response_binary),
    group_levels = c("SD/PD","PR/CR"),
    adjust       = c("plate_barcode","new_group"),
    title        = "Poster 1: IO-containing responders vs non-responders"
  ),
  poster_02_chemo_and_other_response = list(
    df           = analysis1_df %>%
      filter(treatment_class %in% CHEMO_CLASSES) %>%
      mutate(poster_outcome = response_binary),
    group_levels = c("SD/PD","PR/CR"),
    adjust       = c("plate_barcode","new_group"),
    title        = "Poster 2: Chemo/RT responders vs non-responders"
  ),
  poster_03_surgery_recurrence = list(
    df           = analysis1_df %>%
      filter(new_group == "Sx") %>%
      mutate(poster_outcome = recurrence_binary_sensitivity),
    group_levels = c("No recurrence","Recurrence"),
    adjust       = c("plate_barcode"),
    title        = "Poster 3: Surgery — recurrence vs non-recurrence (interval-expanded)"
  )
)

# ── Group counts ──────────────────────────────────────────────
poster_counts <- bind_rows(lapply(names(poster_specs), function(nm) {
  poster_specs[[nm]]$df %>%
    count(poster_outcome, name = "n") %>%
    mutate(comparison = nm, .before = 1)
}))
safe_write_csv(poster_counts,
  file.path(paths$tables, "09_poster_comparisons", "poster_group_counts.csv"))
message("Poster group counts:"); print(poster_counts)

# Count reconciliation table (vs poster-stated n)
poster_count_recon <- tibble(
  comparison = names(poster_specs),
  raw_clinical_n = c(
    sum(clinical$treatment_class %in% IO_CLASSES,    na.rm=TRUE),
    sum(clinical$treatment_class %in% CHEMO_CLASSES, na.rm=TRUE),
    sum(clinical$new_group == "Sx", na.rm=TRUE)
  ),
  matched_baseline_n = c(
    sum(analysis1_df$treatment_class %in% IO_CLASSES,    na.rm=TRUE),
    sum(analysis1_df$treatment_class %in% CHEMO_CLASSES, na.rm=TRUE),
    sum(analysis1_df$new_group == "Sx", na.rm=TRUE)
  ),
  note = c(
    "IO-containing classes only",
    "Chemo and RT classes only",
    "Recurrence = dated recurrence + interval-censored progression"
  )
)
safe_write_csv(poster_count_recon,
  file.path(paths$tables, "00_cohort_audit", "poster_count_reconciliation.csv"))

# ── Limma + volcano per poster ────────────────────────────────
message("\n=== 10 | Poster limma and volcano plots ===")
poster_results <- list()

for (nm in names(poster_specs)) {
  spec <- poster_specs[[nm]]
  message("\n--- ", nm, " ---")
  n_pos <- sum(spec$df$poster_outcome == spec$group_levels[2], na.rm=TRUE)
  n_neg <- sum(spec$df$poster_outcome == spec$group_levels[1], na.rm=TRUE)
  message("  ", spec$group_levels[1], ": n=", n_neg, " | ", spec$group_levels[2], ": n=", n_pos)

  res <- run_binary_limma(
    expr_log2_rep, spec$df, "poster_outcome",
    group_levels  = spec$group_levels,
    adjust_cols   = spec$adjust,
    min_per_group = MIN_PER_GROUP
  )
  if (is.null(res) || nrow(res) == 0) {
    message("  Skipped: insufficient group sizes")
    next
  }
  res <- add_feature_annotation(res, feature_meta)
  poster_results[[nm]] <- res
  safe_write_csv(res,
    file.path(paths$tables, "09_poster_comparisons", paste0(nm, ".csv")))
  make_volcano(res, spec$title,
    file.path(paths$figures, "09_poster_comparisons", nm))
  message("  Proteins: ", nrow(res),
          " | FDR<0.1: ", sum(res$adj_p_value < 0.1, na.rm=TRUE),
          " | nominal p<0.05: ", sum(res$p_value < 0.05, na.rm=TRUE))
}
