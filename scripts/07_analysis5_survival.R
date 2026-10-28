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
load(file.path(paths$rdata, "checkpoint_04.RData"))
message("Checkpoints 01, 03, 04 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
FDR_CUTOFF          <- 0.1
MIN_N               <- 35
MIN_EVENTS          <- 12
DELTA_OS_MIN_N      <- 10
DELTA_OS_MIN_EVENTS <- 5
KM_TOP_N            <- 4
ICENREG_BS_SAMPLES  <- 40
FOREST_N            <- 20

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

max_or_na <- function(x, ...) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  max(x, ...)
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

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op  <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out[tolower(raw_chr) %in% c("dead", "deceased", "event", "yes", "true")]  <- 1
  out[tolower(raw_chr) %in% c("alive", "censored", "no", "false")]          <- 0
  out[!is.finite(out) | !(out %in% c(0, 1))] <- NA_real_
  as.integer(out)
}

summarise_survival_endpoint <- function(sample_df, time_col, event_col) {
  op          <- options(warn = -1)
  time_value  <- as.numeric(sample_df[[time_col]])
  options(op)
  event_value <- coerce_event_indicator(sample_df[[event_col]])
  keep        <- !is.na(time_value) & !is.na(event_value)
  tibble(
    samples  = sum(keep),
    events   = sum(keep & event_value == 1, na.rm = TRUE),
    censored = sum(keep & event_value == 0, na.rm = TRUE)
  )
}

run_univariate_cox <- function(mat_samples_by_features, sample_df, time_col, event_col,
                               min_n = MIN_N, min_events = MIN_EVENTS) {
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

run_icenreg_sensitivity <- function(mat_samples_by_features, sample_df, candidate_ids,
                                    min_n = MIN_N, min_events = MIN_EVENTS,
                                    bs_samples = ICENREG_BS_SAMPLES) {
  candidate_ids <- intersect(candidate_ids, colnames(mat_samples_by_features))
  if (!length(candidate_ids)) return(tibble())
  op <- options(warn = -1)
  df <- sample_df %>%
    filter(sample_id %in% rownames(mat_samples_by_features)) %>%
    mutate(interval_lower = as.numeric(pfs_interval_lower),
           interval_upper = as.numeric(pfs_interval_upper)) %>%
    filter(!is.na(interval_lower), !is.na(interval_upper))
  options(op)
  if (nrow(df) < min_n) return(tibble())
  results <- lapply(candidate_ids, function(feature) {
    x    <- mat_samples_by_features[df$sample_id, feature]
    keep <- is.finite(x) & !is.na(df$interval_lower) & !is.na(df$interval_upper)
    if (sum(keep) < min_n) return(NULL)
    dat  <- tibble(lower = df$interval_lower[keep], upper = df$interval_upper[keep],
                   expr_scaled = as.numeric(scale(x[keep])))
    dat$upper[is.infinite(dat$upper)] <- Inf
    dat$upper[is.finite(dat$upper) & dat$upper <= dat$lower] <-
      dat$lower[is.finite(dat$upper) & dat$upper <= dat$lower] + 0.01
    n_events <- sum(is.finite(dat$upper), na.rm = TRUE)
    if (n_events < min_events) return(NULL)
    fit <- tryCatch(
      icenReg::ic_sp(survival::Surv(lower, upper, type = "interval2") ~ expr_scaled,
                     data = dat, model = "ph", bs_samples = bs_samples),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm   <- tryCatch(summary(fit)$summaryParameters, error = function(e) NULL)
    if (is.null(sm) || !nrow(as.data.frame(sm))) return(NULL)
    sm_df   <- as.data.frame(sm)
    est_col <- intersect(c("Estimate", "estimate", "coef"), names(sm_df))[1]
    se_col  <- intersect(c("Std.Error", "std.error", "se"),   names(sm_df))[1]
    p_col   <- intersect(c("p", "p.value", "Pr(>|z|)"),       names(sm_df))[1]
    z_col   <- intersect(c("z-value", "z", "z.value"),        names(sm_df))[1]
    est <- sm_df[[est_col]][1]; se <- if (!is.na(se_col)) sm_df[[se_col]][1] else NA_real_
    p_val <- if (!is.na(p_col)) sm_df[[p_col]][1] else NA_real_
    z_val <- if (!is.na(z_col)) sm_df[[z_col]][1] else NA_real_
    tibble(nELISA_ID = feature, beta = est, HR = exp(est),
           HR_low  = if (is.finite(se)) exp(est - 1.96*se) else NA_real_,
           HR_high = if (is.finite(se)) exp(est + 1.96*se) else NA_real_,
           z = z_val, p_value = p_val, n = nrow(dat), n_events = n_events,
           model = "icenReg interval-censored PH sensitivity")
  })
  bind_rows(results) %>% mutate(adj_p_value = p.adjust(p_value, method = "BH"))
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
          legend.text = element_text(size = 9), plot.margin = margin(68, 56, 46, 56))
}

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

make_survival_heatmap <- function(mat_samples_by_features, sample_df, result_df,
                                  time_col, event_col, title, path_stub, n = 20) {
  if (nrow(result_df) == 0) return(invisible(NULL))
  op  <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  df <- sample_df %>%
    mutate(time_value  = time_val,
           event_value = coerce_event_indicator(.data[[event_col]])) %>%
    filter(!is.na(time_value), !is.na(event_value),
           sample_id %in% rownames(mat_samples_by_features))
  if (nrow(df) < 10) return(invisible(NULL))
  top_features <- result_df %>%
    filter(nELISA_ID %in% colnames(mat_samples_by_features)) %>%
    arrange(p_value) %>% slice_head(n = min(n, nrow(.))) %>% pull(nELISA_ID)
  if (!length(top_features)) return(invisible(NULL))
  sample_order <- df %>% arrange(desc(event_value), time_value) %>% pull(sample_id)
  heat_mat <- t(mat_samples_by_features[sample_order, top_features, drop=FALSE])
  ann <- df %>%
    mutate(Event  = factor(event_value, levels=c(0,1), labels=c("Censored","Event")),
           Months = round(time_value, 1)) %>%
    select(sample_id, Event, Months) %>% tibble::column_to_rownames("sample_id")
  ann <- ann[sample_order, , drop=FALSE]
  row_labels <- result_df %>% filter(nELISA_ID %in% top_features) %>%
    distinct(nELISA_ID, Protein_Label) %>% tibble::deframe()
  rownames(heat_mat) <- unname(row_labels[rownames(heat_mat)])
  scaled_mat <- scale_rows_na(heat_mat)
  cairo_pdf(paste0(path_stub, ".pdf"), width=12, height=9.5)
  ph <- pheatmap::pheatmap(scaled_mat, annotation_col=ann, show_colnames=FALSE,
    color=colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
    na_col="grey90", main=title, silent=TRUE)
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  dev.off()
  png(paste0(path_stub, ".png"), width=12, height=9.5, units="in", res=300, bg="white")
  ph <- pheatmap::pheatmap(scaled_mat, annotation_col=ann, show_colnames=FALSE,
    color=colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
    na_col="grey90", main=title, silent=TRUE)
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  dev.off()
}

run_km_curves <- function(mat_samples_by_features, sample_df, result_df,
                           time_col, event_col, endpoint_label, figure_dir, table_dir,
                           top_n = KM_TOP_N) {
  if (nrow(result_df) == 0) return(tibble())
  op  <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  df <- sample_df %>%
    mutate(time_value  = time_val,
           event_value = coerce_event_indicator(.data[[event_col]])) %>%
    filter(!is.na(time_value), !is.na(event_value),
           sample_id %in% rownames(mat_samples_by_features))
  if (nrow(df) < 15 || sum(df$event_value == 1, na.rm=TRUE) < 6) return(tibble())
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
    km_df <- df[keep, , drop=FALSE] %>%
      mutate(expr_value = value_i[keep],
             expr_group = factor(ifelse(expr_value > median(expr_value, na.rm=TRUE), "High", "Low"),
                                 levels=c("Low","High")))
    if (any(table(km_df$expr_group) < 4) || length(unique(km_df$event_value)) < 2 ||
        sum(km_df$event_value==1, na.rm=TRUE) < 5) next
    fit <- tryCatch(
      survival::survfit(survival::Surv(time_value, event_value) ~ expr_group, data=km_df),
      error=function(e) NULL)
    diff_fit <- tryCatch(
      survival::survdiff(survival::Surv(time_value, event_value) ~ expr_group, data=km_df),
      error=function(e) NULL)
    if (is.null(fit)) next
    logrank_p <- if (!is.null(diff_fit))
      stats::pchisq(diff_fit$chisq, df=max(length(diff_fit$n)-1,1), lower.tail=FALSE) else NA_real_
    title_stub <- paste0("Analysis 5: ", endpoint_label, " KM — ", candidates$Protein_Label[i])
    km_plot <- survminer::ggsurvplot(
      fit, data=km_df, pval=TRUE, risk.table=FALSE, conf.int=FALSE,
      palette=c("#2166AC","#B2182B"), legend.title="Expression",
      legend.labs=c("Low","High"), ggtheme=ggplot2::theme_bw(base_size=11),
      title=title_stub,
      subtitle=paste0("Median split | HR=", round(candidates$HR[i],2),
                      " | Cox adj.P=", scales::number(candidates$adj_p_value[i], accuracy=0.001)),
      xlab="Months", ylab="Survival probability")
    km_risk <- survminer::ggsurvplot(
      fit, data=km_df, pval=FALSE, risk.table=TRUE, conf.int=FALSE,
      palette=c("#2166AC","#B2182B"), legend="none",
      ggtheme=ggplot2::theme_bw(base_size=11), risk.table.y.text=FALSE,
      xlab="Months", ylab="Survival probability")
    saved <- saved + 1L
    path_stub_km <- file.path(figure_dir,
      paste0("analysis5_km_", tolower(endpoint_label), "_",
             sprintf("%02d_", saved), sanitize_filename(candidates$Protein_Label[i])))
    cairo_pdf(paste0(path_stub_km, ".pdf"), width=7.8, height=5.5)
    print(km_plot$plot); dev.off()
    png(paste0(path_stub_km, ".png"), width=7.8, height=5.5, units="in", res=300, bg="white")
    print(km_plot$plot); dev.off()
    cairo_pdf(paste0(path_stub_km, "_risk.pdf"), width=7.8, height=2.6)
    print(km_risk$table); dev.off()
    png(paste0(path_stub_km, "_risk.png"), width=7.8, height=2.6, units="in", res=300, bg="white")
    print(km_risk$table); dev.off()
    km_rows[[saved]] <- tibble(
      endpoint=endpoint_label, rank=saved,
      nELISA_ID=feature_i, Protein_Label=candidates$Protein_Label[i],
      Gene_Symbol=candidates$Gene_Symbol[i],
      HR=candidates$HR[i], p_value=candidates$p_value[i],
      adj_p_value=candidates$adj_p_value[i], logrank_p_value=logrank_p,
      n=nrow(km_df), n_events=sum(km_df$event_value==1, na.rm=TRUE),
      n_low=sum(km_df$expr_group=="Low"), n_high=sum(km_df$expr_group=="High"))
  }
  km_tbl <- bind_rows(km_rows)
  if (nrow(km_tbl))
    safe_write_csv(km_tbl,
      file.path(table_dir, paste0("analysis5_", tolower(endpoint_label), "_km_summary.csv")))
  km_tbl
}

make_swimmer_plot <- function(clinical_df, path_stub, title, subset_label = NULL) {
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  GROUP_COL <- c(
    "Sx"="#E69F00","CTx"="#D55E00","CTx + IO"="#009E73","CTx + RTx"="#0072B2",
    "IO"="#56B4E9","TKI"="#CC79A7","Other"="#999999")
  swimmer_df <- clinical_df %>%
    mutate(
      followup_months   = pmax(dplyr::coalesce(os_months,0),
                               dplyr::coalesce(efs_months,0),
                               dplyr::coalesce(pfs_months_sensitivity,0), na.rm=TRUE),
      treatment_months  = pmax(dplyr::coalesce(treatment_duration_months,0), 0.15),
      progression_month = ifelse(pfs_event_sensitivity==1, pfs_months_sensitivity, NA_real_),
      death_month       = ifelse(os_event==1, os_months, NA_real_),
      best_response_month = ifelse(response_binary=="PR/CR",
                                   pmax(treatment_months*0.65,0.15), NA_real_),
      ongoing_month     = ifelse(is.na(death_month) & is.finite(followup_months),
                                 followup_months, NA_real_)
    ) %>%
    filter(is.finite(followup_months), followup_months > 0) %>%
    arrange(desc(followup_months), desc(!is.na(death_month)), desc(!is.na(progression_month))) %>%
    mutate(patient_plot = factor(patient_id, levels=rev(patient_id)))
  if (!nrow(swimmer_df)) return(invisible(NULL))
  subtitle_text <- paste0(
    subset_label %||% "Entire cohort", " | n=", nrow(swimmer_df),
    " | progression events=", sum(!is.na(swimmer_df$progression_month)),
    " | deaths=", sum(!is.na(swimmer_df$death_month)))
  p <- ggplot(swimmer_df, aes(y=patient_plot)) +
    geom_segment(aes(x=0, xend=followup_months, yend=patient_plot),
                 color="grey82", linewidth=1.3) +
    geom_segment(aes(x=0, xend=pmin(treatment_months,followup_months),
                     yend=patient_plot, color=new_group), linewidth=3.8, lineend="butt") +
    geom_point(data=swimmer_df %>% filter(!is.na(progression_month)),
               aes(x=progression_month, shape="Progression"),
               size=2.3, stroke=0.8, fill="#8dd3e7", color="#377eb8") +
    geom_point(data=swimmer_df %>% filter(!is.na(death_month)),
               aes(x=death_month, shape="Death"), size=2.6, color="black") +
    geom_point(data=swimmer_df %>% filter(!is.na(best_response_month)),
               aes(x=best_response_month, shape="Best response"),
               size=2.1, stroke=0.7, fill="white", color="grey25") +
    geom_point(data=swimmer_df %>% filter(!is.na(ongoing_month)),
               aes(x=ongoing_month, shape="Ongoing survival"),
               size=2.0, stroke=0.8, fill="white", color="grey30") +
    scale_color_manual(values=GROUP_COL, na.value="grey60", name="Treatment") +
    scale_shape_manual(
      values=c("Progression"=23,"Death"=8,"Best response"=21,"Ongoing survival"=23),
      name="Event") +
    labs(title=title, subtitle=subtitle_text, x="Months from treatment start", y="Patient") +
    theme_bw(base_size=11) +
    theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
          legend.position="top", legend.box="vertical",
          axis.text.y=element_text(size=5.8))
  ht <- min(max(10, 0.12*nrow(swimmer_df)+2.5), 24)
  cairo_pdf(paste0(path_stub,".pdf"), width=12.5, height=ht); print(p); dev.off()
  png(paste0(path_stub,".png"), width=12.5, height=ht, units="in", res=300, bg="white")
  print(p); dev.off()
}

make_upset_plot <- function(set_membership, path_stub, title) {
  set_membership <- set_membership[lengths(set_membership) > 0]
  if (!length(set_membership)) return(invisible(NULL))
  all_ids <- sort(unique(unlist(set_membership)))
  if (!length(all_ids)) return(invisible(NULL))
  upset_df <- tibble(item_id = all_ids)
  for (nm in names(set_membership)) upset_df[[nm]] <- all_ids %in% set_membership[[nm]]
  set_names <- names(set_membership)
  if (!nrow(upset_df) || !length(set_names)) return(invisible(NULL))
  intersection_tbl <- upset_df %>%
    count(across(all_of(set_names)), name="intersection_size") %>%
    filter(if_any(all_of(set_names), identity)) %>%
    arrange(desc(intersection_size)) %>%
    slice_head(n = min(15, nrow(.))) %>%
    mutate(intersection_id = paste0("I", seq_len(n())))
  sets_rev   <- rev(set_names)
  matrix_df  <- intersection_tbl %>%
    select(intersection_id, all_of(set_names)) %>%
    pivot_longer(cols=all_of(set_names), names_to="set_name", values_to="active") %>%
    mutate(set_name=factor(set_name, levels=sets_rev),
           set_num=match(as.character(set_name), sets_rev))
  connector_df <- matrix_df %>% filter(active) %>% group_by(intersection_id) %>%
    summarise(ymin=min(set_num), ymax=max(set_num), active_sets=n(), .groups="drop") %>%
    filter(active_sets > 1)
  bar_plot <- ggplot(intersection_tbl,
                     aes(x=factor(intersection_id, levels=intersection_id), y=intersection_size)) +
    geom_col(fill="grey30", width=0.72) +
    geom_text(aes(label=intersection_size), vjust=-0.35, size=3.6) +
    labs(title=title, x=NULL, y="Intersection size") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=14))
  matrix_plot <- ggplot(matrix_df,
                        aes(x=factor(intersection_id, levels=intersection_tbl$intersection_id),
                            y=set_num)) +
    geom_segment(data=connector_df,
                 aes(x=factor(intersection_id, levels=intersection_tbl$intersection_id),
                     xend=factor(intersection_id, levels=intersection_tbl$intersection_id),
                     y=ymin, yend=ymax),
                 inherit.aes=FALSE, linewidth=0.7, color="grey25") +
    geom_point(aes(fill=active), shape=21, size=3.2, color="black", stroke=0.45) +
    scale_fill_manual(values=c("TRUE"="black","FALSE"="white"), guide="none") +
    scale_y_continuous(breaks=seq_along(sets_rev), labels=sets_rev,
                       expand=expansion(add=c(0.4,0.4))) +
    labs(x=NULL, y=NULL) +
    theme_bw(base_size=11) +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_text(face="bold"))
  cairo_pdf(paste0(path_stub,".pdf"), width=11.5, height=8.4)
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow=2, ncol=1, heights=grid::unit(c(3.8,2.4),"null"))
  grid::pushViewport(grid::viewport(layout=lay))
  print(bar_plot,    vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(matrix_plot, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))
  grid::popViewport()
  dev.off()
  png(paste0(path_stub,".png"), width=11.5, height=8.4, units="in", res=300, bg="white")
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow=2, ncol=1, heights=grid::unit(c(3.8,2.4),"null"))
  grid::pushViewport(grid::viewport(layout=lay))
  print(bar_plot,    vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(matrix_plot, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))
  grid::popViewport()
  dev.off()
}

# ── Survival data frame ───────────────────────────────────────
message("\n=== 07 | Analysis 5 setup ===")
survival_df <- analysis1_df

surv_summary <- bind_rows(
  summarise_survival_endpoint(survival_df, "os_months",            "os_event")          %>% mutate(endpoint="OS",                       .before=1),
  summarise_survival_endpoint(survival_df, "pfs_months_main",      "pfs_event_main")    %>% mutate(endpoint="PFS_main",                  .before=1),
  summarise_survival_endpoint(survival_df, "pfs_months_sensitivity","pfs_event_sensitivity") %>% mutate(endpoint="PFS_sensitivity",       .before=1),
  summarise_survival_endpoint(survival_df, "efs_months",           "efs_event")         %>% mutate(endpoint="EFS",                       .before=1)
)
safe_write_csv(surv_summary,
  file.path(paths$tables, "06_analysis5_survival", "analysis5_survival_summary.csv"))
message("Survival endpoints:"); print(surv_summary)

# ── Univariate Cox screens ────────────────────────────────────
message("\n=== 07 | Univariate Cox ===")
os_res              <- add_feature_annotation(run_univariate_cox(expr_log2_rep, survival_df, "os_months",             "os_event"),            feature_meta)
pfs_res             <- add_feature_annotation(run_univariate_cox(expr_log2_rep, survival_df, "pfs_months_main",       "pfs_event_main"),       feature_meta)
pfs_res_sensitivity <- add_feature_annotation(run_univariate_cox(expr_log2_rep, survival_df, "pfs_months_sensitivity","pfs_event_sensitivity"),feature_meta)
efs_res             <- add_feature_annotation(run_univariate_cox(expr_log2_rep, survival_df, "efs_months",            "efs_event"),            feature_meta)

safe_write_csv(os_res,              file.path(paths$tables, "06_analysis5_survival", "analysis5_os_univariate_cox.csv"))
safe_write_csv(pfs_res,             file.path(paths$tables, "06_analysis5_survival", "analysis5_pfs_univariate_cox_main.csv"))
safe_write_csv(pfs_res_sensitivity, file.path(paths$tables, "06_analysis5_survival", "analysis5_pfs_univariate_cox_sensitivity.csv"))
safe_write_csv(efs_res,             file.path(paths$tables, "06_analysis5_survival", "analysis5_efs_univariate_cox.csv"))

message("OS:  ", nrow(os_res),  " proteins | FDR<0.1: ", sum(os_res$adj_p_value  < FDR_CUTOFF, na.rm=TRUE))
message("PFS: ", nrow(pfs_res), " proteins | FDR<0.1: ", sum(pfs_res$adj_p_value < FDR_CUTOFF, na.rm=TRUE))
message("EFS: ", nrow(efs_res), " proteins | FDR<0.1: ", sum(efs_res$adj_p_value < FDR_CUTOFF, na.rm=TRUE))

# ── Forest + volcano plots ────────────────────────────────────
message("\n=== 07 | Forest and volcano plots ===")
fig6 <- file.path(paths$figures, "06_analysis5_survival")
make_forest_plot(os_res,              "Analysis 5: OS Cox top proteins",           file.path(fig6,"analysis5_os_forest"))
make_forest_plot(pfs_res,             "Analysis 5: PFS Cox top proteins",          file.path(fig6,"analysis5_pfs_forest"))
make_forest_plot(pfs_res_sensitivity, "Analysis 5: PFS sensitivity Cox top proteins", file.path(fig6,"analysis5_pfs_sensitivity_forest"))
make_forest_plot(efs_res,             "Analysis 5: EFS Cox top proteins",          file.path(fig6,"analysis5_efs_forest"))
make_survival_volcano(os_res,              "Analysis 5: OS survival volcano",            file.path(fig6,"analysis5_os_volcano"))
make_survival_volcano(pfs_res,             "Analysis 5: PFS survival volcano",           file.path(fig6,"analysis5_pfs_volcano"))
make_survival_volcano(pfs_res_sensitivity, "Analysis 5: PFS sensitivity volcano",        file.path(fig6,"analysis5_pfs_sensitivity_volcano"))
make_survival_volcano(efs_res,             "Analysis 5: EFS survival volcano",           file.path(fig6,"analysis5_efs_volcano"))
message("Forest and volcano plots saved")

# ── Heatmaps ──────────────────────────────────────────────────
message("\n=== 07 | Survival heatmaps ===")
make_survival_heatmap(expr_log2_rep, survival_df, os_res,  "os_months", "os_event", "Analysis 5: OS top proteins", file.path(fig6,"analysis5_os_heatmap"))
make_survival_heatmap(expr_log2_rep, survival_df, pfs_res, "pfs_months_main","pfs_event_main", "Analysis 5: PFS top proteins", file.path(fig6,"analysis5_pfs_heatmap"))
make_survival_heatmap(expr_log2_rep, survival_df, efs_res, "efs_months","efs_event", "Analysis 5: EFS top proteins", file.path(fig6,"analysis5_efs_heatmap"))
message("Heatmaps saved")

# ── KM curves ─────────────────────────────────────────────────
message("\n=== 07 | KM curves ===")
tbl6 <- file.path(paths$tables, "06_analysis5_survival")
os_km  <- run_km_curves(expr_log2_rep, survival_df, os_res,  "os_months",       "os_event",       "OS",  fig6, tbl6)
pfs_km <- run_km_curves(expr_log2_rep, survival_df, pfs_res, "pfs_months_main", "pfs_event_main", "PFS", fig6, tbl6)
efs_km <- run_km_curves(expr_log2_rep, survival_df, efs_res, "efs_months",      "efs_event",      "EFS", fig6, tbl6)
message("KM curves saved — OS: ", nrow(os_km), " | PFS: ", nrow(pfs_km), " | EFS: ", nrow(efs_km))

# ── icenReg sensitivity (PFS interval-censored) ──────────────
message("\n=== 07 | icenReg PFS sensitivity ===")
pfs_icenreg_candidates <- pfs_res %>%
  filter(is.finite(p_value)) %>% arrange(p_value) %>%
  slice_head(n = min(20, nrow(.))) %>% pull(nELISA_ID)
pfs_icenreg <- add_feature_annotation(
  run_icenreg_sensitivity(expr_log2_rep, survival_df, pfs_icenreg_candidates),
  feature_meta)
safe_write_csv(pfs_icenreg,
  file.path(tbl6, "analysis5_pfs_icenreg_sensitivity_top20.csv"))
if (nrow(pfs_icenreg) > 0) {
  make_forest_plot(pfs_icenreg, "Analysis 5: PFS icenReg sensitivity",
                   file.path(fig6,"analysis5_pfs_icenreg_forest"))
  make_survival_volcano(pfs_icenreg, "Analysis 5: PFS icenReg sensitivity volcano",
                        file.path(fig6,"analysis5_pfs_icenreg_volcano"))
}
message("icenReg: ", nrow(pfs_icenreg), " proteins tested")

# ── Swimmer plots ─────────────────────────────────────────────
message("\n=== 07 | Swimmer plots ===")
make_swimmer_plot(clinical,    file.path(fig6,"analysis5_swimmer_entire_clinical"),
                  "Analysis 5: Swimmer plot — entire clinical cohort", "Whole clinical cohort")
make_swimmer_plot(survival_df, file.path(fig6,"analysis5_swimmer_matched_proteomics"),
                  "Analysis 5: Swimmer plot — matched proteomics cohort", "Matched baseline proteomics")
message("Swimmer plots saved")

# ── Delta (POST-PRE) Cox ──────────────────────────────────────
message("\n=== 07 | Delta Cox ===")
delta_os  <- tibble()
delta_pfs <- tibble()
if (length(paired_ids) >= DELTA_OS_MIN_N) {
  pre_pair  <- paired_pre  %>% arrange(patient_id)
  post_pair <- paired_post %>% arrange(patient_id)
  pair_intersect <- intersect(pre_pair$patient_id, post_pair$patient_id)
  pre_pair  <- pre_pair  %>% filter(patient_id %in% pair_intersect) %>% arrange(patient_id)
  post_pair <- post_pair %>% filter(patient_id %in% pair_intersect) %>% arrange(patient_id)
  delta_mat <- expr_log2_rep[post_pair$sample_id, , drop=FALSE] -
               expr_log2_rep[pre_pair$sample_id,  , drop=FALSE]
  rownames(delta_mat) <- pre_pair$sample_id
  delta_df  <- pre_pair %>% select(sample_id, patient_id) %>%
    left_join(clinical, by="patient_id")
  delta_os  <- add_feature_annotation(
    run_univariate_cox(delta_mat, delta_df, "os_months", "os_event",
                       min_n=DELTA_OS_MIN_N, min_events=DELTA_OS_MIN_EVENTS), feature_meta)
  delta_pfs <- add_feature_annotation(
    run_univariate_cox(delta_mat, delta_df, "pfs_months_main", "pfs_event_main",
                       min_n=20, min_events=6), feature_meta)
  safe_write_csv(delta_os,  file.path(tbl6, "analysis5_delta_os_cox.csv"))
  safe_write_csv(delta_pfs, file.path(tbl6, "analysis5_delta_pfs_cox.csv"))
  if (nrow(delta_os) > 0) {
    make_forest_plot(delta_os, "Analysis 5: Delta OS Cox", file.path(fig6,"analysis5_delta_os_forest"))
    make_survival_volcano(delta_os, "Analysis 5: Delta OS volcano", file.path(fig6,"analysis5_delta_os_volcano"))
  }
  if (nrow(delta_pfs) > 0) {
    make_forest_plot(delta_pfs, "Analysis 5: Delta PFS Cox", file.path(fig6,"analysis5_delta_pfs_forest"))
    make_survival_volcano(delta_pfs, "Analysis 5: Delta PFS volcano", file.path(fig6,"analysis5_delta_pfs_volcano"))
  }
  message("Delta OS: ", nrow(delta_os), " proteins | Delta PFS: ", nrow(delta_pfs))
} else {
  message("Delta Cox skipped: fewer than ", DELTA_OS_MIN_N, " pairs")
}

# ── Survival consistency audit ────────────────────────────────
surv_validation <- tibble(
  endpoint          = c("OS","PFS_main","PFS_sensitivity","EFS"),
  cohort_samples    = surv_summary$samples,
  cohort_events     = surv_summary$events,
  cox_proteins_tested = c(nrow(os_res),nrow(pfs_res),nrow(pfs_res_sensitivity),nrow(efs_res)),
  max_feature_events  = c(max_or_na(os_res$n_events),max_or_na(pfs_res$n_events),
                           max_or_na(pfs_res_sensitivity$n_events),max_or_na(efs_res$n_events))
)
safe_write_csv(surv_validation,
  file.path(paths$tables,"00_cohort_audit","survival_consistency_audit.csv"))

# ── Cohort upset plot ─────────────────────────────────────────
message("\n=== 07 | Cohort upset plot ===")
meaningful_sets <- list(
  IO_baseline    = analysis1_df %>% filter(treatment_class %in% c("Immunotherapy","Chemo + IO","Chemo + IO + RT","RT + IO (sequential)","IO + RT")) %>% pull(patient_id) %>% unique(),
  IO_paired      = paired_clinical %>% filter(treatment_class %in% c("Immunotherapy","Chemo + IO","Chemo + IO + RT","RT + IO (sequential)","IO + RT")) %>% pull(patient_id) %>% unique(),
  Surgery_baseline = analysis1_df %>% filter(treatment_class == "Surgery") %>% pull(patient_id) %>% unique(),
  Surgery_paired   = paired_clinical %>% filter(treatment_class == "Surgery") %>% pull(patient_id) %>% unique()
)
make_upset_plot(meaningful_sets,
  file.path(paths$figures,"00_overview","02_meaningful_cohort_upset"),
  "Meaningful cohort overlap")
safe_write_csv(
  tibble(patient_id=sort(unique(unlist(meaningful_sets))),
    IO_baseline      = as.integer(sort(unique(unlist(meaningful_sets))) %in% meaningful_sets$IO_baseline),
    IO_paired        = as.integer(sort(unique(unlist(meaningful_sets))) %in% meaningful_sets$IO_paired),
    Surgery_baseline = as.integer(sort(unique(unlist(meaningful_sets))) %in% meaningful_sets$Surgery_baseline),
    Surgery_paired   = as.integer(sort(unique(unlist(meaningful_sets))) %in% meaningful_sets$Surgery_paired)),
  file.path(paths$tables,"00_cohort_audit","meaningful_cohort_upset_membership.csv"))
message("Upset plot saved")

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 07 | Saving checkpoint_07.RData ===")
save(
  survival_df,
  os_res, pfs_res, pfs_res_sensitivity, efs_res,
  os_km,  pfs_km,  efs_km,
  pfs_icenreg,
  delta_os, delta_pfs,
  surv_summary,
  file = file.path(paths$rdata, "checkpoint_07.RData")
)
message("07_analysis5_survival.R complete")
