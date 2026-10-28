if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}

sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())
if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
library(broom)

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
MIN_KM_N      <- 10
MIN_KM_EVENTS <- 5
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning", "consensus_km")
OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op  <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out[!is.finite(out) | !(out %in% c(0, 1))] <- NA_real_
  as.integer(out)
}

# ── build_km_curve_data / plot_km_facet_grid: copied verbatim from
#    08f_consensus_km.R so both scripts stay independently runnable ──
build_km_curve_data <- function(sample_df, expr_vec, time_col, event_col,
                                 protein_label, min_n = MIN_KM_N, min_events = MIN_KM_EVENTS) {
  op   <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  event_val <- coerce_event_indicator(sample_df[[event_col]])
  keep <- is.finite(expr_vec) & !is.na(time_val) & !is.na(event_val)
  if (sum(keep) < min_n) { message("    SKIP (n<", min_n, "): ", protein_label); return(NULL) }
  df <- tibble::tibble(
    time_value  = time_val[keep],
    event_value = event_val[keep],
    expr_value  = expr_vec[keep]
  ) %>%
    dplyr::mutate(expr_group = factor(
      ifelse(expr_value > median(expr_value, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")))
  if (any(table(df$expr_group) < 4) || sum(df$event_value == 1) < min_events) {
    message("    SKIP (small stratum/events): ", protein_label); return(NULL)
  }
  fit <- survival::survfit(survival::Surv(time_value, event_value) ~ expr_group, data = df)
  diff_fit <- tryCatch(
    survival::survdiff(survival::Surv(time_value, event_value) ~ expr_group, data = df),
    error = function(e) NULL)
  logrank_p <- if (!is.null(diff_fit))
    stats::pchisq(diff_fit$chisq, df = max(length(diff_fit$n) - 1, 1), lower.tail = FALSE) else NA_real_

  step_df <- broom::tidy(fit) %>%
    dplyr::mutate(
      expr_group = stringr::str_remove(strata, "expr_group="),
      protein    = protein_label,
      p_label    = ifelse(logrank_p < 0.001, "p<0.001",
                           paste0("p=", scales::number(logrank_p, accuracy = 0.001))),
      facet_label = paste0(protein_label, "  (", p_label, ")")
    )
  list(
    curve   = step_df,
    summary = tibble::tibble(protein = protein_label, n = nrow(df),
                             n_events = sum(df$event_value == 1),
                             n_low = sum(df$expr_group == "Low"),
                             n_high = sum(df$expr_group == "High"),
                             logrank_p = logrank_p)
  )
}

plot_km_facet_grid <- function(curve_list, title, path_stub, ncol = 2) {
  curve_df <- dplyr::bind_rows(lapply(curve_list, function(x) x$curve))
  if (!nrow(curve_df)) { message("  No curves to plot for: ", title); return(invisible(NULL)) }
  facet_order <- dplyr::bind_rows(lapply(curve_list, function(x) x$summary)) %>%
    dplyr::arrange(logrank_p) %>% dplyr::pull(protein)
  curve_df <- curve_df %>%
    dplyr::mutate(protein = factor(protein, levels = facet_order)) %>%
    dplyr::arrange(protein)
  facet_labels <- curve_df %>% dplyr::distinct(protein, facet_label) %>%
    dplyr::mutate(facet_label = factor(facet_label, levels = facet_label[match(facet_order, protein)]))
  curve_df <- curve_df %>%
    dplyr::left_join(facet_labels %>% dplyr::select(protein, facet_label_ord = facet_label), by = "protein")

  p <- ggplot(curve_df, aes(time, estimate, colour = expr_group)) +
    geom_step(linewidth = 0.8) +
    geom_point(data = curve_df %>% dplyr::filter(n.censor > 0), shape = 3, size = 1.4, show.legend = FALSE) +
    scale_colour_manual(values = c("Low" = "#2166AC", "High" = "#B2182B")) +
    facet_wrap(~facet_label_ord, ncol = ncol, scales = "free_x") +
    labs(title = title, x = "Months", y = "Survival probability", colour = "Expression") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", strip.text = element_text(face = "bold", size = 9))

  n_rows <- ceiling(length(unique(curve_df$protein)) / ncol)
  ht <- max(5.5, n_rows * 3.2)
  cairo_pdf(paste0(path_stub, ".pdf"), width = 9, height = ht); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 9, height = ht, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# ── Current OS consensus (GDF15, GZMA) — NOT GDF11 ─────────────
message("\n=== 08zd | Current OS consensus protein lookup ===")
consensus_genes <- c("GDF15", "GZMA")

protein_lookup <- feature_meta %>%
  dplyr::filter(Gene_Symbol %in% consensus_genes) %>%
  dplyr::distinct(Gene_Symbol, nELISA_ID, Protein_Label) %>%
  dplyr::group_by(Gene_Symbol) %>% dplyr::slice(1) %>% dplyr::ungroup()

missing_genes <- setdiff(consensus_genes, protein_lookup$Gene_Symbol)
if (length(missing_genes)) {
  stop("Not found in feature_meta: ", paste(missing_genes, collapse = ", "),
       " -- check Gene_Symbol spelling (panel-native names sometimes differ, e.g. 'GDF-15 (MIC-1)').")
}
message("Found ", nrow(protein_lookup), " of ", length(consensus_genes), " proteins")
print(protein_lookup)

# ── Build KM curve data: OS, whole cohort ──────────────────────
message("\n=== 08zd | Building OS KM curve data ===")
os_curves <- list()

for (i in seq_len(nrow(protein_lookup))) {
  id  <- protein_lookup$nELISA_ID[i]
  lbl <- protein_lookup$Gene_Symbol[i]
  if (!id %in% colnames(expr_log2_rep)) {
    message("  ", lbl, " - nELISA_ID not in expression matrix, skipping"); next
  }
  message("  ", lbl)
  res_os <- build_km_curve_data(survival_df, expr_log2_rep[survival_df$sample_id, id],
                                "os_months", "os_event", lbl)
  if (!is.null(res_os)) os_curves[[lbl]] <- res_os
}

# ── Plot ─────────────────────────────────────────────────────
message("\n=== 08zd | Plotting GDF15/GZMA OS KM grid ===")
plot_km_facet_grid(os_curves, "OS consensus biomarkers (Cooper-integrated: GDF15, GZMA) \u2014 median split, whole matched cohort",
                   file.path(OUT_FIG, "gdf15_gzma_OS_facet_grid"), ncol = 2)

# ── Summary table ────────────────────────────────────────────
km_summary <- dplyr::bind_rows(lapply(os_curves, function(x) x$summary)) %>%
  dplyr::arrange(logrank_p)
safe_write_csv(km_summary, file.path(OUT_TBL, "gdf15_gzma_os_km_summary.csv"))
message("\nGDF15/GZMA OS KM summary:")
print(km_summary)

message("\n=== 08zd complete ===")
message("Use gdf15_gzma_OS_facet_grid.png as the replacement Panel B in Sup_Fig_S1_OS ",
        "(the old consensus_OS_facet_grid.png is stale -- GDF11 not GDF15, predates Cooper integration).")
