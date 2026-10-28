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

# ── Local functions ────────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op  <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out[!is.finite(out) | !(out %in% c(0, 1))] <- NA_real_
  as.integer(out)
}

# Median-split helper: classify expression into High/Low groups


# Build KM step-curve data (no plotting yet) - returns a long
# data frame ready for ggplot facetting, plus the logrank p-value.
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
      facet_label = paste0(protein_label, "  (p=", scales::number(logrank_p, accuracy = 0.001), ")")
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

# Faceted combined KM plot: 2 columns x 5 rows (10 proteins)
plot_km_facet_grid <- function(curve_list, title, path_stub, ncol = 2) {
  curve_df <- bind_rows(lapply(curve_list, function(x) x$curve))
  if (!nrow(curve_df)) return(invisible(NULL))
  # Order facets by ascending logrank p (most significant first)
  facet_order <- bind_rows(lapply(curve_list, function(x) x$summary)) %>%
    arrange(logrank_p) %>% pull(protein)
  curve_df <- curve_df %>%
    mutate(protein = factor(protein, levels = facet_order)) %>%
    arrange(protein)
  facet_labels <- curve_df %>% distinct(protein, facet_label) %>%
    mutate(facet_label = factor(facet_label, levels = facet_label[match(facet_order, protein)]))
  curve_df <- curve_df %>%
    left_join(facet_labels %>% select(protein, facet_label_ord = facet_label), by = "protein")

  p <- ggplot(curve_df, aes(time, estimate, colour = expr_group)) +
    geom_step(linewidth = 0.8) +
    geom_point(data = curve_df %>% filter(n.censor > 0), shape = 3, size = 1.4, show.legend = FALSE) +
    scale_colour_manual(values = c("Low" = "#2166AC", "High" = "#B2182B")) +
    facet_wrap(~facet_label_ord, ncol = ncol, scales = "free_x") +
    labs(title = title, x = "Months", y = "Survival probability", colour = "Expression") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", strip.text = element_text(face = "bold", size = 9))

  n_rows <- ceiling(length(unique(curve_df$protein)) / ncol)
  ht <- max(9, n_rows * 2.6)
  cairo_pdf(paste0(path_stub, ".pdf"), width = 9, height = ht); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 9, height = ht, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# ── Consensus protein list ────────────────────────────────────
message("\n=== 08f | Consensus protein lookup ===")
consensus_genes <- c("SHH", "CEACAM5", "GDF11", "POSTN", "MMP7",
                     "AHSG", "SNCG", "NOTUM", "FLT4", "TNFRSF10A")

protein_lookup <- feature_meta %>%
  filter(Gene_Symbol %in% consensus_genes) %>%
  distinct(Gene_Symbol, nELISA_ID, Protein_Label) %>%
  group_by(Gene_Symbol) %>% dplyr::slice(1) %>% ungroup()

missing_genes <- setdiff(consensus_genes, protein_lookup$Gene_Symbol)
if (length(missing_genes)) message("  Not found in feature_meta: ", paste(missing_genes, collapse = ", "))
message("  Found ", nrow(protein_lookup), " of ", length(consensus_genes), " proteins")
print(protein_lookup)

# ── Build KM curve data: OS and PFS, whole cohort ─────────────
message("\n=== 08f | Building KM curve data ===")
pfs_curves <- list(); os_curves <- list()

for (i in seq_len(nrow(protein_lookup))) {
  id  <- protein_lookup$nELISA_ID[i]
  lbl <- protein_lookup$Gene_Symbol[i]
  if (!id %in% colnames(expr_log2_rep)) {
    message("  ", lbl, " - nELISA_ID not in expression matrix, skipping"); next
  }
  message("  ", lbl)
  res_pfs <- build_km_curve_data(survival_df, expr_log2_rep[survival_df$sample_id, id],
                                 "pfs_months_main", "pfs_event_main", lbl)
  if (!is.null(res_pfs)) pfs_curves[[lbl]] <- res_pfs
  res_os <- build_km_curve_data(survival_df, expr_log2_rep[survival_df$sample_id, id],
                                "os_months", "os_event", lbl)
  if (!is.null(res_os)) os_curves[[lbl]] <- res_os
}

# ── Faceted combined plots (2 cols x 5 rows) ──────────────────
message("\n=== 08f | Faceted combined KM plots ===")
plot_km_facet_grid(pfs_curves, "Consensus biomarkers \u2014 PFS (median split, whole matched cohort)",
                   file.path(OUT_FIG, "consensus_PFS_facet_grid"), ncol = 2)
plot_km_facet_grid(os_curves,  "Consensus biomarkers \u2014 OS (median split, whole matched cohort)",
                   file.path(OUT_FIG, "consensus_OS_facet_grid"), ncol = 2)

# ── Summary table ──────────────────────────────────────────────
km_summary <- bind_rows(
  bind_rows(lapply(pfs_curves, function(x) x$summary)) %>% mutate(endpoint = "PFS"),
  bind_rows(lapply(os_curves,  function(x) x$summary)) %>% mutate(endpoint = "OS")
) %>% arrange(endpoint, logrank_p)
safe_write_csv(km_summary, file.path(OUT_TBL, "consensus_km_summary.csv"))
message("\nConsensus KM summary:")
print(km_summary)

message("\n08f_consensus_km.R complete")
