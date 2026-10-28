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
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 03, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
MIN_KM_N        <- 20    # minimum total n across all 4 groups
MIN_GROUP_N     <- 4     # minimum n per individual group
MIN_KM_EVENTS   <- 6     # minimum total events across all 4 groups
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")
OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")

# ── Local functions ────────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op  <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out[!is.finite(out) | !(out %in% c(0, 1))] <- NA_real_
  as.integer(out)
}

# Build 4-group (response x expression) KM curve data for one
# protein. Expression High/Low is a median split computed across
# the protein's full evaluable cohort (same convention as 08f).
build_response_x_expr_km_data <- function(sample_df, expr_vec, time_col, event_col,
                                          protein_label,
                                          min_n = MIN_KM_N, min_group_n = MIN_GROUP_N,
                                          min_events = MIN_KM_EVENTS) {
  op   <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  event_val <- coerce_event_indicator(sample_df[[event_col]])
  keep <- is.finite(expr_vec) & !is.na(time_val) & !is.na(event_val) & !is.na(sample_df$response_binary)
  if (sum(keep) < min_n) { message("    SKIP (n<", min_n, "): ", protein_label); return(NULL) }

  expr_sub <- expr_vec[keep]
  expr_grp <- ifelse(expr_sub > median(expr_sub, na.rm = TRUE), "High", "Low")

  df <- tibble::tibble(
    time_value  = time_val[keep],
    event_value = event_val[keep],
    response    = sample_df$response_binary[keep],
    expr_group  = expr_grp
  ) %>%
    dplyr::mutate(
      group4 = factor(paste0(response, " - ", expr_group),
                      levels = c("SD/PD - Low", "SD/PD - High", "PR/CR - Low", "PR/CR - High"))
    )

  group_counts <- table(df$group4)
  if (any(group_counts < min_group_n) || sum(df$event_value == 1) < min_events) {
    message("    SKIP (small group/events): ", protein_label,
            " - counts: ", paste(names(group_counts), group_counts, collapse = ", "))
    return(NULL)
  }

  fit <- survival::survfit(survival::Surv(time_value, event_value) ~ group4, data = df)
  diff_fit <- tryCatch(
    survival::survdiff(survival::Surv(time_value, event_value) ~ group4, data = df),
    error = function(e) NULL)
  logrank_p <- if (!is.null(diff_fit))
    stats::pchisq(diff_fit$chisq, df = max(length(diff_fit$n) - 1, 1), lower.tail = FALSE) else NA_real_

  step_df <- broom::tidy(fit) %>%
    dplyr::mutate(
      group4      = stringr::str_remove(strata, "group4="),
      protein     = protein_label,
      facet_label = paste0(protein_label, "  (n=", nrow(df), ", global p=",
                           scales::number(logrank_p, accuracy = 0.001), ")")
    )
  list(
    curve   = step_df,
    summary = tibble::tibble(
      protein = protein_label, n = nrow(df), n_events = sum(df$event_value == 1),
      n_sdpd_low  = group_counts["SD/PD - Low"],  n_sdpd_high = group_counts["SD/PD - High"],
      n_prcr_low  = group_counts["PR/CR - Low"],  n_prcr_high = group_counts["PR/CR - High"],
      logrank_p_global = logrank_p
    )
  )
}

# Faceted combined KM plot: 2 columns x 5 rows, 4 curves per panel
plot_4group_km_facet_grid <- function(curve_list, title, path_stub, ncol = 2) {
  curve_df <- bind_rows(lapply(curve_list, function(x) x$curve))
  if (!nrow(curve_df)) return(invisible(NULL))
  facet_order <- bind_rows(lapply(curve_list, function(x) x$summary)) %>%
    arrange(logrank_p_global) %>% pull(protein)
  curve_df <- curve_df %>% mutate(protein = factor(protein, levels = facet_order)) %>% arrange(protein)
  facet_lbl_map <- curve_df %>% distinct(protein, facet_label)
  curve_df <- curve_df %>%
    mutate(facet_label = factor(facet_label,
                                levels = facet_lbl_map$facet_label[match(facet_order, facet_lbl_map$protein)]),
           group4 = factor(group4, levels = c("SD/PD - Low", "SD/PD - High", "PR/CR - Low", "PR/CR - High")))

  group_colours <- c(
    "SD/PD - Low"  = "#90A4AE",  # light grey-blue: worst expected
    "SD/PD - High" = "#B71C1C",  # dark red: non-responder, high expression
    "PR/CR - Low"  = "#1B5E20",  # dark green: best expected
    "PR/CR - High" = "#F9A825"   # amber: responder but high expression - the interesting group
  )

  p <- ggplot(curve_df, aes(time, estimate, colour = group4)) +
    geom_step(linewidth = 0.75) +
    geom_point(data = curve_df %>% filter(n.censor > 0), shape = 3, size = 1.2, show.legend = FALSE) +
    scale_colour_manual(values = group_colours) +
    facet_wrap(~facet_label, ncol = ncol, scales = "free_x") +
    labs(title = title, x = "Months", y = "Survival probability", colour = "Response x Expression") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", legend.text = element_text(size = 8),
          strip.text = element_text(face = "bold", size = 9))

  n_rows <- ceiling(length(unique(curve_df$protein)) / ncol)
  ht <- max(10, n_rows * 2.9)
  cairo_pdf(paste0(path_stub, ".pdf"), width = 9.5, height = ht); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 9.5, height = ht, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# ── Consensus protein list ────────────────────────────────────
message("\n=== 08i | Consensus protein lookup ===")
consensus_genes <- c("SHH", "CEACAM5", "GDF11", "POSTN", "MMP7",
                     "AHSG", "SNCG", "NOTUM", "FLT4", "TNFRSF10A")

protein_lookup <- feature_meta %>%
  filter(Gene_Symbol %in% consensus_genes) %>%
  distinct(Gene_Symbol, nELISA_ID, Protein_Label) %>%
  group_by(Gene_Symbol) %>% dplyr::slice(1) %>% ungroup()
message("  Found ", nrow(protein_lookup), " of ", length(consensus_genes), " proteins")

# ── Build 4-group KM data per protein ─────────────────────────
message("\n=== 08i | Building response x expression KM data ===")
pfs_curves <- list(); os_curves <- list()

for (i in seq_len(nrow(protein_lookup))) {
  id  <- protein_lookup$nELISA_ID[i]
  lbl <- protein_lookup$Gene_Symbol[i]
  if (!id %in% colnames(expr_log2_rep)) {
    message("  ", lbl, " - not in expression matrix, skipping"); next
  }
  message("  ", lbl)
  expr_vec <- expr_log2_rep[survival_df$sample_id, id]

  res_pfs <- build_response_x_expr_km_data(survival_df, expr_vec, "pfs_months_main", "pfs_event_main", lbl)
  if (!is.null(res_pfs)) pfs_curves[[lbl]] <- res_pfs

  res_os <- build_response_x_expr_km_data(survival_df, expr_vec, "os_months", "os_event", lbl)
  if (!is.null(res_os)) os_curves[[lbl]] <- res_os
}

# ── Faceted combined plots ────────────────────────────────────
message("\n=== 08i | Faceted 4-group KM grids ===")
plot_4group_km_facet_grid(
  pfs_curves,
  "Response x Expression (4 groups) \u2014 PFS, per protein evaluable cohort",
  file.path(OUT_FIG, "consensus_response_x_expression_KM_PFS_facet_grid"), ncol = 2)
plot_4group_km_facet_grid(
  os_curves,
  "Response x Expression (4 groups) \u2014 OS, per protein evaluable cohort",
  file.path(OUT_FIG, "consensus_response_x_expression_KM_OS_facet_grid"), ncol = 2)

# ── Summary table ──────────────────────────────────────────────
km_summary <- bind_rows(
  bind_rows(lapply(pfs_curves, function(x) x$summary)) %>% mutate(endpoint = "PFS"),
  bind_rows(lapply(os_curves,  function(x) x$summary)) %>% mutate(endpoint = "OS")
) %>% arrange(endpoint, logrank_p_global)
safe_write_csv(km_summary, file.path(OUT_TBL, "consensus_response_x_expression_km_summary.csv"))
message("\nResponse x Expression 4-group KM summary:")
print(km_summary)

message("\n08i_consensus_response_expression_km.R complete")
