sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 03, 07 loaded")

while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
ALTML_TBL_DIR <- file.path(paths$tables, "07_analysis6_machine_learning", "altml")
OUT_TBL       <- file.path(paths$tables,  "07_analysis6_machine_learning")
KM_FIG_DIR    <- file.path(paths$figures, "07_analysis6_machine_learning", "consensus_km")
RESP_FIG_DIR  <- file.path(paths$figures, "07_analysis6_machine_learning")
dir.create(KM_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_METHODS_PER_COMPARISON <- 2   # of 4: en_selected, splsda_selected, boruta_confirmed, stabl_selected
MIN_COMPARISONS            <- 2   # "recurring" = appears in at least this many of the 5 response contrasts
MIN_KM_N      <- 10
MIN_KM_EVENTS <- 5
MIN_GROUP_N   <- 8

RESPONSE_FILES <- c(
  "altml_01_CR_vs_PD_cross_method_overlap.csv",
  "altml_02_CRPR_vs_PD_cross_method_overlap.csv",
  "altml_03_CRPR_vs_SDPD_cross_method_overlap.csv",
  "altml_04_CRPRSD_vs_PD_cross_method_overlap.csv",
  "altml_05_PR_vs_PD_cross_method_overlap.csv"
)

# ── Local functions ────────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out
}

# ── Step 1: rebuild the true 4-method consensus list ──────────
message("\n=== 08o | Rebuilding consensus gene list from cross_method_overlap tables ===")

per_comparison <- purrr::map_dfr(RESPONSE_FILES, function(f) {
  path <- file.path(ALTML_TBL_DIR, f)
  if (!file.exists(path)) {
    warning("Missing cross_method_overlap file: ", f)
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::mutate(
      comparison    = sub("_cross_method_overlap\\.csv$", "", f),
      n_methods_true = as.integer(en_selected) + as.integer(splsda_selected) +
                        as.integer(boruta_confirmed) + as.integer(stabl_selected)
    ) %>%
    dplyr::filter(n_methods_true >= MIN_METHODS_PER_COMPARISON) %>%
    dplyr::select(comparison, nELISA_ID, Protein_Label, n_methods_true)
})

consensus_summary <- per_comparison %>%
  dplyr::group_by(nELISA_ID, Protein_Label) %>%
  dplyr::summarise(
    n_comparisons  = dplyr::n(),
    comparisons    = paste(comparison, collapse = "; "),
    max_n_methods  = max(n_methods_true),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_comparisons >= MIN_COMPARISONS) %>%
  dplyr::arrange(dplyr::desc(n_comparisons), dplyr::desc(max_n_methods))

safe_write_csv(consensus_summary, file.path(OUT_TBL, "consensus_gene_list_rebuilt.csv"))

message("Rebuilt consensus list (", nrow(consensus_summary), " proteins, all confirmed present in current panel):")
print(consensus_summary)

consensus_ids  <- consensus_summary$nELISA_ID
consensus_syms <- consensus_summary$Protein_Label

if (length(consensus_ids) == 0) stop("No proteins met the consensus criteria — check thresholds.")

# ── Step 2: OS / PFS KM facet grid (rebuilt) ───────────────────
message("\n=== 08o | Building OS/PFS KM facet grid for rebuilt consensus list ===")

build_km_curve <- function(sample_df, expr_vec, protein_label, time_col, event_col) {
  df <- sample_df %>%
    dplyr::mutate(
      expr  = expr_vec[sample_id],
      time  = as.numeric(.data[[time_col]]),
      event = coerce_event_indicator(.data[[event_col]])
    ) %>%
    dplyr::filter(!is.na(expr), !is.na(time), !is.na(event))
  if (nrow(df) < MIN_KM_N || sum(df$event) < MIN_KM_EVENTS) return(NULL)

  med <- median(df$expr, na.rm = TRUE)
  df$group <- ifelse(df$expr >= med, "High", "Low")
  fit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)
  sd  <- survival::survdiff(survival::Surv(time, event) ~ group, data = df)
  p   <- pchisq(sd$chisq, df = 1, lower.tail = FALSE)

  broom::tidy(fit) %>%
    dplyr::mutate(
      group   = gsub("^group=", "", strata),
      protein = paste0(protein_label, "  (p=", signif(p, 3), ")")
    )
}

os_curves  <- list()
pfs_curves <- list()
for (i in seq_along(consensus_ids)) {
  id  <- consensus_ids[i]; lbl <- consensus_syms[i]
  if (!id %in% colnames(expr_log2_rep)) next
  expr_vec <- expr_log2_rep[, id]; names(expr_vec) <- rownames(expr_log2_rep)
  message("  ", lbl)
  c_os  <- build_km_curve(survival_df, expr_vec, lbl, "os_months", "os_event")
  c_pfs <- build_km_curve(survival_df, expr_vec, lbl, "pfs_months_main", "pfs_event_main")
  if (!is.null(c_os))  os_curves[[lbl]]  <- c_os
  if (!is.null(c_pfs)) pfs_curves[[lbl]] <- c_pfs
}

plot_km_grid <- function(curve_list, title, stub) {
  if (length(curve_list) == 0) return(invisible(NULL))
  df <- dplyr::bind_rows(curve_list)
  p <- ggplot(df, aes(x = time, y = estimate, colour = group)) +
    geom_step(linewidth = 0.7) +
    scale_colour_manual(values = c("Low" = "#2f6fb0", "High" = "#c9772b"), name = "Expression") +
    facet_wrap(~ protein, ncol = 3) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title, x = "Months", y = "Survival probability") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.text = element_text(face = "bold", size = 9))
  cairo_pdf(paste0(stub, ".pdf"), width = 11, height = 8); print(p); dev.off()
  png(paste0(stub, ".png"), width = 11, height = 8, units = "in", res = 300, bg = "white"); print(p); dev.off()
}

plot_km_grid(os_curves,  "Consensus biomarkers (rebuilt) \u2014 OS",
             file.path(KM_FIG_DIR, "consensus_OS_facet_grid_rebuilt"))
plot_km_grid(pfs_curves, "Consensus biomarkers (rebuilt) \u2014 PFS",
             file.path(KM_FIG_DIR, "consensus_PFS_facet_grid_rebuilt"))

message("OS/PFS KM facet grids saved (rebuilt consensus list)")

# ── Step 3: response boxplot facet grid (rebuilt) ──────────────
message("\n=== 08o | Building response boxplot facet grid for rebuilt consensus list ===")

response_data <- list()
for (i in seq_along(consensus_ids)) {
  id  <- consensus_ids[i]; lbl <- consensus_syms[i]
  if (!id %in% colnames(expr_log2_rep)) {
    message("  ", lbl, " — skipped (nELISA_ID not found in expr_log2_rep)")
    next
  }
  expr_vec <- expr_log2_rep[, id]; names(expr_vec) <- rownames(expr_log2_rep)

  df <- analysis1_df %>%
    dplyr::mutate(expr = expr_vec[sample_id]) %>%
    dplyr::filter(!is.na(response_binary), !is.na(expr))

  n_sdpd <- sum(df$response_binary == "SD/PD"); n_prcr <- sum(df$response_binary == "PR/CR")
  if (n_sdpd < MIN_GROUP_N || n_prcr < MIN_GROUP_N) {
    message("  ", lbl, " — skipped (n_sdpd=", n_sdpd, ", n_prcr=", n_prcr,
            "; below MIN_GROUP_N=", MIN_GROUP_N, ")")
    next
  }

  wt <- wilcox.test(expr ~ response_binary, data = df)
  message("  ", lbl, " — included (n_sdpd=", n_sdpd, ", n_prcr=", n_prcr,
          ", Wilcoxon p=", signif(wt$p.value, 3), ")")
  response_data[[lbl]] <- df %>%
    dplyr::mutate(protein = paste0(lbl, "  (Wilcoxon p=", signif(wt$p.value, 3), ")"))
}

message("Boxplot panels built: ", length(response_data), " of ", length(consensus_ids), " consensus proteins")

if (length(response_data) > 0) {
  resp_df <- dplyr::bind_rows(response_data)
  p_resp <- ggplot(resp_df, aes(x = response_binary, y = expr, fill = response_binary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.9) +
    scale_fill_manual(values = c("SD/PD" = "#a6bddb", "PR/CR" = "#74c476"), guide = "none") +
    facet_wrap(~ protein, ncol = 3, scales = "free_y") +
    labs(
      title = "Consensus biomarkers (rebuilt) \u2014 baseline expression by response",
      x = NULL, y = "log2 concentration"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 9))

  stub_resp <- file.path(RESP_FIG_DIR, "consensus_response_facet_grid_rebuilt")
  cairo_pdf(paste0(stub_resp, ".pdf"), width = 11, height = 8); print(p_resp); dev.off()
  png(paste0(stub_resp, ".png"), width = 11, height = 8, units = "in", res = 300, bg = "white"); print(p_resp); dev.off()
  message("Response boxplot facet grid saved (rebuilt consensus list)")
}

message("\n08o_consensus_gene_list_rebuild.R complete")
