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

# ── Load checkpoints ───────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
MIN_GROUP_N <- 8
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")
OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")

# ── Local functions ────────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# Wilcoxon rank-sum test + expression data, for one protein
build_response_data <- function(sample_df, expr_vec, protein_label, min_n = MIN_GROUP_N) {
  df <- sample_df %>%
    mutate(expr_value = expr_vec[match(sample_id, names(expr_vec))]) %>%
    filter(!is.na(response_binary), is.finite(expr_value))
  if (nrow(df) < 2 * min_n) {
    message("    SKIP (n<", min_n, " per group): ", protein_label); return(NULL)
  }
  n_sdpd <- sum(df$response_binary == "SD/PD"); n_prcr <- sum(df$response_binary == "PR/CR")
  if (n_sdpd < min_n || n_prcr < min_n) {
    message("    SKIP (small group): ", protein_label, " - SD/PD=", n_sdpd, " PR/CR=", n_prcr)
    return(NULL)
  }
  op <- options(warn = -1)
  wilcox_p <- tryCatch(
    stats::wilcox.test(expr_value ~ response_binary, data = df, exact = FALSE)$p.value,
    error = function(e) NA_real_)
  options(op)
  df %>%
    transmute(
      patient_id, expr_value,
      response_group = factor(response_binary, levels = c("SD/PD", "PR/CR")),
      protein = protein_label,
      facet_label = paste0(protein_label, "  (Wilcoxon p=", scales::number(wilcox_p, accuracy = 0.001), ")"),
      wilcox_p = wilcox_p
    )
}

# Faceted combined boxplot: 2 columns x 5 rows
plot_response_facet_grid <- function(data_list, title, path_stub, ncol = 2) {
  combined <- bind_rows(data_list)
  if (!nrow(combined)) return(invisible(NULL))
  facet_order <- combined %>% distinct(protein, wilcox_p) %>% arrange(wilcox_p) %>% pull(protein)
  combined <- combined %>%
    mutate(protein = factor(protein, levels = facet_order)) %>%
    arrange(protein)
  facet_lbl_map <- combined %>% distinct(protein, facet_label)
  combined <- combined %>%
    mutate(facet_label = factor(facet_label, levels = facet_lbl_map$facet_label[match(facet_order, facet_lbl_map$protein)]))

  p <- ggplot(combined, aes(response_group, expr_value, fill = response_group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.6, aes(colour = response_group), show.legend = FALSE) +
    scale_fill_manual(values  = c("SD/PD" = "#90A4AE", "PR/CR" = "#2E7D32")) +
    scale_colour_manual(values = c("SD/PD" = "#546E7A", "PR/CR" = "#1B5E20")) +
    facet_wrap(~facet_label, ncol = ncol, scales = "free_y") +
    labs(title = title, x = NULL, y = "log2 concentration", fill = "Response") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", strip.text = element_text(face = "bold", size = 9))

  n_rows <- ceiling(length(unique(combined$protein)) / ncol)
  ht <- max(9, n_rows * 2.6)
  cairo_pdf(paste0(path_stub, ".pdf"), width = 9, height = ht); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 9, height = ht, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# ── Consensus protein list ────────────────────────────────────
message("\n=== 08g | Consensus protein lookup ===")
consensus_genes <- c("SHH", "CEACAM5", "GDF11", "POSTN", "MMP7",
                     "AHSG", "SNCG", "NOTUM", "FLT4", "TNFRSF10A")

protein_lookup <- feature_meta %>%
  filter(Gene_Symbol %in% consensus_genes) %>%
  distinct(Gene_Symbol, nELISA_ID, Protein_Label) %>%
  group_by(Gene_Symbol) %>% dplyr::slice(1) %>% ungroup()
message("  Found ", nrow(protein_lookup), " of ", length(consensus_genes), " proteins")

# ── Build response comparison data ────────────────────────────
message("\n=== 08g | Building response comparison data ===")
response_data_list <- list()

for (i in seq_len(nrow(protein_lookup))) {
  id  <- protein_lookup$nELISA_ID[i]
  lbl <- protein_lookup$Gene_Symbol[i]
  if (!id %in% colnames(expr_log2_rep)) {
    message("  ", lbl, " - not in expression matrix, skipping"); next
  }
  message("  ", lbl)
  expr_vec <- expr_log2_rep[, id]
  res <- build_response_data(analysis1_df, expr_vec, lbl)
  if (!is.null(res)) response_data_list[[lbl]] <- res
}

# ── Faceted combined plot ──────────────────────────────────────
message("\n=== 08g | Faceted response boxplot grid ===")
plot_response_facet_grid(
  response_data_list,
  "Consensus biomarkers \u2014 baseline expression by response (PR/CR vs SD/PD)",
  file.path(OUT_FIG, "consensus_response_facet_grid"), ncol = 2
)

# ── Summary table ──────────────────────────────────────────────
response_summary <- bind_rows(lapply(response_data_list, function(d) {
  d %>% distinct(protein, wilcox_p) %>%
    mutate(
      n_sdpd     = sum(d$response_group == "SD/PD"),
      n_prcr     = sum(d$response_group == "PR/CR"),
      median_sdpd = median(d$expr_value[d$response_group == "SD/PD"]),
      median_prcr = median(d$expr_value[d$response_group == "PR/CR"])
    )
})) %>% arrange(wilcox_p)
safe_write_csv(response_summary, file.path(OUT_TBL, "consensus_response_summary.csv"))
message("\nConsensus response summary:")
print(response_summary)

message("\n08g_consensus_response_boxplot.R complete")
