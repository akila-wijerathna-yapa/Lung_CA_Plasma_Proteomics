if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_05.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
load(file.path(paths$rdata, "checkpoint_09.RData"))   # gene_sets, fgsea_log_tbl
message("Checkpoints 01, 03, 05, 07, 09 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
RNG_SEED      <- 42
MIN_GENE_SET  <- 5
MAX_GENE_SET  <- 500
FGSEA_FDR     <- 0.25
NES_SHIFT_TOL <- 0.15   # flag as "moderately dependent" if |NES| shrinks by >15%
if (!exists("ANALYSIS8_DIR")) ANALYSIS8_DIR <- "08_pathway_and_network_analysis"

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

clean_pathway_name <- function(x) {
  x %>%
    stringr::str_replace_all("^HALLMARK_", "") %>%
    stringr::str_replace_all("^REACTOME_", "") %>%
    stringr::str_replace_all("^GOBP_|^GO:BP_", "") %>%
    stringr::str_replace_all("^GOMF_|^GO:MF_", "") %>%
    stringr::str_replace_all("^GOCC_|^GO:CC_", "") %>%
    stringr::str_replace_all("^KEGG_(MEDICUS|LEGACY)_", "") %>%
    stringr::str_replace_all("^KEGG_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish()
}

prepare_rank_vector <- function(result_df, stat_col = "t_stat") {
  rank_df <- result_df %>%
    filter(!is.na(Gene_Symbol), !is.na(.data[[stat_col]]),
           is.finite(.data[[stat_col]])) %>%
    group_by(Gene_Symbol) %>%
    summarise(stat = .data[[stat_col]][which.max(abs(.data[[stat_col]]))], .groups = "drop") %>%
    filter(!is.na(stat), is.finite(stat))
  stats <- rank_df$stat
  names(stats) <- rank_df$Gene_Symbol
  sort(stats, decreasing = TRUE)
}

# ── Reconstruct enrichment_specs exactly as in 09_pathway_network.R ────
enrichment_specs <- list(
  analysis1_ctx_vs_rest    = list(df = analysis1_ovr_results[["analysis1_ctx_vs_rest"]],    stat = "t_stat"),
  analysis1_io_vs_rest     = list(df = analysis1_ovr_results[["analysis1_io_vs_rest"]],     stat = "t_stat"),
  analysis1_ctxio_vs_rest  = list(df = analysis1_ovr_results[["analysis1_ctxio_vs_rest"]],  stat = "t_stat"),
  analysis1_ctxrtx_vs_rest = list(df = analysis1_ovr_results[["analysis1_ctxrtx_vs_rest"]], stat = "t_stat"),
  analysis1_tki_vs_rest    = list(df = analysis1_ovr_results[["analysis1_tki_vs_rest"]],    stat = "t_stat"),
  analysis1_sx_vs_rest     = list(df = analysis1_ovr_results[["analysis1_sx_vs_rest"]],     stat = "t_stat"),
  analysis1_other_vs_rest  = list(df = analysis1_ovr_results[["analysis1_other_vs_rest"]],  stat = "t_stat"),
  analysis3_primary_response    = list(df = analysis3_primary,       stat = "t_stat"),
  analysis3_secondary_response  = list(df = analysis3_secondary,     stat = "t_stat"),
  analysis3_primary_expl        = list(df = analysis3_primary_expl,  stat = "t_stat"),
  analysis3_secondary_expl      = list(df = analysis3_secondary_expl, stat = "t_stat"),
  analysis5_os_survival   = list(df = os_res,  stat = "z"),
  analysis5_pfs_survival  = list(df = pfs_res, stat = "z"),
  analysis5_efs_survival  = list(df = efs_res, stat = "z"),
  analysis5_delta_os      = list(df = delta_os,  stat = "z"),
  analysis5_delta_pfs     = list(df = delta_pfs, stat = "z")
)
enrichment_specs <- enrichment_specs[!sapply(enrichment_specs, function(s) is.null(s$df) || nrow(s$df) == 0)]

# ── Load hit list from 09b ─────────────────────────────────────
base_tbl <- file.path(paths$tables, ANALYSIS8_DIR)
hits_path <- file.path(base_tbl, "target_protein_pathway_hits.csv")
if (!file.exists(hits_path)) stop("target_protein_pathway_hits.csv not found — run 09b_target_protein_enrichment_plots.R first.")
hit_ranked <- readr::read_csv(hits_path, show_col_types = FALSE)
message("Hits loaded for LOO check: ", nrow(hit_ranked))

# ── LOO fgsea rerun, one full-collection rerun per unique (label, gene, collection) ──
message("\n=== 09c | Running leave-one-out fgsea per hit ===")
stats_cache <- list()
loo_rows <- vector("list", nrow(hit_ranked))

for (i in seq_len(nrow(hit_ranked))) {
  gene_i  <- hit_ranked$Gene_Symbol[i]
  label_i <- hit_ranked$label[i]
  coll_i  <- hit_ranked$collection[i]
  pw_i    <- hit_ranked$pathway[i]
  coll_name <- toupper(coll_i)

  spec_i <- enrichment_specs[[label_i]]
  if (is.null(spec_i)) { message("  Skip: no result_df for ", label_i); next }

  cache_key <- label_i
  if (is.null(stats_cache[[cache_key]])) stats_cache[[cache_key]] <- prepare_rank_vector(spec_i$df, stat_col = spec_i$stat)
  stats_full <- stats_cache[[cache_key]]

  if (!(gene_i %in% names(stats_full))) { message("  Skip: ", gene_i, " not in ranked stats for ", label_i); next }
  stats_loo <- stats_full[names(stats_full) != gene_i]

  pathways_use <- lapply(gene_sets[[coll_name]], function(gs) intersect(gs, names(stats_loo)))
  pathways_use <- pathways_use[lengths(pathways_use) >= MIN_GENE_SET]
  if (!(pw_i %in% names(pathways_use))) { message("  Skip: ", pw_i, " below min size after removing ", gene_i); next }

  set.seed(RNG_SEED)
  fg_loo <- tryCatch(
    fgsea::fgsea(pathways = pathways_use, stats = stats_loo, minSize = MIN_GENE_SET, maxSize = MAX_GENE_SET, eps = 0),
    error = function(e) { message("  fgsea LOO error: ", conditionMessage(e)); NULL })
  if (is.null(fg_loo) || !nrow(fg_loo)) next

  row_loo <- fg_loo %>% filter(pathway == pw_i) %>% slice_head(n = 1)
  if (!nrow(row_loo)) next

  NES_orig  <- hit_ranked$NES[i]
  padj_orig <- hit_ranked$padj[i]
  NES_loo   <- row_loo$NES[1]
  padj_loo  <- row_loo$padj[1]
  nes_pct_change <- 100 * (abs(NES_loo) - abs(NES_orig)) / abs(NES_orig)

  verdict <- dplyr::case_when(
    padj_orig <= FGSEA_FDR & (is.na(padj_loo) | padj_loo > FGSEA_FDR) ~ "Collapses without protein (protein-driven)",
    sign(NES_loo) != sign(NES_orig) ~ "Sign flips without protein (protein-driven)",
    abs(nes_pct_change) > (NES_SHIFT_TOL * 100) ~ "Moderately protein-dependent",
    TRUE ~ "Robust to removal"
  )

  loo_rows[[i]] <- tibble(
    Gene_Symbol = gene_i, label = label_i, collection = coll_i, pathway = pw_i,
    pathway_clean = clean_pathway_name(pw_i),
    NES_original = NES_orig, padj_original = padj_orig,
    NES_loo = NES_loo, padj_loo = padj_loo,
    nes_pct_change = round(nes_pct_change, 1),
    verdict = verdict
  )
  message("  ", gene_i, " / ", pw_i, " -> ", verdict)
}

loo_tbl <- bind_rows(loo_rows)
safe_write_csv(loo_tbl, file.path(base_tbl, "target_protein_pathway_hits_LOO.csv"))
message("\nLOO results written: ", nrow(loo_tbl), " rows")
print(loo_tbl %>% select(Gene_Symbol, pathway_clean, NES_original, NES_loo, padj_original, padj_loo, verdict))

# ── Dumbbell plot: NES original vs LOO per hit ──────────────────
if (nrow(loo_tbl)) {
  verdict_colours <- c(
    "Robust to removal"                          = "#2a9d8f",
    "Moderately protein-dependent"                = "#f4a261",
    "Collapses without protein (protein-driven)"  = "#d7301f",
    "Sign flips without protein (protein-driven)" = "#d7301f"
  )
  plot_df <- loo_tbl %>%
    mutate(label_row = paste0(Gene_Symbol, " | ", stringr::str_trunc(pathway_clean, 38), " (", label, ")")) %>%
    arrange(verdict, Gene_Symbol) %>%
    mutate(label_row = factor(label_row, levels = rev(unique(label_row))))

  p_dumbbell <- ggplot(plot_df) +
    geom_segment(aes(x = NES_original, xend = NES_loo, y = label_row, yend = label_row, color = verdict),
                 linewidth = 1.1, alpha = 0.6) +
    geom_point(aes(x = NES_original, y = label_row), shape = 21, fill = "grey30", color = "black", size = 3) +
    geom_point(aes(x = NES_loo, y = label_row, color = verdict), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = verdict_colours, name = "LOO verdict") +
    labs(title = "Leave-one-protein-out sensitivity: NES original (grey) vs NES without protein (colored)",
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 9))

  out_fig <- file.path(paths$figures, ANALYSIS8_DIR, "target_protein_highlights")
  dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    cairo_pdf(file.path(out_fig, "LOO_sensitivity_dumbbell.pdf"), width = 10, height = max(4, 0.5 * nrow(plot_df) + 2))
    print(p_dumbbell); dev.off()
    png(file.path(out_fig, "LOO_sensitivity_dumbbell.png"), width = 10, height = max(4, 0.5 * nrow(plot_df) + 2),
        units = "in", res = 300, bg = "white")
    print(p_dumbbell); dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("WARN: dumbbell plot save failed — ", conditionMessage(e))
  })
  message("Dumbbell plot saved: ", out_fig)
}

message("\n09c_target_protein_leave_one_out_check.R complete")
