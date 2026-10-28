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
TARGET_PROTEINS   <- c("GDF15", "GZMA", "CEACAM5", "SHH", "NOTUM", "PLAUR")
TOP_N_PER_PROTEIN  <- 2      # pathways to plot per protein, ranked by padj
if (!exists("ANALYSIS8_DIR")) ANALYSIS8_DIR <- "08_pathway_and_network_analysis"

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  out <- substr(out, 1, 60)
  ifelse(nzchar(out), out, "item")
}

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

# ── Locate every pathway in which each target protein sits in the leading edge ──
message("\n=== 09b | Scanning leading_edge.csv files for target proteins ===")
base_tbl <- file.path(paths$tables, ANALYSIS8_DIR)
le_files <- list.files(base_tbl, pattern = "_leading_edge\\.csv$", recursive = TRUE, full.names = TRUE)
message("Leading-edge files found: ", length(le_files))

hit_list <- bind_rows(lapply(le_files, function(f) {
  d <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d <- d %>% filter(Gene_Symbol %in% TARGET_PROTEINS)
  if (!nrow(d)) return(NULL)
  parts     <- basename(f)
  coll_dir  <- basename(dirname(f))
  lbl_dir   <- basename(dirname(dirname(f)))
  d %>% mutate(source_file = f, label = lbl_dir, collection = coll_dir)
}))

if (is.null(hit_list) || !nrow(hit_list)) {
  message("No target proteins found in any leading edge across all analyses — stopping.")
  quit(save = "no")
}

hit_ranked <- hit_list %>%
  distinct(Gene_Symbol, label, collection, pathway, .keep_all = TRUE) %>%
  arrange(Gene_Symbol, padj, desc(abs(NES))) %>%
  group_by(Gene_Symbol) %>%
  slice_head(n = TOP_N_PER_PROTEIN) %>%
  ungroup()

safe_write_csv(hit_ranked, file.path(base_tbl, "target_protein_pathway_hits.csv"))
message("Pathway hits selected for plotting: ", nrow(hit_ranked))
missing_proteins <- setdiff(TARGET_PROTEINS, unique(hit_ranked$Gene_Symbol))
if (length(missing_proteins)) {
  message("No pathway-level leading-edge evidence found for: ", paste(missing_proteins, collapse = ", "))
}

# ── Output directory ─────────────────────────────────────────
out_fig <- file.path(paths$figures, ANALYSIS8_DIR, "target_protein_highlights")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# ── Rebuild stats vector per label once, then plot each hit ──
message("\n=== 09b | Rebuilding ranked stats vectors and plotting ===")
stats_cache <- list()

for (i in seq_len(nrow(hit_ranked))) {
  gene_i  <- hit_ranked$Gene_Symbol[i]
  label_i <- hit_ranked$label[i]
  coll_i  <- hit_ranked$collection[i]
  pw_i    <- hit_ranked$pathway[i]
  NES_i   <- hit_ranked$NES[i]
  padj_i  <- hit_ranked$padj[i]

  spec_i <- enrichment_specs[[label_i]]
  if (is.null(spec_i)) { message("  Skip: no result_df for label ", label_i); next }

  if (is.null(stats_cache[[label_i]])) {
    stats_cache[[label_i]] <- prepare_rank_vector(spec_i$df, stat_col = spec_i$stat)
  }
  stats_i <- stats_cache[[label_i]]

  coll_name <- toupper(coll_i)
  pathway_genes <- gene_sets[[coll_name]][[pw_i]]
  if (is.null(pathway_genes)) { message("  Skip: pathway not found — ", pw_i, " (", coll_name, ")"); next }
  pathway_genes <- intersect(pathway_genes, names(stats_i))
  if (length(pathway_genes) < 5 || !(gene_i %in% names(stats_i))) {
    message("  Skip: insufficient overlap for ", gene_i, " / ", pw_i); next
  }

  statsAdj  <- sort(stats_i, decreasing = TRUE)
  rank_gene <- match(gene_i, names(statsAdj))
  pctile    <- round(100 * rank_gene / length(statsAdj), 1)

  p_enrich <- fgsea::plotEnrichment(pathway_genes, stats_i) +
    geom_vline(xintercept = rank_gene, color = "firebrick", linetype = "dashed", linewidth = 0.8) +
    annotate("point", x = rank_gene, y = 0, color = "firebrick", size = 2.6) +
    annotate("text", x = rank_gene, y = Inf, label = gene_i, color = "firebrick",
             fontface = "bold", vjust = 1.3, hjust = -0.15, size = 3.6) +
    labs(title    = paste(gene_i, "-", gsub("_", " ", label_i)),
         subtitle = paste0(clean_pathway_name(pw_i), " (", coll_name, ") | NES=", round(NES_i, 2),
                           ", FDR=", scales::number(padj_i, accuracy = 0.001),
                           " | ", gene_i, " rank ", rank_gene, "/", length(statsAdj),
                           " (top ", pctile, "%)"),
         x = "Ranked proteins", y = "Enrichment score") +
    theme_bw(base_size = 11)

  fname <- paste0(gene_i, "_", sanitize_filename(label_i), "_", coll_i, "_",
                  sanitize_filename(clean_pathway_name(pw_i)))
  tryCatch({
    cairo_pdf(file.path(out_fig, paste0(fname, ".pdf")), width = 8.6, height = 5.4)
    print(p_enrich); dev.off()
    png(file.path(out_fig, paste0(fname, ".png")), width = 8.6, height = 5.4, units = "in", res = 300, bg = "white")
    print(p_enrich); dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("    WARN: save failed for ", fname, " — ", conditionMessage(e))
  })
  message("  Plotted: ", fname)
}

message("\n09b_target_protein_enrichment_plots.R complete. Output: ", out_fig)
