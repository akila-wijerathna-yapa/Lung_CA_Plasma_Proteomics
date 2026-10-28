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
load(file.path(paths$rdata, "checkpoint_05.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 03, 05, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
RNG_SEED       <- 42
FGSEA_FDR      <- 0.25    # fgsea significance threshold
MIN_GENE_SET   <- 5       # minimum genes from ranked list in pathway
MAX_GENE_SET   <- 500
TOP_PLOT       <- 10      # top pathways per collection in lollipop
EMAP_TOP       <- 10      # top pathways in enrichment map
RUN_STRING     <- FALSE   # set TRUE to query STRING API (requires internet)

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

get_msigdbr_table <- function(collection, subcollection = NULL) {
  fn_args <- names(formals(msigdbr::msigdbr))
  args    <- list(species = "Homo sapiens")
  if ("collection" %in% fn_args) {
    args$collection <- collection
    if (!is.null(subcollection)) args$subcollection <- subcollection
  } else {
    args$category <- collection
    if (!is.null(subcollection)) args$subcategory <- subcollection
  }
  do.call(msigdbr::msigdbr, args)
}

build_gene_sets <- function() {
  kegg_medicus <- tryCatch(get_msigdbr_table("C2","CP:KEGG_MEDICUS"), error=function(e) tibble())
  kegg_legacy  <- tryCatch(get_msigdbr_table("C2","CP:KEGG_LEGACY"),  error=function(e) tibble())
  sets <- list(
    HALLMARK = get_msigdbr_table("H"),
    GO_BP    = get_msigdbr_table("C5","GO:BP"),
    GO_MF    = get_msigdbr_table("C5","GO:MF"),
    GO_CC    = get_msigdbr_table("C5","GO:CC"),
    KEGG     = bind_rows(
      kegg_medicus %>% mutate(gs_name = paste0("KEGG_MEDICUS_", gs_name)),
      kegg_legacy  %>% mutate(gs_name = paste0("KEGG_LEGACY_",  gs_name))
    ),
    REACTOME = get_msigdbr_table("C2","CP:REACTOME")
  )
  out <- lapply(sets, function(df) {
    if (is.null(df) || nrow(df) == 0) return(list())
    df %>%
      select(gs_name, gene_symbol) %>%
      filter(!is.na(gs_name), !is.na(gene_symbol), gene_symbol != "") %>%
      distinct() %>%
      split(.$gs_name) %>%
      lapply(function(tbl) unique(tbl$gene_symbol))
  })
  out[lengths(out) > 0]
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

run_fgsea_one_collection <- function(pathways_use, stats, label, collection_name,
                                     coll_tbl_dir, coll_fig_dir) {
  pathways_sub <- lapply(pathways_use, function(gs) intersect(gs, names(stats)))
  pathways_sub <- pathways_sub[lengths(pathways_sub) >= MIN_GENE_SET]
  if (!length(pathways_sub)) return(invisible(NULL))

  set.seed(RNG_SEED)
  fg <- tryCatch(
    fgsea::fgsea(pathways = pathways_sub, stats = stats,
                 minSize = MIN_GENE_SET, maxSize = MAX_GENE_SET, eps = 0),
    error = function(e) { message("  fgsea error (", collection_name, "): ", conditionMessage(e)); NULL })
  if (is.null(fg) || nrow(fg) == 0) return(invisible(NULL))

  fg <- fg %>% as_tibble() %>%
    mutate(
      pathway_clean = clean_pathway_name(pathway),
      leadingEdge   = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    ) %>%
    arrange(padj, desc(abs(NES)))

  sig_fg <- fg %>% filter(!is.na(padj), padj < FGSEA_FDR)

  # Write full results — use subfolder/short filenames to respect Windows MAX_PATH
  lbl_short <- sanitize_filename(label)
  col_short  <- tolower(collection_name)
  safe_write_csv(fg,
    file.path(coll_tbl_dir, paste0(lbl_short, "_", col_short, "_fgsea.csv")))
  safe_write_csv(sig_fg,
    file.path(coll_tbl_dir, paste0(lbl_short, "_", col_short, "_fgsea_sig.csv")))
  safe_write_csv(
    fg %>% select(pathway, pathway_clean, NES, pval, padj, leadingEdge) %>%
      tidyr::separate_rows(leadingEdge, sep = ";") %>%
      rename(Gene_Symbol = leadingEdge) %>%
      filter(!is.na(Gene_Symbol), Gene_Symbol != ""),
    file.path(coll_tbl_dir, paste0(lbl_short, "_", col_short, "_leading_edge.csv")))

  # Lollipop plot — top pathways
  top_plot <- fg %>%
    filter(!is.na(padj), is.finite(NES)) %>%
    arrange(padj, desc(abs(NES))) %>%
    slice_head(n = min(TOP_PLOT, nrow(.))) %>%
    mutate(
      rank             = row_number(),
      abs_nes          = abs(NES),
      neglog_fdr       = -log10(pmax(padj, 1e-300)),
      core_gene_count  = vapply(strsplit(dplyr::coalesce(leadingEdge,""), ";", fixed=TRUE),
                                function(x) sum(nzchar(x)), integer(1)),
      pathway_plot     = paste0(rank, ". ", pathway_clean),
      pathway_plot     = factor(stringr::str_wrap(pathway_plot, width = 48),
                                levels = rev(stringr::str_wrap(pathway_plot, width = 48)))
    )
  if (!nrow(top_plot)) return(invisible(NULL))

  subtitle_str <- if (nrow(sig_fg) >= 1) {
    paste(gsub("_"," ", label), "| ranked by adjusted p")
  } else {
    paste(gsub("_"," ", label), "| top pathways by p; none reached FDR<0.25")
  }
  x_pad <- max(top_plot$abs_nes, na.rm=TRUE) * 0.22
  p <- ggplot(top_plot, aes(abs_nes, pathway_plot)) +
    geom_segment(aes(x=0, xend=abs_nes, yend=pathway_plot), color="grey82", linewidth=0.7) +
    geom_point(aes(size=core_gene_count, fill=neglog_fdr), shape=21, color="#53306d", stroke=0.7) +
    scale_fill_gradientn(
      colors = c("#f3e5f5","#d499eb","#a334d0","#6f1597"),
      name   = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(4, 13), name = "Core genes") +
    labs(title    = paste(collection_name, "top pathways"),
         subtitle = stringr::str_wrap(subtitle_str, width=68),
         x = "Normalized Enrichment Score (|NES|)", y = NULL) +
    coord_cartesian(xlim = c(0, max(top_plot$abs_nes, na.rm=TRUE) + x_pad), clip="off") +
    theme_minimal(base_size = 12) +
    theme(legend.position="right", axis.text.y=element_text(face="bold"),
          panel.grid.major.y=element_blank(),
          plot.title=element_text(size=12.5, face="bold"),
          plot.subtitle=element_text(size=10.2), plot.margin=ggplot2::margin(8,10,8,8))
  tryCatch({
    cairo_pdf(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_fgsea.pdf")),
              width=9, height=6)
    print(p); dev.off()
    png(file.path(coll_fig_dir, paste0(lbl_short, "_", col_short, "_fgsea.png")),
        width=9, height=6, units="in", res=300, bg="white")
    print(p); dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("    WARN: could not save lollipop plot — ", conditionMessage(e),
            "\n    -> Close any open PDF viewers and re-run.")
  })

  # Enrichment plots for top positive + top negative pathway
  enrich_candidates <- bind_rows(
    fg %>% filter(NES > 0) %>% slice_min(order_by=padj, n=1, with_ties=FALSE),
    fg %>% filter(NES < 0) %>% slice_min(order_by=padj, n=1, with_ties=FALSE)
  ) %>% distinct(pathway, .keep_all=TRUE)
  for (i in seq_len(nrow(enrich_candidates))) {
    pw_i <- enrich_candidates$pathway[i]
    if (!pw_i %in% names(pathways_sub)) next
    p_enrich <- fgsea::plotEnrichment(pathways_sub[[pw_i]], stats) +
      labs(title    = paste(label, "-", collection_name),
           subtitle = paste0(clean_pathway_name(pw_i),
                             " | NES=", round(enrich_candidates$NES[i],2),
                             ", FDR=", scales::number(enrich_candidates$padj[i], accuracy=0.001)),
           x = "Ranked proteins", y = "Enrichment score") +
      theme_bw(base_size=11)
    fname_enrich <- sanitize_filename(clean_pathway_name(pw_i))
    tryCatch({
      cairo_pdf(file.path(coll_fig_dir,
                          paste0(lbl_short, "_", col_short, "_enrich_", fname_enrich, ".pdf")),
                width=8.6, height=5.4)
      print(p_enrich); dev.off()
      png(file.path(coll_fig_dir,
                    paste0(lbl_short, "_", col_short, "_enrich_", fname_enrich, ".png")),
          width=8.6, height=5.4, units="in", res=300, bg="white")
      print(p_enrich); dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("    WARN: enrichment plot save failed — ", conditionMessage(e))
    })
  }

  # EMAP (enrichment map) — top pathways by shared leading-edge genes
  emap_nodes <- if (nrow(sig_fg) >= 3) sig_fg else fg %>% arrange(padj, desc(abs(NES)))
  emap_nodes <- emap_nodes %>% slice_head(n = min(EMAP_TOP, nrow(.))) %>%
    mutate(leading_edge_list = strsplit(leadingEdge, ";", fixed=TRUE),
           core_gene_count   = lengths(leading_edge_list))
  if (nrow(emap_nodes) >= 3 && requireNamespace("igraph", quietly=TRUE)) {
    edge_list <- bind_rows(lapply(combn(seq_len(nrow(emap_nodes)), 2, simplify=FALSE),
      function(idx) {
        ga <- unique(emap_nodes$leading_edge_list[[idx[1]]])
        gb <- unique(emap_nodes$leading_edge_list[[idx[2]]])
        shared <- length(intersect(ga, gb))
        union_n <- length(unique(c(ga, gb)))
        jaccard <- ifelse(union_n > 0, shared/union_n, 0)
        if (!is.finite(jaccard) || jaccard < 0.15) return(NULL)
        tibble(from=emap_nodes$pathway_clean[idx[1]], to=emap_nodes$pathway_clean[idx[2]],
               weight=jaccard, shared_genes=shared)
      }))
    vertex_df <- emap_nodes %>%
      transmute(name=pathway_clean, NES, neglog_fdr=-log10(pmax(padj,1e-300)),
                core_gene_count)
    if (is.null(edge_list) || nrow(edge_list) == 0)
      edge_list <- tibble(from=character(), to=character(), weight=numeric(), shared_genes=integer())
    g_emap <- igraph::graph_from_data_frame(edge_list, directed=FALSE, vertices=vertex_df)
    emap_title <- if (nrow(sig_fg) >= 3) {
      paste(label, "-", collection_name, "EMAP (FDR<0.25 pathways)")
    } else {
      paste(label, "-", collection_name, "EMAP top pathways (none FDR<0.25)")
    }
    tryCatch({
      cairo_pdf(file.path(coll_fig_dir, paste0(lbl_short,"_",col_short,"_emap.pdf")),
                width=9.2, height=7.2)
      set.seed(RNG_SEED)
      lay <- if (igraph::ecount(g_emap) > 0) igraph::layout_with_fr(g_emap) else igraph::layout_in_circle(g_emap)
      vcols  <- ifelse(igraph::V(g_emap)$NES >= 0, "#B2182B", "#2166AC")
      vsizes <- scales::rescale(igraph::V(g_emap)$core_gene_count, to=c(14,28))
      ewidth <- if (igraph::ecount(g_emap) > 0) scales::rescale(igraph::E(g_emap)$weight, to=c(1.2,4.5)) else numeric()
      plot(g_emap, layout=lay, vertex.color=vcols, vertex.size=vsizes,
           vertex.label=igraph::V(g_emap)$name, vertex.label.cex=0.72,
           vertex.label.family="sans", vertex.frame.color="grey30",
           edge.width=ewidth, edge.color="grey70", main=emap_title)
      legend("topleft", legend=c("Positive NES","Negative NES"),
             col=c("#B2182B","#2166AC"), pch=16, bty="n", cex=0.9)
      dev.off()
      png(file.path(coll_fig_dir, paste0(lbl_short,"_",col_short,"_emap.png")),
          width=9.2, height=7.2, units="in", res=300, bg="white")
      set.seed(RNG_SEED)
      plot(g_emap, layout=lay, vertex.color=vcols, vertex.size=vsizes,
           vertex.label=igraph::V(g_emap)$name, vertex.label.cex=0.72,
           vertex.label.family="sans", vertex.frame.color="grey30",
           edge.width=ewidth, edge.color="grey70", main=emap_title)
      legend("topleft", legend=c("Positive NES","Negative NES"),
             col=c("#B2182B","#2166AC"), pch=16, bty="n", cex=0.9)
      dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("    WARN: EMAP save failed — ", conditionMessage(e))
    })
  }

  message("    ", collection_name, ": ", nrow(fg), " pathways | FDR<0.25: ", nrow(sig_fg))
  invisible(fg)
}

run_fgsea_suite <- function(result_df, label, base_tbl_dir, base_fig_dir,
                             stat_col = "t_stat") {
  if (!length(gene_sets)) { message("  No gene sets available — skipped"); return(invisible(NULL)) }
  stats <- prepare_rank_vector(result_df, stat_col = stat_col)
  if (length(stats) < 20) { message("  Fewer than 20 ranked genes — skipped"); return(invisible(NULL)) }
  message("  Ranked list: ", length(stats), " genes")
  for (collection_name in names(gene_sets)) {
    # Create label/collection subfolder for Windows MAX_PATH compliance
    coll_tbl_dir <- file.path(base_tbl_dir, sanitize_filename(label), tolower(collection_name))
    coll_fig_dir <- file.path(base_fig_dir, sanitize_filename(label), tolower(collection_name))
    dir.create(coll_tbl_dir, recursive=TRUE, showWarnings=FALSE)
    dir.create(coll_fig_dir, recursive=TRUE, showWarnings=FALSE)
    run_fgsea_one_collection(gene_sets[[collection_name]], stats, label,
                             collection_name, coll_tbl_dir, coll_fig_dir)
  }
}

run_string_ppi <- function(result_df, label, base_tbl_dir, base_fig_dir,
                            max_genes = 30, nominal_fallback = TRUE) {
  lbl_dir <- file.path(base_tbl_dir, sanitize_filename(label))
  dir.create(lbl_dir, recursive=TRUE, showWarnings=FALSE)
  genes_fdr <- result_df %>%
    filter(!is.na(Gene_Symbol), adj_p_value < 0.1) %>%
    arrange(adj_p_value, p_value) %>%
    distinct(Gene_Symbol, .keep_all=TRUE) %>%
    slice_head(n = max_genes)
  genes_tbl       <- genes_fdr
  selection_rule  <- "FDR<0.10"
  if (nominal_fallback && nrow(genes_tbl) < 5 && "p_value" %in% names(result_df)) {
    genes_nom <- result_df %>%
      filter(!is.na(Gene_Symbol), p_value < 0.05) %>%
      arrange(p_value) %>%
      distinct(Gene_Symbol, .keep_all=TRUE) %>%
      slice_head(n = max_genes)
    if (nrow(genes_nom) >= 5) { genes_tbl <- genes_nom; selection_rule <- "Nominal p<0.05" }
  }
  safe_write_csv(
    genes_tbl %>% transmute(Gene_Symbol, selection_rule = selection_rule,
                            adj_p_value, p_value),
    file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_input_genes.csv")))
  ppi_status <- tibble(label=label, genes=nrow(genes_tbl), rule=selection_rule,
                       string_attempted=FALSE, status="Input genes saved")
  if (!RUN_STRING) {
    safe_write_csv(ppi_status, file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_status.csv")))
    return(invisible(NULL))
  }
  if (!requireNamespace("igraph", quietly=TRUE) || nrow(genes_tbl) < 5) {
    ppi_status$status <- "Skipped"
    safe_write_csv(ppi_status, file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_status.csv")))
    return(invisible(NULL))
  }
  genes <- genes_tbl$Gene_Symbol
  string_url <- paste0("https://string-db.org/api/tsv/network?species=9606&identifiers=",
                       utils::URLencode(paste(genes, collapse="\r"), reserved=TRUE))
  raw_tsv <- tryCatch(readLines(string_url, warn=FALSE), error=function(e) character())
  if (length(raw_tsv) <= 1) {
    ppi_status$status <- "STRING API returned no network"
    safe_write_csv(ppi_status, file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_status.csv")))
    return(invisible(NULL))
  }
  api_edges <- tryCatch(readr::read_tsv(I(paste(raw_tsv, collapse="\n")), show_col_types=FALSE),
                        error=function(e) NULL)
  if (is.null(api_edges) || !nrow(api_edges)) {
    ppi_status$status <- "STRING parse failed"
    safe_write_csv(ppi_status, file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_status.csv")))
    return(invisible(NULL))
  }
  edges_use <- api_edges %>%
    transmute(from=preferredName_A, to=preferredName_B, score=score) %>%
    filter(from %in% genes | to %in% genes) %>%
    distinct(from, to, .keep_all=TRUE)
  nodes <- tibble(Gene_Symbol=sort(unique(c(edges_use$from, edges_use$to))))
  safe_write_csv(nodes,     file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_nodes.csv")))
  safe_write_csv(edges_use, file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_edges.csv")))
  ppi_status$status <- "STRING network generated"
  ppi_status$nodes  <- nrow(nodes)
  ppi_status$edges  <- nrow(edges_use)
  safe_write_csv(ppi_status, file.path(lbl_dir, paste0(sanitize_filename(label), "_ppi_status.csv")))
  message("  STRING PPI: ", nrow(nodes), " nodes, ", nrow(edges_use), " edges")
}

# ── Build gene sets ───────────────────────────────────────────
message("\n=== 09 | Building gene sets ===")
gene_sets <- tryCatch(build_gene_sets(), error = function(e) {
  message("Gene set build failed: ", conditionMessage(e)); list()
})
message("Collections loaded: ", paste(names(gene_sets), collapse = ", "))

# ── Output directories ────────────────────────────────────────
base_tbl <- file.path(paths$tables,  ANALYSIS8_DIR)
base_fig <- file.path(paths$figures, ANALYSIS8_DIR)

# ── Define all analyses to enrich ────────────────────────────
# Stat column: t_stat for limma; z for Cox; logFC for descriptive
enrichment_specs <- list(

  # Analysis 1: one-vs-rest (7 treatment groups)
  analysis1_ctx_vs_rest    = list(df = analysis1_ovr_results[["analysis1_ctx_vs_rest"]],    stat = "t_stat"),
  analysis1_io_vs_rest     = list(df = analysis1_ovr_results[["analysis1_io_vs_rest"]],     stat = "t_stat"),
  analysis1_ctxio_vs_rest  = list(df = analysis1_ovr_results[["analysis1_ctxio_vs_rest"]],  stat = "t_stat"),
  analysis1_ctxrtx_vs_rest = list(df = analysis1_ovr_results[["analysis1_ctxrtx_vs_rest"]], stat = "t_stat"),
  analysis1_tki_vs_rest    = list(df = analysis1_ovr_results[["analysis1_tki_vs_rest"]],    stat = "t_stat"),
  analysis1_sx_vs_rest     = list(df = analysis1_ovr_results[["analysis1_sx_vs_rest"]],     stat = "t_stat"),
  analysis1_other_vs_rest  = list(df = analysis1_ovr_results[["analysis1_other_vs_rest"]],  stat = "t_stat"),

  # Analysis 3: response (both main and exploratory)
  analysis3_primary_response    = list(df = analysis3_primary,       stat = "t_stat"),
  analysis3_secondary_response  = list(df = analysis3_secondary,     stat = "t_stat"),
  analysis3_primary_expl        = list(df = analysis3_primary_expl,  stat = "t_stat"),
  analysis3_secondary_expl      = list(df = analysis3_secondary_expl, stat = "t_stat"),

  # Analysis 5: survival
  analysis5_os_survival   = list(df = os_res,  stat = "z"),
  analysis5_pfs_survival  = list(df = pfs_res, stat = "z"),
  analysis5_efs_survival  = list(df = efs_res, stat = "z"),
  analysis5_delta_os      = list(df = delta_os,  stat = "z"),
  analysis5_delta_pfs     = list(df = delta_pfs, stat = "z")
)

# Remove NULL entries (analyses that produced no results)
enrichment_specs <- enrichment_specs[!sapply(enrichment_specs, function(s) {
  is.null(s$df) || nrow(s$df) == 0
})]
message("Enrichment analyses to run: ", length(enrichment_specs))

# ── Run fgsea and STRING for each analysis ─────────────────────
fgsea_run_log <- list()
for (label in names(enrichment_specs)) {
  spec <- enrichment_specs[[label]]
  message("\n--- ", label, " ---")
  run_fgsea_suite(spec$df, label, base_tbl, base_fig, stat_col = spec$stat)
  run_string_ppi(spec$df, label, base_tbl, base_fig)
  fgsea_run_log[[label]] <- list(label = label, n_proteins = nrow(spec$df))
}

fgsea_log_tbl <- bind_rows(lapply(fgsea_run_log, function(x) as_tibble(x)))
safe_write_csv(fgsea_log_tbl, file.path(base_tbl, "09_fgsea_run_log.csv"))

# ── Cross-analysis pathway summary ────────────────────────────
message("\n=== 09 | Cross-analysis pathway summary ===")
sig_pathways_all <- bind_rows(lapply(names(enrichment_specs), function(label) {
  lbl_dir <- file.path(base_tbl, sanitize_filename(label))
  if (!dir.exists(lbl_dir)) return(NULL)
  bind_rows(lapply(names(gene_sets), function(coll) {
    f <- file.path(lbl_dir, tolower(coll),
                   paste0(sanitize_filename(label), "_", tolower(coll), "_fgsea_sig.csv"))
    if (!file.exists(f)) return(NULL)
    tryCatch(readr::read_csv(f, show_col_types=FALSE) %>%
               mutate(across(any_of(c("pval","padj","log2err","ES","NES","size")),
                             ~suppressWarnings(as.numeric(.)))) %>%
               mutate(analysis = label, collection = coll),
             error = function(e) NULL)
  }))
}))
if (!is.null(sig_pathways_all) && nrow(sig_pathways_all) > 0) {
  safe_write_csv(sig_pathways_all %>% arrange(padj, desc(abs(NES))),
    file.path(base_tbl, "09_significant_pathways_all_analyses.csv"))
  message("Significant pathways (FDR<0.25) across all analyses: ", nrow(sig_pathways_all))
} else {
  message("No significant pathways at FDR<0.25 — all results exploratory")
}

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 09 | Saving checkpoint_09.RData ===")
save(
  gene_sets,
  fgsea_log_tbl,
  file = file.path(paths$rdata, "checkpoint_09.RData")
)
message("09_pathway_network.R complete")
