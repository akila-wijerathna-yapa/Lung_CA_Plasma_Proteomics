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

for (pkg in c("UpSetR", "cowplot", "gridExtra")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(UpSetR)

# ── Load checkpoints ───────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")

checkpoint_16_path <- file.path(paths$rdata, "checkpoint_16.RData")
if (!file.exists(checkpoint_16_path))
  stop("checkpoint_16.RData not found. Run 16_site_specific_biomarker_overlap.R first.")
load(checkpoint_16_path)
message("Checkpoint 16 loaded — ", length(site_hits), " sites with hits")
while (dev.cur() > 1) dev.off()

fig_dir <- file.path(paths$figures, "14_site_specific_biomarkers")
tbl_dir <- file.path(paths$tables,  "14_site_specific_biomarkers")
invisible(lapply(c(fig_dir, tbl_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Local constants ────────────────────────────────────────────
NCOL_FACET      <- 4
PANEL_W         <- 3.1
PANEL_H         <- 2.7
MAX_LIST_CHARS  <- 90   # per-cell truncation in the visual table (full list always in CSV)
GROUP_COLORS    <- c("No" = "#2166AC", "Yes" = "#B2182B")

# ── Local utility functions ────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")
sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  ifelse(nzchar(out), out, "item")
}
save_ggplot <- function(p, path_stub, width, height) {
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height)
  print(p); dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# ============================================================
# (1) Derive bins (exact site-combinations) from membership_df
# ============================================================
message("\n=== Deriving intersection bins ===")
site_names <- names(site_hits)

bin_df <- membership_df %>%
  mutate(bin_key = apply(across(all_of(site_names)), 1,
                        function(row) paste(site_names[row == 1], collapse = "+")))

bin_table <- bin_df %>%
  group_by(bin_key) %>%
  summarise(
    n_sites      = n_distinct(unlist(strsplit(bin_key[1], "\\+"))),
    n_proteins   = n(),
    protein_list = paste(sort(Protein_Name), collapse = ", "),
    nELISA_ids   = list(nELISA_ID),
    .groups = "drop"
  ) %>%
  arrange(desc(n_proteins))

safe_write_csv(bin_table %>% select(-nELISA_ids), file.path(tbl_dir, "bin_protein_list.csv"))
message("Bins found: ", nrow(bin_table))
print(bin_table %>% select(bin_key, n_sites, n_proteins))

# ============================================================
# (2) UpSet plot + protein table stacked underneath
# ============================================================
message("\n=== UpSet plot with protein table ===")
upset_input <- as.data.frame(membership_df %>% select(all_of(site_names)))
rownames(upset_input) <- membership_df$Protein_Name

res_upset <- UpSetR::upset(upset_input, sets = site_names, order.by = "freq",
                           main.bar.color = "#2166AC", sets.bar.color = "#B2182B",
                           text.scale = 1.3, keep.order = FALSE)

table_display <- bin_table %>%
  transmute(
    Bin = bin_key,
    N   = n_proteins,
    Proteins = ifelse(nchar(protein_list) > MAX_LIST_CHARS,
                      paste0(substr(protein_list, 1, MAX_LIST_CHARS), " ... (+",
                             n_proteins - lengths(regmatches(substr(protein_list, 1, MAX_LIST_CHARS),
                                                            gregexpr(",", substr(protein_list, 1, MAX_LIST_CHARS)))) - 1,
                             " more, see CSV)"),
                      protein_list)
  )

table_grob <- gridExtra::tableGrob(
  table_display, rows = NULL,
  theme = gridExtra::ttheme_minimal(
    core    = list(fg_params = list(hjust = 0, x = 0.02, fontsize = 7)),
    colhead = list(fg_params = list(fontsize = 8, fontface = "bold"))
  )
)

combined <- cowplot::plot_grid(
  res_upset$Main_bar, res_upset$Matrix, table_grob,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(2.4, 1.4, max(1.2, 0.32 * nrow(table_display)))
)

fig_height <- 6 + max(1.2, 0.32 * nrow(table_display))
save_ggplot(combined, file.path(fig_dir, "site_biomarker_upset_with_table"),
           width = 10, height = fig_height)

# ============================================================
# (3) Faceted boxplots per bin, vs each site in that bin's combination
# ============================================================
message("\n=== Faceted boxplots per bin ===")

build_long_pgml <- function(mat_pgml, sample_df, sig_tbl, group_col, sample_id_col = "sample_id") {
  if (is.null(sig_tbl) || nrow(sig_tbl) == 0) return(tibble())
  feature_ids <- intersect(sig_tbl$nELISA_ID, colnames(mat_pgml))
  if (!length(feature_ids)) return(tibble())
  ids_present <- intersect(sample_df[[sample_id_col]], rownames(mat_pgml))
  sub_df  <- sample_df %>% filter(.data[[sample_id_col]] %in% ids_present)
  sub_mat <- mat_pgml[sub_df[[sample_id_col]], feature_ids, drop = FALSE]
  as_tibble(sub_mat, rownames = sample_id_col) %>%
    tidyr::pivot_longer(cols = -all_of(sample_id_col), names_to = "nELISA_ID", values_to = "pgml_value") %>%
    filter(is.finite(pgml_value)) %>%
    left_join(sub_df, by = sample_id_col) %>%
    left_join(sig_tbl %>% select(nELISA_ID, Protein_Name, p_value, adj_p_value), by = "nELISA_ID") %>%
    mutate(facet_label = paste0(Protein_Name, "\np=", scales::number(p_value, accuracy = 0.0001),
                                " | BH-p=", scales::number(adj_p_value, accuracy = 0.0001)))
}

make_facet_boxplot <- function(long_df, group_col, group_levels, group_colors,
                              title, path_stub, ncol_facet = NCOL_FACET) {
  if (is.null(long_df) || nrow(long_df) == 0) return(invisible(NULL))
  plot_df  <- long_df %>% mutate(group_plot = factor(.data[[group_col]], levels = group_levels))
  n_panels <- dplyr::n_distinct(plot_df$facet_label)
  ncol_use <- min(ncol_facet, n_panels); nrow_use <- ceiling(n_panels / ncol_use)
  p <- ggplot(plot_df, aes(x = group_plot, y = pgml_value)) +
    geom_boxplot(aes(fill = group_plot), alpha = 0.35, outlier.shape = NA, width = 0.55, color = "grey30") +
    geom_jitter(aes(color = group_plot), width = 0.12, size = 1.6, alpha = 0.75) +
    scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors) +
    facet_wrap(~facet_label, ncol = ncol_use, scales = "free_y") +
    labs(title = title, x = NULL, y = "pg/mL", fill = NULL, color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", plot.title = element_text(size = 12, face = "bold"),
          strip.text = element_text(face = "bold", size = 8))
  save_ggplot(p, path_stub, width = max(7, ncol_use * PANEL_W), height = max(4.5, nrow_use * PANEL_H + 1.1))
}

for (i in seq_len(nrow(bin_table))) {
  bin_key      <- bin_table$bin_key[i]
  sites_in_bin <- strsplit(bin_key, "\\+")[[1]]
  ids_in_bin   <- bin_table$nELISA_ids[[i]]
  bin_stub     <- paste0("bin", sprintf("%02d", i), "_", sanitize_filename(bin_key))

  for (s in sites_in_bin) {
    sig_tbl <- site_de_results[[s]] %>% filter(nELISA_ID %in% ids_in_bin)
    if (nrow(sig_tbl) == 0) next
    site_status_df <- cohort_df %>%
      mutate(site_status = factor(ifelse(.data[[s]] == 1, "Yes", "No"), levels = c("No", "Yes")))
    long_df <- build_long_pgml(expr_pgml_rep, site_status_df, sig_tbl, "site_status")
    make_facet_boxplot(
      long_df, group_col = "site_status", group_levels = c("No", "Yes"), group_colors = GROUP_COLORS,
      title = paste0("Bin [", bin_key, "] — vs ", s, " status (n=", nrow(sig_tbl), " proteins)"),
      path_stub = file.path(fig_dir, paste0(bin_stub, "_vs_", s, "_boxplots"))
    )
    message("Bin [", bin_key, "] vs ", s, ": ", nrow(sig_tbl), " proteins plotted")
  }
}

message("\n17_site_bin_protein_tables_boxplots.R complete")
