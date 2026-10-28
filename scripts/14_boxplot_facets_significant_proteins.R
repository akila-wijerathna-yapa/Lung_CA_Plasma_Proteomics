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
message("Checkpoint 01 loaded")

checkpoint_13_path <- file.path(paths$rdata, "checkpoint_13.RData")
if (!file.exists(checkpoint_13_path))
  stop("checkpoint_13.RData not found. Run 13_exploratory_subgroup_analyses.R first.")
load(checkpoint_13_path)
message("Checkpoint 13 loaded")
while (dev.cur() > 1) dev.off()

# ── Output directories ─────────────────────────────────────────
tbl_dir <- file.path(paths$tables,  "12_boxplot_facets")
fig_dir <- file.path(paths$figures, "12_boxplot_facets")
invisible(lapply(c(tbl_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Local constants (must match 13_exploratory_subgroup_analyses.R) ─
PVAL_CUTOFF <- 0.1   # nominal p-value threshold (NOT FDR-adjusted)
LFC_THR     <- 0.5
NCOL_FACET  <- 4
PANEL_W     <- 3.1
PANEL_H     <- 2.7

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
  png(paste0(path_stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# Identify the nominally significant proteins per direction, using the
# SAME rule as the volcano plots (nominal p_value < PVAL_CUTOFF & |logFC| >= LFC_THR).
get_significant_features <- function(result_df, feature_meta,
                                     lfc_thr = LFC_THR, p_thr = PVAL_CUTOFF) {
  if (is.null(result_df) || nrow(result_df) == 0)
    return(list(higher_a = tibble(), higher_b = tibble()))
  df <- result_df
  if (!"Protein_Name" %in% names(df))
    df <- df %>% left_join(feature_meta %>% select(nELISA_ID, Protein_Name), by = "nELISA_ID")
  df <- df %>% mutate(Protein_Name = dplyr::coalesce(Protein_Name, nELISA_ID))
  higher_a <- df %>%
    filter(!is.na(p_value), p_value < p_thr, logFC <= -lfc_thr) %>%
    arrange(p_value)
  higher_b <- df %>%
    filter(!is.na(p_value), p_value < p_thr, logFC >=  lfc_thr) %>%
    arrange(p_value)
  list(higher_a = higher_a, higher_b = higher_b)
}

# Build a long sample-level pg/mL table for a chosen set of features.
build_long_pgml <- function(mat_pgml, sample_df, sig_tbl, group_col,
                           sample_id_col = "sample_id") {
  if (is.null(sig_tbl) || nrow(sig_tbl) == 0) return(tibble())
  feature_ids <- intersect(sig_tbl$nELISA_ID, colnames(mat_pgml))
  if (!length(feature_ids)) return(tibble())
  ids_present <- intersect(sample_df[[sample_id_col]], rownames(mat_pgml))
  sub_df  <- sample_df %>% filter(.data[[sample_id_col]] %in% ids_present)
  sub_mat <- mat_pgml[sub_df[[sample_id_col]], feature_ids, drop = FALSE]
  long_df <- as_tibble(sub_mat, rownames = sample_id_col) %>%
    tidyr::pivot_longer(cols = -all_of(sample_id_col),
                        names_to = "nELISA_ID", values_to = "pgml_value") %>%
    filter(is.finite(pgml_value)) %>%
    left_join(sub_df, by = sample_id_col) %>%
    left_join(sig_tbl %>% select(nELISA_ID, Protein_Name, logFC, p_value, adj_p_value),
             by = "nELISA_ID") %>%
    mutate(
      facet_label = paste0(Protein_Name, "\np=", scales::number(p_value, accuracy = 0.0001),
                           " | BH-p=", scales::number(adj_p_value, accuracy = 0.0001))
    )
  long_df
}

# Non-paired faceted boxplot (groups B, C, D)
make_facet_boxplot <- function(long_df, group_col, group_levels, group_colors,
                              title, path_stub, ncol_facet = NCOL_FACET) {
  if (is.null(long_df) || nrow(long_df) == 0) return(invisible(NULL))
  plot_df <- long_df %>%
    mutate(group_plot = factor(.data[[group_col]], levels = group_levels))
  n_panels <- dplyr::n_distinct(plot_df$facet_label)
  ncol_use <- min(ncol_facet, n_panels)
  nrow_use <- ceiling(n_panels / ncol_use)
  p <- ggplot(plot_df, aes(x = group_plot, y = pgml_value)) +
    geom_boxplot(aes(fill = group_plot), alpha = 0.35, outlier.shape = NA,
                width = 0.55, color = "grey30") +
    geom_jitter(aes(color = group_plot), width = 0.12, size = 1.6, alpha = 0.75) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~facet_label, ncol = ncol_use, scales = "free_y") +
    labs(title = title, x = NULL, y = "pg/mL", fill = NULL, color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top",
          plot.title = element_text(size = 12, face = "bold"),
          strip.text = element_text(face = "bold", size = 8))
  save_ggplot(p, path_stub,
             width  = max(7, ncol_use * PANEL_W),
             height = max(4.5, nrow_use * PANEL_H + 1.1))
}

# Paired faceted boxplot with per-patient connecting lines (surgery only)
make_facet_boxplot_paired <- function(long_df, time_col, time_levels, group_colors,
                                      pair_id_col, title, path_stub, ncol_facet = NCOL_FACET) {
  if (is.null(long_df) || nrow(long_df) == 0) return(invisible(NULL))
  plot_df <- long_df %>%
    mutate(timepoint = factor(.data[[time_col]], levels = time_levels))
  n_panels <- dplyr::n_distinct(plot_df$facet_label)
  ncol_use <- min(ncol_facet, n_panels)
  nrow_use <- ceiling(n_panels / ncol_use)
  p <- ggplot(plot_df, aes(x = timepoint, y = pgml_value)) +
    geom_boxplot(aes(fill = timepoint), alpha = 0.35, outlier.shape = NA,
                width = 0.55, color = "grey30") +
    geom_line(aes(group = .data[[pair_id_col]]), alpha = 0.4, color = "grey65", linewidth = 0.35) +
    geom_point(aes(color = timepoint), size = 1.7) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~facet_label, ncol = ncol_use, scales = "free_y") +
    labs(title = title, x = NULL, y = "pg/mL", fill = NULL, color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top",
          plot.title = element_text(size = 12, face = "bold"),
          strip.text = element_text(face = "bold", size = 8))
  save_ggplot(p, path_stub,
             width  = max(7, ncol_use * PANEL_W),
             height = max(4.5, nrow_use * PANEL_H + 1.1))
}

run_direction_pair <- function(result_df, sample_df, mat_pgml, feature_meta,
                              group_col, group_levels, group_colors,
                              group_a_label, group_b_label,
                              file_stub, title_stub,
                              paired = FALSE, pair_id_col = "patient_id") {
  sig <- get_significant_features(result_df, feature_meta)

  for (side in c("higher_a", "higher_b")) {
    sig_tbl <- sig[[side]]
    side_label <- if (side == "higher_a") group_a_label else group_b_label
    stub_side  <- if (side == "higher_a") "higher_in_a" else "higher_in_b"
    if (nrow(sig_tbl) == 0) {
      message(title_stub, " — Higher in ", side_label, ": no nominally significant proteins — skipped")
      next
    }
    safe_write_csv(
      sig_tbl %>% select(nELISA_ID, Protein_Name, logFC, p_value, adj_p_value),
      file.path(tbl_dir, paste0(file_stub, "_", stub_side, "_protein_list.csv"))
    )
    long_df <- build_long_pgml(mat_pgml, sample_df, sig_tbl, group_col)
    if (paired) {
      make_facet_boxplot_paired(
        long_df, time_col = group_col, time_levels = group_levels,
        group_colors = group_colors, pair_id_col = pair_id_col,
        title = paste0(title_stub, " — Higher in ", side_label, " (n=", nrow(sig_tbl), " proteins)"),
        path_stub = file.path(fig_dir, paste0(file_stub, "_", stub_side, "_boxplots"))
      )
    } else {
      make_facet_boxplot(
        long_df, group_col = group_col, group_levels = group_levels, group_colors = group_colors,
        title = paste0(title_stub, " — Higher in ", side_label, " (n=", nrow(sig_tbl), " proteins)"),
        path_stub = file.path(fig_dir, paste0(file_stub, "_", stub_side, "_boxplots"))
      )
    }
    message(title_stub, " — Higher in ", side_label, ": ", nrow(sig_tbl),
            " proteins plotted")
  }
}

GROUP_COLORS <- c("#2166AC", "#B2182B")

# ============================================================
# (A) Surgery: Post vs Pre (paired) — rebuild sample_df exactly
#     as in 13_exploratory_subgroup_analyses.R
# ============================================================
message("\n=== A | Surgery: Post vs Pre — boxplot facets ===")
surgery_pair_ids  <- paired_clinical %>% filter(new_group == "Sx") %>% pull(patient_id)
paired_pre_sx     <- pre_matched  %>% filter(patient_id %in% surgery_pair_ids)
paired_post_sx    <- post_matched %>% filter(patient_id %in% surgery_pair_ids)
surgery_sample_df <- bind_rows(
  paired_pre_sx  %>% mutate(pre_post = "pre"),
  paired_post_sx %>% mutate(pre_post = "post")
)

run_direction_pair(
  surgery_ppost_result, surgery_sample_df, expr_pgml_rep, feature_meta,
  group_col = "pre_post", group_levels = c("pre", "post"), group_colors = GROUP_COLORS,
  group_a_label = "Pre-surgery", group_b_label = "Post-surgery",
  file_stub = "A_surgery_post_vs_pre", title_stub = "Surgery: Post vs Pre",
  paired = TRUE, pair_id_col = "patient_id"
)

# ============================================================
# (B) Immunotherapy response: PR/CR vs SD/PD (new_group == "IO")
# ============================================================
message("\n=== B | Immunotherapy response — boxplot facets ===")
io_df <- pre_matched %>%
  filter(new_group == "IO") %>%
  mutate(response_binary2 = response_binary)

run_direction_pair(
  io_response_result, io_df, expr_pgml_rep, feature_meta,
  group_col = "response_binary2", group_levels = c("SD/PD", "PR/CR"), group_colors = GROUP_COLORS,
  group_a_label = "SD/PD", group_b_label = "PR/CR",
  file_stub = "B_immunotherapy_response_prcr_vs_sdpd",
  title_stub = "Immunotherapy response: PR/CR vs SD/PD"
)

# ============================================================
# (C) Surgery recurrence: Recurrence vs No recurrence (baseline PRE)
# ============================================================
message("\n=== C | Surgery recurrence — boxplot facets ===")
sx_df <- pre_matched %>%
  filter(new_group == "Sx") %>%
  mutate(recurrence_binary2 = recurrence_binary)

run_direction_pair(
  sx_recurrence_result, sx_df, expr_pgml_rep, feature_meta,
  group_col = "recurrence_binary2", group_levels = c("No recurrence", "Recurrence"),
  group_colors = GROUP_COLORS,
  group_a_label = "No recurrence", group_b_label = "Recurrence",
  file_stub = "C_surgery_recurrence_vs_no_recurrence",
  title_stub = "Surgery: Recurrence vs No recurrence"
)

# ============================================================
# (D) Chemotherapy response: PR/CR vs SD/PD (new_group == "CTx")
# ============================================================
message("\n=== D | Chemotherapy response — boxplot facets ===")
ctx_df <- pre_matched %>%
  filter(new_group == "CTx") %>%
  mutate(response_binary2 = response_binary)

run_direction_pair(
  ctx_response_result, ctx_df, expr_pgml_rep, feature_meta,
  group_col = "response_binary2", group_levels = c("SD/PD", "PR/CR"), group_colors = GROUP_COLORS,
  group_a_label = "SD/PD", group_b_label = "PR/CR",
  file_stub = "D_chemotherapy_response_prcr_vs_sdpd",
  title_stub = "Chemotherapy response: PR/CR vs SD/PD"
)

message("\n14_boxplot_facets_significant_proteins.R complete")
