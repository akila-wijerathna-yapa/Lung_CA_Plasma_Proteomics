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

# ── Load checkpoint ────────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")
while (dev.cur() > 1) dev.off()

# ── Load curated recurrence-sites file ─────────────────────────
recurrence_sites_path <- file.path(PROJECT_DIR, "data", "1L_Sites_Recurrence_curated.xlsx")
if (!file.exists(recurrence_sites_path))
  stop("1L_Sites_Recurrence_curated.xlsx not found at: ", recurrence_sites_path,
       " — place the uploaded file under PROJECT_DIR/data/")

norm_id_local <- function(x) {
  out <- as.character(x)
  out <- trimws(out)
  out <- gsub("\\.0$", "", out)
  out <- gsub("[^A-Za-z0-9]+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  toupper(out)
}

site_cols <- c("brain", "bone", "liver", "lung", "lymph_node", "pleural", "adrenal", "other_site")
recurrence_sites_raw <- readxl::read_excel(recurrence_sites_path)
recurrence_sites <- recurrence_sites_raw %>%
  transmute(
    patient_id            = norm_id_local(Patient_ID),
    recurrence_date_file   = suppressWarnings(as.Date(`1L_Date_Recurrence`)),
    recurrence_site_text   = as.character(`1L_Sites_Recurrence`),
    across(all_of(site_cols), ~ as.integer(coalesce(as.numeric(.x), 0)))
  )
message("Recurrence-sites file: ", nrow(recurrence_sites), " patients")

# ── Output directories ─────────────────────────────────────────
tbl_dir <- file.path(paths$tables,  "13_responder_mets_vs_surgery_recurrence")
fig_dir <- file.path(paths$figures, "13_responder_mets_vs_surgery_recurrence")
invisible(lapply(c(tbl_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Local constants ────────────────────────────────────────────
PVAL_CUTOFF   <- 0.1   # nominal p-value threshold (NOT FDR-adjusted)
LFC_THR       <- 0.5
MIN_PER_GROUP <- 5
NCOL_FACET    <- 4
PANEL_W       <- 3.1
PANEL_H       <- 2.7

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

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

# ============================================================
# Build Group A / Group B
# ============================================================
message("\n=== Build comparison groups ===")

group_a_df <- pre_matched %>%
  filter(
    trimws(toupper(stage_group)) == "IV",
    response_binary == "PR/CR",
    treatment_class != "Surgery"
  ) %>%
  mutate(comparison_group = "Responder_existing_mets")

group_b_ids <- recurrence_sites$patient_id
group_b_df <- pre_matched %>%
  filter(new_group == "Sx", patient_id %in% group_b_ids) %>%
  mutate(comparison_group = "Surgery_recurred_new_mets")

message("Group A (responder, stage IV, non-surgical): n = ", nrow(group_a_df),
        " (", n_distinct(group_a_df$patient_id), " patients)")
message("Group B (surgery, recurred w/ documented site): n = ", nrow(group_b_df),
        " (", n_distinct(group_b_df$patient_id), " patients)")

sx_total   <- pre_matched %>% filter(new_group == "Sx") %>% pull(patient_id) %>% n_distinct()
sx_matched <- n_distinct(group_b_df$patient_id)
message("  -> ", sx_matched, " of ", sx_total,
        " surgery patients matched to the recurrence-sites file")
if (sx_matched == 0)
  warning("Zero surgery patients matched the recurrence-sites file by patient_id — check ID formats.")

if (nrow(group_a_df) < MIN_PER_GROUP || nrow(group_b_df) < MIN_PER_GROUP)
  warning("One or both groups have < ", MIN_PER_GROUP, " samples — results will be unstable.")

comparison_df <- bind_rows(group_a_df, group_b_df)
safe_write_csv(
  comparison_df %>% select(patient_id, sample_id, comparison_group, stage_group,
                           treatment_class, new_group, response_binary,
                           recurrence_binary, sex, histology, age_at_diagnosis),
  file.path(tbl_dir, "00_group_membership.csv")
)

# ============================================================
# (1) Descriptive / clinical comparison
# ============================================================
message("\n=== Descriptive comparison ===")

demo_summary <- comparison_df %>%
  group_by(comparison_group) %>%
  summarise(
    n              = n_distinct(patient_id),
    age_mean       = mean(age_at_diagnosis, na.rm = TRUE),
    age_sd         = sd(age_at_diagnosis,   na.rm = TRUE),
    pct_male       = mean(sex == "MALE", na.rm = TRUE) * 100,
    n_histology    = n_distinct(na.omit(histology)),
    .groups = "drop"
  )
safe_write_csv(demo_summary, file.path(tbl_dir, "01_demographic_summary.csv"))

age_test <- tryCatch(
  wilcox.test(age_at_diagnosis ~ comparison_group, data = comparison_df),
  error = function(e) NULL)
sex_test <- tryCatch(
  chisq.test(table(comparison_df$comparison_group, comparison_df$sex)),
  error = function(e) NULL)
message("Age (Wilcoxon) p = ", if (!is.null(age_test)) signif(age_test$p.value, 3) else NA)
message("Sex (Chi-square) p = ", if (!is.null(sex_test)) signif(sex_test$p.value, 3) else NA)

# Recurrence/metastasis site pattern — available directly for Group B;
# Group A patients are also checked against the same file in case a
# stage-IV responder later progressed with a documented site.
site_pattern <- comparison_df %>%
  distinct(patient_id, comparison_group) %>%
  left_join(recurrence_sites, by = "patient_id") %>%
  group_by(comparison_group) %>%
  summarise(
    n_with_site_data = sum(!is.na(recurrence_site_text)),
    across(all_of(site_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )
safe_write_csv(site_pattern, file.path(tbl_dir, "02_site_pattern_summary.csv"))
message("Group A patients with a documented site in the recurrence file: ",
        site_pattern$n_with_site_data[site_pattern$comparison_group == "Responder_existing_mets"] %||% 0)

site_plot_df <- site_pattern %>%
  tidyr::pivot_longer(all_of(site_cols), names_to = "site", values_to = "n_patients")
p_sites <- ggplot(site_plot_df, aes(x = site, y = n_patients, fill = comparison_group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Responder_existing_mets" = "#2166AC",
                              "Surgery_recurred_new_mets" = "#B2182B")) +
  labs(title = "Recurrence/metastasis site pattern", x = NULL, y = "N patients", fill = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top",
                                   axis.text.x = element_text(angle = 30, hjust = 1))
save_ggplot(p_sites, file.path(fig_dir, "site_pattern_barplot"), width = 7.5, height = 5.5)

# ============================================================
# (2) Proteomic differential expression: Group A vs Group B
# ============================================================
message("\n=== Proteomic DE: Group A vs Group B ===")

run_binary_limma <- function(mat_samples_by_features, sample_df, group_col,
                             group_levels, adjust_cols = character(),
                             min_per_group = MIN_PER_GROUP) {
  keep <- !is.na(sample_df[[group_col]]) & sample_df[[group_col]] %in% group_levels
  df   <- sample_df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  grp  <- factor(df[[group_col]], levels = group_levels)
  if (length(unique(grp)) < 2 || any(table(grp) < min_per_group)) return(NULL)
  mat_sub <- mat_samples_by_features[df$sample_id, , drop = FALSE]
  feature_keep <- colSums(!is.na(mat_sub[grp == group_levels[1], , drop = FALSE])) >= min_per_group &
    colSums(!is.na(mat_sub[grp == group_levels[2], , drop = FALSE])) >= min_per_group
  mat_sub <- mat_sub[, feature_keep, drop = FALSE]
  if (ncol(mat_sub) == 0) return(NULL)
  design_df <- data.frame(primary_group = relevel(grp, ref = group_levels[1]))
  for (col in adjust_cols) {
    vals <- df[[col]]
    if (is.numeric(vals)) {
      if (length(unique(vals[!is.na(vals)])) > 1) design_df[[col]] <- vals
    } else {
      fac <- factor(vals)
      if (length(unique(fac[!is.na(fac)])) > 1) design_df[[col]] <- fac
    }
  }
  design    <- model.matrix(~ ., data = design_df)
  coef_name <- grep("^primary_group", colnames(design), value = TRUE)[1]
  fit <- limma::lmFit(t(mat_sub), design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("nELISA_ID") %>%
    as_tibble() %>%
    transmute(nELISA_ID, logFC, t_stat = t, p_value = P.Value,
              adj_p_value = adj.P.Val, avg_expr = AveExpr)
  means     <- lapply(group_levels, function(gl) colMeans(mat_sub[grp == gl, , drop = FALSE], na.rm = TRUE))
  n_nonmiss <- lapply(group_levels, function(gl) colSums(!is.na(mat_sub[grp == gl, , drop = FALSE])))
  tt %>% mutate(
    group_a      = group_levels[1],
    group_b      = group_levels[2],
    mean_group_a = unname(means[[1]][match(nELISA_ID, names(means[[1]]))]),
    mean_group_b = unname(means[[2]][match(nELISA_ID, names(means[[2]]))]),
    n_group_a    = unname(n_nonmiss[[1]][match(nELISA_ID, names(n_nonmiss[[1]]))]),
    n_group_b    = unname(n_nonmiss[[2]][match(nELISA_ID, names(n_nonmiss[[2]]))])
  )
}

de_result <- run_binary_limma(
  expr_log2_rep, comparison_df, "comparison_group",
  group_levels  = c("Responder_existing_mets", "Surgery_recurred_new_mets"),
  adjust_cols   = c("plate_barcode"),
  min_per_group = MIN_PER_GROUP
)

if (is.null(de_result)) {
  message("DE skipped: insufficient sample size in one or both groups")
} else {
  de_result <- de_result %>%
    left_join(feature_meta %>% select(nELISA_ID, Protein_Name), by = "nELISA_ID") %>%
    mutate(Protein_Name = dplyr::coalesce(Protein_Name, nELISA_ID))
  safe_write_csv(de_result, file.path(tbl_dir, "03_DE_responder_mets_vs_surgery_recurrence.csv"))
  message("DE: ", nrow(de_result), " proteins | nominal p<", PVAL_CUTOFF, ": ",
          sum(de_result$p_value < PVAL_CUTOFF, na.rm = TRUE),
          " | BH-adjusted p<", PVAL_CUTOFF, ": ",
          sum(de_result$adj_p_value < PVAL_CUTOFF, na.rm = TRUE))
}

# ── Volcano + MA (nominal p thresholding, matches 13_/14_ convention) ─
akila_volcano_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      panel.grid.major  = element_blank(), panel.grid.minor = element_blank(),
      axis.text         = element_text(color = "black", size = 12),
      axis.title        = element_text(color = "black", size = 14, face = "bold"),
      legend.position   = "bottom", legend.direction = "horizontal",
      legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
      legend.title      = element_text(face = "bold", size = 10),
      legend.text       = element_text(size = 9),
      plot.margin       = margin(68, 56, 46, 56)
    )
}

make_volcano <- function(df, title, path_stub, lfc_thr = LFC_THR, p_thr = PVAL_CUTOFF,
                         label_col = "Protein_Name", group_a_label = NULL, group_b_label = NULL) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  group_a_label <- group_a_label %||% (unique(na.omit(df$group_a))[1] %||% "Group A")
  group_b_label <- group_b_label %||% (unique(na.omit(df$group_b))[1] %||% "Group B")
  secondary_lbl <- paste0("Secondary (nominal p<=", formatC(p_thr, format = "f", digits = 2), ")")
  left_label    <- paste("Higher in", group_a_label)
  right_label   <- paste("Higher in", group_b_label)
  plot_df <- df %>%
    mutate(
      neg_log10 = -log10(pmax(p_value, 1e-300)),
      status = case_when(
        !is.na(p_value) & p_value < p_thr & logFC >=  lfc_thr ~ right_label,
        !is.na(p_value) & p_value < p_thr & logFC <= -lfc_thr ~ left_label,
        TRUE ~ "Not significant"
      ),
      sig_level = case_when(
        status == "Not significant" ~ "Not significant",
        !is.na(p_value) & p_value <= 0.05 ~ "Primary (nominal p<=0.05)",
        TRUE ~ secondary_lbl
      )
    )
  ann <- plot_df %>% filter(status != "Not significant") %>%
    arrange(p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  if (nrow(ann) == 0) ann <- plot_df %>% arrange(p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  x_rng <- range(plot_df$logFC, na.rm = TRUE)
  y_thr <- -log10(p_thr)
  y_max <- max(plot_df$neg_log10, y_thr, na.rm = TRUE)
  x_pad <- max(0.35, diff(x_rng) * 0.22); y_pad <- max(0.18, y_max * 0.16)
  x_panel_left <- x_rng[1] - 0.45 * x_pad; x_panel_right <- x_rng[2] + 0.45 * x_pad
  y_panel_top  <- y_max + 0.12 * y_pad;    y_label <- y_panel_top + 0.9 * y_pad
  p <- ggplot(plot_df, aes(x = logFC, y = neg_log10)) +
    geom_point(aes(color = status, shape = sig_level), size = 2.4, alpha = 0.82, stroke = 0.7) +
    geom_vline(xintercept = c(-lfc_thr, 0, lfc_thr),
               linetype = c("dashed", "solid", "dashed"), linewidth = c(0.6, 0.4, 0.6), color = "black") +
    geom_hline(yintercept = y_thr, linetype = "dashed", linewidth = 0.6, color = "black") +
    annotate("text", x = x_panel_left - 0.18 * x_pad, y = y_label,
             label = paste0("Higher in\n", group_a_label), hjust = 0, vjust = 0,
             size = 4.8, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_panel_right + 0.18 * x_pad, y = y_label,
             label = paste0("Higher in\n", group_b_label), hjust = 1, vjust = 0,
             size = 4.8, fontface = "bold", color = "#B2182B") +
    annotate("text", x = x_rng[1] + 0.02 * diff(x_rng), y = y_thr,
             label = paste0(" Nominal P = ", formatC(p_thr, format = "f", digits = 2)),
             hjust = 0, vjust = -0.6, size = 3.8, fontface = "italic") +
    scale_color_manual(values = c(setNames("#2166AC", left_label), setNames("#B2182B", right_label),
                                 "Not significant" = "grey85"), drop = FALSE) +
    scale_shape_manual(values = c("Not significant" = 1, "Primary (nominal p<=0.05)" = 19,
                                 setNames(17, secondary_lbl)),
                       breaks = c("Not significant", "Primary (nominal p<=0.05)", secondary_lbl), drop = FALSE) +
    ggrepel::geom_text_repel(data = ann, aes(label = .data[[label_col]]), size = 3.1, fontface = "bold",
                             max.overlaps = 30, box.padding = 0.5, point.padding = 0.2,
                             segment.color = "grey60", segment.size = 0.2, show.legend = FALSE) +
    labs(title = title, x = paste0("Fold change log2(", group_b_label, " / ", group_a_label, ")"),
         y = "-log10(nominal p value)", color = NULL, shape = "Significance") +
    guides(color = "none", shape = guide_legend(override.aes = list(color = "black"))) +
    coord_cartesian(xlim = c(x_panel_left, x_panel_right), ylim = c(0, y_panel_top), clip = "off") +
    akila_volcano_theme()
  save_ggplot(p, path_stub, width = 8.6, height = 8.2)
}

make_ma_plot <- function(df, title, path_stub, lfc_thr = LFC_THR, p_thr = PVAL_CUTOFF,
                         label_col = "Protein_Name", group_a_label = NULL, group_b_label = NULL) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  group_a_label <- group_a_label %||% (unique(na.omit(df$group_a))[1] %||% "Group A")
  group_b_label <- group_b_label %||% (unique(na.omit(df$group_b))[1] %||% "Group B")
  left_label <- paste("Higher in", group_a_label); right_label <- paste("Higher in", group_b_label)
  plot_df <- df %>%
    mutate(status = case_when(
      !is.na(p_value) & p_value < p_thr & logFC >=  lfc_thr ~ right_label,
      !is.na(p_value) & p_value < p_thr & logFC <= -lfc_thr ~ left_label,
      TRUE ~ "Not significant"
    ))
  ann <- plot_df %>% filter(status != "Not significant") %>%
    arrange(p_value, desc(abs(logFC))) %>% slice_head(n = 18)
  p <- ggplot(plot_df, aes(x = avg_expr, y = logFC)) +
    geom_point(aes(color = status), size = 2.2, alpha = 0.8, stroke = 0.6) +
    geom_hline(yintercept = c(-lfc_thr, 0, lfc_thr),
               linetype = c("dashed", "solid", "dashed"), linewidth = c(0.6, 0.4, 0.6), color = "black") +
    scale_color_manual(values = c(setNames("#2166AC", left_label), setNames("#B2182B", right_label),
                                 "Not significant" = "grey80"), drop = FALSE) +
    ggrepel::geom_text_repel(data = ann, aes(label = .data[[label_col]]), size = 3.1, fontface = "bold",
                             max.overlaps = 30, box.padding = 0.5, point.padding = 0.2,
                             segment.color = "grey60", segment.size = 0.2, show.legend = FALSE) +
    labs(title = title, x = "Average log2 expression (A)",
         y = paste0("log2 fold change (M): ", group_b_label, " / ", group_a_label), color = NULL) +
    akila_volcano_theme()
  save_ggplot(p, path_stub, width = 8.6, height = 8.2)
}

group_a_label <- "Responder (existing mets)"
group_b_label <- "Surgery (recurred/new mets)"

if (!is.null(de_result)) {
  make_volcano(de_result, "Responder-existing-mets vs Surgery-recurrence — Volcano",
              file.path(fig_dir, "responder_mets_vs_surgery_recurrence_volcano"),
              group_a_label = group_a_label, group_b_label = group_b_label)
  make_ma_plot(de_result, "Responder-existing-mets vs Surgery-recurrence — MA plot",
              file.path(fig_dir, "responder_mets_vs_surgery_recurrence_ma"),
              group_a_label = group_a_label, group_b_label = group_b_label)
}

# ============================================================
# Boxplot facets for nominally significant proteins (both directions)
# ============================================================
message("\n=== Boxplot facets for nominally significant proteins ===")

get_significant_features <- function(result_df, lfc_thr = LFC_THR, p_thr = PVAL_CUTOFF) {
  if (is.null(result_df) || nrow(result_df) == 0)
    return(list(higher_a = tibble(), higher_b = tibble()))
  list(
    higher_a = result_df %>% filter(!is.na(p_value), p_value < p_thr, logFC <= -lfc_thr) %>% arrange(p_value),
    higher_b = result_df %>% filter(!is.na(p_value), p_value < p_thr, logFC >=  lfc_thr) %>% arrange(p_value)
  )
}

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
    left_join(sig_tbl %>% select(nELISA_ID, Protein_Name, logFC, p_value, adj_p_value), by = "nELISA_ID") %>%
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

GROUP_COLORS <- setNames(c("#2166AC", "#B2182B"),
                         c("Responder_existing_mets", "Surgery_recurred_new_mets"))

if (!is.null(de_result)) {
  sig <- get_significant_features(de_result)
  for (side in c("higher_a", "higher_b")) {
    sig_tbl    <- sig[[side]]
    side_label <- if (side == "higher_a") group_a_label else group_b_label
    stub_side  <- if (side == "higher_a") "higher_in_responder_mets" else "higher_in_surgery_recurrence"
    if (nrow(sig_tbl) == 0) {
      message("Higher in ", side_label, ": no nominally significant proteins — skipped")
      next
    }
    safe_write_csv(sig_tbl %>% select(nELISA_ID, Protein_Name, logFC, p_value, adj_p_value),
                  file.path(tbl_dir, paste0("04_", stub_side, "_protein_list.csv")))
    long_df <- build_long_pgml(expr_pgml_rep, comparison_df, sig_tbl, "comparison_group")
    make_facet_boxplot(
      long_df, group_col = "comparison_group",
      group_levels = c("Responder_existing_mets", "Surgery_recurred_new_mets"),
      group_colors = GROUP_COLORS,
      title = paste0("Higher in ", side_label, " (n=", nrow(sig_tbl), " proteins)"),
      path_stub = file.path(fig_dir, paste0(stub_side, "_boxplots"))
    )
    message("Higher in ", side_label, ": ", nrow(sig_tbl), " proteins plotted")
  }
}

message("\n15_responder_mets_vs_surgery_recurrence.R complete")
