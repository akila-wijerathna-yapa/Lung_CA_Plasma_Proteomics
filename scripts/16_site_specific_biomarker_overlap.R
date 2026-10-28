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

if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
library(UpSetR)

# ── Load checkpoint ────────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")
while (dev.cur() > 1) dev.off()

# ── Load curated recurrence-sites file ─────────────────────────
recurrence_sites_path <- file.path(PROJECT_DIR, "data", "1L_Sites_Recurrence_curated.xlsx")
if (!file.exists(recurrence_sites_path))
  stop("1L_Sites_Recurrence_curated.xlsx not found at: ", recurrence_sites_path)

norm_id_local <- function(x) {
  out <- as.character(x); out <- trimws(out)
  out <- gsub("\\.0$", "", out)
  out <- gsub("[^A-Za-z0-9]+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  toupper(out)
}

SITE_COLS <- c("brain", "bone", "liver", "lung", "lymph_node", "pleural", "adrenal", "other_site")
recurrence_sites_raw <- readxl::read_excel(recurrence_sites_path)
recurrence_sites <- recurrence_sites_raw %>%
  transmute(
    patient_id = norm_id_local(Patient_ID),
    across(all_of(SITE_COLS), ~ as.integer(coalesce(as.numeric(.x), 0)))
  )
message("Recurrence-sites file: ", nrow(recurrence_sites), " patients")

# ── Output directories ─────────────────────────────────────────
tbl_dir <- file.path(paths$tables,  "14_site_specific_biomarkers")
fig_dir <- file.path(paths$figures, "14_site_specific_biomarkers")
invisible(lapply(c(tbl_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Local constants ────────────────────────────────────────────
PVAL_CUTOFF   <- 0.1   # nominal p-value threshold (NOT FDR-adjusted)
LFC_THR       <- 0.5
MIN_PER_GROUP <- 4     # sites with fewer positive cases are skipped with a warning
COHORT_FILTER <- NULL  # e.g. "Sx" to restrict to surgery patients; NULL = whole recurrence cohort

# ── Build analysis cohort: recurrence-sites patients with baseline PRE sample ─
cohort_df <- pre_matched %>%
  inner_join(recurrence_sites, by = "patient_id")
if (!is.null(COHORT_FILTER)) cohort_df <- cohort_df %>% filter(new_group == COHORT_FILTER)
message("Analysis cohort (matched baseline proteomics): ", nrow(cohort_df), " patients")
for (s in SITE_COLS) message("  ", s, ": n_positive = ", sum(cohort_df[[s]] == 1, na.rm = TRUE),
                             " / n_total = ", nrow(cohort_df))

# ── Local utility functions ────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")
save_ggplot <- function(p, path_stub, width, height) {
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height)
  print(p); dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

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
  limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("nELISA_ID") %>%
    as_tibble() %>%
    transmute(nELISA_ID, logFC, t_stat = t, p_value = P.Value,
              adj_p_value = adj.P.Val, avg_expr = AveExpr)
}

# ============================================================
# Per-site DE: site-positive vs site-negative (within cohort_df)
# ============================================================
message("\n=== Per-site DE ===")
site_de_results <- list()
site_hits        <- list()   # nELISA_ID vectors, higher-in-site-positive
site_hits_labeled <- list()  # Protein_Name vectors

for (s in SITE_COLS) {
  n_pos <- sum(cohort_df[[s]] == 1, na.rm = TRUE)
  n_neg <- sum(cohort_df[[s]] == 0, na.rm = TRUE)
  if (n_pos < MIN_PER_GROUP || n_neg < MIN_PER_GROUP) {
    message(s, ": skipped (n_positive=", n_pos, ", n_negative=", n_neg,
            ", need >= ", MIN_PER_GROUP, " each)")
    next
  }
  site_df <- cohort_df %>%
    mutate(site_status = factor(ifelse(.data[[s]] == 1, "Yes", "No"), levels = c("No", "Yes")))
  res <- run_binary_limma(expr_log2_rep, site_df, "site_status", c("No", "Yes"),
                          adjust_cols = c("plate_barcode"), min_per_group = MIN_PER_GROUP)
  if (is.null(res)) { message(s, ": DE failed/insufficient — skipped"); next }
  res <- res %>%
    left_join(feature_meta %>% select(nELISA_ID, Protein_Name), by = "nELISA_ID") %>%
    mutate(Protein_Name = dplyr::coalesce(Protein_Name, nELISA_ID), site = s)
  site_de_results[[s]] <- res
  safe_write_csv(res, file.path(tbl_dir, paste0("site_DE_", s, ".csv")))
  hits <- res %>% filter(!is.na(p_value), p_value < PVAL_CUTOFF, logFC >= LFC_THR) %>% arrange(p_value)
  site_hits[[s]]         <- hits$nELISA_ID
  site_hits_labeled[[s]] <- hits$Protein_Name
  message(s, ": n_pos=", n_pos, " n_neg=", n_neg, " | proteins higher-in-site (nominal p<",
          PVAL_CUTOFF, ", |logFC|>=", LFC_THR, "): ", nrow(hits))
}

site_hits <- site_hits[lengths(site_hits) > 0]
if (length(site_hits) < 2) {
  message("Fewer than 2 sites produced hits — overlap analysis skipped.")
} else {

# ============================================================
# Membership matrix + summary
# ============================================================
message("\n=== Building protein x site membership matrix ===")
all_ids <- unique(unlist(site_hits))
id_to_name <- bind_rows(site_de_results) %>% distinct(nELISA_ID, Protein_Name)

membership_df <- tibble(nELISA_ID = all_ids) %>%
  left_join(id_to_name, by = "nELISA_ID")
for (s in names(site_hits)) {
  membership_df[[s]] <- as.integer(membership_df$nELISA_ID %in% site_hits[[s]])
}
membership_df <- membership_df %>%
  mutate(n_sites = rowSums(across(all_of(names(site_hits))))) %>%
  arrange(desc(n_sites), Protein_Name)
safe_write_csv(membership_df, file.path(tbl_dir, "protein_site_membership_matrix.csv"))

shared_summary <- membership_df %>% count(n_sites, name = "n_proteins") %>% arrange(desc(n_sites))
safe_write_csv(shared_summary, file.path(tbl_dir, "shared_vs_unique_summary.csv"))
message("Shared/unique breakdown (n_sites -> n_proteins):")
print(shared_summary)

# ============================================================
# (1) UpSet plot
# ============================================================
message("\n=== UpSet plot ===")
upset_input <- as.data.frame(membership_df %>% select(all_of(names(site_hits))))
rownames(upset_input) <- membership_df$Protein_Name

pdf(file.path(fig_dir, "site_biomarker_upset.pdf"), width = 9, height = 6, onefile = FALSE)
print(UpSetR::upset(upset_input, sets = names(site_hits), order.by = "freq",
                    main.bar.color = "#2166AC", sets.bar.color = "#B2182B",
                    text.scale = 1.3, keep.order = FALSE))
dev.off()
png(file.path(fig_dir, "site_biomarker_upset.png"), width = 9, height = 6, units = "in", res = 300, bg = "white")
print(UpSetR::upset(upset_input, sets = names(site_hits), order.by = "freq",
                    main.bar.color = "#2166AC", sets.bar.color = "#B2182B",
                    text.scale = 1.3, keep.order = FALSE))
dev.off()

# ============================================================
# (2) Binary membership heatmap (protein x site), clustered
# ============================================================
message("\n=== Membership heatmap ===")
heat_mat <- as.matrix(membership_df %>% select(all_of(names(site_hits))))
rownames(heat_mat) <- membership_df$Protein_Name

ph <- pheatmap::pheatmap(
  heat_mat, color = c("grey95", "#B2182B"), legend_breaks = c(0, 1),
  cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = nrow(heat_mat) <= 60,
  fontsize_row = 6, fontsize_col = 10, main = "Site-associated protein membership",
  silent = TRUE
)
cairo_pdf(file.path(fig_dir, "site_biomarker_membership_heatmap.pdf"),
         width = 7, height = max(5, nrow(heat_mat) * 0.12 + 2))
grid::grid.newpage(); grid::grid.draw(ph$gtable)
dev.off()
png(file.path(fig_dir, "site_biomarker_membership_heatmap.png"),
   width = 7, height = max(5, nrow(heat_mat) * 0.12 + 2), units = "in", res = 300, bg = "white")
grid::grid.newpage(); grid::grid.draw(ph$gtable)
dev.off()

# ============================================================
# (3) Site x site Jaccard similarity heatmap
# ============================================================
message("\n=== Site x site Jaccard similarity ===")
site_names <- names(site_hits)
jaccard_mat <- matrix(NA_real_, length(site_names), length(site_names),
                     dimnames = list(site_names, site_names))
for (i in site_names) for (j in site_names) {
  a <- site_hits[[i]]; b <- site_hits[[j]]
  jaccard_mat[i, j] <- if (length(union(a, b)) == 0) 0 else length(intersect(a, b)) / length(union(a, b))
}
safe_write_csv(as_tibble(jaccard_mat, rownames = "site"),
              file.path(tbl_dir, "site_pairwise_jaccard.csv"))

ph_jac <- pheatmap::pheatmap(
  jaccard_mat, display_numbers = TRUE, number_format = "%.2f",
  color = colorRampPalette(c("white", "#2166AC"))(50),
  main = "Site pairwise Jaccard similarity (higher-in-site protein sets)",
  silent = TRUE
)
cairo_pdf(file.path(fig_dir, "site_pairwise_jaccard_heatmap.pdf"), width = 6, height = 5.5)
grid::grid.newpage(); grid::grid.draw(ph_jac$gtable)
dev.off()
png(file.path(fig_dir, "site_pairwise_jaccard_heatmap.png"), width = 6, height = 5.5,
   units = "in", res = 300, bg = "white")
grid::grid.newpage(); grid::grid.draw(ph_jac$gtable)
dev.off()

# ── Save checkpoint for downstream scripts (e.g. bin protein tables/boxplots) ─
message("\n=== Saving checkpoint_16.RData ===")
save(cohort_df, site_de_results, site_hits, site_hits_labeled, membership_df, id_to_name,
    file = file.path(paths$rdata, "checkpoint_16.RData"))

}

message("\n16_site_specific_biomarker_overlap.R complete")
