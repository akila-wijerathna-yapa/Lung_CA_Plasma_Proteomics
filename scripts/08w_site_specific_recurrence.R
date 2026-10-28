sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")

# ── Local constants ────────────────────────────────────────────
CURATED_SITE_FILE <- file.path(PROJECT_DIR, "data", "1L_Sites_Recurrence_curated.xlsx")
SITE_COLS   <- c("brain", "bone", "liver", "lung", "lymph_node", "pleural", "adrenal", "other_site")
MIN_GROUP_N <- 6   # minimum n per arm (site-positive / no-recurrence) to attempt testing

OUT_DIR <- file.path(paths$tables, "08_site_specific_recurrence")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# ── Load curated site flags and merge onto analysis1_df ────────
site_raw <- readxl::read_excel(CURATED_SITE_FILE, sheet = 1)
site_raw <- site_raw %>%
  dplyr::mutate(patient_id = norm_id(Patient_ID)) %>%
  dplyr::select(patient_id, dplyr::all_of(SITE_COLS))

analysis1_sites <- analysis1_df %>%
  dplyr::left_join(site_raw, by = "patient_id")

for (col in SITE_COLS) {
  analysis1_sites[[col]][is.na(analysis1_sites[[col]])] <- 0L
}

# ── Define comparison groups ────────────────────────────────────
# Control arm: confirmed no recurrence at all (drop patients whose
# recurrence exists but at a site not curated / not in this table).
no_recurrence_ids <- analysis1_sites %>%
  dplyr::filter(recurrence_binary_sensitivity == "No recurrence") %>%
  dplyr::pull(sample_id)

# ── Local Wilcoxon-per-protein helper (self-contained) ──────────
run_wilcoxon_site <- function(expr_mat, site_ids, control_ids, feature_meta_df) {
  ids_use <- c(site_ids, control_ids)
  X <- expr_mat[ids_use, , drop = FALSE]
  grp <- factor(ifelse(ids_use %in% site_ids, "Site", "NoRecurrence"),
                levels = c("NoRecurrence", "Site"))

  res <- lapply(colnames(X), function(feat) {
    x <- X[, feat]
    keep <- !is.na(x)
    if (sum(keep & grp == "Site") < 3 || sum(keep & grp == "NoRecurrence") < 3) {
      return(tibble::tibble(nELISA_ID = feat, p_value = NA_real_,
                             median_site = NA_real_, median_norecurrence = NA_real_,
                             log2_median_diff = NA_real_,
                             n_site = sum(keep & grp == "Site"),
                             n_norecurrence = sum(keep & grp == "NoRecurrence")))
    }
    wt <- suppressWarnings(wilcox.test(x[keep] ~ grp[keep]))
    med_site <- median(x[keep & grp == "Site"], na.rm = TRUE)
    med_no   <- median(x[keep & grp == "NoRecurrence"], na.rm = TRUE)
    tibble::tibble(
      nELISA_ID = feat, p_value = wt$p.value,
      median_site = med_site, median_norecurrence = med_no,
      log2_median_diff = med_site - med_no,
      n_site = sum(keep & grp == "Site"),
      n_norecurrence = sum(keep & grp == "NoRecurrence")
    )
  })

  out <- dplyr::bind_rows(res) %>%
    dplyr::mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(p_value)

  if (exists("add_feature_annotation")) out <- add_feature_annotation(out, feature_meta_df)
  out
}

# ── Group-size audit (written regardless of what gets tested) ──
group_counts <- lapply(SITE_COLS, function(col) {
  site_ids <- analysis1_sites %>% dplyr::filter(.data[[col]] == 1) %>% dplyr::pull(sample_id)
  tibble::tibble(site = col, n_site = length(site_ids), n_no_recurrence = length(no_recurrence_ids))
}) %>% dplyr::bind_rows()
safe_write_csv(group_counts, file.path(OUT_DIR, "site_recurrence_group_counts.csv"))
message("\nGroup sizes:")
print(group_counts)

# ── Run per-site Wilcoxon tests ──────────────────────────────────
message("\n=== 08w | Site-specific recurrence protein tests ===")
site_results <- list()
for (col in SITE_COLS) {
  site_ids <- analysis1_sites %>% dplyr::filter(.data[[col]] == 1) %>% dplyr::pull(sample_id)
  site_ids <- intersect(site_ids, rownames(expr_log2_rep))
  ctrl_ids <- intersect(no_recurrence_ids, rownames(expr_log2_rep))

  if (length(site_ids) < MIN_GROUP_N || length(ctrl_ids) < MIN_GROUP_N) {
    message("  SKIP ", col, ": n_site=", length(site_ids), " n_control=", length(ctrl_ids),
            " (below MIN_GROUP_N=", MIN_GROUP_N, ")")
    next
  }

  message("  Testing ", col, ": n_site=", length(site_ids), " n_control=", length(ctrl_ids))
  res <- run_wilcoxon_site(expr_log2_rep, site_ids, ctrl_ids, feature_meta)
  site_results[[col]] <- res
  safe_write_csv(res, file.path(OUT_DIR, paste0(col, "_wilcoxon.csv")))
  message("    Uncorrected p<0.05: ", sum(res$p_value < 0.05, na.rm = TRUE),
          " | BH q<0.10: ", sum(res$adj_p_value < 0.10, na.rm = TRUE))
}

message("\n08w complete. Outputs written to: ", OUT_DIR)
message("NOTE: 'other_site' group includes several ambiguous/likely-non-progression entries")
message("(second primary cancer, explicit treatment-response text) - interpret with caution.")
