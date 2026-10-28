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
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded (feature_meta)")

ALTML_DIR <- file.path(paths$tables, "07_analysis6_machine_learning", "altml")
STABL_DIR <- file.path(paths$tables, "07_analysis6_machine_learning", "stabl")
OUT_DIR   <- file.path(paths$tables, "07_analysis6_machine_learning", "cooper_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# UPDATE this path to wherever cooper_selected_proteins_structured.csv was saved
COOPER_CSV <- file.path(PROJECT_DIR, "data", "cooper_selected_proteins_structured.csv")

MIN_SELECTION_RATE          <- 0.20
MIN_METHODS_PER_COMPARISON  <- 2
MIN_COMPARISONS             <- 2

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# ── Step 1: Cooper's per-protein selected sets by endpoint ─────
if (!file.exists(COOPER_CSV)) stop("Cooper CSV not found at: ", COOPER_CSV)

cooper <- readr::read_csv(COOPER_CSV, show_col_types = FALSE) %>%
  dplyr::filter(feature_type == "protein", selection_rate >= MIN_SELECTION_RATE) %>%
  dplyr::mutate(uniprot_clean = sub("-\\d+$", "", uniprot_id))

cooper_set <- function(ep) unique(cooper$uniprot_clean[cooper$endpoint == ep])
cooper_os_ids       <- cooper_set("OS")
cooper_pfs_ids      <- cooper_set("PFS")
cooper_response_ids <- cooper_set("Response")

message("Cooper protein sets — OS: ", length(cooper_os_ids),
        " | PFS: ", length(cooper_pfs_ids),
        " | Response: ", length(cooper_response_ids))

# ══════════════════════════════════════════════════════════════
# Response: recompute true 4-method + Cooper 5-method consensus
# ══════════════════════════════════════════════════════════════
RESPONSE_FILES <- c(
  "altml_01_CR_vs_PD_cross_method_overlap.csv",
  "altml_02_CRPR_vs_PD_cross_method_overlap.csv",
  "altml_03_CRPR_vs_SDPD_cross_method_overlap.csv",
  "altml_04_CRPRSD_vs_PD_cross_method_overlap.csv",
  "altml_05_PR_vs_PD_cross_method_overlap.csv"
)

resp_long <- purrr::map_dfr(RESPONSE_FILES, function(f) {
  readr::read_csv(file.path(ALTML_DIR, f), show_col_types = FALSE) %>%
    dplyr::transmute(
      nELISA_ID, Protein_Label,
      comparison = sub("_cross_method_overlap\\.csv$", "", f),
      n_methods_true_4 = as.integer(en_selected) + as.integer(splsda_selected) +
                          as.integer(boruta_confirmed) + as.integer(stabl_selected)
    )
}) %>%
  dplyr::left_join(feature_meta %>% dplyr::select(nELISA_ID, Uniprot_ID), by = "nELISA_ID") %>%
  dplyr::mutate(
    uniprot_clean    = sub("-\\d+$", "", Uniprot_ID),
    cooper_selected  = as.integer(uniprot_clean %in% cooper_response_ids),
    n_methods_true_5 = n_methods_true_4 + cooper_selected
  )

build_consensus <- function(df, n_col, min_methods, min_comparisons) {
  df %>%
    dplyr::filter(.data[[n_col]] >= min_methods) %>%
    dplyr::group_by(nELISA_ID, Protein_Label) %>%
    dplyr::summarise(
      n_comparisons  = dplyr::n(),
      comparisons    = paste(sort(unique(comparison)), collapse = "; "),
      max_n_methods  = max(.data[[n_col]]),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_comparisons >= min_comparisons) %>%
    dplyr::arrange(dplyr::desc(n_comparisons), dplyr::desc(max_n_methods))
}

response_4 <- build_consensus(resp_long, "n_methods_true_4", MIN_METHODS_PER_COMPARISON, MIN_COMPARISONS)
response_5 <- build_consensus(resp_long, "n_methods_true_5", MIN_METHODS_PER_COMPARISON, MIN_COMPARISONS)

message("Response consensus >=2/4 (Akila only): ", nrow(response_4), " proteins — ",
        paste(response_4$Protein_Label, collapse = ", "))
message("Response consensus >=2/5 (+Cooper):     ", nrow(response_5), " proteins — ",
        paste(response_5$Protein_Label, collapse = ", "))

safe_write_csv(response_4, file.path(OUT_DIR, "consensus_response_2of4_akila_only.csv"))
safe_write_csv(response_5, file.path(OUT_DIR, "consensus_response_2of5_with_cooper.csv"))

# ══════════════════════════════════════════════════════════════
# OS / PFS: EN-Cox + STABL (2 flags) -> + Cooper (3 flags)
# ══════════════════════════════════════════════════════════════
build_endpoint <- function(endpoint, altml_file, stabl_file, cooper_ids) {
  altml <- readr::read_csv(file.path(ALTML_DIR, altml_file), show_col_types = FALSE)
  stabl <- readr::read_csv(file.path(STABL_DIR, stabl_file), show_col_types = FALSE)
  stabl_ids <- stabl$nELISA_ID[stabl$stabl_selected]

  altml %>%
    dplyr::transmute(
      nELISA_ID, Gene_Symbol, Uniprot_Clean,
      en_selected      = as.logical(en_selected),
      stabl_selected   = nELISA_ID %in% stabl_ids,
      n_methods_true_2 = as.integer(en_selected) + as.integer(stabl_selected),
      cooper_selected  = as.integer(Uniprot_Clean %in% cooper_ids),
      n_methods_true_3 = n_methods_true_2 + cooper_selected
    ) %>%
    dplyr::filter(n_methods_true_3 > 0) %>%
    dplyr::arrange(dplyr::desc(n_methods_true_3)) %>%
    dplyr::mutate(endpoint = endpoint, .before = 1)
}

os_full  <- build_endpoint("OS",  "altml_14_cox_OS_all_methods.csv",  "stabl_14_cox_OS_stabl_selected_features.csv",  cooper_os_ids)
pfs_full <- build_endpoint("PFS", "altml_15_cox_PFS_all_methods.csv", "stabl_15_cox_PFS_stabl_selected_features.csv", cooper_pfs_ids)

message("OS  EN-Cox selected count: ",  sum(os_full$en_selected),
        " (0 means EN never co-selects with STABL — check before trusting the 2/2 rule)")
message("PFS EN-Cox selected count: ", sum(pfs_full$en_selected))

os_2  <- os_full  %>% dplyr::filter(n_methods_true_2 >= 2)
os_3  <- os_full  %>% dplyr::filter(n_methods_true_3 >= 2)
pfs_2 <- pfs_full %>% dplyr::filter(n_methods_true_2 >= 2)
pfs_3 <- pfs_full %>% dplyr::filter(n_methods_true_3 >= 2)

message("OS  >=2/2 (Akila only): ", paste(os_2$Gene_Symbol, collapse = ", "))
message("OS  >=2/3 (+Cooper):    ", paste(os_3$Gene_Symbol, collapse = ", "))
message("PFS >=2/2 (Akila only): ", paste(pfs_2$Gene_Symbol, collapse = ", "))
message("PFS >=2/3 (+Cooper):    ", paste(pfs_3$Gene_Symbol, collapse = ", "))

safe_write_csv(os_full,  file.path(OUT_DIR, "consensus_OS_full_with_cooper.csv"))
safe_write_csv(pfs_full, file.path(OUT_DIR, "consensus_PFS_full_with_cooper.csv"))

# ── Cross-endpoint flip summary ─────────────────────────────────
flip_row <- function(gene, endpoint, rule_a, rule_b, in_a, in_b) {
  tibble::tibble(
    endpoint = endpoint, gene = gene,
    rule_akila_only = rule_a, rule_with_cooper = rule_b,
    in_akila_only_consensus = in_a, in_with_cooper_consensus = in_b,
    status = ifelse(in_a, "unchanged", "NEW via Cooper")
  )
}

resp_summary <- dplyr::bind_rows(
  purrr::map_dfr(response_5$Protein_Label, ~ flip_row(
    .x, "Response", "2/4", "2/5",
    .x %in% response_4$Protein_Label, TRUE)),
  purrr::map_dfr(setdiff(response_4$Protein_Label, response_5$Protein_Label), ~ flip_row(
    .x, "Response", "2/4", "2/5", TRUE, FALSE))
)

os_summary  <- purrr::map_dfr(os_3$Gene_Symbol,  ~ flip_row(.x, "OS",  "2/2 (EN+STABL)", "2/3 (+Cooper)", .x %in% os_2$Gene_Symbol,  TRUE))
pfs_summary <- purrr::map_dfr(pfs_3$Gene_Symbol, ~ flip_row(.x, "PFS", "2/2 (EN+STABL)", "2/3 (+Cooper)", .x %in% pfs_2$Gene_Symbol, TRUE))

summary_all <- dplyr::bind_rows(resp_summary, os_summary, pfs_summary)
safe_write_csv(summary_all, file.path(OUT_DIR, "consensus_4_vs_5_method_summary.csv"))
print(summary_all)

