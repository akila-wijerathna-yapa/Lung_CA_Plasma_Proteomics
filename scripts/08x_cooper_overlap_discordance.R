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

# ── Local constants / paths ────────────────────────────────────
OUT_TBL  <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_DIR  <- file.path(OUT_TBL,   "cooper_comparison")
FIG_DIR  <- file.path(paths$figures, "07_analysis6_machine_learning", "cooper_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# UPDATE this path to wherever cooper_selected_proteins_structured.csv was saved
COOPER_CSV <- file.path(PROJECT_DIR, "data", "cooper_selected_proteins_structured.csv")

MIN_SELECTION_RATE <- 0.20   # inclusion threshold for Cooper's CV-selected proteins
CV_GAP_FLAG         <- 0.15  # method-level insample-vs-CV gap considered "large"

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# ── Step 1: load Cooper's structured selection data ────────────
if (!file.exists(COOPER_CSV)) {
  stop("Cooper's structured CSV not found at: ", COOPER_CSV,
       " — update COOPER_CSV path at top of script.")
}
cooper <- readr::read_csv(COOPER_CSV, show_col_types = FALSE) %>%
  dplyr::filter(feature_type == "protein", selection_rate >= MIN_SELECTION_RATE) %>%
  dplyr::mutate(uniprot_clean = sub("-\\d+$", "", uniprot_id))

message("Cooper proteins retained (selection_rate >= ", MIN_SELECTION_RATE, "): ", nrow(cooper))

# Map Cooper's endpoint/model onto Akila's three endpoints
cooper_by_endpoint <- split(cooper, cooper$endpoint)

# ── Step 2: load Akila's own selected-protein sets ──────────────
read_if_exists <- function(path) {
  if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else NULL
}

STABL_DIR <- file.path(OUT_TBL, "stabl")
akila_os_stabl  <- read_if_exists(file.path(STABL_DIR, "stabl_14_cox_OS_stabl_selected_features.csv"))
akila_pfs_stabl <- read_if_exists(file.path(STABL_DIR, "stabl_15_cox_PFS_stabl_selected_features.csv"))
akila_consensus <- read_if_exists(file.path(OUT_TBL, "consensus_gene_list_rebuilt.csv"))  # response endpoint

if (is.null(akila_os_stabl) || is.null(akila_pfs_stabl) || is.null(akila_consensus)) {
  stop("One or more Akila selection files missing from ", OUT_TBL,
       " — check filenames match 08b/08o outputs.")
}

akila_by_endpoint <- list(
  OS       = feature_meta %>%
    dplyr::filter(nELISA_ID %in% akila_os_stabl$nELISA_ID) %>%
    dplyr::mutate(uniprot_clean = sub("-\\d+$", "", Uniprot_ID),
                   akila_method_flag = "STABL (freq >= tau*)"),
  PFS      = feature_meta %>%
    dplyr::filter(nELISA_ID %in% akila_pfs_stabl$nELISA_ID) %>%
    dplyr::mutate(uniprot_clean = sub("-\\d+$", "", Uniprot_ID),
                   akila_method_flag = "STABL (freq >= tau*)"),
  Response = feature_meta %>%
    dplyr::filter(nELISA_ID %in% akila_consensus$nELISA_ID) %>%
    dplyr::mutate(uniprot_clean = sub("-\\d+$", "", Uniprot_ID),
                   akila_method_flag = paste0(">=2/4 methods, >=2 comparisons"))
)

# ── Step 3: per-endpoint overlap ─────────────────────────────────
build_overlap <- function(endpoint) {
  ak <- akila_by_endpoint[[endpoint]]
  cp <- cooper_by_endpoint[[endpoint]]
  if (is.null(cp)) cp <- cooper[0, ]

  ak_ids <- unique(ak$uniprot_clean)
  cp_ids <- unique(cp$uniprot_clean)
  all_ids <- union(ak_ids, cp_ids)

  tibble::tibble(uniprot_clean = all_ids) %>%
    dplyr::left_join(ak %>% dplyr::select(uniprot_clean, Gene_Symbol, akila_method_flag) %>%
                        dplyr::distinct(uniprot_clean, .keep_all = TRUE),
                      by = "uniprot_clean") %>%
    dplyr::left_join(cp %>% dplyr::group_by(uniprot_clean) %>%
                        dplyr::summarise(cooper_gene_label = dplyr::first(gene_label_cooper),
                                          cooper_models = paste(unique(model), collapse = "; "),
                                          cooper_max_selection_rate = max(as.numeric(selection_rate)),
                                          .groups = "drop"),
                      by = "uniprot_clean") %>%
    dplyr::mutate(
      endpoint       = endpoint,
      gene_symbol    = dplyr::coalesce(Gene_Symbol, cooper_gene_label),
      in_akila       = !is.na(akila_method_flag),
      in_cooper      = !is.na(cooper_models),
      overlap_status = dplyr::case_when(
        in_akila & in_cooper   ~ "both",
        in_akila & !in_cooper  ~ "akila_only",
        !in_akila & in_cooper  ~ "cooper_only",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(endpoint, gene_symbol, uniprot_clean, overlap_status,
                  akila_method_flag, cooper_models, cooper_max_selection_rate)
}

overlap_all <- purrr::map_dfr(c("OS", "PFS", "Response"), build_overlap) %>%
  dplyr::arrange(endpoint, dplyr::desc(overlap_status == "both"), gene_symbol)

safe_write_csv(overlap_all, file.path(OUT_DIR, "cooper_akila_overlap_by_endpoint.csv"))
message("Overlap table written: ", nrow(overlap_all), " rows")
print(table(overlap_all$endpoint, overlap_all$overlap_status))

# ── Step 4: Venn figure per endpoint ─────────────────────────────
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  message("Installing ggVennDiagram ...")
  install.packages("ggVennDiagram")
}

make_venn <- function(endpoint) {
  df <- overlap_all %>% dplyr::filter(endpoint == !!endpoint)
  set_list <- list(
    Akila  = df$gene_symbol[df$overlap_status %in% c("both", "akila_only")],
    Cooper = df$gene_symbol[df$overlap_status %in% c("both", "cooper_only")]
  )
  p <- ggVennDiagram::ggVennDiagram(set_list, label = "count") +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2166AC") +
    ggplot2::labs(title = paste0("Protein selection overlap — ", endpoint)) +
    ggplot2::theme(legend.position = "none")

  cairo_pdf(file.path(FIG_DIR, paste0("cooper_akila_venn_", endpoint, ".pdf")),
            width = 6, height = 5); print(p); dev.off()
  png(file.path(FIG_DIR, paste0("cooper_akila_venn_", endpoint, ".png")),
      width = 6, height = 5, units = "in", res = 300, bg = "white")
  print(p); dev.off()
  invisible(p)
}
invisible(lapply(c("OS", "PFS", "Response"), make_venn))
message("Venn figures written to ", FIG_DIR)

# ── Step 5: discordance classification ───────────────────────────
# CV robustness proxy: Cooper's own per-protein selection_rate (fold stability)
# joined against Akila's method-level insample-vs-CV gap for the corresponding
# endpoint (method-level, not per-protein — flagged explicitly in output).
cv_gap_binom <- read_if_exists(file.path(OUT_TBL, "insample_vs_cv_binomial_summary.csv"))
cv_gap_cox   <- read_if_exists(file.path(OUT_TBL, "insample_vs_cv_cox_summary.csv"))

endpoint_gap <- dplyr::bind_rows(
  if (!is.null(cv_gap_cox))
    cv_gap_cox %>% dplyr::group_by(endpoint) %>%
      dplyr::summarise(max_method_overfit_gap = max(overfit_gap, na.rm = TRUE), .groups = "drop"),
  if (!is.null(cv_gap_binom))
    tibble::tibble(endpoint = "Response",
                    max_method_overfit_gap = max(cv_gap_binom$overfit_gap, na.rm = TRUE))
)

discordance <- overlap_all %>%
  dplyr::left_join(endpoint_gap, by = "endpoint") %>%
  dplyr::mutate(
    akila_meets_consensus = in_akila_only_col <- overlap_status %in% c("both", "akila_only"),
    classification = dplyr::case_when(
      overlap_status == "both" ~ "confirmed",
      overlap_status %in% c("akila_only", "cooper_only") &
        (dplyr::coalesce(cooper_max_selection_rate, 1) < 0.4 |
           dplyr::coalesce(max_method_overfit_gap, 0) >= CV_GAP_FLAG) ~ "discordant_informative",
      TRUE ~ "hypothesis_generating"
    )
  ) %>%
  dplyr::select(endpoint, gene_symbol, uniprot_clean, overlap_status, classification,
                akila_method_flag, cooper_models, cooper_max_selection_rate,
                max_method_overfit_gap)

safe_write_csv(discordance, file.path(OUT_DIR, "cooper_akila_discordance_classification.csv"))
message("Discordance table written: ", nrow(discordance), " rows")
print(table(discordance$endpoint, discordance$classification))
