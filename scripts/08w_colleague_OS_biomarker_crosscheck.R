sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

library(dplyr)
library(readr)
library(readxl)
library(stringr)

# ---- 1. Colleague's candidate list ------------------------------------------
colleague_path <- file.path(PROJECT_DIR, "data", "OS_feature_Proteins_CJ.xlsx")
stopifnot(file.exists(colleague_path))

colleague <- read_excel(colleague_path, sheet = 1) %>%
  rename(Colleague_Protein_Name = 1, Uniprot_ID_raw = 2) %>%
  mutate(Uniprot_ID_clean = str_trim(str_remove(Uniprot_ID_raw, "-\\d+$")))

# ---- 2. Panel feature metadata (nELISA_ID <-> Uniprot_ID <-> Protein_Label) -
feat_path <- file.path(paths$tables, "00_cohort_audit", "proteomics_feature_annotation.csv")
stopifnot(file.exists(feat_path))
feature_meta <- read_csv(feat_path, show_col_types = FALSE) %>%
  mutate(Uniprot_ID_clean = str_trim(str_remove(Uniprot_ID, "-\\d+$")))

# ---- 3. Match colleague proteins to nELISA_ID -------------------------------
matched <- colleague %>%
  left_join(
    feature_meta %>% select(nELISA_ID, Uniprot_ID, Gene_Symbol, Protein_Label, Uniprot_ID_clean),
    by = "Uniprot_ID_clean"
  )

message("Colleague proteins matched to panel: ",
        sum(!is.na(matched$nELISA_ID)), " / ", nrow(matched))
if (any(is.na(matched$nELISA_ID))) {
  message("Unmatched (check UniProt ID formatting): ",
          paste(matched$Colleague_Protein_Name[is.na(matched$nELISA_ID)], collapse = ", "))
}

# ---- 4. OS univariate Cox results -------------------------------------------
os_path <- file.path(paths$tables, "06_analysis5_survival", "analysis5_os_univariate_cox.csv")
os_res  <- if (file.exists(os_path)) read_csv(os_path, show_col_types = FALSE) else {
  message("MISSING: ", os_path, " — run 07_analysis5_survival.R first.")
  tibble()
}

# ---- 5. ML cross-method selection export ------------------------------------
ml_path <- file.path(paths$tables, "07_analysis6_machine_learning", "biomarker_selection_full.csv")
if (!file.exists(ml_path)) {
  message("MISSING: ", ml_path)
  message("-> 00_biomarker_selection_full_export.R has not been run, or its PROJECT_DIR/OUTPUT_CSV ",
          "does not point at this project's paths$tables. Check that script's hardcoded ",
          "PROJECT_DIR line matches: ", PROJECT_DIR)
}
ml_sel <- if (file.exists(ml_path)) read_csv(ml_path, show_col_types = FALSE) else tibble()

ml_summary <- if (nrow(ml_sel)) {
  ml_sel %>%
    group_by(nELISA_ID) %>%
    summarise(
      n_ML_selections  = n(),
      ML_selected_where = paste(unique(paste0(comparison, ":", method)), collapse = "; "),
      .groups = "drop"
    )
} else tibble(nELISA_ID = character(), n_ML_selections = integer(), ML_selected_where = character())

# ---- 6. Merge everything -----------------------------------------------------
crosscheck <- matched %>%
  left_join(
    if (nrow(os_res)) os_res %>% select(nELISA_ID, HR, HR_low, HR_high, p_value, adj_p_value, n, n_events)
    else tibble(nELISA_ID = character()),
    by = "nELISA_ID"
  ) %>%
  left_join(ml_summary, by = "nELISA_ID") %>%
  mutate(
    in_final_QC_passed_set = nELISA_ID %in% (if (nrow(os_res)) os_res$nELISA_ID else character()),
    n_ML_selections = tidyr::replace_na(n_ML_selections, 0L)
  ) %>%
  arrange(adj_p_value)

dir.create(file.path(paths$tables, "06_analysis5_survival"), recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(paths$tables, "06_analysis5_survival", "colleague_OS_candidates_crosscheck.csv")
write_csv(crosscheck, out_path)

message("\nCrosscheck written: ", out_path)
message("Present in final QC-passed analytic set: ",
        sum(crosscheck$in_final_QC_passed_set, na.rm = TRUE), " / ", nrow(crosscheck))
message("Nominal OS Cox p<0.05: ", sum(crosscheck$p_value < 0.05, na.rm = TRUE))
message("Flagged by >=1 ML method in >=1 comparison: ", sum(crosscheck$n_ML_selections > 0))
print(crosscheck, n = nrow(crosscheck))
