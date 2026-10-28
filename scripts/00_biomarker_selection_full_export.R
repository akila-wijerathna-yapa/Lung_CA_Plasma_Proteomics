library(dplyr)
library(readr)
library(purrr)

# ---- EDIT FOR YOUR MACHINE --------------------------------------------------
PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
ML_DIR      <- file.path(PROJECT_DIR, "data_tables", "07_analysis6_machine_learning")
OUTPUT_CSV  <- file.path(PROJECT_DIR, "data_tables", "07_analysis6_machine_learning", "biomarker_selection_full.csv")
# -----------------------------------------------------------------------------

# ---- 1. EN / sPLS-DA / Boruta / STABL, from cross_method_overlap.csv --------
overlap_files <- list.files(file.path(ML_DIR, "altml"), pattern = "_cross_method_overlap\\.csv$", full.names = TRUE)

overlap_long <- map_dfr(overlap_files, function(f) {
  comparison <- sub("_cross_method_overlap\\.csv$", "", basename(f))
  df <- read_csv(f, show_col_types = FALSE)
  method_cols <- c(en_selected = "Elastic Net", splsda_selected = "sPLS-DA", boruta_confirmed = "Boruta")
  if (!"stabl_selected" %in% names(df)) df$stabl_selected <- FALSE
  method_cols <- c(method_cols, stabl_selected = "STABL")
  bind_rows(lapply(names(method_cols), function(col) {
    if (!col %in% names(df)) return(tibble())
    df %>%
      filter(.data[[col]] == TRUE) %>%
      transmute(comparison = comparison, method = method_cols[[col]],
                nELISA_ID, Protein_Label)
  }))
})

# ---- 2. LASSO, from analysis6a/6b lasso_stability.csv -----------------------
lasso_files <- list.files(file.path(ML_DIR, "data_tables_ml"), pattern = "_lasso_stability\\.csv$", full.names = TRUE)

lasso_long <- map_dfr(lasso_files, function(f) {
  comparison <- sub("_lasso_stability\\.csv$", "", basename(f))
  df <- read_csv(f, show_col_types = FALSE)
  df %>%
    filter(lasso_selected == TRUE) %>%
    transmute(comparison = comparison, method = "LASSO",
              nELISA_ID, Protein_Label = coalesce(Gene_Symbol, nELISA_ID),
              stability_fraction)
})

# ---- 3. STABL final selected features, from stabl_*_selected_features.csv --
stabl_files <- list.files(file.path(ML_DIR, "stabl"), pattern = "_stabl_selected_features\\.csv$", full.names = TRUE)

stabl_long <- map_dfr(stabl_files, function(f) {
  comparison <- sub("_stabl_selected_features\\.csv$", "", basename(f))
  df <- read_csv(f, show_col_types = FALSE)
  df %>%
    transmute(comparison = comparison, method = "STABL (final model)",
              nELISA_ID, Protein_Label = coalesce(Gene_Symbol, nELISA_ID))
})

# ---- combine and write -------------------------------------------------------
full_table <- bind_rows(overlap_long, lasso_long, stabl_long) %>%
  arrange(comparison, method, Protein_Label) %>%
  distinct(comparison, method, Protein_Label, .keep_all = TRUE)

dir.create(dirname(OUTPUT_CSV), recursive = TRUE, showWarnings = FALSE)
write_csv(full_table, OUTPUT_CSV)

message("Rows written: ", nrow(full_table))
message("Full biomarker selection table written to: ", OUTPUT_CSV)
