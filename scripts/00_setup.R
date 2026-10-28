# ── Auto-detect PROJECT_DIR ──────────────────────────────────
# Skipped if PROJECT_DIR already exists in the calling environment.
if (!exists("PROJECT_DIR")) {
  .this_file <- tryCatch({
    # Works when called via source("path/to/00_setup.R")
    frame_files <- lapply(sys.frames(), function(f) f$ofile)
    frame_files <- Filter(Negate(is.null), frame_files)
    if (length(frame_files))
      normalizePath(frame_files[[length(frame_files)]], mustWork = FALSE)
    else
      NA_character_
  }, error = function(e) NA_character_)

  # Fallback: active document in Positron / RStudio
  if (is.na(.this_file) || !nzchar(.this_file)) {
    .this_file <- tryCatch(
      normalizePath(rstudioapi::getActiveDocumentContext()$path, mustWork = FALSE),
      error = function(e) NA_character_
    )
  }

  if (!is.na(.this_file) && nzchar(.this_file) && file.exists(.this_file)) {
    PROJECT_DIR <- normalizePath(file.path(dirname(.this_file), ".."), mustWork = FALSE)
  } else {
    stop(
      "Cannot auto-detect PROJECT_DIR.\n",
      "Set it manually then re-source:\n",
      "  PROJECT_DIR <- 'your/project/path'\n",
      "  source('scripts/00_setup.R')"
    )
  }
  rm(.this_file)
}

message("PROJECT_DIR: ", PROJECT_DIR)

options(
  stringsAsFactors = FALSE,
  scipen            = 999,
  repos             = c(CRAN = "https://cloud.r-project.org"),
  timeout           = max(1200, getOption("timeout"))
)

# ── Package installer ────────────────────────────────────────
ensure_package <- function(pkg, source = c("cran", "bioc"), required = TRUE) {
  source <- match.arg(source)
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  message("Installing: ", pkg, " [", source, "]")
  ok <- tryCatch({
    if (source == "cran") {
      install.packages(pkg, quiet = TRUE)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", quiet = TRUE)
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
    }
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) {
    message("Failed to install '", pkg, "': ", conditionMessage(e))
    FALSE
  })
  if (!ok && required) stop("Required package '", pkg, "' is not available.")
  invisible(ok)
}

# ── Package lists ────────────────────────────────────────────
cran_required <- c(
  "readxl", "dplyr", "tidyr", "tibble", "stringr", "readr",
  "ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "scales",
  "lubridate", "glmnet", "pROC", "msigdbr",
  "survival", "survminer", "icenReg"
)
cran_optional  <- c("igraph")
bioc_required  <- c("limma", "fgsea", "org.Hs.eg.db", "AnnotationDbi")
bioc_optional  <- c("STRINGdb")

for (pkg in cran_required) ensure_package(pkg, "cran", required = TRUE)
for (pkg in cran_optional) ensure_package(pkg, "cran", required = FALSE)
for (pkg in bioc_required) ensure_package(pkg, "bioc", required = TRUE)
for (pkg in bioc_optional) ensure_package(pkg, "bioc", required = FALSE)

# ── Load libraries ───────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(lubridate)
  library(glmnet)
  library(pROC)
  library(msigdbr)
  library(survival)
  library(survminer)
  library(limma)
  library(fgsea)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

if (requireNamespace("igraph",   quietly = TRUE))
  suppressPackageStartupMessages(library(igraph))
if (requireNamespace("STRINGdb", quietly = TRUE))
  suppressPackageStartupMessages(library(STRINGdb))

# ── Namespace guardrails ─────────────────────────────────────
select  <- dplyr::select
rename  <- dplyr::rename
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
count   <- dplyr::count

# ── Shared path constant ─────────────────────────────────────
ANALYSIS8_DIR <- "08_pathway_and_network_analysis"

# ── Paths list ───────────────────────────────────────────────
paths <- list(
  # Input files — all data files live in data/ subdirectory
  clinical   = file.path(PROJECT_DIR, "data", "Clinical_Akila_FINAL_cleaned_20260610_Original.xlsx"),
  proteomics = file.path(PROJECT_DIR, "data", "Protein_pg_ml_Akila_202603271608_Original.xlsx"),
  sample_map = file.path(PROJECT_DIR, "data", "sample_IDs_before_after.xlsx"),
  tables     = file.path(PROJECT_DIR, "data_tables"),
  figures    = file.path(PROJECT_DIR, "figures"),
  reports    = file.path(PROJECT_DIR, "reports"),
  rdata      = file.path(PROJECT_DIR, "Analysis_Output", "rdata"),
  logs       = file.path(PROJECT_DIR, "logs")
)

# ── Create all output directories ───────────────────────────
output_dirs <- c(
  file.path(paths$tables, "00_cohort_audit"),
  file.path(paths$tables, "01_qc"),
  file.path(paths$tables, "02_analysis1_seven_group"),
  file.path(paths$tables, "03_analysis2_before_after"),
  file.path(paths$tables, "04_analysis3_pr_vs_pd"),
  file.path(paths$tables, "05_analysis4_tki_vs_io"),
  file.path(paths$tables, "06_analysis5_survival"),
  file.path(paths$tables, "07_analysis6_machine_learning"),
  file.path(paths$tables, ANALYSIS8_DIR),
  file.path(paths$tables, "09_poster_comparisons"),
  file.path(paths$figures, "00_overview"),
  file.path(paths$figures, "01_qc"),
  file.path(paths$figures, "02_analysis1_seven_group"),
  file.path(paths$figures, "03_analysis2_before_after"),
  file.path(paths$figures, "04_analysis3_pr_vs_pd"),
  file.path(paths$figures, "05_analysis4_tki_vs_io"),
  file.path(paths$figures, "06_analysis5_survival"),
  file.path(paths$figures, "07_analysis6_machine_learning"),
  file.path(paths$figures, ANALYSIS8_DIR),
  file.path(paths$figures, "09_poster_comparisons"),
  paths$reports,
  paths$logs,
  paths$rdata
)

invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

message("00_setup.R complete")
