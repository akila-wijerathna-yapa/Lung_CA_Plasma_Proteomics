## ---- 0. Setup ---------------------------------------------------------- ##

options(warn = -1)
op <- options()
options(op)

PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"

paths <- list(
  clinical = file.path(PROJECT_DIR, "data", "Clinical_Akila_FINAL_cleaned_20260610.xlsx"),
  protein  = file.path(PROJECT_DIR, "data", "Protein_pg_ml_Akila_202603271608.xlsx"),
  out_tab  = file.path(PROJECT_DIR, "data_tables", "09d_stage_sensitivity"),
  out_fig  = file.path(PROJECT_DIR, "figures", "09d_stage_sensitivity")
)
dir.create(paths$out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$out_fig, recursive = TRUE, showWarnings = FALSE)

need <- c("readxl", "dplyr", "survival", "ggplot2")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(readxl); library(dplyr); library(survival); library(ggplot2)

RESTRICT_HISTOLOGY <- TRUE

## ---- 1. Load and clean clinical data ------------------------------------ ##

sheet_names <- readxl::excel_sheets(paths$clinical)
target_sheet <- sheet_names[grepl("clinical|clean", sheet_names, ignore.case = TRUE)][1]
if (is.na(target_sheet)) target_sheet <- sheet_names[1]
clinical_raw <- read_excel(paths$clinical, sheet = target_sheet)

clinical <- clinical_raw %>%
  transmute(
    patient_id = as.integer(Patient_ID),
    sex_female = as.integer(Sex == "F"),
    histology  = Histological_Subtype,
    stage      = Stage_Group,
    tx_raw     = `1L_Treatment_Class`,
    dob        = as.Date(DOB, format = "%d/%m/%Y"),
    tx_start   = as.Date(`1L_Start_Date`, format = "%d/%m/%Y"),
    response_binary = as.integer(Response_Binary == "Responder"),
    pfs_months = as.numeric(PFS_months),
    pfs_event  = as.integer(PFS_Event),
    os_months  = as.numeric(OS_months),
    os_event   = as.integer(OS_Event_1Dead)
  ) %>%
  mutate(
    age = as.numeric(difftime(tx_start, dob, units = "days")) / 365.25,
    tx_chemo   = as.integer(grepl("Chemo",     tx_raw, ignore.case = TRUE)),
    tx_io      = as.integer(grepl("IO",         tx_raw, ignore.case = TRUE)),
    tx_rt      = as.integer(grepl("RT",         tx_raw, ignore.case = TRUE)),
    tx_surgery = as.integer(grepl("Surgery",    tx_raw, ignore.case = TRUE)),
    tx_tki     = as.integer(grepl("TKI",        tx_raw, ignore.case = TRUE)),
    hist_adeno = as.integer(histology == "Adenocarcinoma"),
    stage_collapsed = as.integer(stage %in% c("III", "IV"))   # 0 = I/II, 1 = III/IV
  )

if (RESTRICT_HISTOLOGY) {
  clinical <- clinical %>% filter(histology %in% c("Adenocarcinoma", "Squamous Cell Carcinoma"))
}

CLIN_COVARS       <- c("age", "sex_female", "hist_adeno",
                        "tx_chemo", "tx_io", "tx_rt", "tx_surgery", "tx_tki")
CLIN_COVARS_STAGE <- c(CLIN_COVARS, "stage_collapsed")

cat("CHECKPOINT 1: clinical cohort ->", nrow(clinical), "patients\n")
cat("Missingness, stage-adjusted covariate set:\n")
print(colSums(is.na(clinical[, CLIN_COVARS_STAGE])))

## ---- 2. Load and clean protein data -------------------------------------- ##

read_protein_workbook <- function(path) {
  raw <- read_excel(path, sheet = 1, col_names = FALSE)
  protein_name_row <- 3; ulod_row <- 4; llod_row <- 5
  header_row <- 5; data_start_row <- 6

  protein_cols <- which(!is.na(unlist(raw[protein_name_row, 12:ncol(raw)]))) + 11
  protein_names <- as.character(unlist(raw[protein_name_row, protein_cols], use.names = FALSE))
  llod <- as.numeric(unlist(raw[llod_row, protein_cols]))
  ulod <- as.numeric(unlist(raw[ulod_row, protein_cols]))

  if (anyDuplicated(protein_names)) {
    cat("WARNING:", sum(duplicated(protein_names)),
        "duplicate protein names -- disambiguating with make.unique()\n")
    protein_names <- make.unique(protein_names)
  }

  meta <- raw[data_start_row:nrow(raw), 1:10]
  names(meta) <- unlist(raw[header_row, 1:10])

  vals <- as.data.frame(raw[data_start_row:nrow(raw), protein_cols], check.names = FALSE)
  names(vals) <- protein_names

  list(meta = meta, vals = vals, protein_names = protein_names, llod = llod, ulod = ulod)
}

pw <- read_protein_workbook(paths$protein)

keep <- which(pw$meta$sample_passed_qc == 1 & pw$meta$pre_post == "pre")
meta_pre <- pw$meta[keep, ]
vals_pre <- pw$vals[keep, ]

impute_lod <- function(x, llod_i, ulod_i) {
  x_chr <- as.character(x)
  is_below <- grepl("^<LLOD$", x_chr)
  is_above <- grepl("^>ULOD$", x_chr)
  is_bad   <- grepl("^(NaN|Low bead count)$", x_chr)
  numeric_mask <- grepl("^-?[0-9.]+([eE][-+]?[0-9]+)?$", x_chr)
  out <- rep(NA_real_, length(x))
  out[numeric_mask] <- as.numeric(x_chr[numeric_mask])
  out[is_below] <- 0.5 * llod_i
  out[is_above] <- 2 * ulod_i
  out[is_bad] <- NA_real_
  out
}

expr_pg <- mapply(
  function(col, i) impute_lod(vals_pre[[col]], pw$llod[i], pw$ulod[i]),
  seq_along(pw$protein_names), seq_along(pw$protein_names),
  SIMPLIFY = FALSE
)
expr_pg <- as.data.frame(expr_pg)
names(expr_pg) <- pw$protein_names
expr_pg$patient_id <- as.integer(meta_pre$id_number)

geo_mean <- function(x) exp(mean(log(x[x > 0 & !is.na(x)])))
expr_pg_combined <- expr_pg %>%
  group_by(patient_id) %>%
  summarise(across(all_of(pw$protein_names), geo_mean), .groups = "drop")

protein_cols_actual <- setdiff(names(expr_pg_combined), "patient_id")
expr_log10 <- expr_pg_combined %>%
  mutate(across(all_of(protein_cols_actual), ~ log10(.x)))

cat("CHECKPOINT 2: protein matrix ->", nrow(expr_log10), "patients x",
    length(protein_cols_actual), "proteins\n")

## ---- 3. Merge clinical + protein ------------------------------------------ ##

full_df <- clinical %>% inner_join(expr_log10, by = "patient_id")

protein_names <- protein_cols_actual
miss_frac <- sapply(protein_names, function(p) mean(is.na(full_df[[p]])))
protein_names <- protein_names[miss_frac <= 0.20]
for (p in protein_names) {
  full_df[[p]][is.na(full_df[[p]])] <- median(full_df[[p]], na.rm = TRUE)
}

cat("CHECKPOINT 3: merged cohort ->", nrow(full_df), "patients x",
    length(protein_names), "proteins retained\n")


## ---- 4. Candidate protein check ------------------------------------------- ##

CANDIDATE_PROTEINS <- c(
  GDF15   = "GDF-15 (MIC-1)",
  GZMA    = "GZMA",
  NOTUM   = "NOTUM",
  PLAUR   = "uPA R",
  SHH     = "SHH",
  CEACAM5 = "CEACAM5"
)
matched_idx <- CANDIDATE_PROTEINS %in% protein_names
matched     <- CANDIDATE_PROTEINS[matched_idx]
unmatched   <- CANDIDATE_PROTEINS[!matched_idx]
cat("CHECKPOINT 4: matched candidate proteins ->",
    paste(names(matched), "=", matched, collapse = "; "), "\n")
if (length(unmatched))
  cat("  UNMATCHED (check naming against pw$protein_names) ->",
      paste(names(unmatched), "=", unmatched, collapse = "; "), "\n")
stopifnot(length(matched) > 0)

## ---- 5. Univariate, clinical-adjusted Cox (with / without stage) --------- ##

univariate_adjusted_cox <- function(df, time_col, event_col, protein_names, clin_covars) {
  df <- df[!is.na(df[[time_col]]) & !is.na(df[[event_col]]), ]
  res <- lapply(protein_names, function(p) {
    covar_str <- paste(clin_covars, collapse = " + ")
    f <- as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ `",
                            p, "` + ", covar_str))
    fit <- tryCatch(survival::coxph(f, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    if (!is.null(fit$iter) && fit$iter[1] >= 20) return(NULL)
    if (any(abs(coef(fit)) > 10, na.rm = TRUE)) return(NULL)
    s <- summary(fit)$coefficients
    p_term <- paste0("`", p, "`")
    if (!(p_term %in% rownames(s))) {
      if (p %in% rownames(s)) p_term <- p else return(NULL)
    }
    row <- s[p_term, ]
    data.frame(protein = p, HR = exp(row["coef"]),
               lower = exp(row["coef"] - 1.96 * row["se(coef)"]),
               upper = exp(row["coef"] + 1.96 * row["se(coef)"]),
               p = row["Pr(>|z|)"], n = fit$n, n_events = fit$nevent)
  })
  out <- do.call(rbind, res)
  if (is.null(out) || nrow(out) == 0) return(out)
  out$fdr <- p.adjust(out$p, method = "BH")
  out[order(out$fdr), ]
}

cox_os_base  <- univariate_adjusted_cox(full_df, "os_months", "os_event", unname(matched), CLIN_COVARS)
cox_os_stage <- univariate_adjusted_cox(full_df, "os_months", "os_event", unname(matched), CLIN_COVARS_STAGE)

# relabel panel names back to gene symbols for reporting (single application, post-computation)
label_map <- setNames(names(matched), matched)
cox_os_base$protein  <- label_map[cox_os_base$protein]
cox_os_stage$protein <- label_map[cox_os_stage$protein]

cox_os_base$model  <- "no_stage"
cox_os_stage$model <- "stage_adjusted"

comparison <- rbind(cox_os_base, cox_os_stage) %>%
  arrange(protein, model)

write.csv(comparison, file.path(paths$out_tab, "os_stage_sensitivity_comparison.csv"), row.names = FALSE)
print(comparison)

## ---- 6. Forest plot: stage-adjusted panel only ---------------------------- ##

if (nrow(cox_os_stage) > 0) {
  plot_df <- cox_os_stage
  plot_df$protein <- factor(plot_df$protein, levels = rev(plot_df$protein))
  n_sig <- sum(plot_df$fdr < 0.1, na.rm = TRUE)

  p <- ggplot(plot_df, aes(x = HR, y = protein)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(xmin = lower, xmax = upper)) +
    labs(title = "Univariate OS, clinical + stage-adjusted (candidate panel)",
         subtitle = paste0(n_sig, " of ", nrow(plot_df), " protein(s) reach FDR < 0.1"),
         x = "Hazard ratio (95% CI), clinical + stage-adjusted", y = NULL) +
    theme_bw()

  cairo_pdf(file.path(paths$out_fig, "forest_OS_stage_adjusted_candidates.pdf"),
            width = 7, height = max(4, nrow(plot_df) * 0.3) + 0.5)
  print(p)
  dev.off()
  png(file.path(paths$out_fig, "forest_OS_stage_adjusted_candidates.png"),
      width = 7, height = max(4, nrow(plot_df) * 0.3) + 0.5, units = "in", res = 300, bg = "white")
  print(p)
  dev.off()
}

cat("CHECKPOINT 6: stage sensitivity check complete. Outputs in\n  ",
    paths$out_tab, "\n  ", paths$out_fig, "\n")

