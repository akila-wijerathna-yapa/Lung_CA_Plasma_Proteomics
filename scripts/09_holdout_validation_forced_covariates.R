## ---- 0. Setup ---------------------------------------------------------- ##

options(warn = -1)  # numeric coercion warnings suppressed at option level,
op <- options()      # NOT suppressWarnings() (crashes in loops on R 4.5.2)
options(op)

PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"

paths <- list(
  clinical = file.path(PROJECT_DIR, "data", "Clinical_Akila_FINAL_cleaned_20260610.xlsx"),
  protein  = file.path(PROJECT_DIR, "data", "Protein_pg_ml_Akila_202603271608.xlsx"),
  out_tab  = file.path(PROJECT_DIR, "data_tables", "09_holdout_validation"),
  out_fig  = file.path(PROJECT_DIR, "figures", "09_holdout_validation")
)
dir.create(paths$out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$out_fig, recursive = TRUE, showWarnings = FALSE)

need <- c("readxl", "dplyr", "glmnet", "survival", "ggplot2")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(readxl); library(dplyr); library(glmnet); library(survival); library(ggplot2)

# Toggle: restrict to adenocarcinoma/squamous only, to mirror the other
# scientist's cohort definition. Set FALSE to use the full histology set.
RESTRICT_HISTOLOGY <- TRUE

## ---- 1. Load and clean clinical data ------------------------------------ ##

sheet_names <- readxl::excel_sheets(paths$clinical)
cat("Sheets found in clinical file:", paste(sheet_names, collapse = ", "), "\n")
target_sheet <- sheet_names[grepl("clinical|clean", sheet_names, ignore.case = TRUE)][1]
if (is.na(target_sheet)) target_sheet <- sheet_names[1]
cat("Using sheet:", target_sheet, "\n")
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
    hist_adeno = as.integer(histology == "Adenocarcinoma")
  )

if (RESTRICT_HISTOLOGY) {
  clinical <- clinical %>% filter(histology %in% c("Adenocarcinoma", "Squamous Cell Carcinoma"))
}

CLIN_COVARS <- c("age", "sex_female", "hist_adeno",
                  "tx_chemo", "tx_io", "tx_rt", "tx_surgery", "tx_tki")

cat("CHECKPOINT 1: clinical cohort after cleaning ->", nrow(clinical), "patients\n")
cat("Missingness in clinical covariates:\n")
print(colSums(is.na(clinical[, CLIN_COVARS])))

## ---- 2. Load and clean protein data -------------------------------------- ##

read_protein_workbook <- function(path) {
  raw <- read_excel(path, sheet = 1, col_names = FALSE)

  protein_name_row <- 3
  llod_row <- 5
  ulod_row <- 4
  header_row <- 5      # sample-metadata column names live in row 5, cols 1:10
  data_start_row <- 6

  protein_cols <- which(!is.na(unlist(raw[protein_name_row, 12:ncol(raw)]))) + 11
  protein_names <- as.character(unlist(raw[protein_name_row, protein_cols], use.names = FALSE))
  llod <- as.numeric(unlist(raw[llod_row, protein_cols]))
  ulod <- as.numeric(unlist(raw[ulod_row, protein_cols]))

  if (anyDuplicated(protein_names)) {
    cat("WARNING:", sum(duplicated(protein_names)),
        "duplicate protein names found -- disambiguating with make.unique()\n")
    protein_names <- make.unique(protein_names)
  }

  meta <- raw[data_start_row:nrow(raw), 1:10]
  names(meta) <- unlist(raw[header_row, 1:10])

  # plain data.frame + check.names=FALSE: tibble's `names<-` silently
  # repairs non-unique/non-syntactic names to positional "...N" placeholders;
  # base data.frame with check.names=FALSE does not.
  vals <- as.data.frame(raw[data_start_row:nrow(raw), protein_cols], check.names = FALSE)
  names(vals) <- protein_names

  list(meta = meta, vals = vals, protein_names = protein_names, llod = llod, ulod = ulod)
}

pw <- read_protein_workbook(paths$protein)

# QC filter: passed_qc == 1, pre-treatment (baseline) samples only
keep <- which(pw$meta$sample_passed_qc == 1 & pw$meta$pre_post == "pre")
meta_pre <- pw$meta[keep, ]
vals_pre <- pw$vals[keep, ]

# LOD imputation, applied column-wise (regex pre-filter before as.numeric,
# per R 4.5.2 crash-avoidance convention)
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

# Combine same-patient replicate samples via geometric mean
geo_mean <- function(x) exp(mean(log(x[x > 0 & !is.na(x)])))
expr_pg_combined <- expr_pg %>%
  group_by(patient_id) %>%
  summarise(across(all_of(pw$protein_names), geo_mean), .groups = "drop")

# Log10 scale
# NB: use the column names actually present in expr_pg_combined, not the
# original pw$protein_names vector -- unlist() on a tibble row can silently
# alter whitespace/encoding, causing all_of(pw$protein_names) to miss.
protein_cols_actual <- setdiff(names(expr_pg_combined), "patient_id")
stopifnot(length(protein_cols_actual) > 0)

expr_log10 <- expr_pg_combined %>%
  mutate(across(all_of(protein_cols_actual), ~ log10(.x)))

cat("CHECKPOINT 2: protein matrix ->", nrow(expr_log10), "patients x",
    length(protein_cols_actual), "proteins\n")
cat("Names dropped vs pw$protein_names (diagnostic only):",
    length(setdiff(pw$protein_names, protein_cols_actual)), "\n")

## ---- 3. Merge clinical + protein ------------------------------------------ ##

full_df <- clinical %>%
  inner_join(expr_log10, by = "patient_id")

protein_names <- protein_cols_actual
# drop proteins with excessive missingness in the merged cohort (>20%)
miss_frac <- sapply(protein_names, function(p) mean(is.na(full_df[[p]])))
protein_names <- protein_names[miss_frac <= 0.20]
# median-impute any remaining sparse missingness
for (p in protein_names) {
  full_df[[p]][is.na(full_df[[p]])] <- median(full_df[[p]], na.rm = TRUE)
}

cat("CHECKPOINT 3: merged cohort ->", nrow(full_df), "patients x",
    length(protein_names), "proteins retained after missingness filter\n")

## ---- 4. Forced-covariate design matrices --------------------------------- ##

build_xy <- function(df, outcome_type, protein_names, clin_covars) {
  df <- df[complete.cases(df[, clin_covars]), ]
  if (outcome_type == "response") {
    df <- df[!is.na(df$response_binary), ]
    y <- df$response_binary
  } else if (outcome_type == "os") {
    df <- df[!is.na(df$os_months) & !is.na(df$os_event), ]
    y <- survival::Surv(df$os_months, df$os_event)
  } else if (outcome_type == "pfs") {
    df <- df[!is.na(df$pfs_months) & !is.na(df$pfs_event), ]
    y <- survival::Surv(df$pfs_months, df$pfs_event)
  }
  x_clin <- as.matrix(df[, clin_covars])
  x_prot <- as.matrix(df[, protein_names])
  # drop zero-variance columns FIRST -- scale() divides by SD, and a
  # constant column (SD = 0) produces NaN, which glmnet then rejects as
  # missing data even though there was no NA before scaling.
  col_sd <- apply(x_prot, 2, sd, na.rm = TRUE)
  x_prot <- x_prot[, col_sd > 0, drop = FALSE]
  x_prot <- scale(x_prot)
  protein_names <- colnames(x_prot)
  list(x = cbind(x_clin, x_prot), y = y, df = df,
       clin_covars = clin_covars, protein_names = protein_names)
}

## ---- 5. Stratified split + forced-penalty model fit ---------------------- ##

stratified_split <- function(strata, holdout_frac = 0.25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx_holdout <- unlist(lapply(split(seq_along(strata), strata), function(idx) {
    sample(idx, size = max(1, round(length(idx) * holdout_frac)))
  }))
  list(train = setdiff(seq_along(strata), idx_holdout), holdout = idx_holdout)
}

# Three-way split (train / validation / test) -- see caveat printed at
# CHECKPOINT 4b: a fixed validation partition is statistically LESS
# efficient than cv.glmnet's internal k-fold CV at this sample size
# (single ~20-30 patient draw vs. fold-averaged use of all training
# data). Included for comparison / sensitivity-check purposes only.
stratified_split_3way <- function(strata, val_frac = 0.20, holdout_frac = 0.20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  by_stratum <- split(seq_along(strata), strata)
  idx_holdout <- unlist(lapply(by_stratum, function(idx)
    sample(idx, size = max(1, round(length(idx) * holdout_frac)))))
  remaining <- setdiff(seq_along(strata), idx_holdout)
  by_stratum_rem <- split(remaining, strata[remaining])
  idx_val <- unlist(lapply(by_stratum_rem, function(idx)
    sample(idx, size = max(1, round(length(idx) * val_frac)))))
  idx_train <- setdiff(remaining, idx_val)
  list(train = idx_train, val = idx_val, holdout = idx_holdout)
}

fit_forced_penalty <- function(xy, family, alpha = 1, train_idx) {
  pf <- c(rep(0, length(xy$clin_covars)), rep(1, length(xy$protein_names)))
  y_train <- if (family == "cox") xy$y[train_idx, ] else xy$y[train_idx]
  cv <- glmnet::cv.glmnet(xy$x[train_idx, ], y_train, family = family,
                           alpha = alpha, penalty.factor = pf, nfolds = 5)
  cv
}

# Explicit train -> validation -> test fitting: lambda chosen by
# evaluating a full glmnet path on the validation partition (instead
# of cv.glmnet's internal folds), final coefficients refit on
# train+validation combined at that lambda, then evaluated once on
# the untouched holdout.
fit_forced_penalty_3way <- function(xy, family, alpha, train_idx, val_idx) {
  pf <- c(rep(0, length(xy$clin_covars)), rep(1, length(xy$protein_names)))
  y_train <- if (family == "cox") xy$y[train_idx, ] else xy$y[train_idx]
  path <- glmnet::glmnet(xy$x[train_idx, ], y_train, family = family,
                          alpha = alpha, penalty.factor = pf)
  pred_val <- predict(path, newx = xy$x[val_idx, ], type = "response")
  if (family == "binomial") {
    y_val <- xy$y[val_idx]
    val_perf <- apply(pred_val, 2, function(p) {
      tryCatch(as.numeric(pROC::auc(pROC::roc(y_val, as.numeric(p), quiet = TRUE))),
               error = function(e) NA_real_)
    })
  } else {
    y_val <- xy$y[val_idx, ]
    val_perf <- apply(pred_val, 2, function(p) {
      tryCatch(survival::concordance(y_val ~ as.numeric(p), reverse = TRUE)$concordance,
                error = function(e) NA_real_)
    })
  }
  best <- which.max(val_perf)
  lambda_star <- path$lambda[best]
  trainval_idx <- c(train_idx, val_idx)
  y_trainval <- if (family == "cox") xy$y[trainval_idx, ] else xy$y[trainval_idx]
  final <- glmnet::glmnet(xy$x[trainval_idx, ], y_trainval, family = family,
                           alpha = alpha, penalty.factor = pf, lambda = lambda_star)
  list(model = final, lambda_star = lambda_star, val_perf_at_star = val_perf[best])
}

evaluate_holdout <- function(cv, xy, family, holdout_idx) {
  pred <- predict(cv, newx = xy$x[holdout_idx, ], s = "lambda.min", type = "response")
  if (family == "binomial") {
    y_holdout <- xy$y[holdout_idx]
    auc <- pROC::auc(pROC::roc(y_holdout, as.numeric(pred), quiet = TRUE))
    perf <- as.numeric(auc)
  } else if (family == "cox") {
    y_holdout <- xy$y[holdout_idx, ]
    cindex <- survival::concordance(y_holdout ~ pred, reverse = TRUE)$concordance
    perf <- as.numeric(cindex)
  }
  perf
}

selected_proteins <- function(cv, xy) {
  coefs <- as.matrix(coef(cv, s = "lambda.min"))
  prot_coefs <- coefs[xy$protein_names, , drop = FALSE]
  rownames(prot_coefs)[prot_coefs[, 1] != 0]
}

## ---- 6. Single true holdout run ------------------------------------------- ##

true_holdout_eval <- function(df, outcome_type, protein_names, clin_covars,
                               family, alpha = 1, holdout_frac = 0.25, seed = 1) {
  xy <- build_xy(df, outcome_type, protein_names, clin_covars)
  strata <- if (family == "binomial") xy$y else xy$y[, "status"]
  split_idx <- stratified_split(strata, holdout_frac, seed = seed)
  cv <- fit_forced_penalty(xy, family, alpha, split_idx$train)
  perf <- evaluate_holdout(cv, xy, family, split_idx$holdout)
  sel <- selected_proteins(cv, xy)
  list(perf = perf, n_selected = length(sel), selected = sel,
       n_train = length(split_idx$train), n_holdout = length(split_idx$holdout))
}

res_response <- true_holdout_eval(full_df, "response", protein_names, CLIN_COVARS,
                                   family = "binomial", alpha = 1, seed = 1)
res_os <- true_holdout_eval(full_df, "os", protein_names, CLIN_COVARS,
                             family = "cox", alpha = 0.9, seed = 1)
res_pfs <- true_holdout_eval(full_df, "pfs", protein_names, CLIN_COVARS,
                              family = "cox", alpha = 0.9, seed = 1)

## ---- 6b. Optional: explicit three-way train/validation/test run ---------- ##
# See caveat above fit_forced_penalty_3way: this is a sensitivity check
# against the nested-CV approach above, not a statistical improvement
# on it at this sample size. Toggle to compare.
RUN_THREE_WAY_SPLIT <- TRUE

true_holdout_eval_3way <- function(df, outcome_type, protein_names, clin_covars,
                                    family, alpha = 1, val_frac = 0.20,
                                    holdout_frac = 0.20, seed = 1) {
  xy <- build_xy(df, outcome_type, protein_names, clin_covars)
  strata <- if (family == "binomial") xy$y else xy$y[, "status"]
  split_idx <- stratified_split_3way(strata, val_frac, holdout_frac, seed = seed)
  fit <- fit_forced_penalty_3way(xy, family, alpha, split_idx$train, split_idx$val)
  pred_holdout <- predict(fit$model, newx = xy$x[split_idx$holdout, ], type = "response")
  if (family == "binomial") {
    y_holdout <- xy$y[split_idx$holdout]
    perf <- tryCatch(as.numeric(pROC::auc(pROC::roc(y_holdout, as.numeric(pred_holdout), quiet = TRUE))),
                      error = function(e) NA_real_)
  } else {
    y_holdout <- xy$y[split_idx$holdout, ]
    perf <- tryCatch(survival::concordance(y_holdout ~ as.numeric(pred_holdout), reverse = TRUE)$concordance,
                      error = function(e) NA_real_)
  }
  coefs <- as.matrix(coef(fit$model))
  prot_coefs <- coefs[xy$protein_names, , drop = FALSE]
  sel <- rownames(prot_coefs)[prot_coefs[, 1] != 0]
  list(perf = perf, n_selected = length(sel), selected = sel,
       lambda_star = fit$lambda_star, val_perf_at_star = fit$val_perf_at_star,
       n_train = length(split_idx$train), n_val = length(split_idx$val),
       n_holdout = length(split_idx$holdout))
}

if (RUN_THREE_WAY_SPLIT) {
  res_response_3w <- true_holdout_eval_3way(full_df, "response", protein_names, CLIN_COVARS,
                                             family = "binomial", alpha = 1, seed = 1)
  res_os_3w <- true_holdout_eval_3way(full_df, "os", protein_names, CLIN_COVARS,
                                       family = "cox", alpha = 0.9, seed = 1)
  res_pfs_3w <- true_holdout_eval_3way(full_df, "pfs", protein_names, CLIN_COVARS,
                                        family = "cox", alpha = 0.9, seed = 1)

  cat("\nCHECKPOINT 4b: three-way train/validation/test comparison\n")
  cat("Response  - 2-way holdout AUC:", round(res_response$perf, 3),
      "| 3-way holdout AUC:", round(res_response_3w$perf, 3),
      "(n:", res_response_3w$n_train, "/", res_response_3w$n_val, "/",
      res_response_3w$n_holdout, ")\n")
  cat("OS        - 2-way holdout C-index:", round(res_os$perf, 3),
      "| 3-way holdout C-index:", round(res_os_3w$perf, 3),
      "(n:", res_os_3w$n_train, "/", res_os_3w$n_val, "/", res_os_3w$n_holdout, ")\n")
  cat("PFS       - 2-way holdout C-index:", round(res_pfs$perf, 3),
      "| 3-way holdout C-index:", round(res_pfs_3w$perf, 3),
      "(n:", res_pfs_3w$n_train, "/", res_pfs_3w$n_val, "/", res_pfs_3w$n_holdout, ")\n")
  cat("If 3-way numbers are noticeably less stable across reruns with\n",
      "different seeds than the 2-way nested-CV numbers, that confirms\n",
      "the efficiency caveat above -- prefer the 2-way nested-CV result.\n")
}

cat("CHECKPOINT 4: single true-holdout run\n")
cat("Response  - holdout AUC:", round(res_response$perf, 3),
    "| n selected proteins:", res_response$n_selected,
    "| train/holdout n:", res_response$n_train, "/", res_response$n_holdout, "\n")
cat("OS        - holdout C-index:", round(res_os$perf, 3),
    "| n selected proteins:", res_os$n_selected,
    "| train/holdout n:", res_os$n_train, "/", res_os$n_holdout, "\n")
cat("PFS       - holdout C-index:", round(res_pfs$perf, 3),
    "| n selected proteins:", res_pfs$n_selected,
    "| train/holdout n:", res_pfs$n_train, "/", res_pfs$n_holdout, "\n")

## ---- 7. Repeated / bootstrapped holdout ----------------------------------- ##

repeated_holdout <- function(df, outcome_type, protein_names, clin_covars,
                              family, alpha = 1, holdout_frac = 0.25, n_reps = 100) {
  perf_vec <- numeric(n_reps)
  sel_list <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    out <- tryCatch(
      true_holdout_eval(df, outcome_type, protein_names, clin_covars,
                         family, alpha, holdout_frac, seed = r),
      error = function(e) NULL
    )
    if (is.null(out)) next
    perf_vec[r] <- out$perf
    sel_list[[r]] <- out$selected
  }
  all_sel <- unlist(sel_list)
  selection_rate <- sort(table(all_sel) / n_reps, decreasing = TRUE)
  list(perf = perf_vec, selection_rate = selection_rate)
}

rep_response <- repeated_holdout(full_df, "response", protein_names, CLIN_COVARS,
                                  family = "binomial", alpha = 1, n_reps = 100)
rep_os <- repeated_holdout(full_df, "os", protein_names, CLIN_COVARS,
                            family = "cox", alpha = 0.9, n_reps = 100)
rep_pfs <- repeated_holdout(full_df, "pfs", protein_names, CLIN_COVARS,
                             family = "cox", alpha = 0.9, n_reps = 100)

cat("CHECKPOINT 5: repeated holdout (100 reps) summary\n")
cat("Response  - mean AUC:", round(mean(rep_response$perf), 3),
    "sd:", round(sd(rep_response$perf), 3), "\n")
cat("OS        - mean C-index:", round(mean(rep_os$perf), 3),
    "sd:", round(sd(rep_os$perf), 3), "\n")
cat("PFS       - mean C-index:", round(mean(rep_pfs$perf), 3),
    "sd:", round(sd(rep_pfs$perf), 3), "\n")

write.csv(as.data.frame(rep_response$selection_rate),
          file.path(paths$out_tab, "response_protein_selection_rate.csv"), row.names = FALSE)
write.csv(as.data.frame(rep_os$selection_rate),
          file.path(paths$out_tab, "os_protein_selection_rate.csv"), row.names = FALSE)
write.csv(as.data.frame(rep_pfs$selection_rate),
          file.path(paths$out_tab, "pfs_protein_selection_rate.csv"), row.names = FALSE)

# Histogram of holdout performance across repeats, per outcome
perf_df <- data.frame(
  outcome = rep(c("Response (AUC)", "OS (C-index)", "PFS (C-index)"), each = 100),
  value = c(rep_response$perf, rep_os$perf, rep_pfs$perf)
)
p <- ggplot(perf_df, aes(x = value)) +
  geom_histogram(bins = 20, fill = "#4C72B0", colour = "white") +
  facet_wrap(~ outcome, scales = "free_x") +
  labs(title = "Holdout performance across 100 random 75/25 splits",
       x = "Holdout performance", y = "Count") +
  theme_bw()
cairo_pdf(file.path(paths$out_fig, "repeated_holdout_performance.pdf"), width = 10, height = 4)
print(p)
dev.off()

## ---- 8. Univariate, clinical-adjusted hazard forest plot ----------------- ##

univariate_adjusted_cox <- function(df, time_col, event_col, protein_names, clin_covars) {
  df <- df[!is.na(df[[time_col]]) & !is.na(df[[event_col]]), ]
  res <- lapply(protein_names, function(p) {
    covar_str <- paste(clin_covars, collapse = " + ")
    f <- as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ `",
                            p, "` + ", covar_str))
    fit <- tryCatch(survival::coxph(f, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    if (!is.null(fit$iter) && fit$iter[1] >= 20) return(NULL)
    if (any(abs(coef(fit)) > 10, na.rm = TRUE)) return(NULL)  # likely separation
    s <- summary(fit)$coefficients
    p_term <- paste0("`", p, "`")
    if (!(p_term %in% rownames(s))) {
      if (p %in% rownames(s)) p_term <- p else return(NULL)
    }
    row <- s[p_term, ]
    data.frame(protein = p, HR = exp(row["coef"]),
               lower = exp(row["coef"] - 1.96 * row["se(coef)"]),
               upper = exp(row["coef"] + 1.96 * row["se(coef)"]),
               p = row["Pr(>|z|)"])
  })
  out <- do.call(rbind, res)
  out$fdr <- p.adjust(out$p, method = "BH")
  out[order(out$fdr), ]
}

# Restrict to a candidate list here (e.g. STABL / consensus hits) rather than
# all 973 proteins, to avoid re-incurring the full multiplicity penalty.
# Edit CANDIDATE_PROTEINS to your consensus panel; defaults to full panel.
CANDIDATE_PROTEINS <- protein_names

cox_os <- univariate_adjusted_cox(full_df, "os_months", "os_event",
                                   CANDIDATE_PROTEINS, CLIN_COVARS)
cox_pfs <- univariate_adjusted_cox(full_df, "pfs_months", "pfs_event",
                                    CANDIDATE_PROTEINS, CLIN_COVARS)

write.csv(cox_os,  file.path(paths$out_tab, "univariate_adjusted_cox_OS.csv"), row.names = FALSE)
write.csv(cox_pfs, file.path(paths$out_tab, "univariate_adjusted_cox_PFS.csv"), row.names = FALSE)

sig_os  <- cox_os[cox_os$fdr < 0.1, ]
sig_pfs <- cox_pfs[cox_pfs$fdr < 0.1, ]

cat("CHECKPOINT 6: univariate clinical-adjusted Cox\n")
cat("OS  - significant at FDR<0.1:", nrow(sig_os), "of", nrow(cox_os), "tested\n")
cat("PFS - significant at FDR<0.1:", nrow(sig_pfs), "of", nrow(cox_pfs), "tested\n")

make_forest <- function(res_df, title, path_stub, n = 20, fdr_thresh = 0.1) {
  plot_df <- head(res_df[order(res_df$fdr), ], n)
  plot_df$protein <- factor(plot_df$protein, levels = rev(plot_df$protein))
  n_sig <- sum(res_df$fdr < fdr_thresh, na.rm = TRUE)
  subtitle_txt <- if (n_sig == 0) {
    paste0("No proteins reach FDR < ", fdr_thresh, " (min FDR = ",
           round(min(res_df$fdr, na.rm = TRUE), 3), "); top ", n,
           " shown by rank only, ties broken by input order")
  } else {
    paste0(n_sig, " protein(s) reach FDR < ", fdr_thresh)
  }
  p <- ggplot(plot_df, aes(x = HR, y = protein)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(xmin = lower, xmax = upper)) +
    labs(title = title, subtitle = subtitle_txt,
         x = "Hazard ratio (95% CI), clinical-adjusted", y = NULL) +
    theme_bw()
  cairo_pdf(paste0(path_stub, ".pdf"), width = 7, height = max(4, n * 0.3) + 0.5)
  print(p)
  dev.off()
}

if (nrow(cox_os) > 0)  make_forest(cox_os,  "Univariate OS, clinical-adjusted",  file.path(paths$out_fig, "forest_OS"))
if (nrow(cox_pfs) > 0) make_forest(cox_pfs, "Univariate PFS, clinical-adjusted", file.path(paths$out_fig, "forest_PFS"))

## ---- 9. Train/holdout KM + time-dependent AUC, and per-protein risk-score KM panels ----- ##
# Matches the other researcher's figure style: (a) overall train vs holdout
# KM with a time-dependent ROC AUC curve comparing clinical-only vs
# clinical+protein models; (b) per-candidate-protein risk-score panels,
# median-split on TRAIN, same cutoff applied to holdout, log-rank p-value
# from the holdout stratification (the translational-relevant test).
need2 <- c("timeROC")
for (pkg in need2) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(timeROC)

# Reuse the same seed=1 split used at CHECKPOINT 4 for OS/PFS, so these
# figures correspond exactly to the res_os/res_pfs numbers already reported.
xy_os <- build_xy(full_df, "os", protein_names, CLIN_COVARS)
strata_os <- xy_os$y[, "status"]
split_os <- stratified_split(strata_os, holdout_frac = 0.25, seed = 1)

xy_pfs <- build_xy(full_df, "pfs", protein_names, CLIN_COVARS)
strata_pfs <- xy_pfs$y[, "status"]
split_pfs <- stratified_split(strata_pfs, holdout_frac = 0.25, seed = 1)

fit_km_and_timeroc <- function(xy, split_idx, time_col, event_col,
                                extra_protein = NULL, time_grid = seq(6, 96, by = 6)) {
  train_idx <- split_idx$train; holdout_idx <- split_idx$holdout
  df_train <- xy$df[train_idx, ]; df_holdout <- xy$df[holdout_idx, ]

  # Clinical-only Cox (train), linear predictor on both sets
  f_clin <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                               paste(CLIN_COVARS, collapse = " + ")))
  fit_clin <- survival::coxph(f_clin, data = df_train)
  lp_clin_train   <- predict(fit_clin, newdata = df_train,   type = "lp")
  lp_clin_holdout <- predict(fit_clin, newdata = df_holdout, type = "lp")

  # Clinical + protein Cox (train), linear predictor on both sets
  lp_full_train <- lp_full_holdout <- NULL
  if (!is.null(extra_protein)) {
    f_full <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ `", extra_protein, "` + ",
                                 paste(CLIN_COVARS, collapse = " + ")))
    fit_full <- survival::coxph(f_full, data = df_train)
    lp_full_train   <- predict(fit_full, newdata = df_train,   type = "lp")
    lp_full_holdout <- predict(fit_full, newdata = df_holdout, type = "lp")
  }

  troc <- function(lp_train, lp_holdout) {
    r_train <- timeROC::timeROC(T = df_train[[time_col]], delta = df_train[[event_col]],
                                 marker = lp_train, cause = 1, times = time_grid,
                                 iid = FALSE)
    r_holdout <- timeROC::timeROC(T = df_holdout[[time_col]], delta = df_holdout[[event_col]],
                                   marker = lp_holdout, cause = 1, times = time_grid,
                                   iid = FALSE)
    list(train = r_train$AUC, holdout = r_holdout$AUC)
  }
  auc_clin <- troc(lp_clin_train, lp_clin_holdout)
  auc_full <- if (!is.null(extra_protein)) troc(lp_full_train, lp_full_holdout) else NULL

  list(df_train = df_train, df_holdout = df_holdout,
       auc_clin = auc_clin, auc_full = auc_full, time_grid = time_grid,
       time_col = time_col, event_col = event_col)
}

plot_train_holdout_km_roc <- function(res, path_stub, outcome_label = "OS", protein_label = NULL) {
  f <- as.formula(paste0("Surv(", res$time_col, ", ", res$event_col, ") ~ 1"))
  km_train   <- survival::survfit(f, data = res$df_train)
  km_holdout <- survival::survfit(f, data = res$df_holdout)

  km_df <- rbind(
    data.frame(time = km_train$time,   surv = km_train$surv,   set = "Train"),
    data.frame(time = km_holdout$time, surv = km_holdout$surv, set = "Holdout"))
  roc_df <- data.frame(time = res$time_grid, auc = res$auc_clin$train,  series = "Clinical only: Train")
  roc_df <- rbind(roc_df,
    data.frame(time = res$time_grid, auc = res$auc_clin$holdout, series = "Clinical only: Holdout"))
  if (!is.null(res$auc_full)) {
    roc_df <- rbind(roc_df,
      data.frame(time = res$time_grid, auc = res$auc_full$train,   series = "Clinical and protein: Train"),
      data.frame(time = res$time_grid, auc = res$auc_full$holdout, series = "Clinical and protein: Holdout"))
  }

  p_km <- ggplot(km_df, aes(x = time, y = surv, colour = set)) +
    geom_step() +
    labs(x = "Months since treatment start", y = "Probability of survival",
         colour = NULL, title = paste0("Train vs holdout ", outcome_label,
                                        if (!is.null(protein_label)) paste0(" (", protein_label, ")") else "")) +
    theme_bw()
  p_roc <- ggplot(roc_df, aes(x = time, y = auc, colour = series, linetype = series)) +
    geom_line() + geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    labs(x = "Months since treatment start", y = "Time-dependent ROC AUC", colour = NULL, linetype = NULL) +
    ylim(0, 1) + theme_bw()

  cairo_pdf(paste0(path_stub, ".pdf"), width = 7, height = 8)
  gridExtra::grid.arrange(p_km, p_roc, ncol = 1)
  dev.off()
}

need3 <- c("gridExtra"); for (pkg in need3) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(gridExtra)

top_os_protein <- if (nrow(cox_os) > 0) cox_os$protein[order(cox_os$fdr)][1] else NULL
res_km_roc_os <- fit_km_and_timeroc(xy_os, split_os, "os_months", "os_event", extra_protein = top_os_protein)
plot_train_holdout_km_roc(res_km_roc_os, file.path(paths$out_fig, "train_holdout_KM_ROC_OS"),
                           outcome_label = "OS", protein_label = top_os_protein)

top_pfs_protein <- if (nrow(cox_pfs) > 0) cox_pfs$protein[order(cox_pfs$fdr)][1] else NULL
res_km_roc_pfs <- fit_km_and_timeroc(xy_pfs, split_pfs, "pfs_months", "pfs_event", extra_protein = top_pfs_protein)
plot_train_holdout_km_roc(res_km_roc_pfs, file.path(paths$out_fig, "train_holdout_KM_ROC_PFS"),
                           outcome_label = "PFS", protein_label = top_pfs_protein)

# Per-protein risk-score KM panels: median-split on TRAIN linear predictor,
# same cutoff applied to holdout, log-rank p-value computed on holdout.
plot_protein_risk_km_panel <- function(xy, split_idx, proteins, path_stub,
                                        time_col, event_col, ncol_panel = 3) {
  train_idx <- split_idx$train; holdout_idx <- split_idx$holdout
  df_train <- xy$df[train_idx, ]; df_holdout <- xy$df[holdout_idx, ]
  panels <- list()
  panel_names <- c(proteins, "Clinical only")
  for (nm in panel_names) {
    f <- if (nm == "Clinical only") {
      as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(CLIN_COVARS, collapse = " + ")))
    } else {
      as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ `", nm, "` + ", paste(CLIN_COVARS, collapse = " + ")))
    }
    fit <- tryCatch(survival::coxph(f, data = df_train), error = function(e) NULL)
    if (is.null(fit)) next
    lp_train   <- predict(fit, newdata = df_train,   type = "lp")
    lp_holdout <- predict(fit, newdata = df_holdout, type = "lp")
    cutoff <- median(lp_train)
    grp_train   <- ifelse(lp_train   >= cutoff, "above median", "below median")
    grp_holdout <- ifelse(lp_holdout >= cutoff, "above median", "below median")

    f_diff <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ grp_holdout"))
    sdiff <- tryCatch(
      survival::survdiff(f_diff, data = df_holdout),
      error = function(e) NULL)
    p_val <- if (!is.null(sdiff)) 1 - pchisq(sdiff$chisq, df = 1) else NA_real_

    plot_df <- rbind(
      data.frame(time = df_train[[time_col]],   event = df_train[[event_col]],   grp = grp_train,   set = "Train"),
      data.frame(time = df_holdout[[time_col]], event = df_holdout[[event_col]], grp = grp_holdout, set = "Holdout"))
    panels[[nm]] <- list(plot_df = plot_df, p_val = p_val, label = nm)
  }

  n_panel <- length(panels)
  n_col <- min(ncol_panel, n_panel)
  n_row <- ceiling(n_panel / n_col)
  cairo_pdf(paste0(path_stub, ".pdf"), width = 3.2 * n_col, height = 3.5 * n_row)
  par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))
  for (nm in names(panels)) {
    pd <- panels[[nm]]$plot_df
    fit <- survival::survfit(Surv(time, event) ~ grp + set, data = pd)
    strata_nm <- names(fit$strata)
    is_train  <- grepl("set=Train", strata_nm)
    is_above  <- grepl("grp=above", strata_nm)
    col_vec <- ifelse(is_train, ifelse(is_above, "#FDB863", "#F4A582"),
                                 ifelse(is_above, "#D6604D", "#B2182B"))
    lty_vec <- ifelse(is_train, 2, 1)
    plot(fit, col = col_vec, lty = lty_vec,
         mark.time = TRUE, xlab = "Time in months", ylab = "Probability of survival",
         main = paste0(panels[[nm]]$label, "\np-value: ", round(panels[[nm]]$p_val, 4)))
  }
  plot.new()
  legend("center", legend = c("Train: below median", "Train: above median",
                               "Holdout: below median", "Holdout: above median"),
         col = c("#F4A582", "#FDB863", "#B2182B", "#D6604D"), lty = c(2, 2, 1, 1), bty = "n")
  dev.off()

  cat("KM risk-score panels saved:", paste0(path_stub, ".pdf"), "\n")
  invisible(panels)
}

# Default candidates: top 2 by FDR from the OS / PFS forests, edit as needed
panel_candidates_os <- if (nrow(cox_os) >= 2) head(cox_os$protein[order(cox_os$fdr)], 2) else character(0)
km_panels_os <- plot_protein_risk_km_panel(xy_os, split_os, panel_candidates_os,
                                            file.path(paths$out_fig, "protein_clinical_KM_panels_OS"),
                                            time_col = "os_months", event_col = "os_event")

panel_candidates_pfs <- if (nrow(cox_pfs) >= 2) head(cox_pfs$protein[order(cox_pfs$fdr)], 2) else character(0)
km_panels_pfs <- plot_protein_risk_km_panel(xy_pfs, split_pfs, panel_candidates_pfs,
                                             file.path(paths$out_fig, "protein_clinical_KM_panels_PFS"),
                                             time_col = "pfs_months", event_col = "pfs_event")

cat("Done. Tables in:", paths$out_tab, "| Figures in:", paths$out_fig, "\n")
