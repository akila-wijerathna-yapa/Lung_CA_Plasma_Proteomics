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

# ── Extra packages for this module only ───────────────────────
ensure_package("Boruta",       "cran", required = TRUE)
ensure_package("randomForest", "cran", required = TRUE)
if (!requireNamespace("randomForestSRC", quietly = TRUE))
  ensure_package("randomForestSRC", "cran", required = FALSE)
if (!requireNamespace("mixOmics", quietly = TRUE))
  ensure_package("mixOmics", "bioc", required = FALSE)
suppressPackageStartupMessages({
  library(Boruta)
  library(randomForest)
  if (requireNamespace("randomForestSRC", quietly = TRUE)) library(randomForestSRC)
  if (requireNamespace("mixOmics", quietly = TRUE)) library(mixOmics)
})

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_05.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 03, 05, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
RNG_SEED       <- 42
EN_ALPHA       <- 0.5    # Elastic Net mixing (0.5 = equal L1/L2)
MIN_CLASS_N    <- 5
MIN_COX_N      <- 35
MIN_COX_EVENTS <- 12
SPLSDA_NCOMP   <- 2
SPLSDA_KEEPX   <- 15     # features kept per component
RF_NTREE       <- 1000
BORUTA_MAXRUNS <- 200

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning", "altml")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning", "altml")
dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ── Local functions ────────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# ── Clinical covariate matrix (forced-unpenalized) ────────────
# age, sex, stage, histology, treatment-class dummies, aligned to
# sample_id order. Zero-variance columns dropped automatically --
# necessary for single-new_group subgroups where treatment dummies
# are constant by construction.
build_clinical_covariates <- function(sample_ids, clinical_df = analysis1_df) {
  cd  <- clinical_df[match(sample_ids, clinical_df$sample_id), ]
  age <- as.numeric(difftime(cd$start_date, cd$dob, units = "days")) / 365.25
  mat <- cbind(
    age          = age,
    sex_female   = as.integer(cd$sex == "FEMALE"),
    stage_III_IV = as.integer(cd$stage_collapsed == "III_IV"),
    hist_adeno   = as.integer(cd$histology == "Adenocarcinoma"),
    tx_chemo     = as.integer(grepl("Chemo",   cd$treatment_class, ignore.case = TRUE)),
    tx_io        = as.integer(grepl("IO",      cd$treatment_class, ignore.case = TRUE)),
    tx_rt        = as.integer(grepl("RT",      cd$treatment_class, ignore.case = TRUE)),
    tx_surgery   = as.integer(grepl("Surgery", cd$treatment_class, ignore.case = TRUE)),
    tx_tki       = as.integer(grepl("TKI",     cd$treatment_class, ignore.case = TRUE))
  )
  for (j in seq_len(ncol(mat))) {
    mat[is.na(mat[, j]), j] <- stats::median(mat[, j], na.rm = TRUE)
  }
  col_sd <- apply(mat, 2, sd, na.rm = TRUE)
  mat <- mat[, col_sd > 0, drop = FALSE]
  rownames(mat) <- sample_ids
  mat
}
add_feature_annotation <- function(df, feature_meta) {
  if (is.null(df) || nrow(df) == 0 || !"nELISA_ID" %in% names(df)) return(tibble())
  left_join(df, feature_meta, by = "nELISA_ID")
}
save_ggplot <- function(p, path_stub, width, height) {
  cairo_pdf(paste0(path_stub, ".pdf"), width = width, height = height); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white"); print(p); dev.off()
}

# ── (1) Elastic Net ────────────────────────────────────────────
run_elastic_net_binomial <- function(X, y, label, X_clin = NULL) {
  X_std <- scale(X); X_std[is.nan(X_std)] <- 0
  if (!is.null(X_clin) && nrow(X_clin) == nrow(X_std)) {
    X_design <- cbind(X_clin[rownames(X_std), , drop = FALSE], X_std)
    pf <- c(rep(0, ncol(X_clin)), rep(1, ncol(X_std)))
  } else {
    X_design <- X_std
    pf <- rep(1, ncol(X_std))
  }
  set.seed(RNG_SEED)
  cv_fit <- glmnet::cv.glmnet(X_design, y, family = "binomial", alpha = EN_ALPHA,
                              penalty.factor = pf, nfolds = min(5, length(y)))
  coef_all <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
  names(coef_all) <- colnames(X_design)
  coef_v <- coef_all[colnames(X)]  # protein-only slice, original order
  selected <- colnames(X)[coef_v != 0]
  prob <- as.numeric(predict(cv_fit, newx = X_design, s = "lambda.min", type = "response"))
  op <- options(warn = -1)
  auc_val <- as.numeric(pROC::auc(pROC::roc(y, prob, quiet = TRUE)))
  options(op)
  tibble(nELISA_ID = colnames(X), en_coef = coef_v,
         en_selected = colnames(X) %in% selected, en_model_auc = auc_val)
}

run_elastic_net_cox <- function(X, surv_obj, label, X_clin = NULL) {
  X_std <- scale(X); X_std[is.nan(X_std)] <- 0
  if (!is.null(X_clin) && nrow(X_clin) == nrow(X_std)) {
    X_design <- cbind(X_clin[rownames(X_std), , drop = FALSE], X_std)
    pf <- c(rep(0, ncol(X_clin)), rep(1, ncol(X_std)))
  } else {
    X_design <- X_std
    pf <- rep(1, ncol(X_std))
  }
  set.seed(RNG_SEED)
  cv_fit <- glmnet::cv.glmnet(X_design, surv_obj, family = "cox", alpha = EN_ALPHA,
                              penalty.factor = pf, nfolds = min(5, nrow(X)))
  coef_all <- as.numeric(coef(cv_fit, s = "lambda.min"))
  names(coef_all) <- colnames(X_design)
  coef_v <- coef_all[colnames(X)]  # protein-only slice, original order
  selected <- colnames(X)[coef_v != 0]
  lp <- as.numeric(predict(cv_fit, newx = X_design, s = "lambda.min", type = "link"))
  c_idx <- tryCatch(survival::concordance(surv_obj ~ lp, reverse = TRUE)$concordance,
                     error = function(e) NA_real_)
  tibble(nELISA_ID = colnames(X), en_coef = coef_v,
         en_selected = colnames(X) %in% selected, en_c_index = c_idx)
}

# ── (2) sPLS-DA ────────────────────────────────────────────────
run_splsda <- function(X, y, label) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    message("    mixOmics not available - skipping sPLS-DA"); return(tibble())
  }
  y_fac <- factor(y, labels = c("Class0", "Class1"))
  keepX <- rep(min(SPLSDA_KEEPX, ncol(X)), SPLSDA_NCOMP)
  fit <- tryCatch(
    mixOmics::splsda(X, y_fac, ncomp = SPLSDA_NCOMP, keepX = keepX),
    error = function(e) { message("    sPLS-DA fit failed: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(tibble())
  vip_mat <- mixOmics::vip(fit)
  vip_comp1 <- vip_mat[, 1]
  selected_ids <- unique(unlist(lapply(seq_len(SPLSDA_NCOMP), function(comp) {
    mixOmics::selectVar(fit, comp = comp)$name
  })))
  pred <- tryCatch(predict(fit, X)$class$max.dist[, SPLSDA_NCOMP], error = function(e) NULL)
  acc <- if (!is.null(pred)) mean(pred == y_fac, na.rm = TRUE) else NA_real_
  tibble(nELISA_ID = names(vip_comp1), splsda_vip = vip_comp1,
         splsda_selected = names(vip_comp1) %in% selected_ids,
         splsda_train_accuracy = acc)
}

# ── (3) Boruta + Random Forest ─────────────────────────────────
run_boruta_rf_binomial <- function(X, y, label) {
  y_fac <- factor(y, labels = c("Class0", "Class1"))
  set.seed(RNG_SEED)
  bor <- tryCatch(
    Boruta::Boruta(x = as.data.frame(X), y = y_fac,
                   maxRuns = BORUTA_MAXRUNS, doTrace = 0),
    error = function(e) { message("    Boruta failed: ", conditionMessage(e)); NULL })
  if (is.null(bor)) return(tibble())
  imp_hist  <- Boruta::attStats(bor)
  confirmed <- rownames(imp_hist)[imp_hist$decision == "Confirmed"]
  set.seed(RNG_SEED)
  rf_fit <- randomForest::randomForest(x = X, y = y_fac, ntree = RF_NTREE, importance = TRUE)
  prob <- predict(rf_fit, X, type = "prob")[, "Class1"]
  op <- options(warn = -1)
  auc_val <- as.numeric(pROC::auc(pROC::roc(y, prob, quiet = TRUE)))
  options(op)
  tibble(nELISA_ID = rownames(imp_hist),
         rf_mean_importance = imp_hist$meanImp,
         boruta_decision = imp_hist$decision,
         boruta_confirmed = rownames(imp_hist) %in% confirmed,
         rf_model_auc = auc_val)
}

run_rf_survival <- function(X, sample_df, time_col, event_col, label) {
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    message("    randomForestSRC not available - skipping RF survival"); return(tibble())
  }
  op  <- options(warn = -1)
  tval <- as.numeric(sample_df[[time_col]])
  options(op)
  eval_ <- suppressWarnings(as.integer(sample_df[[event_col]]))
  keep <- !is.na(tval) & !is.na(eval_) & rownames(X) %in% sample_df$sample_id
  df_rf <- as.data.frame(X[keep, , drop = FALSE])
  df_rf$time  <- tval[keep]
  df_rf$event <- eval_[keep]
  set.seed(RNG_SEED)
  fit <- tryCatch(
    randomForestSRC::rfsrc(Surv(time, event) ~ ., data = df_rf,
                           ntree = RF_NTREE, importance = TRUE),
    error = function(e) { message("    RF-survival failed: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(tibble())
  vimp <- randomForestSRC::vimp(fit)$importance
  c_idx <- 1 - fit$err.rate[RF_NTREE]
  tibble(nELISA_ID = names(vimp), rf_surv_vimp = vimp, rf_surv_c_index = c_idx) %>%
    arrange(desc(rf_surv_vimp))
}

# ── Plotting helpers ──────────────────────────────────────────
plot_top_bar <- function(df, value_col, label_col, title, x_label, path_stub, n = 20,
                         fill_col = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
  plot_df <- df %>%
    filter(is.finite(.data[[value_col]])) %>%
    arrange(desc(abs(.data[[value_col]]))) %>%
    slice_head(n = min(n, nrow(.))) %>%
    mutate(label_order = factor(.data[[label_col]], levels = rev(.data[[label_col]])))
  if (!nrow(plot_df)) return(invisible(NULL))
  p <- if (!is.null(fill_col)) {
    ggplot(plot_df, aes(.data[[value_col]], label_order, fill = .data[[fill_col]])) +
      geom_col(width = 0.75) +
      scale_fill_manual(values = c("TRUE" = "#1565C0", "FALSE" = "#90A4AE"))
  } else {
    ggplot(plot_df, aes(.data[[value_col]], label_order)) +
      geom_col(width = 0.75, fill = "#2E7D32")
  }
  p <- p + labs(title = title, x = x_label, y = NULL) + theme_minimal(base_size = 11)
  save_ggplot(p, path_stub, width = 8.5, height = 7)
}

# ── Comparison specifications (same as 08b, rebuilt locally) ──
message("\n=== 08d | Building comparison specifications ===")
complete_cols <- colnames(expr_log2_rep)[
  colSums(is.na(expr_log2_rep[analysis1_df$sample_id, , drop = FALSE])) == 0]
message("Complete features: ", length(complete_cols))

response_specs <- list(
  altml_01_CR_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("CR","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("CR","PD")) %>%
          mutate(y = as.integer(response_4cat == "CR")) %>% pull(y)),
  altml_02_CRPR_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("CR","PR","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("CR","PR","PD")) %>%
          mutate(y = as.integer(response_4cat %in% c("CR","PR"))) %>% pull(y)),
  altml_03_CRPR_vs_SDPD = list(
    ids = analysis1_df %>% filter(!is.na(response_binary)) %>% pull(sample_id),
    y   = analysis1_df %>% filter(!is.na(response_binary)) %>%
          mutate(y = as.integer(response_binary == "PR/CR")) %>% pull(y)),
  altml_04_CRPRSD_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("CR","PR","SD","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("CR","PR","SD","PD")) %>%
          mutate(y = as.integer(response_4cat != "PD")) %>% pull(y)),
  altml_05_PR_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("PR","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("PR","PD")) %>%
          mutate(y = as.integer(response_4cat == "PR")) %>% pull(y))
)

recurrence_base <- analysis1_df %>% filter(!is.na(recurrence_binary_sensitivity))
recurrence_specs <- list(
  altml_06_CTx_recurrence   = recurrence_base %>% filter(new_group == "CTx"),
  altml_07_Sx_recurrence    = recurrence_base %>% filter(new_group == "Sx"),
  altml_08_IO_recurrence    = recurrence_base %>% filter(new_group == "IO"),
  altml_09_CTxIO_recurrence = recurrence_base %>% filter(new_group == "CTx + IO")
)

# ── Run all binomial comparisons ──────────────────────────────
message("\n=== 08d | Response + recurrence comparisons ===")
all_binomial_specs <- c(
  lapply(response_specs, function(s) list(
    X = expr_log2_rep[s$ids, complete_cols, drop = FALSE], y = s$y)),
  lapply(recurrence_specs, function(df) list(
    X = expr_log2_rep[df$sample_id, complete_cols, drop = FALSE],
    y = as.integer(df$recurrence_binary_sensitivity == "Recurrence")))
)

altml_results <- list()

for (nm in names(all_binomial_specs)) {
  spec <- all_binomial_specs[[nm]]
  n0 <- sum(spec$y == 0); n1 <- sum(spec$y == 1)
  message("\n--- ", nm, " --- n0=", n0, " n1=", n1)
  if (n0 < MIN_CLASS_N || n1 < MIN_CLASS_N) { message("  SKIP: small class"); next }

  # (1) Elastic Net
  message("  Running Elastic Net...")
  X_clin <- build_clinical_covariates(rownames(spec$X))
  en_res <- run_elastic_net_binomial(spec$X, spec$y, nm, X_clin = X_clin)
  message("  EN selected: ", sum(en_res$en_selected), " | AUC: ", round(en_res$en_model_auc[1], 3))

  # (2) sPLS-DA
  message("  Running sPLS-DA...")
  spls_res <- run_splsda(spec$X, spec$y, nm)
  if (nrow(spls_res))
    message("  sPLS-DA selected: ", sum(spls_res$splsda_selected),
            " | train accuracy: ", round(spls_res$splsda_train_accuracy[1], 3))

  # (3) Boruta + RF
  message("  Running Boruta + Random Forest...")
  rf_res <- run_boruta_rf_binomial(spec$X, spec$y, nm)
  if (nrow(rf_res))
    message("  RF confirmed: ", sum(rf_res$boruta_confirmed),
            " | AUC: ", round(rf_res$rf_model_auc[1], 3))

  # Merge all three
  merged <- en_res
  if (nrow(spls_res)) merged <- full_join(merged, spls_res, by = "nELISA_ID")
  if (nrow(rf_res))   merged <- full_join(merged, rf_res,   by = "nELISA_ID")
  merged <- merged %>%
    mutate(
      n_methods_selected = rowSums(cbind(
        en_selected, dplyr::coalesce(splsda_selected, FALSE), dplyr::coalesce(boruta_confirmed, FALSE)
      ), na.rm = TRUE)
    ) %>%
    add_feature_annotation(feature_meta) %>%
    arrange(desc(n_methods_selected), desc(abs(dplyr::coalesce(en_coef, 0))))

  safe_write_csv(merged, file.path(OUT_TBL, paste0(nm, "_all_methods.csv")))
  altml_results[[nm]] <- merged

  # Plots
  en_plot_df <- merged %>% filter(en_selected) %>% mutate(abs_coef = abs(en_coef))
  plot_top_bar(en_plot_df, "abs_coef", "Protein_Label",
              paste(nm, "- Elastic Net selected"), "|Coefficient| (standardized)",
              file.path(OUT_FIG, paste0(nm, "_elasticnet")))
  if (nrow(spls_res))
    plot_top_bar(merged, "splsda_vip", "Protein_Label",
                paste(nm, "- sPLS-DA VIP scores"), "VIP score",
                file.path(OUT_FIG, paste0(nm, "_splsda")), fill_col = "splsda_selected")
  if (nrow(rf_res))
    plot_top_bar(merged %>% filter(boruta_confirmed), "rf_mean_importance", "Protein_Label",
                paste(nm, "- Boruta confirmed (RF importance)"), "Mean importance",
                file.path(OUT_FIG, paste0(nm, "_borutaRF")))

  # Cross-method overlap with STABL (if available)
  stabl_path <- file.path(paths$tables, "07_analysis6_machine_learning", "stabl",
                          gsub("altml", "stabl", nm), fsep = .Platform$file.sep)
  stabl_file <- file.path(paths$tables, "07_analysis6_machine_learning", "stabl",
                          paste0(gsub("altml", "stabl", nm), "_stabl_selected_features.csv"))
  if (file.exists(stabl_file)) {
    stabl_ids <- readr::read_csv(stabl_file, show_col_types = FALSE)$nELISA_ID
    overlap_tbl <- merged %>%
      mutate(stabl_selected = nELISA_ID %in% stabl_ids) %>%
      filter(n_methods_selected > 0 | stabl_selected) %>%
      select(nELISA_ID, Protein_Label, en_selected, splsda_selected,
             boruta_confirmed, stabl_selected, n_methods_selected)
    safe_write_csv(overlap_tbl, file.path(OUT_TBL, paste0(nm, "_cross_method_overlap.csv")))
  }
}

# ── Cox: Elastic Net + RF-survival ────────────────────────────
message("\n=== 08d | Cox OS and PFS (Elastic Net + RF-survival) ===")

op <- options(warn = -1)
os_time  <- as.numeric(analysis1_df$os_months)
os_event <- as.integer(analysis1_df$os_event)
options(op)
cox_keep_os <- !is.na(os_time) & !is.na(os_event)

if (sum(cox_keep_os) >= MIN_COX_N) {
  X_cox_os <- expr_log2_rep[analysis1_df$sample_id[cox_keep_os], complete_cols, drop = FALSE]
  surv_os  <- survival::Surv(os_time[cox_keep_os], os_event[cox_keep_os])

  message("  Elastic Net Cox OS...")
  X_clin_os <- build_clinical_covariates(analysis1_df$sample_id[cox_keep_os])
  en_os <- run_elastic_net_cox(X_cox_os, surv_os, "OS", X_clin = X_clin_os)
  message("  EN-Cox OS selected: ", sum(en_os$en_selected), " | C-index: ", round(en_os$en_c_index[1], 3))

  message("  RF-survival OS...")
  rf_os <- run_rf_survival(X_cox_os, analysis1_df[cox_keep_os, ], "os_months", "os_event", "OS")
  if (nrow(rf_os)) message("  RF-survival OS C-index: ", round(rf_os$rf_surv_c_index[1], 3))

  merged_os <- en_os
  if (nrow(rf_os)) merged_os <- full_join(merged_os, rf_os, by = "nELISA_ID")
  merged_os <- add_feature_annotation(merged_os, feature_meta) %>% arrange(desc(abs(en_coef)))
  safe_write_csv(merged_os, file.path(OUT_TBL, "altml_14_cox_OS_all_methods.csv"))
  altml_results[["altml_14_cox_OS"]] <- merged_os

  plot_top_bar(merged_os %>% filter(en_selected), "en_coef", "Protein_Label",
              "Cox OS - Elastic Net selected", "Coefficient (standardized)",
              file.path(OUT_FIG, "altml_14_cox_OS_elasticnet"))
  if (nrow(rf_os))
    plot_top_bar(merged_os, "rf_surv_vimp", "Protein_Label",
                "Cox OS - RF-survival VIMP", "Variable importance",
                file.path(OUT_FIG, "altml_14_cox_OS_rfsurv"))
}

op <- options(warn = -1)
pfs_time  <- as.numeric(analysis1_df$pfs_months_main)
pfs_event <- as.integer(analysis1_df$pfs_event_main)
options(op)
cox_keep_pfs <- !is.na(pfs_time) & !is.na(pfs_event)

if (sum(cox_keep_pfs) >= MIN_COX_N) {
  X_cox_pfs <- expr_log2_rep[analysis1_df$sample_id[cox_keep_pfs], complete_cols, drop = FALSE]
  surv_pfs  <- survival::Surv(pfs_time[cox_keep_pfs], pfs_event[cox_keep_pfs])

  message("  Elastic Net Cox PFS...")
  X_clin_pfs <- build_clinical_covariates(analysis1_df$sample_id[cox_keep_pfs])
  en_pfs <- run_elastic_net_cox(X_cox_pfs, surv_pfs, "PFS", X_clin = X_clin_pfs)
  message("  EN-Cox PFS selected: ", sum(en_pfs$en_selected), " | C-index: ", round(en_pfs$en_c_index[1], 3))

  message("  RF-survival PFS...")
  rf_pfs <- run_rf_survival(X_cox_pfs, analysis1_df[cox_keep_pfs, ], "pfs_months_main", "pfs_event_main", "PFS")
  if (nrow(rf_pfs)) message("  RF-survival PFS C-index: ", round(rf_pfs$rf_surv_c_index[1], 3))

  merged_pfs <- en_pfs
  if (nrow(rf_pfs)) merged_pfs <- full_join(merged_pfs, rf_pfs, by = "nELISA_ID")
  merged_pfs <- add_feature_annotation(merged_pfs, feature_meta) %>% arrange(desc(abs(en_coef)))
  safe_write_csv(merged_pfs, file.path(OUT_TBL, "altml_15_cox_PFS_all_methods.csv"))
  altml_results[["altml_15_cox_PFS"]] <- merged_pfs

  plot_top_bar(merged_pfs %>% filter(en_selected), "en_coef", "Protein_Label",
              "Cox PFS - Elastic Net selected", "Coefficient (standardized)",
              file.path(OUT_FIG, "altml_15_cox_PFS_elasticnet"))
  if (nrow(rf_pfs))
    plot_top_bar(merged_pfs, "rf_surv_vimp", "Protein_Label",
                "Cox PFS - RF-survival VIMP", "Variable importance",
                file.path(OUT_FIG, "altml_15_cox_PFS_rfsurv"))
}

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 08d | Saving checkpoint_08d.RData ===")
save(altml_results, file = file.path(paths$rdata, "checkpoint_08d.RData"))
message("08d_alternative_ml.R complete")
