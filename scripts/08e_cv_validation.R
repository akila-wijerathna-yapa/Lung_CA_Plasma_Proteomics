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

suppressPackageStartupMessages({
  if (!requireNamespace("Boruta", quietly = TRUE)) install.packages("Boruta")
  if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) install.packages("randomForestSRC")
  library(Boruta)
  library(randomForest)
  if (requireNamespace("randomForestSRC", quietly = TRUE)) library(randomForestSRC)
  if (requireNamespace("mixOmics", quietly = TRUE)) library(mixOmics)
})

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))  # provides analysis1_df — previously missing, relied on session state
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
RNG_SEED        <- 42
K_FOLDS         <- 5
N_REPEATS       <- 5
EN_ALPHA        <- 0.5
MIN_CLASS_N     <- 5
MIN_COX_N       <- 35
MIN_COX_EVENTS  <- 12
SPLSDA_NCOMP    <- 2
SPLSDA_KEEPX    <- 15
RF_NTREE_CV     <- 500     # reduced for speed within CV loop
BORUTA_MAXRUNS_CV <- 60    # reduced for speed within CV loop

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning", "cv")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning", "cv")
dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ── Local functions ────────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# Stratified fold assignment: keeps class balance in every fold
make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold_id <- rep(NA_integer_, length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls)
    idx <- sample(idx)  # shuffle within class
    fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold_id
}

scale_train_apply <- function(X_train, X_test) {
  mu <- colMeans(X_train); sig <- apply(X_train, 2, sd)
  sig[sig == 0 | !is.finite(sig)] <- 1
  X_train_s <- sweep(sweep(X_train, 2, mu, "-"), 2, sig, "/")
  X_test_s  <- sweep(sweep(X_test,  2, mu, "-"), 2, sig, "/")
  X_train_s[!is.finite(X_train_s)] <- 0
  X_test_s[!is.finite(X_test_s)]   <- 0
  list(train = X_train_s, test = X_test_s)
}

# ── CV runner: binomial, all three methods ────────────────────
run_cv_binomial <- function(X, y, label, k = K_FOLDS, n_repeats = N_REPEATS) {
  n <- nrow(X)
  oof <- list(en = rep(NA_real_, n), splsda = rep(NA_real_, n), rf = rep(NA_real_, n))
  rep_auc <- list(en = numeric(0), splsda = numeric(0), rf = numeric(0))

  for (rep_i in seq_len(n_repeats)) {
    fold_id <- make_stratified_folds(y, k, seed = RNG_SEED + rep_i)
    rep_pred <- list(en = rep(NA_real_, n), splsda = rep(NA_real_, n), rf = rep(NA_real_, n))

    for (fold_i in seq_len(k)) {
      test_idx  <- which(fold_id == fold_i)
      train_idx <- which(fold_id != fold_i)
      if (length(unique(y[train_idx])) < 2 || min(table(y[train_idx])) < 3) next

      X_tr <- X[train_idx, , drop = FALSE]; y_tr <- y[train_idx]
      X_te <- X[test_idx,  , drop = FALSE]
      sc   <- scale_train_apply(X_tr, X_te)

      # (1) Elastic Net - refit fresh on this fold's training data
      en_pred <- tryCatch({
        cv_fit <- glmnet::cv.glmnet(sc$train, y_tr, family = "binomial",
                                    alpha = EN_ALPHA, nfolds = min(5, length(y_tr)))
        as.numeric(predict(cv_fit, newx = sc$test, s = "lambda.min", type = "response"))
      }, error = function(e) rep(NA_real_, length(test_idx)))
      rep_pred$en[test_idx] <- en_pred

      # (2) sPLS-DA - refit fresh
      splsda_pred <- tryCatch({
        if (!requireNamespace("mixOmics", quietly = TRUE)) stop("no mixOmics")
        y_fac_tr <- factor(y_tr, labels = c("Class0", "Class1"))
        keepX <- rep(min(SPLSDA_KEEPX, ncol(X_tr)), SPLSDA_NCOMP)
        fit <- mixOmics::splsda(X_tr, y_fac_tr, ncomp = SPLSDA_NCOMP, keepX = keepX)
        pr  <- predict(fit, X_te)
        as.numeric(pr$predict[, "Class1", SPLSDA_NCOMP])
      }, error = function(e) rep(NA_real_, length(test_idx)))
      rep_pred$splsda[test_idx] <- splsda_pred

      # (3) Boruta-confirmed features -> Random Forest - refit fresh
      rf_pred <- tryCatch({
        y_fac_tr <- factor(y_tr, labels = c("Class0", "Class1"))
        bor <- Boruta::Boruta(x = as.data.frame(X_tr), y = y_fac_tr,
                              maxRuns = BORUTA_MAXRUNS_CV, doTrace = 0)
        confirmed <- names(bor$finalDecision)[bor$finalDecision == "Confirmed"]
        if (length(confirmed) < 1) confirmed <- colnames(X_tr)  # fallback: use all
        rf_fit <- randomForest::randomForest(
          x = X_tr[, confirmed, drop = FALSE], y = y_fac_tr, ntree = RF_NTREE_CV)
        predict(rf_fit, X_te[, confirmed, drop = FALSE], type = "prob")[, "Class1"]
      }, error = function(e) rep(NA_real_, length(test_idx)))
      rep_pred$rf[test_idx] <- rf_pred
    }

    # AUC for this repeat (pooled across its 5 folds)
    for (m in c("en", "splsda", "rf")) {
      keep <- !is.na(rep_pred[[m]])
      if (sum(keep) >= 10 && length(unique(y[keep])) == 2) {
        op <- options(warn = -1)
        auc_i <- tryCatch(as.numeric(pROC::auc(pROC::roc(y[keep], rep_pred[[m]][keep], quiet = TRUE))),
                          error = function(e) NA_real_)
        options(op)
        rep_auc[[m]] <- c(rep_auc[[m]], auc_i)
      }
    }
    # Store last repeat's out-of-fold preds for pooled AUC across all repeats
    for (m in c("en", "splsda", "rf")) {
      avail <- !is.na(rep_pred[[m]])
      oof[[m]][avail] <- rowMeans(cbind(oof[[m]][avail], rep_pred[[m]][avail]), na.rm = TRUE)
    }
  }

  tibble(
    comparison = label,
    method     = c("Elastic Net", "sPLS-DA", "Boruta+RF"),
    cv_auc_mean = vapply(c("en","splsda","rf"), function(m) mean(rep_auc[[m]], na.rm=TRUE), numeric(1)),
    cv_auc_sd   = vapply(c("en","splsda","rf"), function(m) sd(rep_auc[[m]], na.rm=TRUE), numeric(1)),
    n_repeats_completed = vapply(c("en","splsda","rf"), function(m) length(rep_auc[[m]]), integer(1))
  )
}

# ── CV runner: Cox (Elastic Net + RF-survival) ────────────────
run_cv_cox <- function(X, time_v, event_v, label, k = K_FOLDS, n_repeats = N_REPEATS) {
  n <- nrow(X)
  rep_cidx <- list(en = numeric(0), rfsurv = numeric(0))

  for (rep_i in seq_len(n_repeats)) {
    set.seed(RNG_SEED + rep_i)
    fold_id <- sample(rep(seq_len(k), length.out = n))
    oof_lp_en <- rep(NA_real_, n)
    oof_lp_rf <- rep(NA_real_, n)

    for (fold_i in seq_len(k)) {
      test_idx  <- which(fold_id == fold_i)
      train_idx <- which(fold_id != fold_i)
      if (sum(event_v[train_idx]) < 5) next

      X_tr <- X[train_idx, , drop=FALSE]; X_te <- X[test_idx, , drop=FALSE]
      sc <- scale_train_apply(X_tr, X_te)
      surv_tr <- survival::Surv(time_v[train_idx], event_v[train_idx])

      en_lp <- tryCatch({
        cv_fit <- glmnet::cv.glmnet(sc$train, surv_tr, family = "cox",
                                    alpha = EN_ALPHA, nfolds = min(5, length(train_idx)))
        as.numeric(predict(cv_fit, newx = sc$test, s = "lambda.min", type = "link"))
      }, error = function(e) rep(NA_real_, length(test_idx)))
      oof_lp_en[test_idx] <- en_lp

      if (requireNamespace("randomForestSRC", quietly = TRUE)) {
        rf_lp <- tryCatch({
          df_tr <- as.data.frame(X_tr); df_tr$time <- time_v[train_idx]; df_tr$event <- event_v[train_idx]
          fit <- randomForestSRC::rfsrc(Surv(time, event) ~ ., data = df_tr, ntree = RF_NTREE_CV)
          df_te <- as.data.frame(X_te)
          as.numeric(predict(fit, newdata = df_te)$predicted)
        }, error = function(e) rep(NA_real_, length(test_idx)))
        oof_lp_rf[test_idx] <- rf_lp
      }
    }

    keep_en <- !is.na(oof_lp_en)
    if (sum(keep_en) >= 15) {
      # reverse = TRUE: higher lp = higher hazard = shorter survival (matches 08d convention).
      # Without this, concordance() assumes the opposite direction and returns 1 - true C-index.
      c_en <- tryCatch(survival::concordance(
        survival::Surv(time_v[keep_en], event_v[keep_en]) ~ oof_lp_en[keep_en],
        reverse = TRUE)$concordance,
        error = function(e) NA_real_)
      rep_cidx$en <- c(rep_cidx$en, c_en)
    }
    keep_rf <- !is.na(oof_lp_rf)
    if (sum(keep_rf) >= 15) {
      # Same fix: RF-survival mortality score is also "higher = worse", needs reverse = TRUE.
      c_rf <- tryCatch(survival::concordance(
        survival::Surv(time_v[keep_rf], event_v[keep_rf]) ~ oof_lp_rf[keep_rf],
        reverse = TRUE)$concordance,
        error = function(e) NA_real_)
      rep_cidx$rfsurv <- c(rep_cidx$rfsurv, c_rf)
    }
  }

  tibble(
    comparison = label,
    method     = c("Elastic Net (Cox)", "RF-survival"),
    cv_cindex_mean = c(mean(rep_cidx$en, na.rm=TRUE), mean(rep_cidx$rfsurv, na.rm=TRUE)),
    cv_cindex_sd   = c(sd(rep_cidx$en, na.rm=TRUE), sd(rep_cidx$rfsurv, na.rm=TRUE)),
    n_repeats_completed = c(length(rep_cidx$en), length(rep_cidx$rfsurv))
  )
}

# ── Plot helper ────────────────────────────────────────────────
plot_cv_auc_summary <- function(df, value_col, sd_col, title, ylab, path_stub) {
  p <- ggplot(df, aes(x = comparison, y = .data[[value_col]], fill = method)) +
    geom_col(position = position_dodge(0.75), width = 0.65) +
    geom_errorbar(aes(ymin = .data[[value_col]] - .data[[sd_col]],
                      ymax = .data[[value_col]] + .data[[sd_col]]),
                 position = position_dodge(0.75), width = 0.25) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c("Elastic Net" = "#2166AC", "sPLS-DA" = "#1B9E77",
                                 "Boruta+RF" = "#D95F02", "Elastic Net (Cox)" = "#2166AC",
                                 "RF-survival" = "#D95F02")) +
    labs(title = title, x = NULL, y = ylab, fill = "Method") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
  cairo_pdf(paste0(path_stub, ".pdf"), width = 11, height = 6.5); print(p); dev.off()
  png(paste0(path_stub, ".png"), width = 11, height = 6.5, units = "in", res = 300, bg = "white")
  print(p); dev.off()
}

# ── Build comparison specs (same as 08b/08d) ──────────────────
message("\n=== 08e | Building comparison specifications ===")
complete_cols <- colnames(expr_log2_rep)[
  colSums(is.na(expr_log2_rep[analysis1_df$sample_id, , drop = FALSE])) == 0]

response_specs <- list(
  cv_01_CR_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("CR","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("CR","PD")) %>%
          mutate(y = as.integer(response_4cat == "CR")) %>% pull(y)),
  cv_02_CRPR_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("CR","PR","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("CR","PR","PD")) %>%
          mutate(y = as.integer(response_4cat %in% c("CR","PR"))) %>% pull(y)),
  cv_03_CRPR_vs_SDPD = list(
    ids = analysis1_df %>% filter(!is.na(response_binary)) %>% pull(sample_id),
    y   = analysis1_df %>% filter(!is.na(response_binary)) %>%
          mutate(y = as.integer(response_binary == "PR/CR")) %>% pull(y)),
  cv_04_CRPRSD_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("CR","PR","SD","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("CR","PR","SD","PD")) %>%
          mutate(y = as.integer(response_4cat != "PD")) %>% pull(y)),
  cv_05_PR_vs_PD = list(
    ids = analysis1_df %>% filter(response_4cat %in% c("PR","PD")) %>% pull(sample_id),
    y   = analysis1_df %>% filter(response_4cat %in% c("PR","PD")) %>%
          mutate(y = as.integer(response_4cat == "PR")) %>% pull(y))
)

recurrence_base <- analysis1_df %>% filter(!is.na(recurrence_binary_sensitivity))
recurrence_specs <- list(
  cv_06_CTx_recurrence   = recurrence_base %>% filter(new_group == "CTx"),
  cv_07_Sx_recurrence    = recurrence_base %>% filter(new_group == "Sx"),
  cv_08_IO_recurrence    = recurrence_base %>% filter(new_group == "IO"),
  cv_09_CTxIO_recurrence = recurrence_base %>% filter(new_group == "CTx + IO")
)

all_binomial_specs <- c(
  lapply(response_specs, function(s) list(
    X = expr_log2_rep[s$ids, complete_cols, drop = FALSE], y = s$y)),
  lapply(recurrence_specs, function(df) list(
    X = expr_log2_rep[df$sample_id, complete_cols, drop = FALSE],
    y = as.integer(df$recurrence_binary_sensitivity == "Recurrence")))
)

# ── Run CV for all binomial comparisons ───────────────────────
message("\n=== 08e | Running repeated CV (", K_FOLDS, "-fold x ", N_REPEATS, " repeats) ===")
message("This will take a while - Boruta+RF refit ", K_FOLDS * N_REPEATS, " times per comparison.")

cv_binomial_results <- list()
for (nm in names(all_binomial_specs)) {
  spec <- all_binomial_specs[[nm]]
  n0 <- sum(spec$y == 0); n1 <- sum(spec$y == 1)
  message("\n--- ", nm, " --- n0=", n0, " n1=", n1)
  if (n0 < MIN_CLASS_N || n1 < MIN_CLASS_N) { message("  SKIP: small class"); next }

  res <- run_cv_binomial(spec$X, spec$y, nm)
  print(res)
  cv_binomial_results[[nm]] <- res
  safe_write_csv(res, file.path(OUT_TBL, paste0(nm, "_cv_results.csv")))
}

cv_binomial_all <- bind_rows(cv_binomial_results)
safe_write_csv(cv_binomial_all, file.path(OUT_TBL, "cv_binomial_all_comparisons.csv"))
plot_cv_auc_summary(cv_binomial_all, "cv_auc_mean", "cv_auc_sd",
                    paste0("Cross-validated AUC (", K_FOLDS, "-fold x ", N_REPEATS, " repeats)"),
                    "CV AUC (mean +/- SD)", file.path(OUT_FIG, "cv_binomial_auc_summary"))

# ── Run CV for Cox OS and PFS ──────────────────────────────────
message("\n=== 08e | Running CV for Cox OS and PFS ===")
cv_cox_results <- list()

op <- options(warn = -1)
os_time  <- as.numeric(analysis1_df$os_months)
os_event <- as.integer(analysis1_df$os_event)
options(op)
cox_keep_os <- !is.na(os_time) & !is.na(os_event)
if (sum(cox_keep_os) >= MIN_COX_N) {
  X_cox_os <- expr_log2_rep[analysis1_df$sample_id[cox_keep_os], complete_cols, drop = FALSE]
  res_os <- run_cv_cox(X_cox_os, os_time[cox_keep_os], os_event[cox_keep_os], "cv_14_cox_OS")
  print(res_os)
  cv_cox_results[["OS"]] <- res_os
  safe_write_csv(res_os, file.path(OUT_TBL, "cv_14_cox_OS_results.csv"))
}

op <- options(warn = -1)
pfs_time  <- as.numeric(analysis1_df$pfs_months_main)
pfs_event <- as.integer(analysis1_df$pfs_event_main)
options(op)
cox_keep_pfs <- !is.na(pfs_time) & !is.na(pfs_event)
if (sum(cox_keep_pfs) >= MIN_COX_N) {
  X_cox_pfs <- expr_log2_rep[analysis1_df$sample_id[cox_keep_pfs], complete_cols, drop = FALSE]
  res_pfs <- run_cv_cox(X_cox_pfs, pfs_time[cox_keep_pfs], pfs_event[cox_keep_pfs], "cv_15_cox_PFS")
  print(res_pfs)
  cv_cox_results[["PFS"]] <- res_pfs
  safe_write_csv(res_pfs, file.path(OUT_TBL, "cv_15_cox_PFS_results.csv"))
}

cv_cox_all <- bind_rows(cv_cox_results)
safe_write_csv(cv_cox_all, file.path(OUT_TBL, "cv_cox_all_comparisons.csv"))
if (nrow(cv_cox_all) > 0)
  plot_cv_auc_summary(cv_cox_all, "cv_cindex_mean", "cv_cindex_sd",
                      paste0("Cross-validated C-index (", K_FOLDS, "-fold x ", N_REPEATS, " repeats)"),
                      "CV C-index (mean +/- SD)", file.path(OUT_FIG, "cv_cox_cindex_summary"))

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 08e | Saving checkpoint_08e.RData ===")
save(cv_binomial_all, cv_cox_all, file = file.path(paths$rdata, "checkpoint_08e.RData"))
message("08e_cv_validation.R complete")
