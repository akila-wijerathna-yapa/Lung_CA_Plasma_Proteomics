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

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_04.RData"))
load(file.path(paths$rdata, "checkpoint_05.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 04, 05, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Output directories ────────────────────────────────────────
OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning", "stabl")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning", "stabl")
dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ── STABL constants ───────────────────────────────────────────
N_BOOT         <- 200    # number of subsampling iterations (200 for stable, reproducible results)
SUBSAMPLE_RATE <- 0.60   # proportion of samples per bootstrap
FDP_CUTOFF     <- 0.20   # maximum allowed false discovery proportion
MIN_CLASS_N    <- 5      # minimum samples per class for binomial
MIN_PAIRS      <- 8      # minimum pairs for paired comparison
MIN_COX_N      <- 35     # minimum samples for Cox
MIN_COX_EVENTS <- 12     # minimum events for Cox
RIDGE_ALPHA    <- 0.0    # alpha for initial ridge (step 1 of ALasso)
LASSO_ALPHA    <- 1.0    # alpha for ALasso (step 2)
WEIGHT_FLOOR   <- 1e-4   # floor for ALasso weights (avoids 1/0)
N_LAMBDA       <- 50     # number of lambda values in glmnet path
RNG_SEED       <- 42
SKIP_IF_DONE   <- FALSE  # FALSE = overwrite all existing results for clean re-run

# ── Helper: safe write ────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# ── Helper: sanitize protein names for use in file paths ──────
# Removes characters illegal in Windows filenames (/ \ : * ? " < > |)
sanitize_for_path <- function(x) {
  out <- gsub('[/\\\\:*?"<>|]', "_", as.character(x))
  out <- gsub("\\s+", "_", out)
  out <- gsub("_+", "_", out)
  substr(out, 1, 60)
}

# ─────────────────────────────────────────────────────────────
# STABL FUNCTION LIBRARY (all local, self-contained)
# ─────────────────────────────────────────────────────────────

# ── Step 1: Adaptive LASSO weights ───────────────────────────
# Adaptive LASSO penalises each feature differently.
# Features with small initial coefficients get a LARGE penalty
# (harder to select); features with large initial coefficients
# get a SMALL penalty (easier to keep).
# Weight for feature j = 1 / |β_ridge_j|, floored at WEIGHT_FLOOR
# so no weight explodes to infinity.
#
# Two-step ALasso (Zou 2006):
#   Step 1a: Fit ridge (alpha=0) on full data → β̂_ridge
#   Step 1b: Weight_j = 1 / max(|β̂_ridge_j|, WEIGHT_FLOOR)
get_alasso_weights <- function(X, y, family = "binomial") {
  # Standardise X so ridge coefficients are on the same scale
  X_std <- scale(X)
  # Zero-variance columns produce NaN after scaling (0/0) — set to 0
  X_std[is.nan(X_std)] <- 0
  X_std[is.nan(X_std)] <- 0  # zero-variance columns -> 0 (avoids glmnet NA error)

  # Ridge regression (alpha = 0): shrinks all coefficients but
  # keeps them all non-zero — gives us a starting estimate for
  # each feature's importance.
  op <- options(warn = -1)
  ridge_fit <- glmnet::cv.glmnet(
    x      = X_std,
    y      = y,
    family = family,
    alpha  = RIDGE_ALPHA,
    nfolds = min(5, length(y))
  )
  options(op)

  # Extract ridge coefficients at optimal lambda
  ridge_coef <- as.numeric(
    coef(ridge_fit, s = "lambda.min")[-1]  # drop intercept
  )

  # Adaptive weights: small |β| → large weight → harder to select
  weights <- 1 / pmax(abs(ridge_coef), WEIGHT_FLOOR)
  names(weights) <- colnames(X)
  weights
}

# Same for Cox PH model
get_alasso_weights_cox <- function(X, surv_obj) {
  X_std <- scale(X)
  X_std[is.nan(X_std)] <- 0
  op <- options(warn = -1)
  ridge_fit <- glmnet::cv.glmnet(
    x      = X_std,
    y      = surv_obj,
    family = "cox",
    alpha  = RIDGE_ALPHA,
    nfolds = min(5, nrow(X))
  )
  options(op)
  ridge_coef <- as.numeric(coef(ridge_fit, s = "lambda.min"))
  weights    <- 1 / pmax(abs(ridge_coef), WEIGHT_FLOOR)
  names(weights) <- colnames(X)
  weights
}

# ── Step 2a: Inject ghost features ───────────────────────────
# Ghost features = column-wise random permutation of the real
# feature matrix. They have the same marginal distribution as
# real features but are UNINFORMATIVE (labels shuffled per column).
# By tracking how often ghosts get selected, we estimate the
# false discovery rate of uninformative features.
inject_ghost_features <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  # Permute each column independently
  ghost <- apply(X, 2, function(col) col[sample(n)])
  colnames(ghost) <- paste0("__ghost__", colnames(X))
  cbind(X, ghost)  # [real | ghost]: 2p columns total
}

# ── Step 2b–3: Subsample loop + selection frequencies ─────────
# Core STABL loop (Hédou et al. 2024, Algorithm 1):
# For each bootstrap iteration:
#   1. Draw SUBSAMPLE_RATE fraction of samples
#   2. Inject ghost features into the subsample
#   3. Fit ALasso with per-feature weights on [real|ghost]
#   4. Record which features (real or ghost) are selected at any λ
# Output: selection frequency matrix (features × λ_index)
run_stabl_loop_binomial <- function(X, y, weights, label = "") {
  n     <- nrow(X)
  p     <- ncol(X)
  n_sub <- max(ceiling(SUBSAMPLE_RATE * n), MIN_CLASS_N * 2)

  # Accumulate selection counts: one row per feature (real + ghost)
  feature_names <- c(colnames(X), paste0("__ghost__", colnames(X)))
  sel_count     <- setNames(rep(0, 2 * p), feature_names)

  set.seed(RNG_SEED)
  for (b in seq_len(N_BOOT)) {
    # Stratified subsample: maintain class balance
    class_0 <- which(y == 0); class_1 <- which(y == 1)
    n0 <- max(MIN_CLASS_N, round(n_sub * length(class_0) / n))
    n1 <- max(MIN_CLASS_N, round(n_sub * length(class_1) / n))
    n0 <- min(n0, length(class_0)); n1 <- min(n1, length(class_1))
    idx <- c(sample(class_0, n0), sample(class_1, n1))

    X_sub <- X[idx, , drop = FALSE]
    y_sub <- y[idx]

    # Inject ghost features into this subsample
    X_aug <- inject_ghost_features(X_sub)  # n_sub × 2p

    # Double the weights vector (ghost features get same weights as their real twin)
    w_aug <- c(weights, weights)

    # Standardise augmented matrix
    op <- options(warn = -1)
    X_aug_std <- scale(X_aug)
    # Replace NaN columns (zero variance) with 0
    X_aug_std[is.nan(X_aug_std)] <- 0

    # Fit ALasso: penalty.factor = w_aug overrides the default equal
    # penalisation so each feature has its own penalty strength
    fit <- tryCatch(
      glmnet::glmnet(
        x              = X_aug_std,
        y              = y_sub,
        family         = "binomial",
        alpha          = LASSO_ALPHA,
        penalty.factor = w_aug,
        nlambda        = N_LAMBDA,
        standardize    = FALSE  # already standardised
      ),
      error = function(e) NULL
    )
    options(op)
    if (is.null(fit)) next

    # A feature is "selected" if it has a non-zero coefficient
    # at ANY lambda on the path (we track maximum selection)
    coef_mat  <- as.matrix(coef(fit))[-1, , drop = FALSE]  # drop intercept
    selected  <- rownames(coef_mat)[rowSums(coef_mat != 0) > 0]
    sel_count[intersect(selected, names(sel_count))] <-
      sel_count[intersect(selected, names(sel_count))] + 1

    if (b %% 25 == 0)
      message("    Bootstrap ", b, "/", N_BOOT, " | ", label)
  }

  # Convert counts to frequencies
  sel_freq <- sel_count / N_BOOT
  list(
    real_freq  = sel_freq[seq_len(p)],
    ghost_freq = sel_freq[seq_len(p) + p]
  )
}

# ── Cox PH version of the subsample loop ──────────────────────
run_stabl_loop_cox <- function(X, surv_obj, weights, label = "") {
  n     <- nrow(X)
  p     <- ncol(X)
  n_sub <- max(ceiling(SUBSAMPLE_RATE * n), 15L)

  feature_names <- c(colnames(X), paste0("__ghost__", colnames(X)))
  sel_count     <- setNames(rep(0, 2 * p), feature_names)

  set.seed(RNG_SEED)
  for (b in seq_len(N_BOOT)) {
    idx   <- sort(sample(seq_len(n), n_sub))
    X_sub <- X[idx, , drop = FALSE]
    s_sub <- surv_obj[idx]

    # Only proceed if there are sufficient events in this subsample
    if (sum(s_sub[, "status"]) < 5) next

    X_aug     <- inject_ghost_features(X_sub)
    w_aug     <- c(weights, weights)
    op        <- options(warn = -1)
    X_aug_std <- scale(X_aug)
    X_aug_std[is.nan(X_aug_std)] <- 0

    fit <- tryCatch(
      glmnet::glmnet(
        x              = X_aug_std,
        y              = s_sub,
        family         = "cox",
        alpha          = LASSO_ALPHA,
        penalty.factor = w_aug,
        nlambda        = N_LAMBDA,
        standardize    = FALSE
      ),
      error = function(e) NULL
    )
    options(op)
    if (is.null(fit)) next

    coef_mat <- as.matrix(coef(fit))
    selected <- rownames(coef_mat)[rowSums(coef_mat != 0) > 0]
    sel_count[intersect(selected, names(sel_count))] <-
      sel_count[intersect(selected, names(sel_count))] + 1

    if (b %% 25 == 0)
      message("    Bootstrap ", b, "/", N_BOOT, " | ", label)
  }

  sel_freq <- sel_count / N_BOOT
  list(
    real_freq  = sel_freq[seq_len(p)],
    ghost_freq = sel_freq[seq_len(p) + p]
  )
}

# ── Paired version: subsample PATIENT PAIRS ───────────────────
# For pre/post comparisons, we must subsample whole pairs to
# avoid a pre sample from patient A being compared to a post
# sample from patient B (data leakage on the pairing structure).
run_stabl_loop_paired <- function(X_pre, X_post, weights, label = "") {
  n_pairs <- nrow(X_pre)  # one row per patient
  p       <- ncol(X_pre)
  n_sub   <- max(ceiling(SUBSAMPLE_RATE * n_pairs), MIN_PAIRS)

  feature_names <- c(colnames(X_pre), paste0("__ghost__", colnames(X_pre)))
  sel_count     <- setNames(rep(0, 2 * p), feature_names)

  set.seed(RNG_SEED)
  for (b in seq_len(N_BOOT)) {
    # Sample PAIRS (patient indices), then stack pre + post rows
    pair_idx <- sort(sample(seq_len(n_pairs), min(n_sub, n_pairs)))
    X_sub    <- rbind(X_pre[pair_idx, , drop = FALSE],
                      X_post[pair_idx, , drop = FALSE])
    # Outcome: 0 = pre, 1 = post (for each patient in the pair)
    y_sub    <- c(rep(0L, length(pair_idx)), rep(1L, length(pair_idx)))

    X_aug     <- inject_ghost_features(X_sub)
    w_aug     <- c(weights, weights)
    op        <- options(warn = -1)
    X_aug_std <- scale(X_aug)
    X_aug_std[is.nan(X_aug_std)] <- 0

    fit <- tryCatch(
      glmnet::glmnet(
        x              = X_aug_std,
        y              = y_sub,
        family         = "binomial",
        alpha          = LASSO_ALPHA,
        penalty.factor = w_aug,
        nlambda        = N_LAMBDA,
        standardize    = FALSE
      ),
      error = function(e) NULL
    )
    options(op)
    if (is.null(fit)) next

    coef_mat <- as.matrix(coef(fit))[-1, , drop = FALSE]
    selected <- rownames(coef_mat)[rowSums(coef_mat != 0) > 0]
    sel_count[intersect(selected, names(sel_count))] <-
      sel_count[intersect(selected, names(sel_count))] + 1

    if (b %% 25 == 0)
      message("    Bootstrap ", b, "/", N_BOOT, " | ", label)
  }

  sel_freq <- sel_count / N_BOOT
  list(
    real_freq  = sel_freq[seq_len(p)],
    ghost_freq = sel_freq[seq_len(p) + p]
  )
}

# ── Step 3: FDP curve + data-driven threshold ─────────────────
# For each candidate threshold τ ∈ [0, 1], estimate the
# false discovery proportion (FDP):
#
#   FDP_hat(τ) = #{ghost features with freq > τ}
#                ─────────────────────────────────
#                max(1, #{real features with freq > τ})
#
# This estimates the proportion of "discoveries" (selected features)
# that are actually uninformative (like the ghost features).
#
# We choose τ* = the HIGHEST τ such that FDP_hat(τ) ≤ FDP_CUTOFF.
# Higher τ = more stringent selection (fewer but more reliable features).
compute_fdp_curve <- function(real_freq, ghost_freq) {
  tau_grid <- seq(0.05, 0.99, by = 0.01)

  fdp_tbl <- tibble::tibble(
    tau         = tau_grid,
    n_real_sel  = vapply(tau_grid, function(t) sum(real_freq  > t), integer(1)),
    n_ghost_sel = vapply(tau_grid, function(t) sum(ghost_freq > t), integer(1)),
    fdp_hat     = vapply(tau_grid, function(t) {
      n_ghost <- sum(ghost_freq > t)
      n_real  <- sum(real_freq  > t)
      n_ghost / max(1, n_real)   # FDP estimate
    }, numeric(1))
  )

  # Optimal threshold: highest τ where FDP ≤ allowed cutoff
  feasible <- fdp_tbl %>% dplyr::filter(fdp_hat <= FDP_CUTOFF, n_real_sel > 0)
  tau_star <- if (nrow(feasible)) max(feasible$tau) else NA_real_

  list(fdp_curve = fdp_tbl, tau_star = tau_star)
}

# ── Step 4: Final model + AUC ─────────────────────────────────
# Once features are selected by STABL, we fit a simple logistic
# or Cox regression on those features only (no penalisation) to:
# (a) get interpretable coefficients, and (b) compute AUC/C-index.
fit_final_model_binomial <- function(X_full, y_full, selected_features,
                                     feature_meta) {
  if (!length(selected_features))
    return(list(model = NULL, auc = NA_real_, summary = tibble()))

  X_sel <- scale(X_full[, selected_features, drop = FALSE])
  X_sel[is.nan(X_sel)] <- 0
  fit   <- tryCatch(
    glm(y_full ~ X_sel, family = binomial),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(glm(y_full ~ X_sel, family = binomial))
  )
  if (is.null(fit)) return(list(model = NULL, auc = NA_real_, summary = tibble()))

  op       <- options(warn = -1)
  prob     <- predict(fit, type = "response")
  roc_obj  <- pROC::roc(y_full, prob, quiet = TRUE, direction = "<")
  auc_val  <- as.numeric(pROC::auc(roc_obj))
  options(op)

  coef_df <- broom::tidy(fit) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      nELISA_ID = selected_features[match(
        stringr::str_remove(term, "^X_sel"),
        paste0("X_sel", seq_along(selected_features) - 1)
      )],
      nELISA_ID = selected_features
    ) %>%
    dplyr::left_join(
      feature_meta %>% dplyr::select(nELISA_ID, Protein_Label, Gene_Symbol),
      by = "nELISA_ID"
    ) %>%
    dplyr::mutate(final_model_auc = auc_val) %>%
    dplyr::select(nELISA_ID, Protein_Label, Gene_Symbol,
                  coef = estimate, std_error = std.error,
                  z = statistic, p_value = p.value, final_model_auc)

  list(model = fit, auc = auc_val, summary = coef_df, roc = roc_obj)
}

# ── Save STABL plots ──────────────────────────────────────────
save_stabl_plots <- function(real_freq, fdp_curve, tau_star, selected,
                              label, fig_dir, roc_obj = NULL) {
  # 1. Stability path: bar chart of selection frequencies
  freq_df <- tibble::tibble(
    nELISA_ID = names(real_freq), freq = real_freq,
    selected  = nELISA_ID %in% selected
  ) %>%
    left_join(feature_meta %>% select(nELISA_ID, Protein_Label), by = "nELISA_ID") %>%
    dplyr::arrange(desc(freq)) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::mutate(label_order = factor(Protein_Label, levels = rev(Protein_Label)))

  p_stab <- ggplot(freq_df, aes(freq, label_order, fill = selected)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = dplyr::coalesce(tau_star, 0.6),
               linetype = "dashed", colour = "#B2182B", linewidth = 0.8) +
    scale_fill_manual(values = c("TRUE" = "#1565C0", "FALSE" = "#90A4AE"),
                      labels = c("TRUE" = "STABL selected", "FALSE" = "Not selected")) +
    labs(title    = paste0("STABL — ", gsub("_", " ", label)),
         subtitle = paste0("Data-driven threshold τ* = ",
                           ifelse(is.na(tau_star), "not found (FDP control failed)",
                                  round(tau_star, 3)),
                           " | FDP ≤ ", FDP_CUTOFF,
                           " | N selected = ", length(selected)),
         x = paste0("Selection frequency (", N_BOOT, " bootstraps)"),
         y = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  cairo_pdf(file.path(fig_dir, paste0(label, "_stability_path.pdf")),
            width = 9, height = 7); print(p_stab); dev.off()
  png(file.path(fig_dir, paste0(label, "_stability_path.png")),
      width = 9, height = 7, units = "in", res = 300, bg = "white")
  print(p_stab); dev.off()

  # 2. FDP curve
  p_fdp <- ggplot(fdp_curve, aes(tau, fdp_hat)) +
    geom_line(colour = "#2166AC", linewidth = 1.1) +
    geom_hline(yintercept = FDP_CUTOFF, linetype = "dashed",
               colour = "#B2182B", linewidth = 0.8) +
    { if (!is.na(tau_star))
        geom_vline(xintercept = tau_star, linetype = "dotted",
                   colour = "#4d4d4d", linewidth = 0.7) } +
    annotate("text", x = 0.05, y = FDP_CUTOFF + 0.01,
             label = paste0("FDP cutoff = ", FDP_CUTOFF),
             hjust = 0, size = 3.5, colour = "#B2182B") +
    { if (!is.na(tau_star))
        annotate("text", x = tau_star + 0.01, y = 0.5,
                 label = paste0("τ* = ", round(tau_star, 3)),
                 hjust = 0, size = 3.5, colour = "#4d4d4d") } +
    labs(title    = paste0("STABL FDP curve — ", gsub("_", " ", label)),
         subtitle = "Estimated false discovery proportion vs stability threshold",
         x = "Stability threshold τ", y = "Estimated FDP") +
    ylim(0, 1) + theme_minimal(base_size = 11)

  cairo_pdf(file.path(fig_dir, paste0(label, "_fdp_curve.pdf")),
            width = 7.5, height = 5); print(p_fdp); dev.off()
  png(file.path(fig_dir, paste0(label, "_fdp_curve.png")),
      width = 7.5, height = 5, units = "in", res = 300, bg = "white")
  print(p_fdp); dev.off()

  # 3. ROC curve (binomial only)
  if (!is.null(roc_obj)) {
    roc_df  <- tibble::tibble(
      fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
      geom_path(colour = "#2166AC", linewidth = 1.1) +
      geom_abline(intercept = 0, slope = 1,
                  linetype = "dashed", colour = "grey60") +
      labs(title    = paste0("STABL ROC — ", gsub("_", " ", label)),
           subtitle = paste0("Final model AUC = ", round(auc_val, 3),
                             " | N features = ", length(selected)),
           x = "False positive rate", y = "True positive rate") +
      theme_minimal(base_size = 11)
    cairo_pdf(file.path(fig_dir, paste0(label, "_roc.pdf")),
              width = 6.5, height = 5); print(p_roc); dev.off()
    png(file.path(fig_dir, paste0(label, "_roc.png")),
        width = 6.5, height = 5, units = "in", res = 300, bg = "white")
    print(p_roc); dev.off()
  }
}

# ── Master pipeline: binomial ─────────────────────────────────
run_stabl_binomial <- function(X_full, y_full, label,
                                feature_meta,
                                X_pre = NULL, X_post = NULL,
                                is_paired = FALSE) {
  # Skip-if-done check
  done_path <- file.path(OUT_TBL, paste0(label, "_stabl_selected_features.csv"))
  if (SKIP_IF_DONE && file.exists(done_path)) {
    message("  SKIP (already done): ", label); return(invisible(NULL))
  }
  # Feasibility check
  n0 <- sum(y_full == 0); n1 <- sum(y_full == 1)
  if (!is_paired && (n0 < MIN_CLASS_N || n1 < MIN_CLASS_N)) {
    message("  SKIP (small class): ", label,
            " — n0=", n0, ", n1=", n1); return(invisible(NULL))
  }
  # Remove features with any NA (STABL requires complete data)
  if (!is_paired) {
    X_use <- X_full[, colSums(is.na(X_full)) == 0, drop = FALSE]
  } else {
    X_use <- X_full[, colSums(is.na(X_pre)) == 0 &
                       colSums(is.na(X_post)) == 0, drop = FALSE]
    X_pre_use  <- X_pre[, colnames(X_use), drop = FALSE]
    X_post_use <- X_post[, colnames(X_use), drop = FALSE]
  }
  p_complete <- ncol(if (is_paired) X_use else X_full[, colSums(is.na(X_full)) == 0, drop = FALSE])
  if (p_complete < 20) {
    message("  SKIP (<20 complete features): ", label); return(invisible(NULL))
  }

  message("\n─── STABL binomial | ", label, " ───")
  message("  n = ", if (!is_paired) nrow(X_use) else nrow(X_pre), " | p = ", p_complete,
          if (is_paired) paste0(" | ", nrow(X_pre), " pairs") else paste0(" | n0=", n0, " n1=", n1))

  # Step 1: ALasso weights from full data
  if (!is_paired) {
    message("  Step 1: computing adaptive LASSO weights ...")
    wts <- get_alasso_weights(
      X_full[, colnames(X_use), drop = FALSE], y_full)
  } else {
    # For paired: use stacked pre+post to compute weights
    X_stack <- rbind(X_pre_use, X_post_use)
    y_stack <- c(rep(0L, nrow(X_pre_use)), rep(1L, nrow(X_post_use)))
    message("  Step 1: computing adaptive LASSO weights (paired) ...")
    wts <- get_alasso_weights(X_stack, y_stack)
  }

  # Step 2–3: Bootstrap + selection frequencies + FDP
  message("  Step 2–3: STABL bootstrap loop + FDP curve ...")
  if (!is_paired) {
    freq_res <- run_stabl_loop_binomial(
      X_full[, colnames(X_use), drop = FALSE], y_full, wts, label)
  } else {
    freq_res <- run_stabl_loop_paired(X_pre_use, X_post_use, wts, label)
  }
  fdp_res  <- compute_fdp_curve(freq_res$real_freq, freq_res$ghost_freq)

  tau_star  <- fdp_res$tau_star
  selected  <- if (!is.na(tau_star))
    names(freq_res$real_freq)[freq_res$real_freq > tau_star] else character(0)

  message("  τ* = ", ifelse(is.na(tau_star), "not found", round(tau_star, 3)),
          " | selected = ", length(selected))

  # Step 4: Final logistic model
  if (!is_paired) {
    message("  Step 4: fitting final logistic model ...")
    final <- fit_final_model_binomial(
      X_full[, colnames(X_use), drop = FALSE], y_full, selected, feature_meta)
    message("  Final model AUC = ", round(dplyr::coalesce(final$auc, NA_real_), 3))
    roc_for_plot <- final$roc
  } else {
    X_stack <- rbind(X_pre_use, X_post_use)
    y_stack <- c(rep(0L, nrow(X_pre_use)), rep(1L, nrow(X_post_use)))
    final   <- fit_final_model_binomial(X_stack, y_stack, selected, feature_meta)
    message("  Final model AUC (pre vs post) = ",
            round(dplyr::coalesce(final$auc, NA_real_), 3))
    roc_for_plot <- final$roc
  }

  # Build selection frequency table with annotations
  freq_tbl <- tibble::tibble(
    nELISA_ID         = names(freq_res$real_freq),
    real_freq         = freq_res$real_freq,
    ghost_freq        = freq_res$ghost_freq,
    tau_star          = tau_star,
    stabl_selected    = nELISA_ID %in% selected,
    final_model_auc   = dplyr::coalesce(final$auc, NA_real_)
  ) %>%
    dplyr::left_join(feature_meta %>% dplyr::select(nELISA_ID, Protein_Label,
                                                     Gene_Symbol, Protein_Name),
                     by = "nELISA_ID") %>%
    dplyr::arrange(desc(real_freq))

  # Save tables
  safe_write_csv(freq_tbl,
    file.path(OUT_TBL, paste0(label, "_stabl_selection_freq.csv")))
  safe_write_csv(fdp_res$fdp_curve,
    file.path(OUT_TBL, paste0(label, "_stabl_fdp_curve.csv")))
  safe_write_csv(freq_tbl %>% dplyr::filter(stabl_selected),
    file.path(OUT_TBL, paste0(label, "_stabl_selected_features.csv")))
  if (nrow(final$summary))
    safe_write_csv(final$summary,
      file.path(OUT_TBL, paste0(label, "_stabl_final_model.csv")))

  # Save plots
  save_stabl_plots(freq_res$real_freq, fdp_res$fdp_curve,
                   tau_star, selected, label, OUT_FIG, roc_for_plot)

  invisible(freq_tbl)
}

# ── Master pipeline: Cox PH ───────────────────────────────────
run_stabl_cox <- function(X_full, surv_obj, time_col, event_col,
                           sample_df, label, feature_meta) {
  done_path <- file.path(OUT_TBL, paste0(label, "_stabl_selected_features.csv"))
  if (SKIP_IF_DONE && file.exists(done_path)) {
    message("  SKIP (already done): ", label); return(invisible(NULL))
  }
  # Event check
  n_events <- sum(surv_obj[, "status"] == 1, na.rm = TRUE)
  if (nrow(X_full) < MIN_COX_N || n_events < MIN_COX_EVENTS) {
    message("  SKIP (Cox threshold): ", label,
            " — n=", nrow(X_full), " events=", n_events)
    return(invisible(NULL))
  }
  X_use <- X_full[, colSums(is.na(X_full)) == 0, drop = FALSE]
  if (ncol(X_use) < 20) {
    message("  SKIP (<20 complete features): ", label); return(invisible(NULL))
  }

  message("\n─── STABL Cox | ", label, " ───")
  message("  n = ", nrow(X_use), " | events = ", n_events, " | p = ", ncol(X_use))

  message("  Step 1: ALasso weights ...")
  wts <- get_alasso_weights_cox(X_use, surv_obj)

  message("  Step 2–3: bootstrap + FDP curve ...")
  freq_res <- run_stabl_loop_cox(X_use, surv_obj, wts, label)
  fdp_res  <- compute_fdp_curve(freq_res$real_freq, freq_res$ghost_freq)

  tau_star <- fdp_res$tau_star
  selected <- if (!is.na(tau_star))
    names(freq_res$real_freq)[freq_res$real_freq > tau_star] else character(0)
  message("  τ* = ", ifelse(is.na(tau_star), "not found", round(tau_star, 3)),
          " | selected = ", length(selected))

  # Final Cox model on selected features
  c_index <- NA_real_
  final_cox_tbl <- tibble()
  if (length(selected)) {
    X_sel <- scale(X_use[, selected, drop = FALSE])
    X_sel[is.nan(X_sel)] <- 0
    fit_cox <- tryCatch(
      survival::coxph(surv_obj ~ X_sel),
      error = function(e) NULL)
    if (!is.null(fit_cox)) {
      c_index <- as.numeric(summary(fit_cox)$concordance["C"])
      message("  C-index = ", round(c_index, 3))
      coef_df <- broom::tidy(fit_cox) %>%
        dplyr::mutate(
          nELISA_ID = selected,
          HR        = exp(estimate),
          c_index   = c_index
        ) %>%
        dplyr::left_join(feature_meta %>% dplyr::select(nELISA_ID, Protein_Label, Gene_Symbol),
                         by = "nELISA_ID") %>%
        dplyr::select(nELISA_ID, Protein_Label, Gene_Symbol,
                      coef = estimate, HR,
                      std_error = std.error, z = statistic,
                      p_value = p.value, c_index)
      final_cox_tbl <- coef_df
    }
  }

  freq_tbl <- tibble::tibble(
    nELISA_ID       = names(freq_res$real_freq),
    real_freq       = freq_res$real_freq,
    ghost_freq      = freq_res$ghost_freq,
    tau_star        = tau_star,
    stabl_selected  = nELISA_ID %in% selected,
    c_index         = c_index
  ) %>%
    dplyr::left_join(feature_meta %>% dplyr::select(nELISA_ID, Protein_Label,
                                                     Gene_Symbol, Protein_Name),
                     by = "nELISA_ID") %>%
    dplyr::arrange(desc(real_freq))

  safe_write_csv(freq_tbl,
    file.path(OUT_TBL, paste0(label, "_stabl_selection_freq.csv")))
  safe_write_csv(fdp_res$fdp_curve,
    file.path(OUT_TBL, paste0(label, "_stabl_fdp_curve.csv")))
  safe_write_csv(freq_tbl %>% dplyr::filter(stabl_selected),
    file.path(OUT_TBL, paste0(label, "_stabl_selected_features.csv")))
  if (nrow(final_cox_tbl))
    safe_write_csv(final_cox_tbl,
      file.path(OUT_TBL, paste0(label, "_stabl_final_cox_model.csv")))

  save_stabl_plots(freq_res$real_freq, fdp_res$fdp_curve,
                   tau_star, selected, label, OUT_FIG, roc_obj = NULL)
  invisible(freq_tbl)
}

# ─────────────────────────────────────────────────────────────
# BUILD COMPARISON SPECIFICATIONS
# ─────────────────────────────────────────────────────────────
message("\n=== 08b | Building comparison specifications ===")

# Convenience: complete-feature matrix for all matched baseline samples
# Only columns with zero NA are used (STABL requires complete data)
complete_cols <- colnames(expr_log2_rep)[
  colSums(is.na(expr_log2_rep[analysis1_df$sample_id, , drop = FALSE])) == 0]
message("Complete features available for STABL: ", length(complete_cols))

# ── Response comparisons (comparisons 1–5) ────────────────────
response_specs <- list(
  stabl_01_CR_vs_PD = list(
    ids    = analysis1_df %>% dplyr::filter(response_4cat %in% c("CR","PD")) %>% pull(sample_id),
    y_vals = analysis1_df %>% dplyr::filter(response_4cat %in% c("CR","PD")) %>%
             dplyr::mutate(y = as.integer(response_4cat == "CR")) %>% pull(y),
    note   = "CR=1 vs PD=0"
  ),
  stabl_02_CRPR_vs_PD = list(
    ids    = analysis1_df %>% dplyr::filter(response_4cat %in% c("CR","PR","PD")) %>% pull(sample_id),
    y_vals = analysis1_df %>% dplyr::filter(response_4cat %in% c("CR","PR","PD")) %>%
             dplyr::mutate(y = as.integer(response_4cat %in% c("CR","PR"))) %>% pull(y),
    note   = "(CR+PR)=1 vs PD=0"
  ),
  stabl_03_CRPR_vs_SDPD = list(
    ids    = analysis1_df %>% dplyr::filter(!is.na(response_binary)) %>% pull(sample_id),
    y_vals = analysis1_df %>% dplyr::filter(!is.na(response_binary)) %>%
             dplyr::mutate(y = as.integer(response_binary == "PR/CR")) %>% pull(y),
    note   = "(CR+PR)=1 vs (SD+PD)=0"
  ),
  stabl_04_CRPRSD_vs_PD = list(
    ids    = analysis1_df %>% dplyr::filter(response_4cat %in% c("CR","PR","SD","PD")) %>% pull(sample_id),
    y_vals = analysis1_df %>% dplyr::filter(response_4cat %in% c("CR","PR","SD","PD")) %>%
             dplyr::mutate(y = as.integer(response_4cat != "PD")) %>% pull(y),
    note   = "(CR+PR+SD)=1 vs PD=0"
  ),
  stabl_05_PR_vs_PD = list(
    ids    = analysis1_df %>% dplyr::filter(response_4cat %in% c("PR","PD")) %>% pull(sample_id),
    y_vals = analysis1_df %>% dplyr::filter(response_4cat %in% c("PR","PD")) %>%
             dplyr::mutate(y = as.integer(response_4cat == "PR")) %>% pull(y),
    note   = "PR=1 vs PD=0"
  )
)

# ── Recurrence comparisons (comparisons 6–11) ─────────────────
recurrence_base <- analysis1_df %>%
  dplyr::filter(!is.na(recurrence_binary_sensitivity))

recurrence_specs <- list(
  stabl_06_CTx_recurrence = list(
    df   = recurrence_base %>% dplyr::filter(new_group == "CTx"),
    note = "CTx: Recurrence=1 vs No recurrence=0"
  ),
  stabl_07_Sx_recurrence = list(
    df   = recurrence_base %>% dplyr::filter(new_group == "Sx"),
    note = "Sx: Recurrence=1 vs No recurrence=0"
  ),
  stabl_08_IO_recurrence = list(
    df   = recurrence_base %>% dplyr::filter(new_group == "IO"),
    note = "IO: Recurrence=1 vs No recurrence=0"
  ),
  stabl_09_CTxIO_recurrence = list(
    df   = recurrence_base %>% dplyr::filter(new_group == "CTx + IO"),
    note = "CTx+IO: Recurrence=1 vs No recurrence=0"
  )
  # CTx+RTx (#10) and TKI (#11) skipped: n_no_recurrence < MIN_CLASS_N
)

# ── Paired comparisons (comparison 12) ───────────────────────
# Pre vs post Surgery (n=13 pairs; comparison 13 CTx+IO skipped: only 4 pairs)
sx_pre_ids  <- paired_pre  %>% dplyr::filter(patient_id %in% surgery_pair_ids) %>%
               dplyr::arrange(patient_id) %>% pull(sample_id)
sx_post_ids <- paired_post %>% dplyr::filter(patient_id %in% surgery_pair_ids) %>%
               dplyr::arrange(patient_id) %>% pull(sample_id)

# ─────────────────────────────────────────────────────────────
# RUN ALL COMPARISONS
# ─────────────────────────────────────────────────────────────

# ── Response: comparisons 1–5 ─────────────────────────────────
message("\n=== 08b | Response comparisons (1–5) ===")
for (nm in names(response_specs)) {
  spec  <- response_specs[[nm]]
  X_use <- expr_log2_rep[spec$ids, complete_cols, drop = FALSE]
  run_stabl_binomial(X_use, spec$y_vals, nm, feature_meta)
}

# ── Recurrence: comparisons 6–9 ──────────────────────────────
message("\n=== 08b | Recurrence comparisons (6–9) ===")
for (nm in names(recurrence_specs)) {
  spec   <- recurrence_specs[[nm]]
  df_sub <- spec$df
  y_sub  <- as.integer(df_sub$recurrence_binary_sensitivity == "Recurrence")
  X_use  <- expr_log2_rep[df_sub$sample_id, complete_cols, drop = FALSE]
  run_stabl_binomial(X_use, y_sub, nm, feature_meta)
}

# ── Paired pre/post Surgery: comparison 12 ───────────────────
message("\n=== 08b | Pre vs post Surgery (comparison 12) ===")
X_pre_sx  <- expr_log2_rep[sx_pre_ids,  complete_cols, drop = FALSE]
X_post_sx <- expr_log2_rep[sx_post_ids, complete_cols, drop = FALSE]
# Stack for weights + feasibility check; pass paired flag
X_stack_sx <- rbind(X_pre_sx, X_post_sx)
y_stack_sx <- c(rep(0L, nrow(X_pre_sx)), rep(1L, nrow(X_post_sx)))
run_stabl_binomial(X_stack_sx, y_stack_sx, "stabl_12_Sx_pre_vs_post", feature_meta,
                   X_pre = X_pre_sx, X_post = X_post_sx, is_paired = TRUE)

# ── Cox OS: comparison 14 ─────────────────────────────────────
message("\n=== 08b | STABL Cox OS (comparison 14) ===")
op <- options(warn = -1)
os_time  <- as.numeric(analysis1_df$os_months)
os_event <- as.integer(analysis1_df$os_event)
options(op)
cox_keep_os <- !is.na(os_time) & !is.na(os_event)
if (sum(cox_keep_os) >= MIN_COX_N) {
  X_cox_os <- expr_log2_rep[analysis1_df$sample_id[cox_keep_os], complete_cols, drop = FALSE]
  surv_os  <- survival::Surv(os_time[cox_keep_os], os_event[cox_keep_os])
  run_stabl_cox(X_cox_os, surv_os, "os_months", "os_event",
                analysis1_df[cox_keep_os, ], "stabl_14_cox_OS", feature_meta)
}

# ── Cox PFS: comparison 15 ────────────────────────────────────
message("\n=== 08b | STABL Cox PFS (comparison 15) ===")
op <- options(warn = -1)
pfs_time  <- as.numeric(analysis1_df$pfs_months_main)
pfs_event <- as.integer(analysis1_df$pfs_event_main)
options(op)
cox_keep_pfs <- !is.na(pfs_time) & !is.na(pfs_event)
if (sum(cox_keep_pfs) >= MIN_COX_N) {
  X_cox_pfs <- expr_log2_rep[analysis1_df$sample_id[cox_keep_pfs], complete_cols, drop = FALSE]
  surv_pfs  <- survival::Surv(pfs_time[cox_keep_pfs], pfs_event[cox_keep_pfs])
  run_stabl_cox(X_cox_pfs, surv_pfs, "pfs_months_main", "pfs_event_main",
                analysis1_df[cox_keep_pfs, ], "stabl_15_cox_PFS", feature_meta)
}

# ─────────────────────────────────────────────────────────────
# SUMMARY TABLE
# ─────────────────────────────────────────────────────────────
message("\n=== 08b | Building STABL summary ===")
sel_files <- list.files(OUT_TBL, pattern = "_stabl_selected_features\\.csv$",
                        full.names = TRUE)
stabl_summary <- bind_rows(lapply(sel_files, function(f) {
  tryCatch({
    df  <- readr::read_csv(f, show_col_types = FALSE)
    lbl <- stringr::str_remove(basename(f), "_stabl_selected_features\\.csv$")
    freq_file <- file.path(OUT_TBL, paste0(lbl, "_stabl_selection_freq.csv"))
    tau_star  <- if (file.exists(freq_file))
      dplyr::first(readr::read_csv(freq_file, show_col_types=FALSE, n_max=1)$tau_star)
      else NA_real_
    tibble(comparison = lbl, n_selected = nrow(df), tau_star = tau_star,
           top_protein = dplyr::first(df$Protein_Label))
  }, error = function(e) NULL)
}))
safe_write_csv(stabl_summary, file.path(OUT_TBL, "stabl_run_summary.csv"))
message("STABL summary:")
print(stabl_summary)

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 08b | Saving checkpoint_08b.RData ===")
save(stabl_summary, file = file.path(paths$rdata, "checkpoint_08b.RData"))
message("08b_stabl_biomarker.R complete")
message("Note: broom package required for tidy() — install if missing:")
message("  install.packages('broom')")
