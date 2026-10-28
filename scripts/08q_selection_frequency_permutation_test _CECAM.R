if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")

while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
RNG_SEED       <- 42
N_BOOT_MAIN    <- 200   # bootstraps for the REAL-label selection frequency
N_BOOT_PERM    <- 30    # bootstraps per permutation (reduced for speed — this is a null-calibration
                        # sample, not a final discovery estimate, so fewer boots per replicate is standard)
N_PERM         <- 200   # number of label permutations
EN_ALPHA       <- 0.5
CV_NFOLDS      <- 5     # reduced glmnet internal CV folds, for speed within the permutation loop
TARGET_PROTEIN <- "SHH"

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

message("\nEstimated fits: ", N_BOOT_MAIN + N_PERM * N_BOOT_PERM,
        " glmnet CV fits total. Reduce N_PERM / N_BOOT_PERM below if this is too slow.")

# ── Build data: response_binary (PR/CR vs SD/PD) ───────────────
message("\n=== 08q | Building response data (PR/CR vs SD/PD) ===")

target_id <- feature_meta %>%
  dplyr::filter(Gene_Symbol == TARGET_PROTEIN) %>%
  dplyr::pull(nELISA_ID) %>% unique()
if (length(target_id) != 1) stop("Target protein not found uniquely on panel: ", TARGET_PROTEIN)

resp_df <- analysis1_df %>%
  dplyr::filter(!is.na(response_binary))

complete_cols <- colnames(expr_log2_rep)[
  colSums(is.na(expr_log2_rep[resp_df$sample_id, , drop = FALSE])) == 0
]
if (!target_id %in% complete_cols) stop(TARGET_PROTEIN, " (", target_id, ") has missing values \u2014 cannot proceed.")

X <- expr_log2_rep[resp_df$sample_id, complete_cols, drop = FALSE]
y <- as.integer(resp_df$response_binary == "PR/CR")   # 1 = PR/CR, 0 = SD/PD
target_col_idx <- which(colnames(X) == target_id)

message("n = ", nrow(X), " (", sum(y == 1), " PR/CR, ", sum(y == 0), " SD/PD) | p = ", ncol(X), " proteins")

# ── Local function: bootstrap selection frequency for one label vector ──
# For each bootstrap resample (stratified by class to preserve balance),
# fit an Elastic Net logistic regression (cv.glmnet) and record whether
# the target protein has a non-zero coefficient at lambda.min.
bootstrap_selection_freq <- function(X, y, target_col_idx, n_boot, alpha, nfolds, seed_offset = 0) {
  n <- nrow(X)
  idx_pos <- which(y == 1); idx_neg <- which(y == 0)
  selected <- logical(n_boot)
  for (b in seq_len(n_boot)) {
    set.seed(RNG_SEED * 1000 + seed_offset + b)
    boot_idx <- c(sample(idx_pos, length(idx_pos), replace = TRUE),
                  sample(idx_neg, length(idx_neg), replace = TRUE))
    X_b <- X[boot_idx, , drop = FALSE]; y_b <- y[boot_idx]
    X_std <- scale(X_b); X_std[is.nan(X_std)] <- 0
    fit <- tryCatch(
      glmnet::cv.glmnet(X_std, y_b, family = "binomial", alpha = alpha, nfolds = nfolds),
      error = function(e) NULL
    )
    if (is.null(fit)) { selected[b] <- FALSE; next }
    coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]  # drop intercept
    selected[b] <- coefs[target_col_idx] != 0
  }
  mean(selected)
}

# ── Observed selection frequency (real labels) ──────────────────
message("\n=== 08q | Computing observed selection frequency (real labels, ", N_BOOT_MAIN, " bootstraps) ===")

observed_freq <- bootstrap_selection_freq(X, y, target_col_idx, N_BOOT_MAIN, EN_ALPHA, CV_NFOLDS, seed_offset = 0)
message("Observed ", TARGET_PROTEIN, " selection frequency (real labels): ", round(observed_freq, 3))

# ── Null distribution (permuted labels) ─────────────────────────
message("\n=== 08q | Building null distribution (", N_PERM, " permutations x ", N_BOOT_PERM, " bootstraps) ===")

null_freqs <- numeric(N_PERM)
for (p in seq_len(N_PERM)) {
  set.seed(RNG_SEED + p)
  y_perm <- sample(y)
  null_freqs[p] <- bootstrap_selection_freq(X, y_perm, target_col_idx, N_BOOT_PERM, EN_ALPHA, CV_NFOLDS,
                                             seed_offset = p * 10000)
  if (p %% 25 == 0) message("  permutation ", p, "/", N_PERM)
}

empirical_p <- (1 + sum(null_freqs >= observed_freq)) / (1 + N_PERM)

message("\nNull distribution summary: mean=", round(mean(null_freqs), 3),
        ", max=", round(max(null_freqs), 3), ", 95th pct=", round(quantile(null_freqs, 0.95), 3))
message("Empirical p-value (observed vs. permutation null): ", signif(empirical_p, 3))

result_summary <- tibble::tibble(
  protein = TARGET_PROTEIN, nELISA_ID = target_id,
  n_patients = nrow(X), n_prcr = sum(y == 1), n_sdpd = sum(y == 0),
  observed_selection_freq = observed_freq,
  null_mean = mean(null_freqs), null_95th_pct = quantile(null_freqs, 0.95), null_max = max(null_freqs),
  n_permutations = N_PERM, n_boot_main = N_BOOT_MAIN, n_boot_perm = N_BOOT_PERM,
  empirical_p_value = empirical_p
)
safe_write_csv(result_summary, file.path(OUT_TBL, "shh_selection_freq_permutation_test.csv"))
safe_write_csv(tibble::tibble(permutation = seq_len(N_PERM), null_selection_freq = null_freqs),
               file.path(OUT_TBL, "shh_selection_freq_permutation_null_distribution.csv"))

# ── Plot: null distribution with observed value marked ──────────
message("\n=== 08q | Building permutation test figure ===")

p_perm <- ggplot(tibble::tibble(null_selection_freq = null_freqs), aes(x = null_selection_freq)) +
  geom_histogram(bins = 25, fill = "#a6bddb", colour = "white") +
  geom_vline(xintercept = observed_freq, colour = "#c9772b", linewidth = 1.1, linetype = "dashed") +
  annotate("text", x = observed_freq, y = Inf, label = paste0("Observed = ", round(observed_freq, 3)),
           colour = "#c9772b", vjust = 1.5, hjust = -0.05, fontface = "bold", size = 3.6) +
  labs(
    title = paste0(TARGET_PROTEIN, " Selection Frequency: Permutation Null Test"),
    subtitle = paste0("Empirical p = ", signif(empirical_p, 3), " (", N_PERM, " label permutations, ",
                       N_BOOT_PERM, " bootstraps each) | Observed: ", N_BOOT_MAIN, " bootstraps, real labels"),
    x = "Bootstrap selection frequency (Elastic Net, non-zero coefficient)", y = "Count (permutations)"
  ) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10.5, colour = "grey30"))

stub <- file.path(OUT_FIG, "shh_selection_freq_permutation_test")
cairo_pdf(paste0(stub, ".pdf"), width = 8.5, height = 5.5); print(p_perm); dev.off()
png(paste0(stub, ".png"), width = 8.5, height = 5.5, units = "in", res = 300, bg = "white"); print(p_perm); dev.off()

message("Permutation test figure saved: ", stub, ".pdf / .png")
message("\n08q_selection_frequency_permutation_test.R complete")
