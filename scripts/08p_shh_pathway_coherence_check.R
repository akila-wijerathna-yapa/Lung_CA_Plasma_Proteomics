if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")

while (dev.cur() > 1) dev.off()

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

akila_volcano_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title        = element_text(face = "bold", size = 13.5),
      plot.subtitle     = element_text(size = 10.5, colour = "grey30")
    )
}

# ── Confirm which Hedgehog pathway members are actually on the panel ──
message("\n=== 08p | Checking Hedgehog pathway member coverage ===")

hh_candidates <- c("SHH", "GLI1", "GLI2", "GLI3", "PTCH1", "PTCH2", "SMO",
                   "HHIP", "SUFU", "DISP1", "CDON", "GAS1", "BOC")

hh_lookup <- feature_meta %>%
  dplyr::filter(Gene_Symbol %in% hh_candidates) %>%
  dplyr::distinct(Gene_Symbol, nELISA_ID)

message("Hedgehog-pathway-related proteins found on panel: ",
        paste(hh_lookup$Gene_Symbol, collapse = ", "))
message("Not on panel (expected \u2014 intracellular/membrane components not ",
        "typically plasma-detectable): ",
        paste(setdiff(hh_candidates, hh_lookup$Gene_Symbol), collapse = ", "))

shh_id  <- hh_lookup$nELISA_ID[hh_lookup$Gene_Symbol == "SHH"]
cdon_id <- hh_lookup$nELISA_ID[hh_lookup$Gene_Symbol == "CDON"]
gas1_id <- hh_lookup$nELISA_ID[hh_lookup$Gene_Symbol == "GAS1"]

if (length(shh_id) != 1) stop("SHH not found on panel — cannot proceed.")

# ── Build expression + response data frame ─────────────────────
message("\n=== 08p | Building SHH / CDON / GAS1 data frame ===")

get_expr <- function(id) {
  v <- expr_log2_rep[, id]; names(v) <- rownames(expr_log2_rep); v
}

shh_expr  <- get_expr(shh_id)
cdon_expr <- if (length(cdon_id) == 1) get_expr(cdon_id) else NULL
gas1_expr <- if (length(gas1_id) == 1) get_expr(gas1_id) else NULL

df <- analysis1_df %>%
  dplyr::mutate(
    SHH  = shh_expr[sample_id],
    CDON = if (!is.null(cdon_expr)) cdon_expr[sample_id] else NA_real_,
    GAS1 = if (!is.null(gas1_expr)) gas1_expr[sample_id] else NA_real_
  )

# ── Correlation: SHH vs. CDON / GAS1 (whole matched cohort) ────
message("\n=== 08p | SHH correlation with CDON / GAS1 ===")

run_corr <- function(x, y, label) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 10) return(tibble::tibble(pair = label, n = sum(keep), rho = NA_real_, p_value = NA_real_))
  ct <- cor.test(x[keep], y[keep], method = "spearman")
  tibble::tibble(pair = label, n = sum(keep), rho = unname(ct$estimate), p_value = ct$p.value)
}

corr_results <- dplyr::bind_rows(
  run_corr(df$SHH, df$CDON, "SHH vs CDON"),
  run_corr(df$SHH, df$GAS1, "SHH vs GAS1")
)
message("Correlation results:"); print(corr_results)

# ── Response association: CDON / GAS1 (same test as SHH/consensus) ──
message("\n=== 08p | CDON / GAS1 response association (Wilcoxon) ===")

run_response_wilcox <- function(expr_vec, label) {
  d <- analysis1_df %>%
    dplyr::mutate(expr = expr_vec[sample_id]) %>%
    dplyr::filter(!is.na(response_binary), !is.na(expr))
  n_sdpd <- sum(d$response_binary == "SD/PD"); n_prcr <- sum(d$response_binary == "PR/CR")
  if (n_sdpd < 8 || n_prcr < 8) {
    return(tibble::tibble(protein = label, n_sdpd = n_sdpd, n_prcr = n_prcr, wilcox_p = NA_real_))
  }
  wt <- wilcox.test(expr ~ response_binary, data = d)
  tibble::tibble(protein = label, n_sdpd = n_sdpd, n_prcr = n_prcr, wilcox_p = wt$p.value,
                 median_sdpd = median(d$expr[d$response_binary == "SD/PD"]),
                 median_prcr = median(d$expr[d$response_binary == "PR/CR"]))
}

response_results <- dplyr::bind_rows(
  run_response_wilcox(shh_expr,  "SHH"),
  if (!is.null(cdon_expr)) run_response_wilcox(cdon_expr, "CDON") else NULL,
  if (!is.null(gas1_expr)) run_response_wilcox(gas1_expr, "GAS1") else NULL
)
message("Response association results:"); print(response_results)

safe_write_csv(
  list(hedgehog_panel_coverage = hh_lookup, correlation = corr_results, response_association = response_results) %>%
    purrr::imap_dfr(~ dplyr::mutate(.x, table = .y)),
  file.path(OUT_TBL, "shh_pathway_coherence_check.csv")
)

# ── Combined figure: correlation scatter + response boxplots ───
message("\n=== 08p | Building coherence check figure ===")

p_scatter_cdon <- if (!is.null(cdon_expr)) {
  ggplot(df, aes(x = SHH, y = CDON)) +
    geom_point(alpha = 0.55, size = 1.6, colour = "#2f6fb0") +
    geom_smooth(method = "lm", se = TRUE, colour = "#c9772b", linewidth = 0.7) +
    labs(title = paste0("SHH vs CDON (Spearman \u03c1=",
                         signif(corr_results$rho[corr_results$pair == "SHH vs CDON"], 3),
                         ", p=", signif(corr_results$p_value[corr_results$pair == "SHH vs CDON"], 3), ")"),
         x = "SHH (log2)", y = "CDON (log2)") +
    akila_volcano_theme()
} else NULL

p_scatter_gas1 <- if (!is.null(gas1_expr)) {
  ggplot(df, aes(x = SHH, y = GAS1)) +
    geom_point(alpha = 0.55, size = 1.6, colour = "#2f9e5b") +
    geom_smooth(method = "lm", se = TRUE, colour = "#c9772b", linewidth = 0.7) +
    labs(title = paste0("SHH vs GAS1 (Spearman \u03c1=",
                         signif(corr_results$rho[corr_results$pair == "SHH vs GAS1"], 3),
                         ", p=", signif(corr_results$p_value[corr_results$pair == "SHH vs GAS1"], 3), ")"),
         x = "SHH (log2)", y = "GAS1 (log2)") +
    akila_volcano_theme()
} else NULL

resp_box_df <- dplyr::bind_rows(
  analysis1_df %>% dplyr::mutate(expr = shh_expr[sample_id], protein = "SHH"),
  if (!is.null(cdon_expr)) analysis1_df %>% dplyr::mutate(expr = cdon_expr[sample_id], protein = "CDON") else NULL,
  if (!is.null(gas1_expr)) analysis1_df %>% dplyr::mutate(expr = gas1_expr[sample_id], protein = "GAS1") else NULL
) %>% dplyr::filter(!is.na(response_binary), !is.na(expr))

p_box <- ggplot(resp_box_df, aes(x = response_binary, y = expr, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
  scale_fill_manual(values = c("SD/PD" = "#a6bddb", "PR/CR" = "#74c476"), guide = "none") +
  facet_wrap(~ protein, scales = "free_y") +
  labs(title = "SHH pathway members by treatment response", x = NULL, y = "log2 concentration") +
  akila_volcano_theme()

if (!is.null(p_scatter_cdon) && !is.null(p_scatter_gas1)) {
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  top_row <- patchwork::wrap_plots(p_scatter_cdon, p_scatter_gas1, ncol = 2)
  combined <- top_row / p_box + patchwork::plot_layout(heights = c(1, 1.1))
} else {
  combined <- p_box
}

stub <- file.path(OUT_FIG, "shh_pathway_coherence_check")
cairo_pdf(paste0(stub, ".pdf"), width = 11, height = 9); print(combined); dev.off()
png(paste0(stub, ".png"), width = 11, height = 9, units = "in", res = 300, bg = "white"); print(combined); dev.off()

message("Coherence check figure saved: ", stub, ".pdf / .png")
message("\n08p_shh_pathway_coherence_check.R complete")
