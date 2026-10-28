if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 03, 07 loaded")

while (dev.cur() > 1) dev.off()

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")
akila_theme <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13.5),
          plot.subtitle = element_text(size = 10.5, colour = "grey30"))
}
coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x)); op <- options(warn = -1)
  out <- as.numeric(raw_chr); options(op); out
}
get_expr <- function(id) { v <- expr_log2_rep[, id]; names(v) <- rownames(expr_log2_rep); v }

# ── CV percentile lookup for all three headline proteins ───────
message("\n=== 08r | Technical CV percentile check (GDF15, CEACAM5, SHH) ===")

cv_tbl <- readr::read_csv(file.path(paths$tables, "01_qc", "protein_cv_distribution.csv"), show_col_types = FALSE)
cv_percentile <- function(gene) {
  target_cv <- cv_tbl$cv[cv_tbl$Gene_Symbol == gene]
  if (length(target_cv) != 1 || is.na(target_cv)) return(NA_real_)
  mean(cv_tbl$cv < target_cv, na.rm = TRUE)
}
cv_summary <- tibble::tibble(
  Gene_Symbol = c("GDF15", "CEACAM5", "SHH"),
  cv          = sapply(c("GDF15", "CEACAM5", "SHH"), function(g) cv_tbl$cv[cv_tbl$Gene_Symbol == g][1]),
  cv_percentile = sapply(c("GDF15", "CEACAM5", "SHH"), cv_percentile)
)
message("Note: interpretation of this CV depends on whether it reflects technical (assay) noise ",
        "or cross-sample biological variability — confirm the exact formula in 02_qc_plots.R before ",
        "drawing conclusions from it.")
print(cv_summary)

# ── GDF15 coherence: IL-6 axis + Growth Hormone ─────────────────
message("\n=== 08r | GDF15 coherence check (IL-6 axis, Growth Hormone) ===")

gdf15_id <- feature_meta %>% dplyr::filter(Gene_Symbol == "GDF15") %>% dplyr::pull(nELISA_ID) %>% unique()
gdf15_secondary <- feature_meta %>%
  dplyr::filter(Gene_Symbol %in% c("IL6", "IL6R", "IL6ST")) %>%
  dplyr::distinct(Gene_Symbol, nELISA_ID, Protein_Label)
gh_id <- feature_meta %>% dplyr::filter(Protein_Label == "Growth Hormone (Somatotropin)") %>%
  dplyr::pull(nELISA_ID) %>% unique()

message("GDF15 secondary proteins found: ", paste(gdf15_secondary$Gene_Symbol, collapse = ", "),
        if (length(gh_id) == 1) ", Growth Hormone (Somatotropin)" else "")

gdf15_expr <- get_expr(gdf15_id)

run_corr <- function(x, y, label) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 10) return(tibble::tibble(pair = label, n = sum(keep), rho = NA_real_, p_value = NA_real_))
  ct <- cor.test(x[keep], y[keep], method = "spearman")
  tibble::tibble(pair = label, n = sum(keep), rho = unname(ct$estimate), p_value = ct$p.value)
}
run_os_km <- function(expr_vec, label) {
  d <- survival_df %>%
    dplyr::mutate(expr = expr_vec[sample_id], time = as.numeric(os_months), event = coerce_event_indicator(os_event)) %>%
    dplyr::filter(!is.na(expr), !is.na(time), !is.na(event))
  if (nrow(d) < 20 || sum(d$event) < 8) return(tibble::tibble(protein = label, n = nrow(d), logrank_p = NA_real_))
  med <- median(d$expr, na.rm = TRUE); d$grp <- ifelse(d$expr >= med, "High", "Low")
  sd_ <- survival::survdiff(survival::Surv(time, event) ~ grp, data = d)
  tibble::tibble(protein = label, n = nrow(d), logrank_p = pchisq(sd_$chisq, df = 1, lower.tail = FALSE))
}

gdf15_corrs <- purrr::pmap_dfr(gdf15_secondary, function(Gene_Symbol, nELISA_ID, Protein_Label) {
  run_corr(gdf15_expr, get_expr(nELISA_ID), paste0("GDF15 vs ", Gene_Symbol))
})
if (length(gh_id) == 1) {
  gdf15_corrs <- dplyr::bind_rows(gdf15_corrs, run_corr(gdf15_expr, get_expr(gh_id), "GDF15 vs Growth Hormone"))
}
gdf15_os_checks <- dplyr::bind_rows(
  run_os_km(gdf15_expr, "GDF15"),
  purrr::pmap_dfr(gdf15_secondary, function(Gene_Symbol, nELISA_ID, Protein_Label) run_os_km(get_expr(nELISA_ID), Gene_Symbol)),
  if (length(gh_id) == 1) run_os_km(get_expr(gh_id), "Growth Hormone") else NULL
)
message("GDF15 correlations:"); print(gdf15_corrs)
message("GDF15 axis OS associations (median-split log-rank):"); print(gdf15_os_checks)

# ── CEACAM5 coherence: CEACAM1, CEACAM6 ─────────────────────────
message("\n=== 08r | CEACAM5 coherence check (CEACAM1, CEACAM6) ===")

ceacam5_id <- feature_meta %>% dplyr::filter(Gene_Symbol == "CEACAM5") %>% dplyr::pull(nELISA_ID) %>% unique()
ceacam_secondary <- feature_meta %>%
  dplyr::filter(Gene_Symbol %in% c("CEACAM1", "CEACAM6")) %>%
  dplyr::distinct(Gene_Symbol, nELISA_ID, Protein_Label)

message("CEACAM5 paralogs found: ", paste(ceacam_secondary$Gene_Symbol, collapse = ", "))

ceacam5_expr <- get_expr(ceacam5_id)
ceacam_corrs <- purrr::pmap_dfr(ceacam_secondary, function(Gene_Symbol, nELISA_ID, Protein_Label) {
  run_corr(ceacam5_expr, get_expr(nELISA_ID), paste0("CEACAM5 vs ", Gene_Symbol))
})

run_response_wilcox <- function(expr_vec, label) {
  d <- analysis1_df %>% dplyr::mutate(expr = expr_vec[sample_id]) %>%
    dplyr::filter(!is.na(response_binary), !is.na(expr))
  n_sdpd <- sum(d$response_binary == "SD/PD"); n_prcr <- sum(d$response_binary == "PR/CR")
  if (n_sdpd < 8 || n_prcr < 8) return(tibble::tibble(protein = label, wilcox_p = NA_real_))
  wt <- wilcox.test(expr ~ response_binary, data = d)
  tibble::tibble(protein = label, n_sdpd = n_sdpd, n_prcr = n_prcr, wilcox_p = wt$p.value)
}
ceacam_response <- dplyr::bind_rows(
  run_response_wilcox(ceacam5_expr, "CEACAM5"),
  purrr::pmap_dfr(ceacam_secondary, function(Gene_Symbol, nELISA_ID, Protein_Label) run_response_wilcox(get_expr(nELISA_ID), Gene_Symbol))
)
message("CEACAM5 correlations:"); print(ceacam_corrs)
message("CEACAM family response associations:"); print(ceacam_response)

# ── Save combined summary ───────────────────────────────────────
safe_write_csv(
  dplyr::bind_rows(
    cv_summary %>% dplyr::mutate(table = "cv_percentile"),
    gdf15_corrs %>% dplyr::mutate(table = "gdf15_correlation"),
    gdf15_os_checks %>% dplyr::mutate(table = "gdf15_os_association"),
    ceacam_corrs %>% dplyr::mutate(table = "ceacam5_correlation"),
    ceacam_response %>% dplyr::mutate(table = "ceacam5_response_association")
  ),
  file.path(OUT_TBL, "gdf15_ceacam5_coherence_check.csv")
)

# ── Figures ──────────────────────────────────────────────────────
message("\n=== 08r | Building coherence check figures ===")

make_scatter <- function(x, y, xlab, ylab, corr_row) {
  d <- tibble::tibble(x = x, y = y) %>% dplyr::filter(!is.na(x), !is.na(y))
  ggplot(d, aes(x, y)) +
    geom_point(alpha = 0.55, size = 1.5, colour = "#2f6fb0") +
    geom_smooth(method = "lm", se = TRUE, colour = "#c9772b", linewidth = 0.6) +
    labs(title = paste0(xlab, " vs ", ylab, " (\u03c1=", signif(corr_row$rho, 3), ", p=", signif(corr_row$p_value, 3), ")"),
         x = paste0(xlab, " (log2)"), y = paste0(ylab, " (log2)")) +
    akila_theme()
}

gdf15_plots <- purrr::pmap(gdf15_secondary, function(Gene_Symbol, nELISA_ID, Protein_Label) {
  row <- gdf15_corrs[gdf15_corrs$pair == paste0("GDF15 vs ", Gene_Symbol), ]
  make_scatter(gdf15_expr, get_expr(nELISA_ID), "GDF15", Gene_Symbol, row)
})
if (length(gh_id) == 1) {
  row <- gdf15_corrs[gdf15_corrs$pair == "GDF15 vs Growth Hormone", ]
  gdf15_plots <- c(gdf15_plots, list(make_scatter(gdf15_expr, get_expr(gh_id), "GDF15", "Growth Hormone", row)))
}
gdf15_combined <- patchwork::wrap_plots(gdf15_plots, ncol = 2) +
  patchwork::plot_annotation(title = "GDF15 coherence check: systemic inflammation / cachexia axis")

stub_g <- file.path(OUT_FIG, "gdf15_coherence_check")
cairo_pdf(paste0(stub_g, ".pdf"), width = 10, height = 8); print(gdf15_combined); dev.off()
png(paste0(stub_g, ".png"), width = 10, height = 8, units = "in", res = 300, bg = "white"); print(gdf15_combined); dev.off()

ceacam_plots <- purrr::pmap(ceacam_secondary, function(Gene_Symbol, nELISA_ID, Protein_Label) {
  row <- ceacam_corrs[ceacam_corrs$pair == paste0("CEACAM5 vs ", Gene_Symbol), ]
  make_scatter(ceacam5_expr, get_expr(nELISA_ID), "CEACAM5", Gene_Symbol, row)
})
ceacam_combined <- patchwork::wrap_plots(ceacam_plots, ncol = 2) +
  patchwork::plot_annotation(title = "CEACAM5 coherence check: paralogous family members")

stub_c <- file.path(OUT_FIG, "ceacam5_coherence_check")
cairo_pdf(paste0(stub_c, ".pdf"), width = 10, height = 5); print(ceacam_combined); dev.off()
png(paste0(stub_c, ".png"), width = 10, height = 5, units = "in", res = 300, bg = "white"); print(ceacam_combined); dev.off()

message("Figures saved: ", stub_g, ".pdf/.png and ", stub_c, ".pdf/.png")
message("\n08r_gdf15_ceacam5_coherence_check.R complete")
