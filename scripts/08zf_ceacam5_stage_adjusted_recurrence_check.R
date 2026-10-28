if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")

while (dev.cur() > 1) dev.off()

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

akila_volcano_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 15),
      plot.subtitle     = element_text(size = 11, colour = "grey30"),
      axis.title        = element_text(face = "bold")
    )
}

get_expr <- function(id) { v <- expr_log2_rep[, id]; names(v) <- rownames(expr_log2_rep); v }

# ── Locate CEACAM5 ────────────────────────────────────────────
message("\n=== 08zf | Locating CEACAM5 ===")
ceacam5_id <- feature_meta %>%
  dplyr::filter(Gene_Symbol == "CEACAM5") %>%
  dplyr::pull(nELISA_ID) %>% unique()
if (length(ceacam5_id) != 1) stop("Expected exactly one CEACAM5 nELISA_ID, found: ", length(ceacam5_id))
message("CEACAM5 nELISA_ID: ", ceacam5_id)
ceacam5_expr <- get_expr(ceacam5_id)

# ── Build analysis frame: recurrence + stage_collapsed + CEACAM5 ──
message("\n=== 08zf | Building analysis frame ===")
df <- analysis1_df %>%
  dplyr::mutate(expr = ceacam5_expr[sample_id]) %>%
  dplyr::filter(!is.na(recurrence_binary_sensitivity), !is.na(stage_collapsed), !is.na(expr)) %>%
  dplyr::mutate(
    recurrence_binary_sensitivity = factor(recurrence_binary_sensitivity,
                                            levels = c("No recurrence", "Recurrence")),
    stage_collapsed = factor(stage_collapsed, levels = c("I_II", "III_IV"))
  )

message("n = ", nrow(df), " | Stage I-II: ", sum(df$stage_collapsed == "I_II"),
        " | Stage III-IV: ", sum(df$stage_collapsed == "III_IV"))
print(table(df$stage_collapsed, df$recurrence_binary_sensitivity))

# ── Check 1: Stage-stratified Wilcoxon sensitivity ─────────────
message("\n=== 08zf | Stage-stratified Wilcoxon sensitivity check ===")

run_stratum_wilcox <- function(sub_df, stratum_label) {
  n_norec <- sum(sub_df$recurrence_binary_sensitivity == "No recurrence")
  n_rec   <- sum(sub_df$recurrence_binary_sensitivity == "Recurrence")
  if (n_norec < 3 || n_rec < 3) {
    return(tibble::tibble(
      stratum = stratum_label, n_no_recurrence = n_norec, n_recurrence = n_rec,
      median_no_recurrence = NA_real_, median_recurrence = NA_real_, wilcox_p = NA_real_
    ))
  }
  wt <- wilcox.test(expr ~ recurrence_binary_sensitivity, data = sub_df)
  tibble::tibble(
    stratum              = stratum_label,
    n_no_recurrence       = n_norec,
    n_recurrence          = n_rec,
    median_no_recurrence  = median(sub_df$expr[sub_df$recurrence_binary_sensitivity == "No recurrence"]),
    median_recurrence     = median(sub_df$expr[sub_df$recurrence_binary_sensitivity == "Recurrence"]),
    wilcox_p              = wt$p.value
  )
}

strat_I_II   <- run_stratum_wilcox(df %>% dplyr::filter(stage_collapsed == "I_II"),   "Stage I-II")
strat_III_IV <- run_stratum_wilcox(df %>% dplyr::filter(stage_collapsed == "III_IV"), "Stage III-IV")
strat_all    <- run_stratum_wilcox(df, "All stages combined")

stratified_summary <- dplyr::bind_rows(strat_I_II, strat_III_IV, strat_all)
message("Stage-stratified CEACAM5 vs recurrence:")
print(stratified_summary)

# ── Check 2: Multivariable logistic model, CEACAM5 + stage_collapsed ──
message("\n=== 08zf | Multivariable logistic model: recurrence ~ CEACAM5 + stage ===")

df_std <- df %>% dplyr::mutate(expr_z = as.numeric(scale(expr)))

fit_adj <- glm(recurrence_binary_sensitivity ~ expr_z + stage_collapsed,
               data = df_std, family = binomial())
fit_unadj <- glm(recurrence_binary_sensitivity ~ expr_z,
                  data = df_std, family = binomial())

adj_summary <- broom::tidy(fit_adj, conf.int = TRUE, exponentiate = TRUE) %>%
  dplyr::filter(term == "expr_z") %>%
  dplyr::mutate(model = "CEACAM5 + stage (adjusted)")
unadj_summary <- broom::tidy(fit_unadj, conf.int = TRUE, exponentiate = TRUE) %>%
  dplyr::filter(term == "expr_z") %>%
  dplyr::mutate(model = "CEACAM5 only (unadjusted)")

logistic_summary <- dplyr::bind_rows(unadj_summary, adj_summary) %>%
  dplyr::select(model, term, OR = estimate, OR_low = conf.low, OR_high = conf.high,
                std.error, p.value)

message("Logistic model OR for CEACAM5 (per 1-SD increase in log2 expression):")
print(logistic_summary)

# ── Combine and save table ──────────────────────────────────────
safe_write_csv(stratified_summary, file.path(OUT_TBL, "ceacam5_stage_stratified_recurrence_check.csv"))
safe_write_csv(logistic_summary,   file.path(OUT_TBL, "ceacam5_stage_adjusted_recurrence_logistic.csv"))

# ── Interpretation message ──────────────────────────────────────
p_I_II   <- strat_I_II$wilcox_p
p_III_IV <- strat_III_IV$wilcox_p
p_adj    <- adj_summary$p.value

message("\n=== 08zf | Interpretation ===")
if (!is.na(p_I_II) && !is.na(p_III_IV) && p_I_II <= 0.10 && p_III_IV <= 0.10) {
  message(">> CEACAM5-recurrence association holds within BOTH stage strata ",
          "(Stage I-II p=", signif(p_I_II, 3), ", Stage III-IV p=", signif(p_III_IV, 3),
          ") — not purely a stage proxy.")
} else {
  message(">> CEACAM5-recurrence association does NOT hold consistently across both stage strata ",
          "(Stage I-II p=", signif(p_I_II, 3), ", Stage III-IV p=", signif(p_III_IV, 3),
          ") — interpret with caution, may be partly stage-driven.")
}
if (!is.na(p_adj)) {
  if (p_adj <= 0.10) {
    message(">> Stage-adjusted multivariable model: CEACAM5 remains associated with recurrence ",
            "after adjusting for stage (adjusted p=", signif(p_adj, 3), ").")
  } else {
    message(">> Stage-adjusted multivariable model: CEACAM5's association with recurrence ",
            "attenuates after adjusting for stage (adjusted p=", signif(p_adj, 3), ") — ",
            "flag as a stage-confounded finding, do not report as stage-independent without this caveat.")
  }
}

# ── Figure: CEACAM5 by recurrence, faceted by stage ─────────────
message("\n=== 08zf | Building stage-faceted boxplot ===")

p_stage <- ggplot(df, aes(x = recurrence_binary_sensitivity, y = expr, fill = recurrence_binary_sensitivity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
  scale_fill_manual(values = c("No recurrence" = "#a6bddb", "Recurrence" = "#fb8072"), guide = "none") +
  facet_wrap(~ stage_collapsed, labeller = as_labeller(c("I_II" = "Stage I-II", "III_IV" = "Stage III-IV"))) +
  labs(
    title    = "CEACAM5 \u2014 Recurrence Association, Stage-Stratified Sensitivity Check",
    subtitle = paste0("Stage I-II: Wilcoxon p=", signif(p_I_II, 3),
                       "  |  Stage III-IV: Wilcoxon p=", signif(p_III_IV, 3),
                       "  |  Stage-adjusted multivariable p=", signif(p_adj, 3)),
    x = NULL, y = "log2 concentration"
  ) +
  akila_volcano_theme()

stub_stage <- file.path(OUT_FIG, "ceacam5_recurrence_by_stage_boxplot")
cairo_pdf(paste0(stub_stage, ".pdf"), width = 9, height = 5.5); print(p_stage); dev.off()
png(paste0(stub_stage, ".png"), width = 9, height = 5.5, units = "in", res = 300, bg = "white"); print(p_stage); dev.off()

