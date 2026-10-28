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

# Runs the Wilcoxon test + returns a tidy one-row summary, for a given
# subset of analysis1_df (all patients, or TKI-excluded), on the
# recurrence_binary_sensitivity outcome.
run_ceacam5_recurrence_wilcox <- function(sample_df, expr_vec, subset_label) {
  df <- sample_df %>%
    dplyr::mutate(expr = expr_vec[sample_id]) %>%
    dplyr::filter(!is.na(recurrence_binary_sensitivity), !is.na(expr))

  n_norec <- sum(df$recurrence_binary_sensitivity == "No recurrence")
  n_rec   <- sum(df$recurrence_binary_sensitivity == "Recurrence")

  if (n_norec < 3 || n_rec < 3) {
    return(tibble::tibble(
      subset = subset_label, n_no_recurrence = n_norec, n_recurrence = n_rec,
      median_no_recurrence = NA_real_, median_recurrence = NA_real_, wilcox_p = NA_real_
    ))
  }

  wt <- wilcox.test(expr ~ recurrence_binary_sensitivity, data = df)

  tibble::tibble(
    subset                = subset_label,
    n_no_recurrence        = n_norec,
    n_recurrence           = n_rec,
    median_no_recurrence   = median(df$expr[df$recurrence_binary_sensitivity == "No recurrence"]),
    median_recurrence      = median(df$expr[df$recurrence_binary_sensitivity == "Recurrence"]),
    wilcox_p               = wt$p.value
  )
}

# ── Locate CEACAM5 ────────────────────────────────────────────
message("\n=== 08zc | Locating CEACAM5 ===")

ceacam5_id <- feature_meta %>%
  dplyr::filter(Gene_Symbol == "CEACAM5") %>%
  dplyr::pull(nELISA_ID) %>%
  unique()

if (length(ceacam5_id) != 1) {
  stop("Expected exactly one CEACAM5 nELISA_ID in feature_meta, found: ",
       length(ceacam5_id))
}
message("CEACAM5 nELISA_ID: ", ceacam5_id)

ceacam5_expr <- expr_log2_rep[, ceacam5_id]
names(ceacam5_expr) <- rownames(expr_log2_rep)

# ── Confirm TKI subgroup size in the recurrence-evaluable set ─
message("\n=== 08zc | TKI subgroup check ===")

tki_n_total     <- sum(analysis1_df$new_group == "TKI")
tki_n_evaluable <- analysis1_df %>%
  dplyr::filter(new_group == "TKI", !is.na(recurrence_binary_sensitivity)) %>%
  nrow()

message("TKI patients (total, baseline PRE): ", tki_n_total)
message("TKI patients evaluable for recurrence_binary_sensitivity: ", tki_n_evaluable)

# ── Run Wilcoxon test: all patients vs. TKI-excluded ──────────
message("\n=== 08zc | Running sensitivity comparison ===")

result_all    <- run_ceacam5_recurrence_wilcox(analysis1_df, ceacam5_expr, "All patients")
result_no_tki <- run_ceacam5_recurrence_wilcox(
  analysis1_df %>% dplyr::filter(new_group != "TKI"),
  ceacam5_expr, "TKI excluded"
)

sensitivity_summary <- dplyr::bind_rows(result_all, result_no_tki)
safe_write_csv(sensitivity_summary, file.path(OUT_TBL, "ceacam5_recurrence_tki_sensitivity_check.csv"))

message("CEACAM5 recurrence association — all patients vs. TKI-excluded:")
print(sensitivity_summary)

p_all    <- result_all$wilcox_p
p_no_tki <- result_no_tki$wilcox_p
if (!is.na(p_all) && !is.na(p_no_tki)) {
  if (p_no_tki <= 0.10) {
    message("\n>> CEACAM5 recurrence association HOLDS with TKI patients excluded ",
            "(p = ", signif(p_no_tki, 3), " vs. p = ", signif(p_all, 3), " with all patients).")
  } else {
    message("\n>> CEACAM5 recurrence association WEAKENS/DISAPPEARS with TKI patients excluded ",
            "(p = ", signif(p_no_tki, 3), " vs. p = ", signif(p_all, 3), " with all patients) — ",
            "flag as a TKI-subgroup-driven confound in the manuscript.")
  }
} else {
  message("\n>> Wilcoxon test could not be computed for one or both subsets (insufficient n per arm). ",
          "Check sensitivity_summary above before using this result.")
}

# ── Boxplot: side-by-side, all patients vs. TKI-excluded ──────
message("\n=== 08zc | Building sensitivity boxplot ===")

plot_df <- dplyr::bind_rows(
  analysis1_df %>%
    dplyr::mutate(expr = ceacam5_expr[sample_id], subset = "All patients"),
  analysis1_df %>%
    dplyr::filter(new_group != "TKI") %>%
    dplyr::mutate(expr = ceacam5_expr[sample_id], subset = "TKI excluded")
) %>%
  dplyr::filter(!is.na(recurrence_binary_sensitivity), !is.na(expr)) %>%
  dplyr::mutate(
    subset = factor(subset, levels = c("All patients", "TKI excluded")),
    is_tki = ifelse(new_group == "TKI", "TKI patient", "Other treatment")
  )

p_sens <- ggplot(plot_df, aes(x = recurrence_binary_sensitivity, y = expr, fill = recurrence_binary_sensitivity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(shape = is_tki, size = is_tki), width = 0.15, alpha = 0.7) +
  scale_shape_manual(values = c("TKI patient" = 17, "Other treatment" = 16), name = NULL) +
  scale_size_manual(values = c("TKI patient" = 2.6, "Other treatment" = 1.6), name = NULL) +
  scale_fill_manual(values = c("No recurrence" = "#a6bddb", "Recurrence" = "#fb8072"), guide = "none") +
  facet_wrap(~ subset) +
  labs(
    title    = "CEACAM5 \u2014 Recurrence/Metastasis Association, TKI Sensitivity Check",
    subtitle = paste0("All patients: Wilcoxon p=", signif(p_all, 3),
                       "  |  TKI excluded: Wilcoxon p=", signif(p_no_tki, 3),
                       "  (TKI patients highlighted as triangles)"),
    x = NULL, y = "log2 concentration"
  ) +
  akila_volcano_theme()

stub_sens <- file.path(OUT_FIG, "ceacam5_recurrence_tki_sensitivity_boxplot")

cairo_pdf(paste0(stub_sens, ".pdf"), width = 9, height = 5.5)
print(p_sens)
dev.off()

png(paste0(stub_sens, ".png"), width = 9, height = 5.5, units = "in", res = 300, bg = "white")
print(p_sens)
dev.off()

message("\n=== 08zc complete ===")
message("NOTE: if 'No recurrence'/'Recurrence' arms are too small for a reliable Wilcoxon ",
        "(check n_no_recurrence/n_recurrence in the summary), treat this as descriptive only.")
