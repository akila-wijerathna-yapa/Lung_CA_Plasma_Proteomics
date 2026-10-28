sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

ALTML_TBL_DIR <- file.path(paths$tables, "07_analysis6_machine_learning", "altml")
CV_TBL_DIR    <- file.path(paths$tables, "07_analysis6_machine_learning", "cv")
OUT_TBL       <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG       <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

ENDPOINT_FILES <- tibble::tribble(
  ~endpoint, ~altml_file,
  "OS",  "altml_14_cox_OS_all_methods.csv",
  "PFS", "altml_15_cox_PFS_all_methods.csv"
)

# ── Load in-sample Cox results ──────────────────────────────────
message("\n=== 08n | Loading in-sample Cox (altml) results ===")

insample_df <- purrr::pmap_dfr(ENDPOINT_FILES, function(endpoint, altml_file) {
  f <- file.path(ALTML_TBL_DIR, altml_file)
  if (!file.exists(f)) {
    warning("Missing altml Cox file: ", f)
    return(NULL)
  }
  row1 <- readr::read_csv(f, show_col_types = FALSE, n_max = 1)
  tibble::tibble(
    endpoint = endpoint,
    method   = c("Elastic Net (Cox)", "RF-survival"),
    c_index_insample = c(row1$en_c_index, row1$rf_surv_c_index)
  )
})

message("In-sample Cox rows loaded: ", nrow(insample_df))

# ── Load CV Cox results ─────────────────────────────────────────
message("\n=== 08n | Loading cross-validated Cox results ===")

cv_csv_path <- file.path(CV_TBL_DIR, "cv_cox_all_comparisons.csv")
if (!file.exists(cv_csv_path)) {
  stop("cv_cox_all_comparisons.csv not found — run 08e_cv_validation.R first ",
       "(with the reverse=TRUE and checkpoint_03 fixes).")
}

cv_cox <- readr::read_csv(cv_csv_path, show_col_types = FALSE) %>%
  dplyr::mutate(endpoint = dplyr::case_when(
    grepl("cox_OS",  comparison) ~ "OS",
    grepl("cox_PFS", comparison) ~ "PFS",
    TRUE ~ NA_character_
  ))

message("CV Cox rows loaded: ", nrow(cv_cox))

# ── Join and reshape ─────────────────────────────────────────────
message("\n=== 08n | Joining in-sample + CV ===")

joined <- insample_df %>%
  dplyr::inner_join(
    cv_cox %>% dplyr::select(endpoint, method, c_index_cv = cv_cindex_mean, cv_cindex_sd),
    by = c("endpoint", "method")
  ) %>%
  dplyr::mutate(
    row_label   = paste0(endpoint, " - ", method),
    overfit_gap = c_index_insample - c_index_cv
  )

missing_pairs <- nrow(ENDPOINT_FILES) * 2 - nrow(joined)
if (missing_pairs != 0) {
  warning(missing_pairs, " expected endpoint x method pairs did not join — ",
          "check comparison/method naming consistency.")
}

safe_write_csv(
  joined %>% dplyr::select(row_label, endpoint, method, c_index_insample, c_index_cv, cv_cindex_sd, overfit_gap),
  file.path(OUT_TBL, "insample_vs_cv_cox_summary.csv")
)

message("In-sample vs. CV Cox summary:")
print(joined %>% dplyr::select(row_label, c_index_insample, c_index_cv, overfit_gap))

# ── Dumbbell plot ────────────────────────────────────────────────
message("\n=== 08n | Building in-sample vs CV Cox dumbbell plot ===")

row_order <- c("PFS - RF-survival", "PFS - Elastic Net (Cox)",
               "OS - RF-survival",  "OS - Elastic Net (Cox)")

joined <- joined %>%
  dplyr::mutate(row_label = factor(row_label, levels = row_order))

long_df <- dplyr::bind_rows(
  joined %>% dplyr::transmute(row_label, phase = "In-sample", c_index = c_index_insample),
  joined %>% dplyr::transmute(row_label, phase = "Cross-validated", c_index = c_index_cv)
) %>%
  dplyr::mutate(phase = factor(phase, levels = c("In-sample", "Cross-validated")))

p_cox_dumbbell <- ggplot() +
  geom_segment(
    data = joined,
    aes(x = row_label, xend = row_label, y = c_index_insample, yend = c_index_cv),
    colour = "grey60", linewidth = 0.7
  ) +
  geom_point(
    data = long_df,
    aes(x = row_label, y = c_index, colour = phase), size = 3.4
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  coord_flip(ylim = c(0.4, 0.8)) +
  scale_colour_manual(values = c("In-sample" = "#c9772b", "Cross-validated" = "#2f6fb0"), name = NULL) +
  labs(
    title    = "Figure 6(b): In-Sample vs. Cross-Validated C-index, Cox Models",
    subtitle = "OS stable under CV (both methods); PFS collapses for Elastic Net specifically",
    x = NULL, y = "C-index"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(size = 10.5, colour = "grey30"),
    legend.position   = "bottom"
  )

stub_cox_dumbbell <- file.path(OUT_FIG, "insample_vs_cv_cox_dumbbell")

cairo_pdf(paste0(stub_cox_dumbbell, ".pdf"), width = 8, height = 5.5)
print(p_cox_dumbbell)
dev.off()

png(paste0(stub_cox_dumbbell, ".png"), width = 8, height = 5.5, units = "in", res = 300, bg = "white")
print(p_cox_dumbbell)
dev.off()

message("Cox dumbbell plot saved: ", stub_cox_dumbbell, ".pdf / .png")
message("\n08n_insample_vs_cv_cox_dumbbell.R complete")
