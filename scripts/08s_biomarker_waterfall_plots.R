if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoints 01 and 03 loaded")

while (dev.cur() > 1) dev.off()

OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

# ── Local constants: palettes verified against actual project scripts ──
# GROUP_COLOURS matches 02_qc_plots.R exactly (not guessed).
GROUP_COLOURS <- c(
  "Sx"        = "#E69F00",
  "CTx"       = "#D55E00",
  "CTx + IO"  = "#009E73",
  "CTx + RTx" = "#0072B2",
  "IO"        = "#56B4E9",
  "TKI"       = "#CC79A7",
  "Other"     = "#999999",
  "Unmatched" = "#DDDDDD"
)
# pal_response matches 01_clinical_proteomics_landscape.R.
PAL_RESPONSE <- c("CR" = "#7cb342", "PR" = "#2a9d8f", "SD" = "#f4a261", "PD" = "#d7301f", "NE" = "#d0d0d0")

TARGET_PROTEINS <- c("SHH", "CEACAM5", "GDF15")

# ── Local function: build one protein's waterfall + response strip ──
make_waterfall <- function(gene_symbol) {
  target_id <- feature_meta %>% dplyr::filter(Gene_Symbol == gene_symbol) %>% dplyr::pull(nELISA_ID) %>% unique()
  if (length(target_id) != 1) { warning("Skipping ", gene_symbol, " — not found uniquely on panel."); return(NULL) }

  expr_vec <- expr_log2_rep[, target_id]; names(expr_vec) <- rownames(expr_log2_rep)

  df <- analysis1_df %>%
    dplyr::mutate(expr = expr_vec[sample_id]) %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::mutate(
      response_display = dplyr::case_when(
        response_4cat %in% c("CR", "PR", "SD", "PD") ~ response_4cat,
        TRUE ~ "NE"
      ),
      new_group_display = factor(new_group, levels = names(GROUP_COLOURS)),
      response_display  = factor(response_display, levels = names(PAL_RESPONSE))
    ) %>%
    dplyr::arrange(expr) %>%
    dplyr::mutate(rank = dplyr::row_number())

  med_expr <- median(df$expr, na.rm = TRUE)

  bar_plot <- ggplot(df, aes(x = rank, y = expr, fill = new_group_display)) +
    geom_col(width = 1) +
    geom_hline(yintercept = med_expr, linetype = "dashed", colour = "black", linewidth = 0.4) +
    scale_fill_manual(values = GROUP_COLOURS, name = "Treatment", drop = FALSE, na.value = "#DDDDDD") +
    labs(
      title = paste0(gene_symbol, " \u2014 Waterfall Plot Across Cohort (n=", nrow(df), ")"),
      subtitle = paste0("Sorted ascending by log2 concentration; dashed line = cohort median (",
                         round(med_expr, 2), ")"),
      x = NULL, y = paste0(gene_symbol, " (log2)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30"),
      legend.position = "right"
    )

  strip_plot <- ggplot(df, aes(x = rank, y = "Response", fill = response_display)) +
    geom_tile(height = 1) +
    scale_fill_manual(values = PAL_RESPONSE, name = "Response", drop = FALSE) +
    labs(x = "Patients (sorted by expression, ascending)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  combined <- bar_plot / strip_plot + patchwork::plot_layout(heights = c(5, 1))

  n_by_group <- df %>% dplyr::count(new_group_display, name = "n") %>% tidyr::drop_na()
  message("  ", gene_symbol, ": n=", nrow(df), " | groups: ",
          paste(paste0(n_by_group$new_group_display, "=", n_by_group$n), collapse = ", "))

  combined
}

# ── Build and save one waterfall plot per protein ───────────────
message("\n=== 08s | Building waterfall plots ===")

for (gene in TARGET_PROTEINS) {
  p <- make_waterfall(gene)
  if (is.null(p)) next
  stub <- file.path(OUT_FIG, paste0("waterfall_", gene))
  cairo_pdf(paste0(stub, ".pdf"), width = 11, height = 6.5); print(p); dev.off()
  png(paste0(stub, ".png"), width = 11, height = 6.5, units = "in", res = 300, bg = "white"); print(p); dev.off()
  message("  Saved: ", stub, ".pdf / .png")
}

message("\n08s_biomarker_waterfall_plots.R complete")
