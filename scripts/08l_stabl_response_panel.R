sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
STABL_TBL_DIR <- file.path(paths$tables,  "07_analysis6_machine_learning", "stabl")
STABL_FIG_DIR <- file.path(paths$figures, "07_analysis6_machine_learning", "stabl")

TOP_N_PER_PANEL <- 25   # proteins shown per mini-panel (matches the individual plots)

# Response comparisons in run order, with a clean display label for each
RESPONSE_COMPARISONS <- tibble::tribble(
  ~file_prefix,               ~display_label,
  "stabl_01_CR_vs_PD",        "CR vs PD",
  "stabl_02_CRPR_vs_PD",      "CR/PR vs PD",
  "stabl_03_CRPR_vs_SDPD",    "CR/PR vs SD/PD",
  "stabl_04_CRPRSD_vs_PD",    "CR/PR/SD vs PD",
  "stabl_05_PR_vs_PD",        "PR vs PD"
)

# ── Local function: one mini selection-frequency panel ─────────
make_stabl_mini_panel <- function(file_prefix, display_label, top_n = TOP_N_PER_PANEL) {
  freq_path <- file.path(STABL_TBL_DIR, paste0(file_prefix, "_stabl_selection_freq.csv"))
  if (!file.exists(freq_path)) {
    warning("Missing selection_freq file for ", file_prefix, " — skipping panel.")
    return(NULL)
  }

  df <- readr::read_csv(freq_path, show_col_types = FALSE)
  tau_star <- unique(df$tau_star)[1]

  plot_df <- df %>%
    dplyr::arrange(dplyr::desc(real_freq)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(Gene_Symbol = factor(Gene_Symbol, levels = rev(Gene_Symbol)))

  n_selected <- sum(df$stabl_selected, na.rm = TRUE)

  ggplot(plot_df, aes(x = Gene_Symbol, y = real_freq, fill = stabl_selected)) +
    geom_col(width = 0.75) +
    geom_hline(yintercept = tau_star, linetype = "dashed", colour = "#d62728", linewidth = 0.4) +
    scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#a6bddb"), guide = "none") +
    coord_flip() +
    labs(
      title = paste0(display_label, "  (\u03c4*=", round(tau_star, 2), ", n selected=", n_selected, ")"),
      x = NULL, y = "Selection frequency"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title        = element_text(face = "bold", size = 9.5),
      axis.text.y       = element_text(size = 6.3),
      axis.title.x      = element_text(size = 8)
    )
}

# ── Build all 5 panels ──────────────────────────────────────────
message("\n=== 08l | Building STABL response selection-frequency panels ===")

panels <- purrr::map2(
  RESPONSE_COMPARISONS$file_prefix, RESPONSE_COMPARISONS$display_label,
  make_stabl_mini_panel
)
panels <- Filter(Negate(is.null), panels)

message("Panels built: ", length(panels), " of ", nrow(RESPONSE_COMPARISONS), " response comparisons")

combined <- patchwork::wrap_plots(panels, ncol = 3) +
  patchwork::plot_annotation(
    title    = "Figure 5(b): STABL Selection Frequency Across Response Contrasts",
    subtitle = "Top 25 proteins per contrast; dashed line = data-driven stability threshold (\u03c4*); dark bars = STABL-selected",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, colour = "grey30"))
  )

stub_panel <- file.path(STABL_FIG_DIR, "stabl_response_panel")

cairo_pdf(paste0(stub_panel, ".pdf"), width = 13, height = 8)
print(combined)
dev.off()

png(paste0(stub_panel, ".png"), width = 13, height = 8, units = "in", res = 300, bg = "white")
print(combined)
dev.off()

message("STABL response panel saved: ", stub_panel, ".pdf / .png")
message("\n08l_stabl_response_panel.R complete")
