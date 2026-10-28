if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoint 03 loaded")

while (dev.cur() > 1) dev.off()

OUT_TBL <- file.path(paths$tables,  "07_analysis6_machine_learning")
OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")
norm_id <- function(x) {
  out <- as.character(x); out <- trimws(out); out <- gsub("\\.0$", "", out)
  out <- gsub("[^A-Za-z0-9]+", "_", out); out <- gsub("^_+|_+$", "", out); toupper(out)
}
simplify_histology <- function(x) {
  x_low <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    is.na(x) ~ "Unknown",
    stringr::str_detect(x_low, "squamous") & !stringr::str_detect(x_low, "adeno") ~ "Squamous Cell Carcinoma",
    stringr::str_detect(x_low, "small cell|neuroendocrine|nec\\b") ~ "Small Cell / Neuroendocrine",
    stringr::str_detect(x_low, "sarcomatoid|spindle|adenosquamous") ~ "Sarcomatoid / Mixed / Other",
    stringr::str_detect(x_low, "adenocarcinoma") ~ "Adenocarcinoma",
    TRUE ~ "Sarcomatoid / Mixed / Other"
  )
}

GROUP_COLOURS <- c("Sx" = "#E69F00", "CTx" = "#D55E00", "CTx + IO" = "#009E73", "CTx + RTx" = "#0072B2",
                   "IO" = "#56B4E9", "TKI" = "#CC79A7", "Other" = "#999999")
PAL_RESPONSE  <- c("CR" = "#7cb342", "PR" = "#2a9d8f", "SD" = "#f4a261", "PD" = "#d7301f")
PAL_STAGE     <- c("I" = "#fbb4ae", "II" = "#fdcdac", "III" = "#b3cde3", "IV" = "#80cdc1")
PAL_HIST      <- c("Adenocarcinoma" = "#7cb342", "Squamous Cell Carcinoma" = "#2a9d8f",
                    "Small Cell / Neuroendocrine" = "#e76f51", "Sarcomatoid / Mixed / Other" = "#8d99ae")

# ── Defensive read of clinical workbook (same checks as 08t/08u) ──
message("\n=== 08v | Reading raw clinical workbook ===")

if (!file.exists(paths$clinical)) stop("Clinical file not found at: ", paths$clinical)
if (file.size(paths$clinical) < 5000) {
  stop("Clinical file at ", paths$clinical, " is suspiciously small (",
       file.size(paths$clinical), " bytes) \u2014 likely a OneDrive sync placeholder. ",
       "Right-click the file in File Explorer > 'Always keep on this device', wait for sync, retry.")
}
sheet_names <- tryCatch(readxl::excel_sheets(paths$clinical), error = function(e) {
  stop("Could not open ", paths$clinical, " (", conditionMessage(e), ") \u2014 check OneDrive sync/file integrity.")
})
target_sheet <- sheet_names[grepl("clinical|clean", sheet_names, ignore.case = TRUE)][1]
if (is.na(target_sheet)) target_sheet <- sheet_names[1]

clin_raw <- readxl::read_excel(paths$clinical, sheet = target_sheet) %>%
  dplyr::transmute(
    patient_id           = norm_id(Patient_ID),
    stage_group_raw      = as.character(Stage_Group),
    histological_subtype = as.character(Histological_Subtype),
    pdl1_tps_pct         = suppressWarnings(as.numeric(PDL1_TPS_pct))
  ) %>%
  dplyr::mutate(
    stage_group_display = ifelse(stage_group_raw %in% c("I","II","III","IV"), stage_group_raw, NA_character_),
    histology_display    = simplify_histology(histological_subtype)
  )

# ── Join + filter ────────────────────────────────────────────────
df <- analysis1_df %>%
  dplyr::distinct(patient_id, sample_id, new_group, response_binary, response_4cat) %>%
  dplyr::inner_join(clin_raw, by = "patient_id") %>%
  dplyr::filter(!is.na(pdl1_tps_pct)) %>%
  dplyr::mutate(response_display = dplyr::case_when(response_4cat %in% names(PAL_RESPONSE) ~ response_4cat, TRUE ~ NA_character_))

message("n with non-missing PDL1_TPS_pct: ", nrow(df))

# ── Local function: one grouped-boxplot panel + test ────────────
make_panel <- function(group_col, palette, title, test_type = c("kruskal", "wilcox")) {
  test_type <- match.arg(test_type)
  d <- df %>% dplyr::filter(!is.na(.data[[group_col]]), .data[[group_col]] %in% names(palette)) %>%
    dplyr::mutate(grp = factor(.data[[group_col]], levels = names(palette)))
  if (nrow(d) < 6 || dplyr::n_distinct(d$grp) < 2) return(list(plot = NULL, stat = NA_real_, n = nrow(d)))

  if (test_type == "kruskal" ) {
    kt <- kruskal.test(pdl1_tps_pct ~ grp, data = d)
    p_val <- kt$p.value; test_label <- "Kruskal-Wallis"
  } else {
    wt <- wilcox.test(pdl1_tps_pct ~ grp, data = d)
    p_val <- wt$p.value; test_label <- "Wilcoxon"
  }

  n_by_grp <- d %>% dplyr::count(grp)
  plt <- ggplot(d, aes(x = grp, y = pdl1_tps_pct, fill = grp)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 1.4, alpha = 0.5) +
    scale_fill_manual(values = palette, guide = "none") +
    labs(title = paste0(title, " (", test_label, " p=", signif(p_val, 3), ")"),
         x = NULL, y = "PD-L1 TPS (%)") +
    theme_bw(base_size = 11.5) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11.5),
          axis.text.x = element_text(angle = 20, hjust = 1))

  list(plot = plt, stat = p_val, n = nrow(d), test = test_label)
}

message("\n=== 08v | Building subgroup panels ===")

p_treat <- make_panel("new_group", GROUP_COLOURS, "PD-L1 by Treatment", "kruskal")
p_resp  <- make_panel("response_display", PAL_RESPONSE, "PD-L1 by Response (whole cohort)", "kruskal")
p_stage <- make_panel("stage_group_display", PAL_STAGE, "PD-L1 by Stage", "kruskal")
p_hist  <- make_panel("histology_display", PAL_HIST, "PD-L1 by Histology", "wilcox")

results_summary <- tibble::tibble(
  variable = c("Treatment", "Response (whole cohort)", "Stage", "Histology"),
  n     = c(p_treat$n, p_resp$n, p_stage$n, p_hist$n),
  test  = c(p_treat$test, p_resp$test, p_stage$test, p_hist$test),
  p_value = c(p_treat$stat, p_resp$stat, p_stage$stat, p_hist$stat)
)
safe_write_csv(results_summary, file.path(OUT_TBL, "pdl1_by_subgroup_tests.csv"))
message("Subgroup test results:"); print(results_summary)
message("NOTE: see 08u_pdl1_io_response_check.R for the mechanistically correct ",
        "IO-restricted response test \u2014 report that result, not the whole-cohort response panel here.")

combined <- (p_treat$plot | p_resp$plot) / (p_stage$plot | p_hist$plot) +
  patchwork::plot_annotation(
    title = "PD-L1 TPS% by Clinical Subgroup",
    subtitle = paste0("n=", nrow(df), " with non-missing PD-L1; boxplot + jitter, one panel per variable"),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10.5, colour = "grey30"))
  )

stub <- file.path(OUT_FIG, "pdl1_by_subgroup_boxplots")
cairo_pdf(paste0(stub, ".pdf"), width = 11, height = 9); print(combined); dev.off()
png(paste0(stub, ".png"), width = 11, height = 9, units = "in", res = 300, bg = "white"); print(combined); dev.off()

message("Saved: ", stub, ".pdf / .png")
message("\n08v_pdl1_by_subgroup_boxplots.R complete")
