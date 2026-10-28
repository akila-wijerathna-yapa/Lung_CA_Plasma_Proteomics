if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# ── Load checkpoint ────────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_03.RData"))
message("Checkpoint 03 loaded")

while (dev.cur() > 1) dev.off()

OUT_FIG <- file.path(paths$figures, "07_analysis6_machine_learning")

# ── Local constants: palettes verified against actual project scripts ──
GROUP_COLOURS <- c(
  "Sx"        = "#E69F00", "CTx"       = "#D55E00", "CTx + IO"  = "#009E73",
  "CTx + RTx" = "#0072B2", "IO"        = "#56B4E9", "TKI"       = "#CC79A7",
  "Other"     = "#999999", "Unmatched" = "#DDDDDD"
)
PAL_RESPONSE <- c("CR" = "#7cb342", "PR" = "#2a9d8f", "SD" = "#f4a261", "PD" = "#d7301f", "NE" = "#d0d0d0")
PAL_STAGE <- c("I" = "#fbb4ae", "II" = "#fdcdac", "III" = "#b3cde3", "IV" = "#80cdc1", "Unknown" = "#d9d9d9")
PAL_HIST  <- c("Adenocarcinoma" = "#7cb342", "Squamous Cell Carcinoma" = "#2a9d8f",
               "Small Cell / Neuroendocrine" = "#e76f51", "Sarcomatoid / Mixed / Other" = "#8d99ae",
               "Unknown" = "#d9d9d9")

# ── Local function: normalize patient IDs (matches project convention) ──
norm_id <- function(x) {
  out <- as.character(x); out <- trimws(out); out <- gsub("\\.0$", "", out)
  out <- gsub("[^A-Za-z0-9]+", "_", out); out <- gsub("^_+|_+$", "", out); toupper(out)
}

# Collapses the ~12 raw free-text histology strings into 4 broad, readable
# categories for the annotation strip. This grouping is a visualization
# simplification only — the underlying raw string is preserved in the
# joined data frame if finer detail is needed later.
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

# ── Read PDL1_TPS_pct, Stage_Group, Histological_Subtype from raw clinical file ──
message("\n=== 08t | Reading raw clinical workbook for PDL1/Stage/Histology ===")

if (!file.exists(paths$clinical)) {
  stop("Clinical file not found at: ", paths$clinical)
}
if (file.size(paths$clinical) < 5000) {
  stop("Clinical file at ", paths$clinical, " is suspiciously small (",
       file.size(paths$clinical), " bytes) for an xlsx workbook \u2014 ",
       "this usually means it's a OneDrive cloud-only placeholder that hasn't ",
       "fully synced. In File Explorer, right-click the file and choose ",
       "'Always keep on this device', wait for the sync icon to turn into a ",
       "solid checkmark, then re-run this script.")
}

sheet_names <- tryCatch(
  readxl::excel_sheets(paths$clinical),
  error = function(e) {
    stop("Could not open ", paths$clinical, " as an Excel workbook (readxl error: ",
         conditionMessage(e), "). File exists and reports a plausible size, but the ",
         "content isn't readable \u2014 most likely a OneDrive sync placeholder ",
         "(check the file's cloud/sync icon in File Explorer) or a corrupted/partially ",
         "downloaded copy. Re-sync or re-download the file and try again before assuming ",
         "this is a script problem.")
  }
)
target_sheet <- sheet_names[grepl("clinical|clean", sheet_names, ignore.case = TRUE)][1]
if (is.na(target_sheet)) target_sheet <- sheet_names[1]
message("Using sheet: ", target_sheet)

clin_raw <- readxl::read_excel(paths$clinical, sheet = target_sheet) %>%
  dplyr::transmute(
    patient_id             = norm_id(Patient_ID),
    stage_group_raw        = as.character(Stage_Group),
    histological_subtype   = as.character(Histological_Subtype),
    pdl1_tps_pct           = suppressWarnings(as.numeric(PDL1_TPS_pct))
  ) %>%
  dplyr::mutate(
    stage_group_display = dplyr::case_when(
      stage_group_raw %in% c("I", "II", "III", "IV") ~ stage_group_raw,
      TRUE ~ "Unknown"
    ),
    histology_display = simplify_histology(histological_subtype)
  )

message("Clinical rows read: ", nrow(clin_raw),
        " | non-missing PDL1_TPS_pct: ", sum(!is.na(clin_raw$pdl1_tps_pct)))

# ── Join to analytic cohort and filter to non-missing PDL1 ──────
message("\n=== 08t | Joining to analytic cohort (144-patient set) and filtering ===")

df <- analysis1_df %>%
  dplyr::distinct(patient_id, sample_id, new_group, response_4cat) %>%
  dplyr::inner_join(clin_raw, by = "patient_id") %>%
  dplyr::filter(!is.na(pdl1_tps_pct)) %>%
  dplyr::mutate(
    response_display   = dplyr::case_when(response_4cat %in% c("CR","PR","SD","PD") ~ response_4cat, TRUE ~ "NE"),
    new_group_display   = factor(new_group, levels = names(GROUP_COLOURS)),
    response_display     = factor(response_display, levels = names(PAL_RESPONSE)),
    stage_group_display = factor(stage_group_display, levels = names(PAL_STAGE)),
    histology_display   = factor(histology_display, levels = names(PAL_HIST))
  ) %>%
  dplyr::arrange(pdl1_tps_pct) %>%
  dplyr::mutate(rank = dplyr::row_number())

message("Patients with non-missing PDL1_TPS_pct in analytic cohort: ", nrow(df),
        " (of ", dplyr::n_distinct(analysis1_df$patient_id), " total)")
message("Excluded for missing PDL1_TPS_pct: ", dplyr::n_distinct(analysis1_df$patient_id) - nrow(df))

# ── Main bar plot: fill = treatment; response shown as a separate strip below ──
message("\n=== 08t | Building annotated bar plot ===")

bar_plot <- ggplot(df, aes(x = rank, y = pdl1_tps_pct, fill = new_group_display)) +
  geom_col(width = 1, colour = NA) +
  scale_fill_manual(values = GROUP_COLOURS, name = "Treatment") +
  labs(
    title = paste0("PD-L1 TPS% Across Cohort (n=", nrow(df), " with non-missing PDL1_TPS_pct)"),
    subtitle = "Sorted ascending by PD-L1 TPS%; fill = treatment group (see Response/Stage/Histology strips below)",
    x = NULL, y = "PD-L1 TPS (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey30")
  )

# ── Annotation strips: Response, Stage_Group, Histological_Subtype ──
response_strip <- ggplot(df, aes(x = rank, y = "Response", fill = response_display)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = PAL_RESPONSE, name = "Response") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

stage_strip <- ggplot(df, aes(x = rank, y = "Stage", fill = stage_group_display)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = PAL_STAGE, name = "Stage") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

histology_strip <- ggplot(df, aes(x = rank, y = "Histology", fill = histology_display)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = PAL_HIST, name = "Histology") +
  labs(x = "Patients (sorted by PD-L1 TPS%, ascending)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

# ── Combine, legend at bottom ────────────────────────────────────
combined <- (bar_plot / response_strip / stage_strip / histology_strip) +
  patchwork::plot_layout(heights = c(5, 1, 1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical")

stub <- file.path(OUT_FIG, "pdl1_tps_annotated_barplot")
cairo_pdf(paste0(stub, ".pdf"), width = 13, height = 8.5); print(combined); dev.off()
png(paste0(stub, ".png"), width = 13, height = 8.5, units = "in", res = 300, bg = "white"); print(combined); dev.off()

message("Saved: ", stub, ".pdf / .png")
message("\n08t_pdl1_tps_annotated_barplot.R complete")
