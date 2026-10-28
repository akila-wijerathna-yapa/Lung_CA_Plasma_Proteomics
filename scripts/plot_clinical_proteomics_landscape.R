suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(lubridate)
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork", quiet = TRUE)
  if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale", quiet = TRUE)
  if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot", quiet = TRUE)
  library(patchwork)
  library(ggnewscale)
  library(cowplot)
})

CLIN_PATH <- r"(C:\Users\uqamudiy\OneDrive - The University of Queensland\Documents\Nomic_Lung_CA\nELISA_20260716\data\Clinical_Akila_MERGED_subset.xlsx)"
MAP_PATH  <- r"(C:\Users\uqamudiy\OneDrive - The University of Queensland\Documents\Nomic_Lung_CA\nELISA_20260716\data\sample_IDs_before_after.xlsx)"
OUT_STUB  <- r"(C:\Users\uqamudiy\OneDrive - The University of Queensland\Documents\Nomic_Lung_CA\figures\01_clinical_proteomics_landscape)"
dir.create(dirname(OUT_STUB), recursive = TRUE, showWarnings = FALSE)

# ---- helpers (verbatim from 01_run_nELISA_clinical_analysis.R) ----
norm_id <- function(x) {
  out <- as.character(x)
  out <- trimws(out)
  out <- gsub("\\.0$", "", out)
  out <- gsub("[^A-Za-z0-9]+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  toupper(out)
}

parse_date_mixed <- function(x) {
  parse_one <- function(v) {
    if (length(v) == 0 || is.na(v) || trimws(as.character(v)) %in% c("", "NA", "NaN", "Not applicable")) {
      return(as.Date(NA))
    }
    if (inherits(v, "Date")) return(as.Date(v))
    if (inherits(v, "POSIXt")) return(as.Date(v))
    if (is.numeric(v)) {
      if (v > 59 && v < 80000) return(as.Date(floor(v), origin = "1899-12-30"))
      return(as.Date(NA))
    }
    txt <- trimws(as.character(v))
    parsed <- suppressWarnings(lubridate::parse_date_time(
      txt, orders = c("dmy", "ymd", "mdy", "d-b-y", "d-b-Y", "d/m/Y", "m/d/Y", "Y-m-d")
    ))
    if (length(parsed) == 0 || is.na(parsed[1])) return(as.Date(NA))
    as.Date(parsed[1])
  }
  as.Date(vapply(x, parse_one, as.Date(NA)))
}

# ---- build_clinical (verbatim logic from 01_run_nELISA_clinical_analysis.R) ----
build_clinical <- function(path) {
  raw <- as_tibble(readxl::read_excel(path, sheet = 1))

  group_map <- c(
    "Chemotherapy" = "CTx",
    "Surgery" = "Sx",
    "Immunotherapy" = "IO",
    "Chemo + IO" = "CTx + IO",
    "Chemo + RT" = "CTx + RTx",
    "TKI" = "TKI",
    "Chemo + IO + RT" = "Other",
    "RT only" = "Other",
    "RT + IO (sequential)" = "Other",
    "IO + RT" = "Other"
  )

  clin <- raw %>%
    transmute(
      patient_id = norm_id(Patient_ID),
      patient_id_raw = as.character(Patient_ID),
      sex = case_when(
        toupper(trimws(as.character(Sex))) %in% c("M", "MALE") ~ "MALE",
        toupper(trimws(as.character(Sex))) %in% c("F", "FEMALE") ~ "FEMALE",
        TRUE ~ NA_character_
      ),
      dob = parse_date_mixed(DOB),
      stage_group = as.character(Stage_Group),
      histology_raw = as.character(Histological_Subtype),
      diagnosis_date = parse_date_mixed(Date_Diagnosis),
      treatment_class = as.character(`1L_Treatment_Class`),
      new_group = unname(group_map[as.character(`1L_Treatment_Class`)]),
      response_4cat = case_when(
        trimws(as.character(`1L_Response_Std`)) %in% c("CR", "PR", "SD", "PD") ~ trimws(as.character(`1L_Response_Std`)),
        TRUE ~ NA_character_
      ),
      pfs_months_raw = suppressWarnings(as.numeric(`PFS_months`)),
      pfs_censoring_type = as.character(`PFS_Censoring_Type`),
      os_months = suppressWarnings(as.numeric(OS_months)),
      alive_code = suppressWarnings(as.integer(`Alive_1_Dead_0`)),
      os_event = suppressWarnings(as.integer(`OS_Event_1Dead`))
    ) %>%
    mutate(
      age_at_diagnosis = ifelse(!is.na(dob) & !is.na(diagnosis_date), floor(as.numeric(diagnosis_date - dob) / 365.25), NA_real_),
      pfs_censor_class = case_when(
        str_detect(pfs_censoring_type, "^Exact") ~ "exact",
        str_detect(pfs_censoring_type, "^Right-censored") ~ "right_censored",
        str_detect(pfs_censoring_type, "^Interval-censored") ~ "interval_censored",
        str_detect(pfs_censoring_type, "^Not applicable") ~ "not_applicable",
        TRUE ~ "unknown"
      )
    )

  pfs_event_main <- case_when(
    clin$pfs_censor_class == "exact" ~ 1L,
    clin$pfs_censor_class == "right_censored" ~ 0L,
    TRUE ~ NA_integer_
  )

  clin %>%
    mutate(
      pfs_months_main = pfs_months_raw,
      pfs_event_main = pfs_event_main,
      alive_status = case_when(
        alive_code == 1 ~ "Alive",
        alive_code == 0 ~ "Dead",
        TRUE ~ NA_character_
      )
    )
}

# ---- pre/post sample map -> real pair_status ----
read_sample_map <- function(path) {
  raw <- readxl::read_excel(path, sheet = 1, col_names = FALSE, skip = 3)
  pre_ids  <- norm_id(raw[[1]][!is.na(raw[[1]])])
  post_ids <- norm_id(raw[[4]][!is.na(raw[[4]])])
  list(pre = pre_ids, post = post_ids)
}

clinical <- build_clinical(CLIN_PATH)
sample_map <- read_sample_map(MAP_PATH)

clinical <- clinical %>%
  mutate(
    pair_status = case_when(
      patient_id %in% sample_map$pre & patient_id %in% sample_map$post ~ "both",
      patient_id %in% sample_map$pre ~ "only_pre",
      patient_id %in% sample_map$post ~ "only_post",
      TRUE ~ "no_proteomics"
    )
  )

message("Loaded ", nrow(clinical), " patients.")
message("pair_status counts:"); print(table(clinical$pair_status, useNA = "ifany"))

# ---- landscape figure ----
make_clinical_landscape_figure <- function(clinical, path_stub) {
  group_levels    <- c("Sx", "CTx", "CTx + IO", "CTx + RTx", "IO", "TKI", "Other", "Unknown")
  response_levels <- c("CR", "PR", "SD", "PD", "NE")
  stage_levels    <- c("I-II", "III-IV", "Unknown")
  histology_levels <- c("Adenocarcinoma", "Squamous Cell Carcinoma", "Other/Mixed", "Unknown")
  pair_levels     <- c("both", "only_pre", "only_post", "no_proteomics")
  strip_levels    <- rev(c("Group", "Sex", "Response", "Stage", "Histology", "Vital", "PFS", "OS", "Pair"))

  plot_df <- clinical %>%
    mutate(
      group_display = coalesce(new_group, "Unknown"),
      sex_display = case_when(
        toupper(coalesce(sex, "")) %in% c("F", "FEMALE") ~ "Female",
        toupper(coalesce(sex, "")) %in% c("M", "MALE") ~ "Male",
        TRUE ~ "Unknown"
      ),
      response_display = case_when(
        response_4cat %in% c("CR", "PR", "SD", "PD") ~ response_4cat,
        TRUE ~ "NE"
      ),
      stage_display = case_when(
        stage_group %in% c("I", "II", "I-II", "I/II", "Stage I-II", "Stage I/II") ~ "I-II",
        stage_group %in% c("III", "IV", "III-IV", "III/IV", "Stage III-IV", "Stage III/IV") ~ "III-IV",
        TRUE ~ "Unknown"
      ),
      histology_display = case_when(
        is.na(histology_raw) ~ "Unknown",
        str_detect(histology_raw, regex("squamous", ignore_case = TRUE)) &
          !str_detect(histology_raw, regex("adeno", ignore_case = TRUE)) ~ "Squamous Cell Carcinoma",
        str_detect(histology_raw, regex("adeno", ignore_case = TRUE)) &
          !str_detect(histology_raw, regex("squamous|sarcomatoid|neuroendocrine|large cell|small cell", ignore_case = TRUE)) ~ "Adenocarcinoma",
        TRUE ~ "Other/Mixed"
      ),
      vital_display = case_when(
        alive_status %in% c("Alive", "Dead") ~ alive_status,
        TRUE ~ "Unknown"
      ),
      pfs_display = case_when(
        is.na(pfs_event_main) | is.na(pfs_months_main) ~ "PFS unavailable",
        pfs_event_main == 1 ~ "PFS event",
        TRUE ~ "PFS censored"
      ),
      os_display = case_when(
        is.na(os_event) | is.na(os_months) ~ "OS unavailable",
        os_event == 1 ~ "OS death",
        TRUE ~ "OS censored"
      )
    ) %>%
    arrange(
      factor(pair_status, levels = pair_levels),
      factor(group_display, levels = group_levels),
      factor(response_display, levels = response_levels),
      factor(stage_display, levels = stage_levels),
      age_at_diagnosis
    ) %>%
    mutate(
      idx = row_number(),
      group_display    = factor(group_display, levels = group_levels),
      sex_display      = factor(sex_display, levels = c("Female", "Male", "Unknown")),
      response_display = factor(response_display, levels = response_levels),
      stage_display     = factor(stage_display, levels = stage_levels),
      histology_display = factor(histology_display, levels = histology_levels),
      vital_display     = factor(vital_display, levels = c("Alive", "Dead", "Unknown")),
      pfs_display       = factor(pfs_display, levels = c("PFS event", "PFS censored", "PFS unavailable")),
      os_display        = factor(os_display, levels = c("OS death", "OS censored", "OS unavailable")),
      pair_status_f     = factor(pair_status, levels = pair_levels)
    )

  if (nrow(plot_df) == 0) {
    warning("plot_df has 0 rows - nothing to plot.")
    return(invisible(NULL))
  }

  pal_group    <- c("Sx" = "#f4a261", "CTx" = "#e76f51", "CTx + IO" = "#2a9d8f", "CTx + RTx" = "#457b9d",
                     "IO" = "#8d99ae", "TKI" = "#6a4c93", "Other" = "#7f7f7f", "Unknown" = "#d9d9d9")
  pal_sex      <- c("Female" = "#f4a6b7", "Male" = "#355070", "Unknown" = "#d9d9d9")
  pal_response <- c("CR" = "#7cb342", "PR" = "#2a9d8f", "SD" = "#f4a261", "PD" = "#d7301f", "NE" = "#d0d0d0")
  pal_stage    <- c("I-II" = "#fbb4ae", "III-IV" = "#80cdc1", "Unknown" = "#d9d9d9")
  pal_histology <- c("Adenocarcinoma" = "#1b9e77", "Squamous Cell Carcinoma" = "#d95f02",
                      "Other/Mixed" = "#7570b3", "Unknown" = "#d9d9d9")
  pal_vital    <- c("Alive" = "#66c2a5", "Dead" = "#d73027", "Unknown" = "#d9d9d9")
  pal_pfs      <- c("PFS event" = "#d73027", "PFS censored" = "#d9d9d9", "PFS unavailable" = "#f0f0f0")
  pal_os       <- c("OS death" = "#b2182b", "OS censored" = "#d9d9d9", "OS unavailable" = "#f0f0f0")
  pal_pair     <- c("both" = "#4d4d4d", "only_pre" = "#fdae61", "only_post" = "#2ca25f", "no_proteomics" = "#ececec")

  id_breaks <- plot_df$idx
  id_labels <- plot_df$patient_id_raw

  age_plot <- ggplot(plot_df, aes(idx, age_at_diagnosis)) +
    geom_point(size = 1.05, color = "black", alpha = 0.85) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.45, color = "#7f7f7f", span = 0.3) +
    labs(
      title = paste0("Clinical and Proteomics Cohort Landscape (n=", nrow(plot_df), ")"),
      subtitle = paste0(
        "Pair status: both=", sum(plot_df$pair_status == "both"),
        ", only_pre=", sum(plot_df$pair_status == "only_pre"),
        ", only_post=", sum(plot_df$pair_status == "only_post"),
        ", no_proteomics=", sum(plot_df$pair_status == "no_proteomics")
      ),
      x = NULL, y = "Age"
    ) +
    scale_x_continuous(breaks = id_breaks, labels = id_labels, expand = expansion(mult = c(0.002, 0.002))) +
    theme_bw(base_size = 10.5) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5, color = "grey30")
    )

  # data panel: no legends here, they're extracted separately below
  tile_plot <- ggplot(plot_df, aes(x = idx)) +
    geom_tile(aes(y = "Group", fill = group_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_group, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "Sex", fill = sex_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_sex, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "Response", fill = response_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_response, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "Stage", fill = stage_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_stage, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "Histology", fill = histology_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_histology, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "Vital", fill = vital_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_vital, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "PFS", fill = pfs_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_pfs, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "OS", fill = os_display), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_os, drop = FALSE, guide = "none") +

    new_scale_fill() +
    geom_tile(aes(y = "Pair", fill = pair_status_f), color = "white", linewidth = 0.12, width = 1) +
    scale_fill_manual(values = pal_pair, drop = FALSE, guide = "none") +

    scale_y_discrete(limits = strip_levels) +
    scale_x_continuous(breaks = id_breaks, labels = id_labels, expand = expansion(mult = c(0.002, 0.002))) +
    labs(x = "Patient ID", y = NULL) +
    theme_minimal(base_size = 10.5) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      axis.ticks.x = element_line(linewidth = 0.2),
      panel.grid = element_blank(),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none"
    )

  # ---- build one standalone legend grob per panel ----
  make_legend <- function(values_vec, title) {
    df <- data.frame(x = 1, y = 1, val = factor(names(values_vec), levels = names(values_vec)))
    p <- ggplot(df, aes(x, y, fill = val)) +
      geom_tile() +
      scale_fill_manual(
        values = values_vec, name = title, drop = FALSE,
        guide = guide_legend(
          title.position = "top", title.theme = element_text(face = "bold", size = 8),
          direction = "horizontal", nrow = 1, byrow = TRUE, override.aes = list(color = NA)
        )
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.key.height = grid::unit(0.3, "cm"),
        legend.key.width = grid::unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        legend.margin = ggplot2::margin(0, 0, 0, 0)
      )
    cowplot::get_legend(p)
  }

  legend_group    <- make_legend(pal_group, "Group")
  legend_sex      <- make_legend(pal_sex, "Sex")
  legend_response <- make_legend(pal_response, "Response")
  legend_stage    <- make_legend(pal_stage, "Stage")
  legend_histology <- make_legend(pal_histology, "Histology")
  legend_vital    <- make_legend(pal_vital, "Vital")
  legend_pfs      <- make_legend(pal_pfs, "PFS")
  legend_os       <- make_legend(pal_os, "OS")
  legend_pair     <- make_legend(pal_pair, "Pair")

  # ---- arrange legends into the 4 requested rows, spread horizontally ----
  row1 <- cowplot::plot_grid(legend_group, legend_response, nrow = 1)
  row2 <- cowplot::plot_grid(legend_vital, legend_pair, legend_sex, nrow = 1)
  row3 <- cowplot::plot_grid(legend_os, legend_stage, legend_pfs, nrow = 1)
  row4 <- cowplot::plot_grid(legend_histology, nrow = 1)
  legend_block <- cowplot::plot_grid(row1, row2, row3, row4, ncol = 1)

  data_block <- age_plot / tile_plot + plot_layout(heights = c(2.1, 4.2))
  combined <- wrap_elements(full = data_block) / legend_block + plot_layout(heights = c(6.3, 3.4))

  ggsave(paste0(path_stub, ".png"), combined, width = 13.2, height = 10.2, dpi = 300, bg = "white")
  ggsave(paste0(path_stub, ".pdf"), combined, width = 13.2, height = 10.2, device = cairo_pdf, bg = "white")

  message("Saved: ", path_stub, ".png / .pdf")
  invisible(plot_df)
}

plot_df <- make_clinical_landscape_figure(clinical, OUT_STUB)
