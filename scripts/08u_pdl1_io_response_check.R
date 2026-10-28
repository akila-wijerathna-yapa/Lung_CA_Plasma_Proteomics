if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716"
}
sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

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

# ── Read PDL1_TPS_pct from raw clinical workbook (defensive) ───
message("\n=== 08u | Reading PDL1_TPS_pct from raw clinical workbook ===")

if (!file.exists(paths$clinical)) stop("Clinical file not found at: ", paths$clinical)
if (file.size(paths$clinical) < 5000) {
  stop("Clinical file at ", paths$clinical, " is suspiciously small (",
       file.size(paths$clinical), " bytes) \u2014 likely a OneDrive sync placeholder. ",
       "Right-click the file in File Explorer > 'Always keep on this device', wait for sync, retry.")
}
sheet_names <- tryCatch(readxl::excel_sheets(paths$clinical), error = function(e) {
  stop("Could not open ", paths$clinical, " as an Excel workbook (", conditionMessage(e),
       ") \u2014 check OneDrive sync status / file integrity before assuming a script problem.")
})
target_sheet <- sheet_names[grepl("clinical|clean", sheet_names, ignore.case = TRUE)][1]
if (is.na(target_sheet)) target_sheet <- sheet_names[1]

clin_raw <- readxl::read_excel(paths$clinical, sheet = target_sheet) %>%
  dplyr::transmute(
    patient_id   = norm_id(Patient_ID),
    pdl1_tps_pct = suppressWarnings(as.numeric(PDL1_TPS_pct))
  )

message("Clinical rows read: ", nrow(clin_raw),
        " | non-missing PDL1_TPS_pct: ", sum(!is.na(clin_raw$pdl1_tps_pct)))

# ── Restrict to IO / CTx+IO, non-missing PDL1 and response ─────
message("\n=== 08u | Restricting to IO-containing regimens ===")

df <- analysis1_df %>%
  dplyr::distinct(patient_id, sample_id, new_group, response_binary, response_4cat) %>%
  dplyr::inner_join(clin_raw, by = "patient_id") %>%
  dplyr::filter(new_group %in% c("IO", "CTx + IO")) %>%
  dplyr::filter(!is.na(pdl1_tps_pct), !is.na(response_binary))

message("IO-containing patients with non-missing PDL1 AND response: ", nrow(df))
message("Breakdown by treatment group:")
print(df %>% dplyr::count(new_group))
message("Breakdown by response:")
print(df %>% dplyr::count(response_binary))

if (nrow(df) < 10) {
  warning("n=", nrow(df), " is very small for this subgroup test \u2014 interpret any result ",
          "as hypothesis-generating only, not confirmatory.")
}

# ── Test 1: PDL1_TPS_pct (continuous) by response_binary ────────
message("\n=== 08u | Wilcoxon test: PDL1_TPS_pct by response (IO-containing only) ===")

wt <- if (nrow(df) >= 6 && dplyr::n_distinct(df$response_binary) == 2) {
  wilcox.test(pdl1_tps_pct ~ response_binary, data = df)
} else NULL

wilcox_summary <- tibble::tibble(
  n = nrow(df),
  n_sdpd = sum(df$response_binary == "SD/PD"),
  n_prcr = sum(df$response_binary == "PR/CR"),
  median_sdpd = median(df$pdl1_tps_pct[df$response_binary == "SD/PD"]),
  median_prcr = median(df$pdl1_tps_pct[df$response_binary == "PR/CR"]),
  wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_
)
message("Result:"); print(wilcox_summary)

# ── Test 2: standard clinical TPS categories vs response (Fisher's exact) ──
message("\n=== 08u | Standard PD-L1 TPS categories vs response (Fisher's exact) ===")

df <- df %>%
  dplyr::mutate(
    pdl1_category = dplyr::case_when(
      pdl1_tps_pct < 1  ~ "<1%",
      pdl1_tps_pct < 50 ~ "1-49%",
      TRUE              ~ "\u226550%"
    ),
    pdl1_category = factor(pdl1_category, levels = c("<1%", "1-49%", "\u226550%"))
  )

cat_table <- table(df$pdl1_category, df$response_binary)
message("Cross-tabulation (PD-L1 category x response):")
print(cat_table)

fisher_result <- tryCatch(fisher.test(cat_table), error = function(e) NULL)
fisher_p <- if (!is.null(fisher_result)) fisher_result$p.value else NA_real_
message("Fisher's exact test p-value: ", signif(fisher_p, 3))

# ── Save combined summary ───────────────────────────────────────
safe_write_csv(
  dplyr::bind_rows(
    wilcox_summary %>% dplyr::mutate(test = "continuous_wilcoxon"),
    tibble::tibble(n = nrow(df), fisher_p = fisher_p, test = "categorical_fisher")
  ),
  file.path(OUT_TBL, "pdl1_io_response_check.csv")
)
safe_write_csv(as.data.frame.matrix(cat_table) %>% tibble::rownames_to_column("pdl1_category"),
               file.path(OUT_TBL, "pdl1_io_response_crosstab.csv"))

# ── Figure: boxplot + jitter, colored by treatment subgroup ─────
message("\n=== 08u | Building figure ===")

p <- ggplot(df, aes(x = response_binary, y = pdl1_tps_pct)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "grey85") +
  geom_jitter(aes(colour = new_group), width = 0.15, size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = c("IO" = "#56B4E9", "CTx + IO" = "#009E73"), name = "Treatment") +
  labs(
    title = "PD-L1 TPS% vs. Response \u2014 IO-Containing Regimens Only",
    subtitle = paste0("n=", nrow(df), " (IO + CTx+IO) | Wilcoxon p=",
                       if (!is.null(wt)) signif(wt$p.value, 3) else "NA",
                       " | Fisher's exact (TPS categories) p=", signif(fisher_p, 3)),
    x = NULL, y = "PD-L1 TPS (%)"
  ) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, colour = "grey30"))

stub <- file.path(OUT_FIG, "pdl1_io_response_check")
cairo_pdf(paste0(stub, ".pdf"), width = 7.5, height = 5.5); print(p); dev.off()
png(paste0(stub, ".png"), width = 7.5, height = 5.5, units = "in", res = 300, bg = "white"); print(p); dev.off()

message("Saved: ", stub, ".pdf / .png")
message("\n08u_pdl1_io_response_check.R complete")
