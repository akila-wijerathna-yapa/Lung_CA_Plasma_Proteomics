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
GROUP_COLOURS <- c("Sx" = "#E69F00", "CTx" = "#D55E00", "CTx + IO" = "#009E73", "CTx + RTx" = "#0072B2",
                   "IO" = "#56B4E9", "TKI" = "#CC79A7", "Other" = "#999999")
PAL_RESP2 <- c("SD/PD" = "#a6bddb", "PR/CR" = "#74c476")

akila_theme <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          strip.text = element_text(face = "bold"))
}

# ── Pull CTSA / CTSL / LAMP1 from the already-QC'd analytic panel ──
message("\n=== 08x | Extracting CTSA, CTSL, LAMP1 from expr_log2_rep ===")

panel_ids <- feature_meta %>%
  dplyr::filter(Gene_Symbol %in% c("CTSA", "CTSL", "LAMP1")) %>%
  dplyr::distinct(Gene_Symbol, nELISA_ID)
message("Found in analytic panel: ", paste(panel_ids$Gene_Symbol, collapse = ", "))

get_expr <- function(id) { v <- expr_log2_rep[, id]; names(v) <- rownames(expr_log2_rep); v }

qc_long <- purrr::pmap_dfr(panel_ids, function(Gene_Symbol, nELISA_ID) {
  v <- get_expr(nELISA_ID)
  tibble::tibble(sample_id = names(v), protein = Gene_Symbol, log2_value = as.numeric(v))
})

# ── Re-extract MAP1LC3A from the raw original proteomics workbook ──
message("\n=== 08x | Re-extracting MAP1LC3A from raw proteomics workbook ===")
message("NOTE: MAP1LC3A was filtered out during QC and is NOT in expr_log2_rep or any ",
        "checkpoint. It is pulled fresh here and log2-transformed independently — it did ",
        "NOT go through the same duplicate-resolution/QC process as CTSA/CTSL/LAMP1. ",
        "Treat any MAP1LC3A result with more caution than the other three proteins.")

if (!file.exists(paths$proteomics)) stop("Proteomics file not found at: ", paths$proteomics)
if (file.size(paths$proteomics) < 5000) {
  stop("Proteomics file at ", paths$proteomics, " is suspiciously small (",
       file.size(paths$proteomics), " bytes) \u2014 likely a OneDrive sync placeholder. ",
       "Right-click in File Explorer > 'Always keep on this device', wait for sync, retry.")
}

raw <- tryCatch(
  readxl::read_excel(paths$proteomics, sheet = 1, col_names = FALSE),
  error = function(e) stop("Could not open ", paths$proteomics, " (", conditionMessage(e), ") \u2014 ",
                            "check OneDrive sync/file integrity before assuming a script problem.")
)

# Header structure (matches the pipeline's own raw proteomics layout):
# row1=nELISA ID, row2=Uniprot ID, row3=Protein Name, row4=ULOD, row5=LLOD, row6+=data
nelisa_row <- as.character(raw[1, ])
name_row   <- as.character(raw[3, ])
maplc3a_col <- which(name_row == "MAP1LC3A")
if (length(maplc3a_col) != 1) stop("MAP1LC3A not found at the expected header row \u2014 raw file layout may differ from expected.")

data_rows <- raw[6:nrow(raw), ]
sample_id_col <- which(as.character(raw[5, ]) == "sample_id_")
# Fallback: sample_id_ is a data column header, not in row 5 (LLOD row) — it's in the metadata
# header, which per the pipeline's own layout is row 5 for column names of the first 10 metadata cols.
if (length(sample_id_col) == 0) {
  meta_header <- as.character(raw[5, 1:10])
  sample_id_col <- which(meta_header == "sample_id_")
}

maplc3a_raw <- tibble::tibble(
  sample_id = as.character(data_rows[[sample_id_col]]),
  pg_ml     = suppressWarnings(as.numeric(unlist(data_rows[[maplc3a_col]])))
) %>%
  dplyr::filter(!is.na(sample_id), sample_id != "")

n_match <- sum(maplc3a_raw$sample_id %in% rownames(expr_log2_rep))
message("MAP1LC3A raw rows read: ", nrow(maplc3a_raw),
        " | matching current analytic sample_id set: ", n_match, " of ", length(rownames(expr_log2_rep)))
if (n_match < 0.5 * length(rownames(expr_log2_rep))) {
  warning("Fewer than half of analytic samples matched by sample_id_ \u2014 check that the raw ",
          "file's sample identifier format matches rownames(expr_log2_rep) before trusting this result.")
}

maplc3a_long <- maplc3a_raw %>%
  dplyr::filter(sample_id %in% rownames(expr_log2_rep)) %>%
  dplyr::mutate(protein = "MAP1LC3A", log2_value = log2(pmax(pg_ml, 1e-3))) %>%
  dplyr::select(sample_id, protein, log2_value)

# ── Combine + join clinical annotation ───────────────────────────
message("\n=== 08x | Combining and joining clinical annotation ===")

expr_long <- dplyr::bind_rows(qc_long, maplc3a_long) %>%
  dplyr::mutate(protein = factor(protein, levels = c("CTSA", "CTSL", "LAMP1", "MAP1LC3A")))

plot_df <- expr_long %>%
  dplyr::inner_join(
    analysis1_df %>% dplyr::distinct(sample_id, new_group, response_binary),
    by = "sample_id"
  )

safe_write_csv(plot_df, file.path(OUT_TBL, "autophagy_proteins_data.csv"))
message("Final combined n rows: ", nrow(plot_df), " across ", dplyr::n_distinct(plot_df$sample_id), " patients")

# ── Figure 1: by response ────────────────────────────────────────
message("\n=== 08x | Building response boxplot ===")

df_resp <- plot_df %>% dplyr::filter(!is.na(response_binary))

resp_stats <- df_resp %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(p = tryCatch(wilcox.test(log2_value ~ response_binary)$p.value, error = function(e) NA_real_), .groups = "drop")

df_resp <- df_resp %>% dplyr::left_join(resp_stats, by = "protein") %>%
  dplyr::mutate(facet_label = paste0(protein, " (Wilcoxon p=", signif(p, 3), ")"))

p_resp <- ggplot(df_resp, aes(x = response_binary, y = log2_value, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_fill_manual(values = PAL_RESP2, guide = "none") +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
  labs(title = "CTSA, CTSL, LAMP1, MAP1LC3A by Treatment Response",
       subtitle = "MAP1LC3A re-extracted from raw panel, independently log2-transformed (not QC-matched to the other three)",
       x = NULL, y = "log2 concentration") +
  akila_theme()

stub_resp <- file.path(OUT_FIG, "autophagy_proteins_by_response")
cairo_pdf(paste0(stub_resp, ".pdf"), width = 9, height = 7.5); print(p_resp); dev.off()
png(paste0(stub_resp, ".png"), width = 9, height = 7.5, units = "in", res = 300, bg = "white"); print(p_resp); dev.off()
message("Saved: ", stub_resp, ".pdf / .png")

# ── Figure 2: by treatment ───────────────────────────────────────
message("\n=== 08x | Building treatment boxplot ===")

df_tx <- plot_df %>% dplyr::filter(new_group %in% names(GROUP_COLOURS)) %>%
  dplyr::mutate(new_group = factor(new_group, levels = names(GROUP_COLOURS)))

tx_stats <- df_tx %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(p = tryCatch(kruskal.test(log2_value ~ new_group)$p.value, error = function(e) NA_real_), .groups = "drop")

df_tx <- df_tx %>% dplyr::left_join(tx_stats, by = "protein") %>%
  dplyr::mutate(facet_label = paste0(protein, " (Kruskal-Wallis p=", signif(p, 3), ")"))

p_tx <- ggplot(df_tx, aes(x = new_group, y = log2_value, fill = new_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
  labs(title = "CTSA, CTSL, LAMP1, MAP1LC3A by Treatment Group",
       subtitle = "MAP1LC3A re-extracted from raw panel, independently log2-transformed (not QC-matched to the other three)",
       x = NULL, y = "log2 concentration") +
  akila_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

stub_tx <- file.path(OUT_FIG, "autophagy_proteins_by_treatment")
cairo_pdf(paste0(stub_tx, ".pdf"), width = 10, height = 8); print(p_tx); dev.off()
png(paste0(stub_tx, ".png"), width = 10, height = 8, units = "in", res = 300, bg = "white"); print(p_tx); dev.off()
message("Saved: ", stub_tx, ".pdf / .png")


