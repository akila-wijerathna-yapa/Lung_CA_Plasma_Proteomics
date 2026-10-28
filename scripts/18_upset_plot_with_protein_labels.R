# ── Preamble ─────────────────────────────────────────────────
PROJECT_DIR <- normalizePath(file.path(dirname(
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    rstudioapi::getActiveDocumentContext()$path
  else {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- grep("^--file=", args, value = TRUE)
    if (length(f)) sub("^--file=", "", f[1]) else stop("Cannot detect script path.")
  }
), ".."), mustWork = FALSE)

sys.source(file.path(PROJECT_DIR, "scripts", "00_setup.R"), envir = globalenv())

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# ── Load checkpoints ───────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
message("Checkpoint 01 loaded")

checkpoint_16_path <- file.path(paths$rdata, "checkpoint_16.RData")
if (!file.exists(checkpoint_16_path))
  stop("checkpoint_16.RData not found. Run 16_site_specific_biomarker_overlap.R first.")
load(checkpoint_16_path)
message("Checkpoint 16 loaded — ", length(site_hits), " sites with hits")
while (dev.cur() > 1) dev.off()

fig_dir <- file.path(paths$figures, "14_site_specific_biomarkers")
tbl_dir <- file.path(paths$tables,  "14_site_specific_biomarkers")
invisible(lapply(c(fig_dir, tbl_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Local constants ────────────────────────────────────────────
LABEL_TRUNCATE_N <- Inf   # set e.g. 20 to cap names drawn per box; Inf = show all
LABEL_FONT_SIZE  <- 2.2
ROW_HEIGHT_IN    <- 0.14  # inches per stacked name row, used for figure sizing

# ── Local utility functions ────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

# ============================================================
# Derive bins (exact site-combinations) from membership_df
# ============================================================
message("\n=== Deriving intersection bins ===")
site_names <- names(site_hits)

bin_df <- membership_df %>%
  mutate(bin_key = apply(across(all_of(site_names)), 1,
                        function(row) paste(site_names[row == 1], collapse = "+")))

# Per-protein sort p-value = min nominal p-value across the site(s) in its own bin
all_site_de <- bind_rows(site_de_results) %>% select(nELISA_ID, site, p_value)
bin_df$sort_p <- vapply(seq_len(nrow(bin_df)), function(i) {
  row_sites <- site_names[as.logical(bin_df[i, site_names])]
  pid <- bin_df$nELISA_ID[i]
  ps  <- all_site_de %>% filter(nELISA_ID == pid, site %in% row_sites) %>% pull(p_value)
  if (length(ps)) min(ps, na.rm = TRUE) else NA_real_
}, numeric(1))

bin_table <- bin_df %>%
  arrange(sort_p) %>%
  group_by(bin_key) %>%
  summarise(
    n_proteins   = n(),
    protein_list = paste(Protein_Name, collapse = ", "),   # already ascending p-value order
    proteins     = list(Protein_Name),                      # ascending p-value order
    p_values     = list(sort_p),
    .groups = "drop"
  ) %>%
  arrange(desc(n_proteins))

safe_write_csv(bin_table %>% select(bin_key, n_proteins, protein_list),
              file.path(tbl_dir, "bin_protein_list.csv"))
message("Bins found: ", nrow(bin_table))
print(bin_table %>% select(bin_key, n_proteins))

bin_levels <- bin_table$bin_key
bin_table  <- bin_table %>% mutate(bin_key = factor(bin_key, levels = bin_levels))

# ============================================================
# Top panel: bars + size number + stacked protein-name labels
# ============================================================
message("\n=== Building bar + label panel ===")

label_df <- bin_table %>%
  mutate(proteins = purrr::map(proteins, ~ head(.x, if (is.infinite(LABEL_TRUNCATE_N)) length(.x) else LABEL_TRUNCATE_N))) %>%
  tidyr::unnest_longer(proteins, values_to = "label") %>%
  group_by(bin_key) %>%
  mutate(rank = row_number(),                     # 1 = lowest p-value (already pre-sorted)
        y    = n_proteins - rank + 0.5) %>%       # top of box = lowest p-value
  ungroup()

p_top <- ggplot() +
  geom_col(data = bin_table, aes(x = bin_key, y = n_proteins),
          fill = "grey92", color = "black", linewidth = 0.5, width = 0.75) +
  geom_text(data = bin_table, aes(x = bin_key, y = n_proteins, label = n_proteins),
           vjust = -0.4, fontface = "bold", size = 3.2) +
  geom_text(data = label_df, aes(x = bin_key, y = y, label = label),
           size = LABEL_FONT_SIZE, vjust = 0.5) +
  scale_x_discrete(limits = bin_levels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(y = "Intersection size (n proteins)", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

# ============================================================
# Bottom panel: standard UpSet dot-matrix (which sites make each bin)
#   — each site gets its own color (filled dots + row label), matching
#     the ComplexUpset-style query coloring requested.
#   NOTE: row labels are drawn with geom_text (not ggtext::element_markdown
#   on the axis) — the markdown-in-theme approach did not render reliably
#   (raw "<span>" tags leaked through), so plain colored geom_text is used
#   instead, which is robust across devices/environments.
#   Layer order matters: grey (unfilled) dots first, then the black
#   connecting line, then colored (filled) dots on top — so the line is
#   never hidden behind grey circles, and colored dots sit cleanly on
#   top of the line at their own position. `site` is forced to the SAME
#   explicit factor levels everywhere (matrix, line, labels) so rows
#   cannot drift out of alignment between panels/layers.
# ============================================================
message("=== Building dot-matrix panel ===")

SITE_COLORS <- setNames(scales::hue_pal()(length(site_names)), site_names)
site_levels <- site_names   # single source of truth for row order, used everywhere below

site_n_pos <- vapply(site_names, function(s) sum(cohort_df[[s]] == 1, na.rm = TRUE), integer(1))
site_label_df <- tibble(
  site  = factor(site_names, levels = site_levels),
  label = paste0(site_names, " (N=", site_n_pos, ")")
)

matrix_df <- tidyr::expand_grid(bin_key = factor(bin_levels, levels = bin_levels), site = site_names) %>%
  mutate(
    site   = factor(site, levels = site_levels),
    filled = purrr::map2_lgl(as.character(bin_key), as.character(site),
                             ~ .y %in% strsplit(.x, "\\+")[[1]])
  )

grey_df   <- matrix_df %>% filter(!filled)
filled_df <- matrix_df %>% filter(filled)
connect_df <- filled_df %>% group_by(bin_key) %>% filter(n() > 1) %>% ungroup()

p_matrix <- ggplot() +
  geom_point(data = grey_df, aes(x = bin_key, y = site), color = "grey85", size = 3.2) +
  geom_line(data = connect_df, aes(x = bin_key, y = site, group = bin_key),
           color = "black", linewidth = 0.7) +
  geom_point(data = filled_df, aes(x = bin_key, y = site, color = site), size = 3.2) +
  scale_color_manual(values = SITE_COLORS, guide = "none") +
  scale_x_discrete(limits = bin_levels) +
  labs(x = "Bin (site combination)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

# Row-label panel — kept as its OWN plot (not text overflowing outside p_matrix's
# panel via clip="off") specifically so patchwork can align it against a matching
# spacer above, and so p_matrix's x-scale expansion is never touched/mismatched
# from p_top's. This is what actually fixes the bar/matrix column misalignment.
p_row_labels <- ggplot(site_label_df, aes(x = 1, y = site, label = label, color = site)) +
  geom_text(hjust = 1, size = 3, fontface = "bold") +
  scale_color_manual(values = SITE_COLORS, guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(plot.margin = margin(0, 2, 0, 2))

# ============================================================
# Combine + save
# ============================================================
message("=== Combining panels ===")
max_n       <- max(bin_table$n_proteins)
n_rows_show <- min(LABEL_TRUNCATE_N, max_n)
top_height  <- 2.5 + n_rows_show * ROW_HEIGHT_IN
fig_height  <- top_height + 3.2
fig_width   <- max(9, length(bin_levels) * 1.3 + 1.6)

# 2x2 layout: spacer(A) above row-labels(C), top-box panel(B) above matrix(D).
# B and D share column 2 -> patchwork forces them to the same panel width, which
# is what actually guarantees the bins line up; A and C share column 1 so the
# row-label column width matches the spacer exactly.
design <- "
AB
CD
"
combined <- plot_spacer() + p_top + p_row_labels + p_matrix +
  plot_layout(design = design, widths = c(1.3, 6), heights = c(top_height, 3))

cairo_pdf(file.path(fig_dir, "site_biomarker_upset_labeled.pdf"), width = fig_width, height = fig_height)
print(combined); dev.off()
png(file.path(fig_dir, "site_biomarker_upset_labeled.png"), width = fig_width, height = fig_height,
   units = "in", res = 300, bg = "white")
print(combined); dev.off()

message("\n18_upset_plot_with_protein_labels.R complete")
