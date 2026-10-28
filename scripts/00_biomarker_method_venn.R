library(VennDiagram)
library(dplyr)
library(readr)

# ---- EDIT FOR YOUR MACHINE --------------------------------------------------
PROJECT_DIR <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA"
INPUT_CSV   <- file.path(PROJECT_DIR, "data_tables", "07_analysis6_machine_learning", "biomarker_selection_full.csv")
OUTPUT_PDF  <- file.path(PROJECT_DIR, "figures", "07_analysis6_machine_learning", "biomarker_method_venn.pdf")
# -----------------------------------------------------------------------------

dir.create(dirname(OUTPUT_PDF), recursive = TRUE, showWarnings = FALSE)

df <- read_csv(INPUT_CSV, show_col_types = FALSE)
df$method[df$method %in% c("STABL", "STABL (final model)")] <- "STABL"

sets <- lapply(split(df$Protein_Label, df$method), unique)
sets <- sets[c("LASSO", "STABL", "Elastic Net", "sPLS-DA", "Boruta")]

venn_grob <- venn.diagram(
  x = sets,
  filename = NULL,
  fill = c("#dbe9f6", "#d7eedb", "#fde7c7", "#e2dcf4", "#f8d9de"),
  col = "grey30",
  cat.col = "grey20",
  cat.cex = 1.1,
  cex = 1,
  margin = 0.08,
  main = "Unique proteins selected per ML method (pooled across comparisons)",
  main.cex = 1.1
)

cairo_pdf(OUTPUT_PDF, width = 10, height = 10)
grid::grid.newpage()
grid::grid.draw(venn_grob)
dev.off()

message("Venn diagram written to: ", OUTPUT_PDF)
