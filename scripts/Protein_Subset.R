library(readxl)

dat <- read_excel(
  "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716/data/Protein_pg_ml_Akila_202603271608.xlsx",
  col_names = FALSE
)

dim(dat)

dat[1:6, 1:15]

ids <- read_excel(
  "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716/data/Protein_IDs.xlsx"
)

names(ids)

head(ids)

#############

library(readxl)
library(openxlsx)

#-----------------------------
# File paths
#-----------------------------
protein_file <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716/data/Protein_pg_ml_Akila_202603271608.xlsx"

id_file <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716/data/Protein_IDs.xlsx"

output_file <- "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716/data/Protein_pg_ml_Akila_202603271608_subset.xlsx"

#-----------------------------
# Read data
#-----------------------------
dat <- read_excel(protein_file, col_names = FALSE)

ids <- read_excel(id_file)

#-----------------------------
# Create lookup IDs
#-----------------------------
ids$key <- paste(ids$`Protein Name:`,
                 ids$`Uniprot ID:`,
                 sep = "|")

#-----------------------------
# Protein columns begin at K
#-----------------------------
protein_cols <- 11:ncol(dat)

protein_names <- as.character(unlist(dat[3, protein_cols]))
uniprot_ids  <- as.character(unlist(dat[2, protein_cols]))

protein_key <- paste(protein_names,
                     uniprot_ids,
                     sep = "|")

# Columns to retain
keep_proteins <- protein_cols[protein_key %in% ids$key]

#-----------------------------
# Keep metadata (A:J)
#-----------------------------
metadata_cols <- 1:10

subset_dat <- dat[, c(metadata_cols, keep_proteins)]

#-----------------------------
# Save
#-----------------------------
write.xlsx(subset_dat,
           output_file,
           overwrite = TRUE)

cat("Number of proteins retained:", length(keep_proteins), "\n")
cat("Total columns in output:", ncol(subset_dat), "\n")



missing <- ids[!(ids$key %in% protein_key), ]

missing

write.xlsx(
  missing,
  "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/nELISA_20260716/data/Missing_Proteins.xlsx",
  overwrite = TRUE
)
