library(readr)
library(readxl)
library(writexl)

# Load the CSV with the Patient IDs to keep
csv_data <- read_csv("C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/data/2_1_clinical_cooper_results_samples_used.csv")

# Extract the Patient IDs
patient_ids <- csv_data$Patient_ID

# Load the Excel file
excel_data <- read_excel("C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/data/Clinical_Akila_MERGED_20260714.xlsx")

# Subset: keep only rows where Patient_ID matches those in the CSV
subset_data <- excel_data[excel_data$Patient_ID %in% patient_ids, ]

# Check how many patients were retained
cat("Patients in CSV:", length(unique(patient_ids)), "\n")
cat("Rows in original Excel:", nrow(excel_data), "\n")
cat("Rows after subsetting:", nrow(subset_data), "\n")

# Save the result
write_xlsx(subset_data, "C:/Users/uqamudiy/OneDrive - The University of Queensland/Documents/Nomic_Lung_CA/data/Clinical_Akila_MERGED_subset.xlsx")

