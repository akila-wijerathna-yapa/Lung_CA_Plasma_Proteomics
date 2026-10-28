suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(tidyr)
})

CLIN_PATH <- r"(C:\Users\uqamudiy\OneDrive - The University of Queensland\Documents\Nomic_Lung_CA\data\Clinical_Akila_MERGED_subset.xlsx)"

# ---- helpers (verbatim from 01_run_nELISA_clinical_analysis.R) ----
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

raw <- as_tibble(read_excel(CLIN_PATH, sheet = 1))
message("Rows loaded: ", nrow(raw))

group_map <- c(
  "Chemotherapy" = "CTx", "Surgery" = "Sx", "Immunotherapy" = "IO",
  "Chemo + IO" = "CTx + IO", "Chemo + RT" = "CTx + RTx", "TKI" = "TKI",
  "Chemo + IO + RT" = "Other", "RT only" = "Other",
  "RT + IO (sequential)" = "Other", "IO + RT" = "Other"
)
group_order <- c("Sx", "CTx", "IO", "CTx + IO", "CTx + RTx", "TKI", "Other")

df <- raw %>%
  mutate(
    new_group = factor(unname(group_map[`1L_Treatment_Class`]), levels = group_order),
    dob = parse_date_mixed(DOB),
    dx_date = parse_date_mixed(Date_Diagnosis),
    age_dx = floor(as.numeric(dx_date - dob) / 365.25),
    stage_bin = case_when(
      Stage_Group %in% c("I", "II") ~ "I-II",
      Stage_Group %in% c("III", "IV") ~ "III-IV",
      TRUE ~ NA_character_
    ),
    vital = case_when(
      Alive_1_Dead_0 == 1 ~ "Alive",
      Alive_1_Dead_0 == 0 ~ "Dead",
      TRUE ~ NA_character_
    ),
    pdl1_cat = case_when(
      is.na(PDL1_TPS_pct) ~ "Unknown",
      PDL1_TPS_pct < 1 ~ "<1%",
      PDL1_TPS_pct < 50 ~ "1-49%",
      TRUE ~ ">=50%"
    )
  )

message("N per group:")
print(table(df$new_group, useNA = "ifany"))

# ---- summary builder: n (%) for categorical, median (IQR) for continuous ----
n_pct <- function(x, level, denom) {
  n <- sum(x == level, na.rm = TRUE)
  sprintf("%d (%.0f%%)", n, 100 * n / denom)
}

med_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.0f (%.0f\u2013%.0f)", q[2], q[1], q[3])
}

summarise_group <- function(sub) {
  n <- nrow(sub)
  tibble(
    n = n,
    age = med_iqr(sub$age_dx),
    sex_F = n_pct(sub$Sex, "F", n),
    sex_M = n_pct(sub$Sex, "M", n),
    stage_early = n_pct(sub$stage_bin, "I-II", n),
    stage_late = n_pct(sub$stage_bin, "III-IV", n),
    hist_adeno = n_pct(sub$Histological_Subtype, "Adenocarcinoma", n),
    hist_squam = n_pct(sub$Histological_Subtype, "Squamous Cell Carcinoma", n),
    pdl1_low = n_pct(sub$pdl1_cat, "<1%", n),
    pdl1_mid = n_pct(sub$pdl1_cat, "1-49%", n),
    pdl1_high = n_pct(sub$pdl1_cat, ">=50%", n),
    pdl1_na = n_pct(sub$pdl1_cat, "Unknown", n),
    resp_CR = n_pct(sub$`1L_Response_Std`, "CR", n),
    resp_PR = n_pct(sub$`1L_Response_Std`, "PR", n),
    resp_SD = n_pct(sub$`1L_Response_Std`, "SD", n),
    resp_PD = n_pct(sub$`1L_Response_Std`, "PD", n),
    responder = n_pct(sub$Response_Binary, "Responder", n),
    nonresponder = n_pct(sub$Response_Binary, "Non-Responder", n),
    vital_alive = n_pct(sub$vital, "Alive", n),
    vital_dead = n_pct(sub$vital, "Dead", n),
    pfs = med_iqr(sub$PFS_months),
    os = med_iqr(sub$OS_months)
  )
}

table1 <- bind_rows(
  lapply(group_order, function(g) summarise_group(df %>% filter(new_group == g)) %>% mutate(group = g)),
  summarise_group(df) %>% mutate(group = "Total")
) %>%
  select(group, everything()) %>%
  mutate(group = factor(group, levels = c(group_order, "Total"))) %>%
  arrange(group)

print(as.data.frame(table1), row.names = FALSE)

# transpose to characteristic-by-column layout matching the docx table
table1_t <- table1 %>%
  mutate(n = as.character(n)) %>%
  pivot_longer(-group, names_to = "characteristic", values_to = "value") %>%
  pivot_wider(names_from = group, values_from = value)

message("\nTable 1 (characteristic x group):")
print(as.data.frame(table1_t), row.names = FALSE)

# ---- save ----
OUT_STUB <- r"(C:\Users\uqamudiy\OneDrive - The University of Queensland\Documents\Nomic_Lung_CA\tables\table1_characteristics)"
dir.create(dirname(OUT_STUB), recursive = TRUE, showWarnings = FALSE)

write.csv(table1_t, paste0(OUT_STUB, ".csv"), row.names = FALSE)

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
openxlsx::write.xlsx(table1_t, paste0(OUT_STUB, ".xlsx"))

message("Saved: ", OUT_STUB, ".csv / .xlsx")

