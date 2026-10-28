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

# ── Local utility functions ───────────────────────────────────
# All functions are defined here; none imported from other modules.

clean_names_simple <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | y == ""] <- paste0("col_", seq_along(y))[is.na(y) | y == ""]
  y <- tolower(gsub("[^A-Za-z0-9]+", "_", y))
  y <- gsub("^_+|_+$", "", y)
  make.unique(y, sep = "_dup")
}

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
    if (length(v) == 0 || is.na(v) ||
        trimws(as.character(v)) %in% c("", "NA", "NaN", "Not applicable"))
      return(as.Date(NA))
    if (inherits(v, "Date"))   return(as.Date(v))
    if (inherits(v, "POSIXt")) return(as.Date(v))
    if (is.numeric(v)) {
      if (v > 59 && v < 80000) return(as.Date(floor(v), origin = "1899-12-30"))
      return(as.Date(NA))
    }
    txt    <- trimws(as.character(v))
    op     <- options(warn = -1)
    parsed <- lubridate::parse_date_time(
      txt,
      orders = c("dmy", "ymd", "mdy", "d-b-y", "d-b-Y", "d/m/Y", "m/d/Y", "Y-m-d")
    )
    options(op)
    if (length(parsed) == 0 || is.na(parsed[1])) return(as.Date(NA))
    as.Date(parsed[1])
  }
  as.Date(vapply(x, parse_one, as.Date(NA)))
}

months_between <- function(start_date, end_date) {
  out <- as.numeric(difftime(end_date, start_date, units = "days")) / 30.4375
  out[!is.finite(out) | out < 0] <- NA_real_
  out
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

annotate_features <- function(analyte_meta) {
  keys <- sub("-.*$", "", analyte_meta$Uniprot_ID)
  op   <- options(warn = -1)
  gene_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys, keytype = "UNIPROT",
    column = "SYMBOL",   multiVals = "first"
  )
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys, keytype = "UNIPROT",
    column = "ENTREZID", multiVals = "first"
  )
  options(op)
  analyte_meta %>%
    mutate(
      Uniprot_Clean = keys,
      Gene_Symbol   = unname(gene_symbol[match(keys, names(gene_symbol))]),
      Entrez_ID     = unname(entrez[match(keys, names(entrez))]),
      Protein_Label = dplyr::coalesce(Gene_Symbol, Protein_Name, nELISA_ID)
    )
}

read_proteomics_workbook <- function(path) {
  raw             <- readxl::read_excel(path, sheet = 1, col_names = FALSE)
  first_row       <- as.character(unlist(raw[1, ]))
  first_analyte_col <- which(grepl("^nEL\\d+", first_row))[1]
  if (is.na(first_analyte_col))
    stop("Could not locate analyte columns in proteomics workbook.")

  meta_headers <- clean_names_simple(unlist(raw[5, seq_len(first_analyte_col - 1)]))
  sample_info  <- raw[6:nrow(raw), seq_len(first_analyte_col - 1), drop = FALSE]
  colnames(sample_info) <- meta_headers

  pick_col <- function(nms, exact = character(), regex = character()) {
    for (nm in exact)  { if (nm %in% nms) return(nm) }
    for (rx in regex)  { hit <- grep(rx, nms, value = TRUE); if (length(hit)) return(hit[1]) }
    NA_character_
  }

  analyte_ids  <- as.character(unlist(raw[1, first_analyte_col:ncol(raw)]))
  op <- options(warn = -1)
  ulod_vals <- as.numeric(unlist(raw[4, first_analyte_col:ncol(raw)]))
  llod_vals <- as.numeric(unlist(raw[5, first_analyte_col:ncol(raw)]))
  options(op)
  analyte_meta <- tibble(
    nELISA_ID    = analyte_ids,
    Uniprot_ID   = as.character(unlist(raw[2, first_analyte_col:ncol(raw)])),
    Protein_Name = as.character(unlist(raw[3, first_analyte_col:ncol(raw)])),
    ULOD_pg_ml   = ulod_vals,
    LLOD_pg_ml   = llod_vals
  )

  expr_raw <- raw[6:nrow(raw), first_analyte_col:ncol(raw), drop = FALSE]
  colnames(expr_raw) <- analyte_ids

  sample_cols   <- names(sample_info)
  sample_id_col <- pick_col(sample_cols, exact = c("sample_id_", "sample_id"), regex = c("^sample_id"))
  pre_post_col  <- pick_col(sample_cols, exact = c("pre_post"),    regex = c("^pre.*post$"))
  comments_col  <- pick_col(sample_cols, exact = c("sample_comments"), regex = c("comment"))
  qc_col        <- pick_col(sample_cols, exact = c("sample_passed_qc"), regex = c("passed.*qc"))
  id_col        <- pick_col(sample_cols, exact = c("id_number"),   regex = c("^id"))
  plate_col     <- pick_col(sample_cols, exact = c("plate_barcode"), regex = c("plate"))
  well_col      <- pick_col(sample_cols, exact = c("well_id"),     regex = c("well"))

  if (is.na(sample_id_col) || is.na(pre_post_col) || is.na(id_col))
    stop("Could not detect required sample metadata columns in proteomics workbook.")

  sample_info <- sample_info %>%
    mutate(
      sample_id       = as.character(.data[[sample_id_col]]),
      pre_post        = tolower(trimws(as.character(.data[[pre_post_col]]))),
      sample_comments = if (!is.na(comments_col)) as.character(.data[[comments_col]]) else NA_character_,
      sample_passed_qc = if (!is.na(qc_col))      as.character(.data[[qc_col]])       else NA_character_,
      id_number       = norm_id(.data[[id_col]]),
      plate_barcode   = if (!is.na(plate_col))    as.character(.data[[plate_col]])    else NA_character_,
      well_id         = if (!is.na(well_col))     as.character(.data[[well_col]])     else NA_character_
    )

  expr_num <- matrix(
    NA_real_, nrow = nrow(expr_raw), ncol = ncol(expr_raw),
    dimnames = list(sample_info$sample_id, analyte_ids)
  )
  for (j in seq_along(analyte_ids)) {
    raw_chr   <- trimws(as.character(expr_raw[[j]]))
    lower_chr <- tolower(raw_chr)
    num_ready     <- gsub(",", "", raw_chr)
    is_num        <- grepl("^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$", num_ready)
    vals          <- rep(NA_real_, length(raw_chr))
    vals[is_num]  <- as.numeric(num_ready[is_num])
    vals[grepl("^<",  raw_chr)]                    <- analyte_meta$LLOD_pg_ml[j] / 2
    vals[grepl("^>",  raw_chr)]                    <- analyte_meta$ULOD_pg_ml[j] * 2
    vals[grepl("low\\s*bead\\s*count", lower_chr)] <- NA_real_
    vals[lower_chr %in% c("", "na", "nan", "null")] <- NA_real_
    vals[!is.finite(vals) | vals <= 0]             <- NA_real_
    expr_num[, j] <- vals
  }

  log2_mat <- log2(expr_num)
  log2_mat[!is.finite(log2_mat)] <- NA_real_

  list(
    sample_info  = sample_info,
    analyte_meta = annotate_features(analyte_meta),
    expr_pgml    = expr_num,
    expr_log2    = log2_mat
  )
}

choose_representatives <- function(sample_df, expr_log2) {
  sample_df <- sample_df %>%
    mutate(
      source_row        = seq_len(nrow(sample_df)),
      sample_id         = as.character(sample_id),
      missing_rate      = rowMeans(is.na(expr_log2[sample_id, , drop = FALSE])),
      observed_features = rowSums(!is.na(expr_log2[sample_id, , drop = FALSE]))
    )

  split_idx  <- split(seq_len(nrow(sample_df)), sample_df$sample_id)
  chosen_idx <- unlist(lapply(split_idx, function(idx) {
    ord <- order(
      sample_df$missing_rate[idx],
      -sample_df$observed_features[idx],
      sample_df$plate_barcode[idx],
      idx
    )
    idx[ord[1]]
  }))

  duplicate_ids <- names(split_idx)[lengths(split_idx) > 1]
  dup_audit <- bind_rows(lapply(duplicate_ids, function(id) {
    idx         <- split_idx[[id]]
    chosen_i    <- intersect(idx, chosen_idx)[1]
    discarded_i <- setdiff(idx, chosen_i)
    sub_mat     <- expr_log2[sample_df$sample_id[idx], , drop = FALSE]
    corr_val    <- if (nrow(sub_mat) >= 2) {
      op <- options(warn = -1)
      cv <- cor(as.numeric(sub_mat[1, ]), as.numeric(sub_mat[2, ]),
                use = "pairwise.complete.obs")
      options(op); cv
    } else NA_real_
    chosen_plate    <- trimws(as.character(sample_df$plate_barcode[chosen_i]))
    discarded_plates <- trimws(as.character(sample_df$plate_barcode[discarded_i]))
    tibble(
      sample_id                          = id,
      n_rows                             = length(idx),
      rows                               = paste(sample_df$source_row[idx], collapse = ";"),
      plates                             = paste(sample_df$plate_barcode[idx], collapse = ";"),
      chosen_source_row                  = sample_df$source_row[chosen_i],
      chosen_plate_barcode               = sample_df$plate_barcode[chosen_i],
      discarded_source_rows              = paste(sample_df$source_row[discarded_i], collapse = ";"),
      discarded_plate_barcodes           = paste(sample_df$plate_barcode[discarded_i], collapse = ";"),
      chosen_missing_rate                = sample_df$missing_rate[chosen_i],
      chosen_observed_features           = sample_df$observed_features[chosen_i],
      cross_plate_duplicate              = length(unique(trimws(as.character(sample_df$plate_barcode[idx])))) > 1,
      chosen_plate_differs_from_discarded = any(discarded_plates != chosen_plate, na.rm = TRUE),
      pearson_correlation                = corr_val
    )
  }))

  if (nrow(dup_audit) && any(dup_audit$chosen_plate_differs_from_discarded, na.rm = TRUE)) {
    message(
      "Duplicate audit: representative from a different plate than discarded for ",
      sum(dup_audit$chosen_plate_differs_from_discarded, na.rm = TRUE), " sample_id(s)."
    )
  }

  representative <- sample_df[chosen_idx, , drop = FALSE] %>%
    mutate(representative_flag = TRUE) %>%
    arrange(pre_post, id_number, sample_id)

  choice_table <- sample_df %>%
    mutate(
      representative_flag   = seq_len(nrow(sample_df)) %in% chosen_idx,
      duplicate_sample_flag = duplicated(sample_id) | duplicated(sample_id, fromLast = TRUE)
    ) %>%
    left_join(
      dup_audit %>% select(sample_id, cross_plate_duplicate, chosen_plate_differs_from_discarded),
      by = "sample_id"
    ) %>%
    arrange(sample_id, desc(representative_flag), missing_rate, plate_barcode)

  list(representative = representative, duplicate_audit = dup_audit, choice_table = choice_table)
}

read_sample_map <- function(path) {
  raw      <- readxl::read_excel(path, sheet = 1, col_names = FALSE)
  pre_ids  <- norm_id(raw[4:nrow(raw), 1, drop = TRUE])
  post_ids <- norm_id(raw[4:nrow(raw), 4, drop = TRUE])
  list(
    pre  = pre_ids[ !is.na(pre_ids)  & pre_ids  != ""],
    post = post_ids[!is.na(post_ids) & post_ids != ""]
  )
}

build_clinical <- function(path) {
  raw <- as_tibble(readxl::read_excel(path, sheet = 1))

  group_map <- c(
    "Chemotherapy"           = "CTx",
    "Surgery"                = "Sx",
    "Immunotherapy"          = "IO",
    "Chemo + IO"             = "CTx + IO",
    "Chemo + RT"             = "CTx + RTx",
    "TKI"                    = "TKI",
    "Chemo + IO + RT"        = "Other",
    "RT only"                = "Other",
    "RT + IO (sequential)"   = "Other",
    "IO + RT"                = "Other"
  )

  clin <- raw %>%
    transmute(
      patient_id         = norm_id(Patient_ID),
      patient_id_raw     = as.character(Patient_ID),
      sex = case_when(
        toupper(trimws(as.character(Sex))) %in% c("M", "MALE")   ~ "MALE",
        toupper(trimws(as.character(Sex))) %in% c("F", "FEMALE") ~ "FEMALE",
        TRUE ~ NA_character_
      ),
      dob                = parse_date_mixed(DOB),
      stage_group        = as.character(Stage_Group),
      stage_collapsed    = case_when(
        str_detect(toupper(as.character(Stage_Group)), "III|IV") ~ "III_IV",
        str_detect(toupper(as.character(Stage_Group)), "I|II")   ~ "I_II",
        TRUE ~ NA_character_
      ),
      histology          = as.character(Histological_Subtype),
      pdl1_tps = {
        txt <- trimws(as.character(PDL1_TPS_pct))
        op  <- options(warn = -1); val <- as.numeric(txt); options(op)
        val[tolower(txt) %in% c("<1", "< 1")] <- 0.5
        val
      },
      pdl1_category = case_when(
        is.na(pdl1_tps)   ~ NA_character_,
        pdl1_tps < 1      ~ "Negative",
        pdl1_tps < 50     ~ "Low",
        TRUE              ~ "High"
      ),
      diagnosis_date     = parse_date_mixed(Date_Diagnosis),
      treatment_class    = as.character(`1L_Treatment_Class`),
      new_group          = unname(group_map[as.character(`1L_Treatment_Class`)]),
      start_date         = parse_date_mixed(`1L_Start_Date`),
      end_date           = parse_date_mixed(`1L_End_Date`),
      response_4cat = case_when(
        trimws(as.character(`1L_Response_Std`)) %in% c("CR", "PR", "SD", "PD") ~
          trimws(as.character(`1L_Response_Std`)),
        TRUE ~ NA_character_
      ),
      response_binary = case_when(
        trimws(as.character(Response_Binary)) == "Responder"     ~ "PR/CR",
        trimws(as.character(Response_Binary)) == "Non-Responder" ~ "SD/PD",
        TRUE ~ NA_character_
      ),
      recurrence_date    = parse_date_mixed(`1L_Date_Recurrence`),
      recurrence_site    = as.character(`1L_Sites_Recurrence`),
      pfs_months_raw     = { op <- options(warn=-1); v <- as.numeric(`PFS_months`); options(op); v },
      pfs_censoring_type = as.character(`PFS_Censoring_Type`),
      os_months          = { op <- options(warn=-1); v <- as.numeric(OS_months); options(op); v },
      alive_code         = { op <- options(warn=-1); v <- as.integer(`Alive_1_Dead_0`); options(op); v },
      os_event           = { op <- options(warn=-1); v <- as.integer(`OS_Event_1Dead`); options(op); v },
      censored_date      = parse_date_mixed(Censored_Date)
    ) %>%
    mutate(
      age_at_diagnosis = ifelse(
        !is.na(dob) & !is.na(diagnosis_date),
        floor(as.numeric(diagnosis_date - dob) / 365.25),
        NA_real_
      ),
      treatment_duration_months = months_between(start_date, end_date),
      pfs_censor_class = case_when(
        str_detect(pfs_censoring_type, "^Exact")             ~ "exact",
        str_detect(pfs_censoring_type, "^Right-censored")    ~ "right_censored",
        str_detect(pfs_censoring_type, "^Interval-censored") ~ "interval_censored",
        str_detect(pfs_censoring_type, "^Not applicable")    ~ "not_applicable",
        TRUE ~ "unknown"
      )
    )

  # PFS time: use spreadsheet value, fall back to date arithmetic
  derived_pfs_time <- clin$pfs_months_raw
  is_na_pfs <- is.na(derived_pfs_time)
  derived_pfs_time[is_na_pfs & !is.na(clin$diagnosis_date) & !is.na(clin$recurrence_date)] <-
    months_between(
      clin$diagnosis_date[is_na_pfs & !is.na(clin$diagnosis_date) & !is.na(clin$recurrence_date)],
      clin$recurrence_date[is_na_pfs & !is.na(clin$diagnosis_date) & !is.na(clin$recurrence_date)]
    )
  is_na_pfs <- is.na(derived_pfs_time)
  derived_pfs_time[is_na_pfs & !is.na(clin$diagnosis_date) & !is.na(clin$censored_date)] <-
    months_between(
      clin$diagnosis_date[is_na_pfs & !is.na(clin$diagnosis_date) & !is.na(clin$censored_date)],
      clin$censored_date[is_na_pfs & !is.na(clin$diagnosis_date) & !is.na(clin$censored_date)]
    )

  pfs_event_main <- case_when(
    clin$pfs_censor_class == "exact"         ~ 1L,
    clin$pfs_censor_class == "right_censored" ~ 0L,
    TRUE ~ NA_integer_
  )

  pfs_event_sensitivity <- case_when(
    clin$pfs_censor_class == "exact"             ~ 1L,
    clin$pfs_censor_class == "right_censored"    ~ 0L,
    clin$pfs_censor_class == "interval_censored" ~ 1L,
    TRUE ~ NA_integer_
  )

  pfs_interval_lower <- case_when(
    clin$pfs_censor_class == "exact"             ~ derived_pfs_time,
    clin$pfs_censor_class == "right_censored"    ~ derived_pfs_time,
    clin$pfs_censor_class == "interval_censored" ~ dplyr::coalesce(
      months_between(clin$diagnosis_date, clin$end_date),
      months_between(clin$diagnosis_date, clin$start_date),
      0
    ),
    TRUE ~ NA_real_
  )

  pfs_interval_upper <- case_when(
    clin$pfs_censor_class == "exact"             ~ derived_pfs_time,
    clin$pfs_censor_class == "right_censored"    ~ Inf,
    clin$pfs_censor_class == "interval_censored" ~ dplyr::coalesce(
      months_between(clin$diagnosis_date, clin$censored_date),
      clin$os_months,
      derived_pfs_time
    ),
    TRUE ~ NA_real_
  )
  needs_epsilon <- is.finite(pfs_interval_upper) & is.finite(pfs_interval_lower) &
    pfs_interval_upper <= pfs_interval_lower
  pfs_interval_upper[needs_epsilon] <- pfs_interval_lower[needs_epsilon] + 0.01

  followup_months <- dplyr::coalesce(
    derived_pfs_time,
    clin$os_months,
    months_between(clin$diagnosis_date, clin$censored_date)
  )

  pfs_months_sensitivity <- case_when(
    clin$pfs_censor_class == "interval_censored" ~ pfs_interval_upper,
    TRUE ~ derived_pfs_time
  )

  efs_event <- ifelse(
    !is.na(pfs_event_sensitivity) & pfs_event_sensitivity == 1L &
      !is.na(pfs_months_sensitivity) &
      (is.na(clin$os_months) | pfs_months_sensitivity <= clin$os_months),
    1L,
    ifelse(!is.na(clin$os_event) & clin$os_event == 1L & !is.na(clin$os_months), 1L,
      ifelse(!is.na(followup_months), 0L, NA_integer_)
    )
  )

  efs_months <- ifelse(
    efs_event == 1L & !is.na(pfs_months_sensitivity) &
      (is.na(clin$os_months) | pfs_months_sensitivity <= clin$os_months),
    pfs_months_sensitivity,
    ifelse(efs_event == 1L & !is.na(clin$os_months), clin$os_months, followup_months)
  )

  clin %>%
    mutate(
      pfs_months_main          = derived_pfs_time,
      pfs_event_main           = pfs_event_main,
      pfs_interval_lower       = pfs_interval_lower,
      pfs_interval_upper       = pfs_interval_upper,
      pfs_months_sensitivity   = pfs_months_sensitivity,
      pfs_event_sensitivity    = pfs_event_sensitivity,
      efs_months               = efs_months,
      efs_event                = efs_event,
      recurrence_binary_strict      = ifelse(!is.na(recurrence_date), "Recurrence", "No recurrence"),
      recurrence_binary_sensitivity = ifelse(
        !is.na(recurrence_date) | pfs_censor_class == "interval_censored",
        "Recurrence", "No recurrence"
      ),
      recurrence_basis = case_when(
        !is.na(recurrence_date)              ~ "dated_recurrence",
        pfs_censor_class == "interval_censored" ~ "interval_censored_progression",
        TRUE ~ "no_documented_recurrence"
      ),
      recurrence_binary = recurrence_binary_sensitivity,
      alive_status = case_when(
        alive_code == 1 ~ "Alive",
        alive_code == 0 ~ "Dead",
        TRUE ~ NA_character_
      )
    )
}

# ── Load inputs ───────────────────────────────────────────────
message("\n=== 01 | Load inputs ===")
clinical   <- build_clinical(paths$clinical)
sample_map <- read_sample_map(paths$sample_map)
proteomics <- read_proteomics_workbook(paths$proteomics)

sample_info_raw <- proteomics$sample_info
expr_pgml_raw   <- proteomics$expr_pgml
expr_log2_raw   <- proteomics$expr_log2
feature_meta    <- proteomics$analyte_meta

message("Clinical patients: ", nrow(clinical))
message("Proteomics rows (raw): ", nrow(sample_info_raw))
message("Proteins: ", ncol(expr_pgml_raw))

# ── QC-fail filter ────────────────────────────────────────────
message("\n=== 01 | QC-fail filter ===")
qc_fail_mask <- str_detect(
  tolower(dplyr::coalesce(sample_info_raw$sample_comments, "")),
  "qc\\s*_?fail"
)
sample_info_qc <- sample_info_raw[!qc_fail_mask, , drop = FALSE]
expr_pgml_qc   <- expr_pgml_raw[!qc_fail_mask, , drop = FALSE]
expr_log2_qc   <- expr_log2_raw[!qc_fail_mask, , drop = FALSE]
rownames(expr_pgml_qc) <- sample_info_qc$sample_id
rownames(expr_log2_qc) <- sample_info_qc$sample_id
message("Rows removed (qc_fail): ", sum(qc_fail_mask))
message("Rows after filter: ", nrow(sample_info_qc))

# ── Duplicate resolution ──────────────────────────────────────
message("\n=== 01 | Duplicate resolution ===")
rep_selection  <- choose_representatives(sample_info_qc, expr_log2_qc)
sample_rep     <- rep_selection$representative
expr_pgml_rep  <- expr_pgml_qc[sample_rep$sample_id, , drop = FALSE]
expr_log2_rep  <- expr_log2_qc[sample_rep$sample_id, , drop = FALSE]
rownames(expr_pgml_rep) <- sample_rep$sample_id
rownames(expr_log2_rep) <- sample_rep$sample_id
message("Representative samples: ", nrow(sample_rep))

# ── Split pre / post ──────────────────────────────────────────
pre_all  <- sample_rep %>% filter(pre_post == "pre")
post_all <- sample_rep %>% filter(pre_post == "post")
message("PRE unique: ", nrow(pre_all), " | POST unique: ", nrow(post_all))

# ── Clinical matching ─────────────────────────────────────────
message("\n=== 01 | Clinical matching ===")
pre_matched  <- pre_all  %>% mutate(patient_id = id_number) %>% inner_join(clinical, by = "patient_id")
post_matched <- post_all %>% mutate(patient_id = id_number) %>% inner_join(clinical, by = "patient_id")

paired_ids      <- intersect(pre_matched$patient_id, post_matched$patient_id)
paired_clinical <- clinical %>% filter(patient_id %in% paired_ids)

pre_rep_ids  <- sort(unique(pre_all$id_number))
post_rep_ids <- sort(unique(post_all$id_number))

message("PRE matched: ",  n_distinct(pre_matched$patient_id))
message("POST matched: ", n_distinct(post_matched$patient_id))
message("Paired: ",       length(paired_ids))

# ── Cohort summary ────────────────────────────────────────────
cohort_summary <- tibble(
  metric = c(
    "Clinical patients",
    "Proteomics rows (raw)",
    "Proteomics rows after qc_fail filter",
    "Baseline PRE rows after qc_fail filter",
    "POST rows after qc_fail filter",
    "Baseline PRE unique samples after duplicate resolution",
    "POST unique samples after duplicate resolution",
    "Baseline PRE clinically matched patients",
    "POST clinically matched patients",
    "Matched paired patients",
    "Matched surgery pairs"
  ),
  value = c(
    nrow(clinical),
    nrow(sample_info_raw),
    nrow(sample_info_qc),
    sum(sample_info_qc$pre_post == "pre"),
    sum(sample_info_qc$pre_post == "post"),
    nrow(pre_all),
    nrow(post_all),
    n_distinct(pre_matched$patient_id),
    n_distinct(post_matched$patient_id),
    length(paired_ids),
    sum(paired_clinical$new_group == "Sx", na.rm = TRUE)
  )
)

# ── Audit mismatches ──────────────────────────────────────────
clinical_ids <- sort(unique(clinical$patient_id))
audit_mismatch <- list(
  clinical_pre_missing_from_proteomics      = tibble(patient_id = setdiff(clinical_ids, pre_rep_ids)),
  proteomics_pre_not_in_clinical            = tibble(patient_id = setdiff(pre_rep_ids, clinical_ids)),
  sample_map_pre_missing_from_proteomics    = tibble(patient_id = setdiff(sample_map$pre,  pre_rep_ids)),
  sample_map_post_missing_from_proteomics   = tibble(patient_id = setdiff(sample_map$post, post_rep_ids))
)

# ── Write audit CSVs ──────────────────────────────────────────
message("\n=== 01 | Writing audit tables ===")
safe_write_csv(clinical,      file.path(paths$tables, "00_cohort_audit", "clinical_cleaned.csv"))
safe_write_csv(feature_meta,  file.path(paths$tables, "00_cohort_audit", "proteomics_feature_annotation.csv"))
safe_write_csv(cohort_summary, file.path(paths$tables, "00_cohort_audit", "cohort_summary.csv"))
safe_write_csv(rep_selection$choice_table,   file.path(paths$tables, "00_cohort_audit", "duplicate_sample_choice_table.csv"))
safe_write_csv(rep_selection$duplicate_audit, file.path(paths$tables, "00_cohort_audit", "duplicate_sample_audit.csv"))
safe_write_csv(sample_info_raw[qc_fail_mask, , drop = FALSE],
               file.path(paths$tables, "00_cohort_audit", "qc_fail_samples_removed.csv"))
invisible(lapply(names(audit_mismatch), function(nm) {
  safe_write_csv(audit_mismatch[[nm]],
                 file.path(paths$tables, "00_cohort_audit", paste0(nm, ".csv")))
}))

writeLines(c(
  "nELISA clinical analysis — key assumptions",
  "-------------------------------------------",
  paste("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "1. Nine original treatment classes mapped to 7 New_Group classes.",
  "2. 'IO + RT' merged into 'Other' (not in the instruction table).",
  "3. Special tokens: <LLOD -> LLOD/2; >ULOD -> ULOD*2; 'Low bead count' -> NA.",
  "4. No imputation for ordinary missing values.",
  "5. Technical duplicate PRE rows resolved by lowest missing-rate per sample_id.",
  "6. Clinical/proteomics mismatches saved as audit tables.",
  "7. Main PFS Cox uses exact and right-censored cases only.",
  "8. Machine learning applied to PR vs PD and PR/CR vs SD/PD.",
  "9. Surgery poster uses interval-expanded recurrence definition."
),
file.path(paths$tables, "00_cohort_audit", "assumptions_and_caveats.txt"))

# ── Save checkpoint ───────────────────────────────────────────
message("\n=== 01 | Saving checkpoint_01.RData ===")
save(
  clinical, feature_meta,
  sample_info_raw, sample_info_qc,
  expr_pgml_raw,  expr_log2_raw,
  expr_pgml_rep,  expr_log2_rep,
  sample_rep, rep_selection,
  pre_all,  post_all,
  pre_matched, post_matched,
  pre_rep_ids, post_rep_ids,
  paired_ids, paired_clinical,
  cohort_summary, sample_map,
  file = file.path(paths$rdata, "checkpoint_01.RData")
)

message("01_data_loading.R complete")
print(cohort_summary)
