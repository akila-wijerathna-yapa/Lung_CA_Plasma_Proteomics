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

# ── Load checkpoints ──────────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_04.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
message("Checkpoints 01, 03, 04, 07 loaded")
while (dev.cur() > 1) dev.off()

# ── Local constants ───────────────────────────────────────────
STABL_TBL_DIR <- file.path(paths$tables, "07_analysis6_machine_learning", "stabl")
OUT_FIG       <- file.path(paths$figures, "07_analysis6_machine_learning", "stabl", "km")
OUT_TBL       <- file.path(paths$tables,  "07_analysis6_machine_learning", "stabl", "km")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)

# Set TRUE to skip plots already saved (avoids Windows file-lock errors
# when PDFs are open in a viewer from a previous run).
SKIP_IF_DONE  <- FALSE

MIN_KM_N      <- 10
MIN_KM_EVENTS <- 5

# ── Local functions ────────────────────────────────────────────
sanitize_filename <- function(x) {
  # Remove characters illegal in Windows filenames: / \ : * ? " < > |
  out <- gsub('[/\\\\:*?"<>|]', "_", as.character(x))
  out <- gsub("\\s+", "_", out)
  out <- gsub("_+", "_", out)
  substr(out, 1, 60)
}

coerce_event_indicator <- function(x) {
  raw_chr <- trimws(as.character(x))
  op  <- options(warn = -1); out <- as.numeric(raw_chr); options(op)
  out[!is.finite(out) | !(out %in% c(0, 1))] <- NA_real_
  as.integer(out)
}

# Generic median-split KM plot — reused for all comparisons
plot_km_median_split <- function(sample_df, expr_vec, time_col, event_col,
                                 protein_label, subtitle_extra, path_stub,
                                 min_n = MIN_KM_N, min_events = MIN_KM_EVENTS) {
  op   <- options(warn = -1)
  time_val <- as.numeric(sample_df[[time_col]])
  options(op)
  event_val <- coerce_event_indicator(sample_df[[event_col]])
  keep <- is.finite(expr_vec) & !is.na(time_val) & !is.na(event_val)
  if (sum(keep) < min_n) {
    message("    SKIP KM (n<", min_n, "): ", protein_label); return(invisible(NULL))
  }
  df <- tibble::tibble(
    time_value  = time_val[keep],
    event_value = event_val[keep],
    expr_value  = expr_vec[keep]
  ) %>%
    dplyr::mutate(expr_group = factor(
      ifelse(expr_value > median(expr_value, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")))
  if (any(table(df$expr_group) < 4) || sum(df$event_value == 1) < min_events) {
    message("    SKIP KM (small stratum/events): ", protein_label); return(invisible(NULL))
  }
  fit <- survival::survfit(survival::Surv(time_value, event_value) ~ expr_group, data = df)
  diff_fit <- tryCatch(
    survival::survdiff(survival::Surv(time_value, event_value) ~ expr_group, data = df),
    error = function(e) NULL)
  logrank_p <- if (!is.null(diff_fit))
    stats::pchisq(diff_fit$chisq, df = max(length(diff_fit$n) - 1, 1), lower.tail = FALSE) else NA_real_

  km_plot <- survminer::ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    palette = c("#2166AC", "#B2182B"),
    legend.title = "Expression", legend.labs = c("Low", "High"),
    title    = paste0("KM \u2014 ", protein_label),
    subtitle = paste0(subtitle_extra, " | logrank p=",
                      scales::number(logrank_p, accuracy = 0.001)),
    xlab = "Months", ylab = "Survival probability",
    ggtheme = ggplot2::theme_bw(base_size = 11)
  )
  pdf_path <- paste0(path_stub, ".pdf")
  png_path <- paste0(path_stub, ".png")
  if (SKIP_IF_DONE && file.exists(pdf_path)) {
    message("    SKIP (already exists): ", basename(pdf_path))
  } else {
    tryCatch({
      cairo_pdf(pdf_path, width = 7.8, height = 6.5)
      print(km_plot); dev.off()
      png(png_path, width = 7.8, height = 6.5, units = "in", res = 300, bg = "white")
      print(km_plot); dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("    ERROR saving ", basename(pdf_path), ": ", conditionMessage(e),
              "\n    -> Close any open PDF viewers and re-run.")
    })
  }

  tibble::tibble(
    protein = protein_label, n = nrow(df), n_events = sum(df$event_value == 1),
    n_low = sum(df$expr_group == "Low"), n_high = sum(df$expr_group == "High"),
    logrank_p = logrank_p
  )
}

# ── Paired dot-line plot (for AZU1, not KM) ───────────────────
plot_paired_dotline <- function(pre_vec, post_vec, patient_ids, protein_label,
                                path_stub, p_value = NA_real_) {
  df <- tibble::tibble(
    patient_id = rep(patient_ids, 2),
    timepoint  = rep(c("PRE", "POST"), each = length(patient_ids)),
    value      = c(pre_vec, post_vec)
  ) %>% dplyr::mutate(timepoint = factor(timepoint, levels = c("PRE", "POST")))

  p <- ggplot(df, aes(timepoint, value)) +
    geom_boxplot(aes(fill = timepoint), alpha = 0.35, outlier.shape = NA,
                width = 0.5, color = "grey30") +
    geom_line(aes(group = patient_id), alpha = 0.45, color = "grey55", linewidth = 0.4) +
    geom_point(aes(color = timepoint), size = 2) +
    scale_fill_manual(values  = c(PRE = "#f4b6b2", POST = "#8dd3e7")) +
    scale_color_manual(values = c(PRE = "#d95f5f", POST = "#1f9fb5")) +
    labs(title    = paste0("Pre vs Post Surgery — ", protein_label),
         subtitle = paste0("Paired Wilcoxon p=", scales::number(p_value, accuracy = 0.0001)),
         x = NULL, y = "log2 concentration") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  pdf_path2 <- paste0(path_stub, ".pdf")
  if (SKIP_IF_DONE && file.exists(pdf_path2)) {
    message("    SKIP (already exists): ", basename(pdf_path2))
  } else {
    tryCatch({
      cairo_pdf(pdf_path2, width = 6, height = 5.5); print(p); dev.off()
      png(paste0(path_stub, ".png"), width = 6, height = 5.5,
          units = "in", res = 300, bg = "white"); print(p); dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("    ERROR saving ", basename(pdf_path2), ": ", conditionMessage(e))
    })
  }
}

# ── Helper: get nELISA_ID for a Gene_Symbol from STABL output ─
get_nelisa_id <- function(stabl_file) {
  if (!file.exists(stabl_file)) return(NULL)
  df <- readr::read_csv(stabl_file, show_col_types = FALSE)
  if (!nrow(df)) return(NULL)
  df$nELISA_ID[1]
}

# ─────────────────────────────────────────────────────────────
# PART 1: Cox-discovered proteins — KM on their own endpoint
# ─────────────────────────────────────────────────────────────
message("\n=== 08c | Part 1: Cox-discovered proteins ===")

cox_os_id  <- readr::read_csv(file.path(STABL_TBL_DIR, "stabl_14_cox_OS_stabl_selected_features.csv"),
                              show_col_types = FALSE)
cox_pfs_id <- readr::read_csv(file.path(STABL_TBL_DIR, "stabl_15_cox_PFS_stabl_selected_features.csv"),
                              show_col_types = FALSE)

km_summary_list <- list()

# OS: GDF15, NOTUM
for (i in seq_len(nrow(cox_os_id))) {
  id  <- cox_os_id$nELISA_ID[i]
  lbl <- cox_os_id$Protein_Label[i]
  message("  KM OS: ", lbl)
  res <- plot_km_median_split(
    survival_df, expr_log2_rep[survival_df$sample_id, id],
    "os_months", "os_event", lbl,
    "STABL Cox OS discovery | whole matched cohort",
    file.path(OUT_FIG, paste0("cox_os_", lbl))
  )
  if (!is.null(res)) km_summary_list[[paste0("OS_", lbl)]] <- res %>% dplyr::mutate(endpoint = "OS")
}

# PFS: CEACAM5
for (i in seq_len(nrow(cox_pfs_id))) {
  id  <- cox_pfs_id$nELISA_ID[i]
  lbl <- cox_pfs_id$Protein_Label[i]
  message("  KM PFS: ", lbl)
  res <- plot_km_median_split(
    survival_df, expr_log2_rep[survival_df$sample_id, id],
    "pfs_months_main", "pfs_event_main", lbl,
    "STABL Cox PFS discovery | whole matched cohort",
    file.path(OUT_FIG, paste0("cox_pfs_", lbl))
  )
  if (!is.null(res)) km_summary_list[[paste0("PFS_", lbl)]] <- res %>% dplyr::mutate(endpoint = "PFS")
}

# ─────────────────────────────────────────────────────────────
# PART 2: Recurrence-discovered proteins — KM for PFS + OS
#         within their discovery subgroup (robustness check)
# ─────────────────────────────────────────────────────────────
message("\n=== 08c | Part 2: Recurrence-discovered proteins (subgroup KM) ===")

recurrence_checks <- list(
  list(file = "stabl_07_Sx_recurrence_stabl_selected_features.csv",
       group = "Sx",       tag = "Sx_recurrence"),
  list(file = "stabl_08_IO_recurrence_stabl_selected_features.csv",
       group = "IO",       tag = "IO_recurrence"),
  list(file = "stabl_09_CTxIO_recurrence_stabl_selected_features.csv",
       group = "CTx + IO", tag = "CTxIO_recurrence")
)

for (check in recurrence_checks) {
  feat_path <- file.path(STABL_TBL_DIR, check$file)
  if (!file.exists(feat_path)) { message("  Skipping (no file): ", check$tag); next }
  feat_df <- readr::read_csv(feat_path, show_col_types = FALSE)
  if (!nrow(feat_df)) next
  id  <- feat_df$nELISA_ID[1]
  lbl <- feat_df$Protein_Label[1]
  sub_df <- analysis1_df %>% dplyr::filter(new_group == check$group)

  message("  ", check$tag, " | ", lbl, " — testing PFS and OS within subgroup (n=", nrow(sub_df), ")")

  res_pfs <- plot_km_median_split(
    sub_df, expr_log2_rep[sub_df$sample_id, id],
    "pfs_months_main", "pfs_event_main", lbl,
    paste0("STABL recurrence discovery (", check$group, ") | PFS within subgroup"),
    file.path(OUT_FIG, paste0(check$tag, "_PFS_", sanitize_filename(lbl))),
    min_n = 10, min_events = 4
  )
  if (!is.null(res_pfs))
    km_summary_list[[paste0(check$tag, "_PFS")]] <-
      res_pfs %>% dplyr::mutate(endpoint = paste0("PFS (", check$group, " subgroup)"))

  res_os <- plot_km_median_split(
    sub_df, expr_log2_rep[sub_df$sample_id, id],
    "os_months", "os_event", lbl,
    paste0("STABL recurrence discovery (", check$group, ") | OS within subgroup"),
    file.path(OUT_FIG, paste0(check$tag, "_OS_", sanitize_filename(lbl))),
    min_n = 10, min_events = 4
  )
  if (!is.null(res_os))
    km_summary_list[[paste0(check$tag, "_OS")]] <-
      res_os %>% dplyr::mutate(endpoint = paste0("OS (", check$group, " subgroup)"))
}

# ─────────────────────────────────────────────────────────────
# PART 3: Pre/post-discovered protein (AZU1) — paired plot
# ─────────────────────────────────────────────────────────────
message("\n=== 08c | Part 3: AZU1 pre/post paired plot ===")

azu1_path <- file.path(STABL_TBL_DIR, "stabl_12_Sx_pre_vs_post_stabl_selected_features.csv")
if (file.exists(azu1_path)) {
  azu1_df <- readr::read_csv(azu1_path, show_col_types = FALSE)
  if (nrow(azu1_df)) {
    id  <- azu1_df$nELISA_ID[1]
    lbl <- azu1_df$Protein_Label[1]

    sx_pre_ids  <- paired_pre  %>% dplyr::filter(patient_id %in% surgery_pair_ids) %>%
                   dplyr::arrange(patient_id)
    sx_post_ids <- paired_post %>% dplyr::filter(patient_id %in% surgery_pair_ids) %>%
                   dplyr::arrange(patient_id)

    pre_vals  <- expr_log2_rep[sx_pre_ids$sample_id,  id]
    post_vals <- expr_log2_rep[sx_post_ids$sample_id, id]

    op <- options(warn = -1)
    wilcox_p <- tryCatch(
      stats::wilcox.test(pre_vals, post_vals, paired = TRUE)$p.value,
      error = function(e) NA_real_)
    options(op)

    message("  Plotting AZU1 paired pre/post (n=", length(pre_vals), " pairs)")
    plot_paired_dotline(pre_vals, post_vals, sx_pre_ids$patient_id, lbl,
                        file.path(OUT_FIG, paste0("paired_", lbl)), wilcox_p)
  }
}

# ─────────────────────────────────────────────────────────────
# PART 4: Response-discovered proteins (CEACAM5, SHH) — KM for
#         PFS + OS on whole cohort (robustness check, same logic
#         as Part 2 for recurrence proteins). Also a sanity-check
#         KM of response category itself vs PFS/OS.
# ─────────────────────────────────────────────────────────────
message("\n=== 08c | Part 4: Response-discovered proteins + response sanity check ===")

# 4a. Response-discovered proteins: do they ALSO track survival?
response_proteins <- tibble::tribble(
  ~Gene_Symbol, ~nELISA_ID,
  "CEACAM5",    "nEL07541",
  "SHH",        "nEL08721"
)

for (i in seq_len(nrow(response_proteins))) {
  id  <- response_proteins$nELISA_ID[i]
  lbl <- response_proteins$Gene_Symbol[i]
  message("  ", lbl, " — testing PFS and OS on whole cohort")

  res_pfs <- plot_km_median_split(
    survival_df, expr_log2_rep[survival_df$sample_id, id],
    "pfs_months_main", "pfs_event_main", lbl,
    "STABL response discovery (CR/PR vs PD) | PFS whole cohort",
    file.path(OUT_FIG, paste0("response_disc_PFS_", lbl))
  )
  if (!is.null(res_pfs))
    km_summary_list[[paste0("response_disc_PFS_", lbl)]] <-
      res_pfs %>% dplyr::mutate(endpoint = "PFS (response-discovered)")

  res_os <- plot_km_median_split(
    survival_df, expr_log2_rep[survival_df$sample_id, id],
    "os_months", "os_event", lbl,
    "STABL response discovery (CR/PR vs PD) | OS whole cohort",
    file.path(OUT_FIG, paste0("response_disc_OS_", lbl))
  )
  if (!is.null(res_os))
    km_summary_list[[paste0("response_disc_OS_", lbl)]] <-
      res_os %>% dplyr::mutate(endpoint = "OS (response-discovered)")
}

# 4b. Sanity check: response category itself vs PFS/OS
# (confirms responders survive longer — expected, but validates
# that response_4cat is behaving as a real clinical signal here)
message("  Sanity check: response category (PR/CR vs SD/PD) vs PFS/OS")

resp_df <- survival_df %>% dplyr::filter(!is.na(response_binary))

for (endpoint in list(
  list(time = "pfs_months_main", event = "pfs_event_main", tag = "PFS"),
  list(time = "os_months",       event = "os_event",       tag = "OS")
)) {
  op <- options(warn = -1)
  time_val  <- as.numeric(resp_df[[endpoint$time]])
  options(op)
  event_val <- coerce_event_indicator(resp_df[[endpoint$event]])
  keep <- !is.na(time_val) & !is.na(event_val)
  df_r <- resp_df[keep, ] %>%
    dplyr::mutate(time_value = time_val[keep], event_value = event_val[keep],
                  response_group = factor(response_binary, levels = c("SD/PD", "PR/CR")))

  fit_r <- survival::survfit(survival::Surv(time_value, event_value) ~ response_group, data = df_r)
  diff_r <- tryCatch(
    survival::survdiff(survival::Surv(time_value, event_value) ~ response_group, data = df_r),
    error = function(e) NULL)
  p_r <- if (!is.null(diff_r))
    stats::pchisq(diff_r$chisq, df = max(length(diff_r$n) - 1, 1), lower.tail = FALSE) else NA_real_

  km_resp <- survminer::ggsurvplot(
    fit_r, data = df_r, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    palette = c("#2166AC", "#B2182B"),
    legend.title = "Response", legend.labs = c("SD/PD", "PR/CR"),
    title    = paste0("KM — Response category sanity check (", endpoint$tag, ")"),
    subtitle = paste0("PR/CR vs SD/PD | logrank p=", scales::number(p_r, accuracy = 0.001)),
    xlab = "Months", ylab = "Survival probability",
    ggtheme = ggplot2::theme_bw(base_size = 11)
  )
  stub_r <- file.path(OUT_FIG, paste0("response_sanity_check_", endpoint$tag))
  if (SKIP_IF_DONE && file.exists(paste0(stub_r, ".pdf"))) {
    message("    SKIP (already exists): response_sanity_check_", endpoint$tag)
  } else {
    tryCatch({
      cairo_pdf(paste0(stub_r, ".pdf"), width = 7.8, height = 6.5); print(km_resp); dev.off()
      png(paste0(stub_r, ".png"), width = 7.8, height = 6.5,
          units = "in", res = 300, bg = "white"); print(km_resp); dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("    ERROR: ", conditionMessage(e), " -> Close PDF viewers and re-run.")
    })
  }

  km_summary_list[[paste0("response_sanity_", endpoint$tag)]] <- tibble::tibble(
    protein = "Response category (PR/CR vs SD/PD)", n = nrow(df_r),
    n_events = sum(df_r$event_value == 1),
    n_low = sum(df_r$response_group == "SD/PD"), n_high = sum(df_r$response_group == "PR/CR"),
    logrank_p = p_r, endpoint = paste0(endpoint$tag, " (response sanity check)")
  )
}
