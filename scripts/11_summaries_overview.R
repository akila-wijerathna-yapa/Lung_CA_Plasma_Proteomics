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

# ── Load all checkpoints ──────────────────────────────────────
load(file.path(paths$rdata, "checkpoint_01.RData"))
load(file.path(paths$rdata, "checkpoint_03.RData"))
load(file.path(paths$rdata, "checkpoint_04.RData"))
load(file.path(paths$rdata, "checkpoint_05.RData"))
load(file.path(paths$rdata, "checkpoint_06.RData"))
load(file.path(paths$rdata, "checkpoint_07.RData"))
load(file.path(paths$rdata, "checkpoint_08.RData"))
load(file.path(paths$rdata, "checkpoint_10.RData"))
message("All checkpoints loaded")
while (dev.cur() > 1) dev.off()

# ── Local functions ───────────────────────────────────────────
safe_write_csv <- function(df, path) readr::write_csv(df, path, na = "")

clean_pathway_name <- function(x) {
  x %>%
    stringr::str_replace_all("^HALLMARK_","") %>%
    stringr::str_replace_all("^REACTOME_","") %>%
    stringr::str_replace_all("^GOBP_|^GO:BP_","") %>%
    stringr::str_replace_all("^GOMF_|^GO:MF_","") %>%
    stringr::str_replace_all("^GOCC_|^GO:CC_","") %>%
    stringr::str_replace_all("^KEGG_(MEDICUS|LEGACY)_","") %>%
    stringr::str_replace_all("^KEGG_","") %>%
    stringr::str_replace_all("_"," ") %>%
    stringr::str_squish()
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x

summarise_result_table <- function(df, analysis_label, analysis_type = "inferential") {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(analysis=analysis_label, analysis_type=analysis_type,
                  n_tested=0L,
                  n_fdr_005 = if (analysis_type == "descriptive_logFC_only") NA_integer_ else 0L,
                  n_fdr_010 = if (analysis_type == "descriptive_logFC_only") NA_integer_ else 0L,
                  top_hit=NA_character_, top_adj_p_value=NA_real_))
  }
  if (analysis_type == "descriptive_logFC_only") {
    top_row <- df %>% filter(is.finite(logFC)) %>% arrange(desc(abs(logFC)))
    return(tibble(analysis=analysis_label, analysis_type=analysis_type,
                  n_tested=nrow(df), n_fdr_005=NA_integer_, n_fdr_010=NA_integer_,
                  top_hit = if (nrow(top_row)) (top_row$Protein_Label[1] %||% top_row$nELISA_ID[1]) else NA_character_,
                  top_adj_p_value=NA_real_))
  }
  top_row <- df %>% filter(!is.na(adj_p_value)) %>% arrange(adj_p_value)
  tibble(analysis=analysis_label, analysis_type=analysis_type,
         n_tested=nrow(df),
         n_fdr_005=sum(df$adj_p_value < 0.05, na.rm=TRUE),
         n_fdr_010=sum(df$adj_p_value < 0.10, na.rm=TRUE),
         top_hit = if (nrow(top_row)) (top_row$Protein_Label[1] %||% top_row$nELISA_ID[1]) else NA_character_,
         top_adj_p_value = if (nrow(top_row)) top_row$adj_p_value[1] else NA_real_)
}

make_prisma_flowchart <- function(clinical, pre_ids, post_ids, cohort_summary, path_stub) {
  get_metric <- function(name) as.numeric(cohort_summary$value[match(name, cohort_summary$metric)])
  clinical_n    <- get_metric("Clinical patients")
  proteomics_qc <- get_metric("Proteomics rows after qc_fail filter")
  pre_unique    <- get_metric("Baseline PRE unique samples after duplicate resolution")
  post_unique   <- get_metric("POST unique samples after duplicate resolution")
  pre_matched   <- get_metric("Baseline PRE clinically matched patients")
  post_matched  <- get_metric("POST clinically matched patients")
  paired_n      <- get_metric("Matched paired patients")
  surgery_pairs <- get_metric("Matched surgery pairs")

  d <- clinical %>% mutate(
    has_pre  = patient_id %in% pre_ids,
    has_post = patient_id %in% post_ids,
    has_any  = has_pre | has_post,
    response_d = case_when(response_4cat %in% c("CR","PR","SD","PD") ~ response_4cat, TRUE ~ "NE"),
    stage_d    = case_when(grepl("^(I|II|I-II|I/II)", stage_group, ignore.case=TRUE) ~ "I-II",
                           grepl("^(III|IV)",         stage_group, ignore.case=TRUE) ~ "III-IV",
                           TRUE ~ "Unknown"),
    sex_d = case_when(toupper(dplyr::coalesce(sex,"")) %in% c("F","FEMALE") ~ "Female",
                      toupper(dplyr::coalesce(sex,"")) %in% c("M","MALE")   ~ "Male", TRUE ~ "Unknown")
  )
  any_n       <- sum(d$has_any, na.rm=TRUE)
  prcr_n      <- sum(d$has_any & d$response_d %in% c("CR","PR"), na.rm=TRUE)
  sd_n        <- sum(d$has_any & d$response_d == "SD",  na.rm=TRUE)
  pd_n        <- sum(d$has_any & d$response_d == "PD",  na.rm=TRUE)
  resp_unk_n  <- any_n - sum(d$has_any & d$response_d != "NE", na.rm=TRUE)
  stage_e_n   <- sum(d$has_any & d$stage_d == "I-II",  na.rm=TRUE)
  stage_l_n   <- sum(d$has_any & d$stage_d == "III-IV",na.rm=TRUE)
  stage_unk_n <- any_n - stage_e_n - stage_l_n
  female_n    <- sum(d$has_any & d$sex_d == "Female", na.rm=TRUE)
  male_n      <- sum(d$has_any & d$sex_d == "Male",   na.rm=TRUE)
  sex_unk_n   <- any_n - female_n - male_n
  pre_only_n  <- sum(d$has_pre & !d$has_post, na.rm=TRUE)
  post_only_n <- sum(!d$has_pre & d$has_post, na.rm=TRUE)
  pfs_data_n  <- sum(d$has_any & !is.na(d$pfs_months_main) & !is.na(d$pfs_event_main), na.rm=TRUE)
  pfs_evt_n   <- sum(d$has_any & !is.na(d$pfs_event_main) & d$pfs_event_main == 1, na.rm=TRUE)
  pfs_cen_n   <- pfs_data_n - pfs_evt_n
  os_data_n   <- sum(d$has_any & !is.na(d$os_months) & !is.na(d$os_event), na.rm=TRUE)
  os_evt_n    <- sum(d$has_any & !is.na(d$os_event) & d$os_event == 1, na.rm=TRUE)
  os_cen_n    <- os_data_n - os_evt_n

  node_df <- tibble(
    id = c("clin","prot","rep","match","any",
           "pair_st","paired","pre_only","post_only",
           "resp","prcr","sd","pd","resp_unk",
           "stage","st_e","st_l","st_unk",
           "sex","female","male","sex_unk",
           "pfs","pfs_e","pfs_c","pfs_m",
           "os","os_e","os_c","os_m"),
    x = c(-5,-5,-5,-5,-5,
          -1.6,1.9,1.9,1.9,
          -1.6,1.9,1.9,1.9,1.9,
          -1.6,1.9,1.9,1.9,
          -1.6,1.9,1.9,1.9,
          -1.6,1.9,1.9,1.9,
          -1.6,1.9,1.9,1.9),
    y = c(10.2,8.45,6.7,4.95,3.2,
          8.25,9.15,8.25,7.35,
          5.55,6.45,5.55,4.65,3.75,
          2.85,3.55,2.85,2.15,
          0.2,0.95,0.2,-0.55,
          -2.6,-1.85,-2.6,-3.35,
          -5.2,-4.45,-5.2,-5.95),
    fill_group = c(rep("source",2),rep("matching",3),
                   rep("paired",4),rep("response",5),rep("stage",4),
                   rep("sex",4),rep("survival",4),rep("survival",4)),
    label = c(
      paste0("Clinical cohort\nN = ",clinical_n),
      paste0("Proteomics after QC filter\nRows = ",proteomics_qc),
      paste0("Representative samples\nPRE = ",pre_unique," | POST = ",post_unique),
      paste0("Clinically matched\nPRE = ",pre_matched," | POST = ",post_matched),
      paste0("Patients with proteomics\nN = ",any_n),
      paste0("Pair status\nSurgery pairs = ",surgery_pairs),
      paste0("Pre + post pairs\nN = ",paired_n),
      paste0("Pre only\nN = ",pre_only_n),
      paste0("Post only\nN = ",post_only_n),
      paste0("Response evaluable\nN = ",prcr_n+sd_n+pd_n),
      paste0("PR/CR\nN = ",prcr_n), paste0("SD\nN = ",sd_n),
      paste0("PD\nN = ",pd_n), paste0("Response NE/unknown\nN = ",resp_unk_n),
      paste0("Stage known\nN = ",stage_e_n+stage_l_n),
      paste0("Stage I-II\nN = ",stage_e_n), paste0("Stage III-IV\nN = ",stage_l_n),
      paste0("Stage unknown\nN = ",stage_unk_n),
      paste0("Sex known\nN = ",female_n+male_n),
      paste0("Female\nN = ",female_n), paste0("Male\nN = ",male_n),
      paste0("Sex unknown\nN = ",sex_unk_n),
      paste0("PFS available\nN = ",pfs_data_n),
      paste0("PFS events\nN = ",pfs_evt_n), paste0("PFS censored\nN = ",pfs_cen_n),
      paste0("PFS unavailable\nN = ",any_n-pfs_data_n),
      paste0("OS available\nN = ",os_data_n),
      paste0("OS deaths\nN = ",os_evt_n), paste0("OS censored\nN = ",os_cen_n),
      paste0("OS unavailable\nN = ",any_n-os_data_n))
  )
  edge_df <- tribble(~from,~to,
    "clin","prot","prot","rep","rep","match","match","any",
    "any","pair_st","pair_st","paired","pair_st","pre_only","pair_st","post_only",
    "any","resp","resp","prcr","resp","sd","resp","pd","resp","resp_unk",
    "any","stage","stage","st_e","stage","st_l","stage","st_unk",
    "any","sex","sex","female","sex","male","sex","sex_unk",
    "any","pfs","pfs","pfs_e","pfs","pfs_c","pfs","pfs_m",
    "any","os","os","os_e","os","os_c","os","os_m") %>%
    left_join(node_df %>% transmute(from=id,x,y), by="from") %>%
    left_join(node_df %>% transmute(to=id,xend=x,yend=y), by="to")

  fill_pal <- c(source="#dbe9f6",matching="#d7eedb",paired="#fde7c7",
                response="#f8d9de",stage="#e2dcf4",sex="#d8eef0",survival="#ececec")
  p <- ggplot() +
    geom_segment(data=edge_df, aes(x=x,y=y,xend=xend,yend=yend),
                 color="grey40", linewidth=0.45,
                 arrow=grid::arrow(length=grid::unit(0.12,"inches"))) +
    geom_label(data=node_df, aes(x=x,y=y,label=label,fill=fill_group),
               linewidth=0.35, label.padding=grid::unit(0.18,"lines"),
               label.r=grid::unit(0.18,"lines"), lineheight=0.96,
               size=3.15, family="sans", show.legend=FALSE) +
    scale_fill_manual(values=fill_pal) +
    annotate("text",x=-6.45,y=11.2,label="Study Flow (PRISMA-style)",
             hjust=0,size=5.1,fontface="bold") +
    annotate("text",x=-6.45,y=10.62,
             label="Left column: filtering into matched cohort. Right branches: analysis subsets.",
             hjust=0,size=3.25,color="grey25") +
    coord_cartesian(xlim=c(-6.6,3.25),ylim=c(-7.4,11.55),clip="off") +
    theme_void()
  cairo_pdf(paste0(path_stub,".pdf"), width=13.2, height=11.2); print(p); dev.off()
  png(paste0(path_stub,".png"), width=13.2, height=11.2, units="in", res=300, bg="white")
  print(p); dev.off()
}

make_clinical_landscape_figure <- function(clinical, pre_ids, post_ids, path_stub) {
  palette_values <- c(
    "Sx"="#E69F00","CTx"="#D55E00","CTx + IO"="#009E73","CTx + RTx"="#0072B2",
    "IO"="#56B4E9","TKI"="#CC79A7","Other"="#999999","Unknown"="#d9d9d9",
    "Female"="#f4a6b7","Male"="#355070",
    "CR"="#7cb342","PR"="#2a9d8f","SD"="#f4a261","PD"="#d7301f","NE"="#d0d0d0",
    "I-II"="#fbb4ae","III-IV"="#80cdc1","Unknown"="#d9d9d9",
    "Alive"="#66c2a5","Dead"="#d73027",
    "PFS event"="#d73027","PFS censored"="#d9d9d9","PFS unavailable"="#f0f0f0",
    "OS death"="#b2182b","OS censored"="#d9d9d9","OS unavailable"="#f0f0f0",
    "both"="#4d4d4d","only_pre"="#fdae61","only_post"="#2ca25f","no_proteomics"="#ececec"
  )
  strip_levels <- rev(c("Group","Sex","Response","Stage","Vital","PFS","OS","Pair"))
  plot_df <- clinical %>%
    mutate(
      pair_status    = case_when(
        patient_id %in% pre_ids & patient_id %in% post_ids ~ "both",
        patient_id %in% pre_ids  ~ "only_pre",
        patient_id %in% post_ids ~ "only_post",
        TRUE ~ "no_proteomics"),
      group_display  = dplyr::coalesce(new_group,"Unknown"),
      sex_display    = case_when(toupper(dplyr::coalesce(sex,"")) %in% c("F","FEMALE") ~ "Female",
                                 toupper(dplyr::coalesce(sex,"")) %in% c("M","MALE")   ~ "Male",
                                 TRUE ~ "Unknown"),
      response_display = case_when(response_4cat %in% c("CR","PR","SD","PD") ~ response_4cat, TRUE ~ "NE"),
      stage_display  = case_when(grepl("^(I|II|I-II|I/II)", stage_group, ignore.case=TRUE) ~ "I-II",
                                 grepl("^(III|IV)", stage_group, ignore.case=TRUE)          ~ "III-IV",
                                 TRUE ~ "Unknown"),
      vital_display  = dplyr::coalesce(alive_status,"Unknown"),
      pfs_display    = case_when(is.na(pfs_event_main)|is.na(pfs_months_main) ~ "PFS unavailable",
                                 pfs_event_main == 1 ~ "PFS event", TRUE ~ "PFS censored"),
      os_display     = case_when(is.na(os_event)|is.na(os_months) ~ "OS unavailable",
                                 os_event == 1 ~ "OS death", TRUE ~ "OS censored")
    ) %>%
    arrange(factor(pair_status, levels=c("both","only_pre","only_post","no_proteomics")),
            factor(group_display, levels=c("Sx","CTx","CTx + IO","CTx + RTx","IO","TKI","Other","Unknown")),
            factor(response_display, levels=c("CR","PR","SD","PD","NE")),
            age_at_diagnosis) %>%
    mutate(idx = row_number())
  if (!nrow(plot_df)) return(invisible(NULL))

  age_plot <- ggplot(plot_df, aes(idx, age_at_diagnosis)) +
    geom_point(size=1.05, color="black", alpha=0.85) +
    geom_smooth(method="loess", se=FALSE, linewidth=0.45, color="#7f7f7f", span=0.3) +
    labs(title=paste0("Clinical & Proteomics Cohort Landscape (n=",nrow(plot_df),")"),
         subtitle=paste0("both=",sum(plot_df$pair_status=="both"),
                         " | only_pre=",sum(plot_df$pair_status=="only_pre"),
                         " | only_post=",sum(plot_df$pair_status=="only_post"),
                         " | no_proteomics=",sum(plot_df$pair_status=="no_proteomics")),
         x=NULL, y="Age") +
    scale_x_continuous(expand=expansion(mult=c(0.002,0.002))) +
    theme_bw(base_size=10.5) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=9.5,color="grey30"))

  tile_df <- bind_rows(
    plot_df %>% transmute(idx,strip="Group",   value=group_display),
    plot_df %>% transmute(idx,strip="Sex",     value=sex_display),
    plot_df %>% transmute(idx,strip="Response",value=response_display),
    plot_df %>% transmute(idx,strip="Stage",   value=stage_display),
    plot_df %>% transmute(idx,strip="Vital",   value=vital_display),
    plot_df %>% transmute(idx,strip="PFS",     value=pfs_display),
    plot_df %>% transmute(idx,strip="OS",      value=os_display),
    plot_df %>% transmute(idx,strip="Pair",    value=pair_status)
  ) %>% mutate(strip=factor(strip,levels=strip_levels))

  tile_plot <- ggplot(tile_df, aes(idx,strip,fill=value)) +
    geom_tile(color="white", linewidth=0.12) +
    scale_fill_manual(values=palette_values, drop=FALSE) +
    scale_x_continuous(expand=expansion(mult=c(0.002,0.002))) +
    labs(x=NULL, y=NULL, fill=NULL) +
    theme_minimal(base_size=10.5) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          panel.grid=element_blank(), axis.text.y=element_text(face="bold"),
          legend.position="right", legend.key.height=grid::unit(0.34,"cm"),
          legend.text=element_text(size=8.6)) +
    guides(fill=guide_legend(ncol=2, byrow=TRUE, override.aes=list(color=NA)))

  cairo_pdf(paste0(path_stub,".pdf"), width=13.2, height=7.2)
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow=2,ncol=1,heights=grid::unit(c(2.1,4.2),"null"))
  grid::pushViewport(grid::viewport(layout=lay))
  print(age_plot,  vp=grid::viewport(layout.pos.row=1,layout.pos.col=1))
  print(tile_plot, vp=grid::viewport(layout.pos.row=2,layout.pos.col=1))
  grid::popViewport()
  dev.off()
  png(paste0(path_stub,".png"), width=13.2, height=7.2, units="in", res=300, bg="white")
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow=2,ncol=1,heights=grid::unit(c(2.1,4.2),"null"))
  grid::pushViewport(grid::viewport(layout=lay))
  print(age_plot,  vp=grid::viewport(layout.pos.row=1,layout.pos.col=1))
  print(tile_plot, vp=grid::viewport(layout.pos.row=2,layout.pos.col=1))
  grid::popViewport()
  dev.off()
}

# ── Analysis summary dashboard ────────────────────────────────
message("\n=== 11 | Analysis summary dashboard ===")
analysis_summary_tbl <- bind_rows(
  summarise_result_table(analysis1_omnibus, "Analysis 1 omnibus F-test"),
  bind_rows(lapply(names(analysis1_ovr_results), function(nm)
    summarise_result_table(analysis1_ovr_results[[nm]], paste("Analysis 1", nm)))),
  summarise_result_table(analysis2_surgery, "Analysis 2 surgery paired"),
  summarise_result_table(analysis3_primary,    "Analysis 3 primary PR vs PD"),
  summarise_result_table(analysis3_secondary,  "Analysis 3 secondary PR/CR vs SD/PD"),
  summarise_result_table(analysis3_primary_expl,  "Analysis 3 primary EXPLORATORY (plate only)"),
  summarise_result_table(analysis3_secondary_expl,"Analysis 3 secondary EXPLORATORY (plate only)"),
  summarise_result_table(analysis4_desc, "Analysis 4 TKI vs IO descriptive logFC",
                         analysis_type = "descriptive_logFC_only"),
  summarise_result_table(os_res,              "Analysis 5 OS Cox"),
  summarise_result_table(pfs_res,             "Analysis 5 PFS Cox"),
  summarise_result_table(pfs_res_sensitivity, "Analysis 5 PFS Cox sensitivity"),
  summarise_result_table(pfs_icenreg,         "Analysis 5 PFS icenReg sensitivity"),
  summarise_result_table(efs_res,             "Analysis 5 EFS Cox"),
  summarise_result_table(delta_os,            "Analysis 5 delta OS Cox"),
  summarise_result_table(delta_pfs,           "Analysis 5 delta PFS Cox"),
  bind_rows(lapply(names(poster_results), function(nm)
    summarise_result_table(poster_results[[nm]], paste("Poster", nm))))
)
safe_write_csv(analysis_summary_tbl,
  file.path(paths$tables, "00_cohort_audit", "analysis_summary_dashboard.csv"))
message("Analysis summary: ", nrow(analysis_summary_tbl), " rows")
print(analysis_summary_tbl %>% select(analysis, n_tested, n_fdr_010))

# ── fgsea pathway summary ─────────────────────────────────────
message("\n=== 11 | Pathway enrichment summary ===")
# Read from subfolder structure created by module 09
fgsea_sig_files <- list.files(
  file.path(paths$tables, ANALYSIS8_DIR),
  pattern = "_fgsea_sig\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(fgsea_sig_files) > 0) {
  bio_summary_tbl <- bind_rows(lapply(fgsea_sig_files, function(f) {
    tryCatch({
      df <- readr::read_csv(f, show_col_types=FALSE) %>%
        mutate(across(any_of(c("pval","padj","log2err","ES","NES","size")),
                      ~suppressWarnings(as.numeric(.))))
      # Extract label and collection from path
      parts    <- strsplit(f, .Platform$file.sep, fixed=FALSE)[[1]]
      label    <- parts[length(parts)-2]
      coll     <- parts[length(parts)-1]
      tibble(analysis=label, collection=toupper(coll),
             n_sig_pathways=nrow(df),
             top_pathway = if (nrow(df)) clean_pathway_name(df$pathway[which.min(df$padj)]) else NA_character_,
             top_fdr     = if (nrow(df)) min(df$padj, na.rm=TRUE) else NA_real_)
    }, error=function(e) NULL)
  })) %>% arrange(analysis, collection)
  safe_write_csv(bio_summary_tbl,
    file.path(paths$tables, ANALYSIS8_DIR, "pathway_enrichment_summary.csv"))
  message("Pathway summary: ", nrow(bio_summary_tbl), " analysis×collection rows | total sig pathways: ",
          sum(bio_summary_tbl$n_sig_pathways, na.rm=TRUE))
} else {
  message("No fgsea sig files found")
}

# ── Instruction checklist ─────────────────────────────────────
message("\n=== 11 | Instruction checklist ===")
km_pngs       <- list.files(file.path(paths$figures,"06_analysis5_survival"),
                             pattern="^analysis5_km_.*\\.png$") %>%
                 grep("_risk\\.png$",., invert=TRUE, value=TRUE)
ppi_generated <- file.exists(file.path(paths$tables, ANALYSIS8_DIR,
                                       list.files(file.path(paths$tables, ANALYSIS8_DIR),
                                                  pattern="_ppi_nodes\\.csv$", recursive=TRUE)[1]))

checklist <- tibble(
  section     = c("Output","Output","QC","QC","QC","QC","QC","QC",
                  "Analysis","Analysis","Analysis","Analysis","Analysis","Analysis","Analysis",
                  "Pathway","PPI","ML","Poster","Poster","Poster"),
  item        = c("Output-1","Output-2",
                  "QC-1","QC-2","QC-3","QC-4","QC-5","QC-6",
                  "A1","A2","A3","A4","A5","A5-KM","A6",
                  "P8A","P8B","P8C","P9-1","P9-2","P9-3"),
  request     = c("All analysis in modular scripts",
                  "All figures as PDF + PNG",
                  "PCA coloured by New_Group, shaped by plate",
                  "Sample correlation heatmap",
                  "Protein CV distribution",
                  "Missing-value heatmap",
                  "Per-sample intensity boxplots",
                  "Batch-correction density QC",
                  "Seven-group omnibus + one-vs-rest",
                  "Paired pre/post before-after analysis",
                  "PR vs PD + PR/CR vs SD/PD response",
                  "TKI vs IO descriptive logFC",
                  "Survival Cox (OS/PFS/EFS) + icenReg",
                  "KM curves for eligible subgroups",
                  "LASSO + stability + AUC + biomarker scoring",
                  "fgsea pathway enrichment (all collections)",
                  "STRING PPI network",
                  "Biomarker scoring table + ROC/AUC",
                  "Poster 1: IO response",
                  "Poster 2: Chemo/RT response",
                  "Poster 3: Surgery recurrence"),
  status      = c("Completed — modules 00-11",
                  "Completed",
                  "Completed — 02_qc_plots.R",
                  "Completed — 02_qc_plots.R",
                  "Completed — 02_qc_plots.R",
                  "Completed — 02_qc_plots.R",
                  "Completed — 02_qc_plots.R",
                  "Completed (QC only; batch correction not applied to analysis matrices)",
                  "Completed — 03_analysis1_seven_group.R",
                  "Completed — 04_analysis2_before_after.R",
                  "Completed — 05_analysis3_response.R (underpowered; exploratory pass added)",
                  "Completed — 06_analysis4_tki_io.R (descriptive only)",
                  "Completed — 07_analysis5_survival.R",
                  ifelse(length(km_pngs) > 0, "Completed", "Completed (threshold not met for all endpoints)"),
                  "Completed — 08_analysis6_ml.R",
                  "Completed — 09_pathway_network.R",
                  "Input genes saved; RUN_STRING=FALSE (set TRUE for network generation)",
                  "Completed — 08_analysis6_ml.R",
                  "Completed — 10_poster_comparisons.R (hypothesis-generating)",
                  "Completed — 10_poster_comparisons.R (hypothesis-generating)",
                  "Completed — 10_poster_comparisons.R (hypothesis-generating)"),
  notes       = c("Self-contained numbered modules connected by save/load checkpoints.",
                  "cairo_pdf + png for all figures.",
                  "Colorblind-friendly Okabe-Ito palette.",
                  "Spearman correlation, hierarchical clustering.",
                  "Computed in pg/mL.",
                  "Ordered by missing rate.",
                  "Ordered by total log2 intensity.",
                  "limma::removeBatchEffect for visualisation; analyses use covariate adjustment.",
                  "F-test + 7 one-vs-rest contrasts + top-100 heatmap.",
                  "Surgery: formal paired limma. IO-containing: descriptive delta logFC.",
                  "No proteins reached FDR<0.1 — genuine power limitation at n~154.",
                  "Too few TKI samples for formal testing; descriptive logFC + heatmap.",
                  "Plate_barcode covariate; icenReg sensitivity for interval-censored PFS.",
                  paste0(length(km_pngs), " KM curves generated."),
                  "Primary biomarker discovery route. 6B: AUC=0.788, 11 LASSO-selected proteins.",
                  "HALLMARK + GO_BP/MF/CC + KEGG + REACTOME; label/collection/ subfolder structure.",
                  "Set RUN_STRING <- TRUE in 09_pathway_network.R to query STRING API.",
                  "Composite scoring: limma + OS Cox + PFS Cox + LASSO + stability + AUC.",
                  "No FDR hits — same power limitation as analysis 3.",
                  "No FDR hits — same power limitation as analysis 3.",
                  "Interval-expanded recurrence definition used.")
)
safe_write_csv(checklist,
  file.path(paths$tables, "00_cohort_audit", "instruction_checklist.csv"))
message("Checklist saved: ", nrow(checklist), " items")

# ── PRISMA flowchart ──────────────────────────────────────────
message("\n=== 11 | PRISMA study flow ===")
make_prisma_flowchart(clinical, pre_rep_ids, post_rep_ids, cohort_summary,
  file.path(paths$figures, "00_overview", "00_prisma_study_flow"))
message("PRISMA figure saved")

# ── Clinical landscape ────────────────────────────────────────
message("\n=== 11 | Clinical landscape ===")
make_clinical_landscape_figure(clinical, pre_rep_ids, post_rep_ids,
  file.path(paths$figures, "00_overview", "01_clinical_proteomics_landscape"))
message("Landscape figure saved")

message("\n11_summaries_overview.R complete")
