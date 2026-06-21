#!/usr/bin/env Rscript

# Project: CRC MMR spatial immune ecology
# Purpose: Score GSE39582 and GSE17536 and perform KM survival analysis.

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(GSVA)
  library(survival)
  library(survminer)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    data_dir = "data/external/GEO",
    gpl570 = "data/external/GPL570.csv",
    marker_list = "signatures/marker_list.Rdata",
    output_dir = "results/07_external_validation/survival",
    download_if_missing = "TRUE"
  )
  for (arg in args) {
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2 && kv[1] %in% names(out)) out[[kv[1]]] <- kv[2]
  }
  out$download_if_missing <- toupper(out$download_if_missing) %in% c("TRUE", "T", "1", "YES")
  out
}

args <- parse_args()
dir.create(args$data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

load(args$marker_list)  # marker_list
if (!exists("marker_list") || !is.list(marker_list)) stop("marker_list.Rdata must contain marker_list")
# Keep marker_list as the final signature source, but remove generic TLS entry
# and rename the internal Stroma4_1 score to a manuscript-facing name.
if ("TLS" %in% names(marker_list)) marker_list[["TLS"]] <- NULL
if ("Stroma4_1" %in% names(marker_list)) names(marker_list)[names(marker_list) == "Stroma4_1"] <- "dMMR_TLS_program"


download_gse <- function(gse, data_dir) {
  expr_file <- file.path(data_dir, paste0(gse, "_exprSet.csv"))
  meta_file <- file.path(data_dir, paste0(gse, "_metadata.csv"))
  if (file.exists(expr_file) && file.exists(meta_file)) return(invisible(TRUE))
  eset <- GEOquery::getGEO(gse, destdir = data_dir, getGPL = FALSE)
  expr_set <- Biobase::exprs(eset[[1]])
  pdata <- Biobase::pData(eset[[1]])
  write.csv(expr_set, expr_file)
  write.csv(pdata, meta_file)
}

read_metadata <- function(gse, data_dir) {
  meta <- read.csv(file.path(data_dir, paste0(gse, "_metadata.csv")), check.names = FALSE)
  rownames(meta) <- meta[[1]]
  if (gse == "GSE39582") {
    # Original workflow used columns 23:24 for OS event/time.
    tmp <- meta[, c(23, 24), drop = FALSE]
    colnames(tmp) <- c("Event", "Time")
    clean <- data.frame(
      Event = sub("os.event: ", "", tmp$Event),
      Time = sub("os.delay \\(months\\): ", "", tmp$Time),
      row.names = rownames(tmp)
    )
  } else if (gse == "GSE17536") {
    # Original workflow used columns 15 and 18 for OS event/time.
    tmp <- meta[, c(15, 18), drop = FALSE]
    colnames(tmp) <- c("Event", "Time")
    clean <- data.frame(
      Event = sub("overall_event \\(death from any cause\\): ", "", tmp$Event),
      Time = sub("overall survival follow-up time: ", "", tmp$Time),
      row.names = rownames(tmp)
    )
    clean$Event <- ifelse(clean$Event == "no death", 0, 1)
  } else {
    stop("Unsupported cohort: ", gse)
  }
  clean <- clean %>% filter(across(everything(), ~ .x != "N/A"))
  clean$Time <- as.numeric(clean$Time)
  clean$Event <- as.numeric(clean$Event)
  clean <- clean[!is.na(clean$Time) & !is.na(clean$Event), , drop = FALSE]
  clean
}

map_probes_to_symbols <- function(expr_file, gpl_file) {
  expr <- data.table::fread(expr_file, data.table = FALSE)
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL
  gpl <- read.csv(gpl_file, check.names = FALSE)
  gpl <- gpl[, c(1, 11)]
  colnames(gpl) <- c("probe_id", "symbol")
  expr %>%
    rownames_to_column("probe_id") %>%
    inner_join(gpl, by = "probe_id") %>%
    filter(!is.na(symbol), symbol != "") %>%
    select(-probe_id) %>%
    select(symbol, everything()) %>%
    mutate(rowMean = rowMeans(across(-symbol), na.rm = TRUE)) %>%
    arrange(desc(rowMean)) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    select(-rowMean) %>%
    column_to_rownames("symbol") %>%
    as.matrix()
}

score_cohort <- function(gse) {
  if (args$download_if_missing) download_gse(gse, args$data_dir)
  meta <- read_metadata(gse, args$data_dir)
  expr <- map_probes_to_symbols(file.path(args$data_dir, paste0(gse, "_exprSet.csv")), args$gpl570)
  common <- intersect(colnames(expr), rownames(meta))
  expr <- expr[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
  gsva_scores <- GSVA::gsva(expr, marker_list, method = "ssgsea", parallel.sz = 1)
  if (requireNamespace("xCell", quietly = TRUE)) {
    xcell_scores <- xCell::xCellAnalysis(expr, cell.types.use = c("B-cells", "Plasma cells", "Tregs"))
    gsva_scores <- rbind(gsva_scores, xcell_scores[intersect(c("B-cells", "Plasma cells", "Tregs"), rownames(xcell_scores)), , drop = FALSE])
  } else {
    warning("xCell is not installed; xCell B/plasma/Treg scores are skipped for ", gse)
  }
  out <- cbind(meta, as.data.frame(t(gsva_scores), check.names = FALSE))
  write.csv(out, file.path(args$output_dir, paste0(gse, "_clinical_signature_scores.csv")), row.names = TRUE, quote = FALSE)
  saveRDS(gsva_scores, file.path(args$output_dir, paste0(gse, "_signature_score_matrix.rds")))
  out
}

run_km <- function(df, cohort, variable) {
  if (!variable %in% colnames(df)) return(NULL)
  dat <- df[, c("Time", "Event", variable), drop = FALSE]
  colnames(dat)[3] <- "score"
  dat <- dat[complete.cases(dat), ]
  if (nrow(dat) < 20 || length(unique(dat$Event)) < 2) return(NULL)
  cut <- tryCatch(survminer::surv_cutpoint(dat, time = "Time", event = "Event", variables = "score"), error = function(e) NULL)
  cutoff <- if (!is.null(cut)) summary(cut)$cutpoint[1, "cutpoint"] else median(dat$score, na.rm = TRUE)
  dat$Group <- ifelse(dat$score > cutoff, paste(variable, "High"), paste(variable, "Low"))
  dat$Group <- factor(dat$Group, levels = c(paste(variable, "Low"), paste(variable, "High")))
  fit <- survfit(Surv(Time, Event) ~ Group, data = dat)
  p <- ggsurvplot(
    fit, data = dat, title = cohort, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    xlab = "Months", ylab = "Overall survival", legend.title = variable,
    palette = c("#084583", "#F94040")
  )
  ggsave(file.path(args$output_dir, paste0(cohort, "_", make.names(variable), "_KM.pdf")), p$plot, width = 5.5, height = 5)
  write.csv(data.frame(cohort = cohort, variable = variable, cutoff = cutoff),
            file.path(args$output_dir, paste0(cohort, "_", make.names(variable), "_cutoff.csv")), row.names = FALSE)
  invisible(list(fit = fit, cutoff = cutoff))
}

cohorts <- c("GSE39582", "GSE17536")
score_tables <- lapply(cohorts, score_cohort)
names(score_tables) <- cohorts

for (cohort in cohorts) {
  df <- score_tables[[cohort]]
  vars <- intersect(c("dMMR_TLS_program", "Plasma_IgG", "B-cells", "Plasma cells"), colnames(df))
  invisible(lapply(vars, function(v) run_km(df, cohort, v)))
}

message("Done: survival scores and KM plots written to ", args$output_dir)
