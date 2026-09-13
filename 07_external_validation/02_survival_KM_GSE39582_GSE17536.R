#!/usr/bin/env Rscript

# ==============================================================================
# Project: CRC MMR spatial immune ecology
# Purpose: Final Figure 5e-h Kaplan-Meier analyses.
#
# Final panel mapping:
#   Fig5e = GSE39582 dMMR_TLS_signature
#   Fig5f = GSE17536 B_cells
#   Fig5g = GSE39582 Plasma_IgG
#   Fig5h = GSE17536 Plasma
#
# Key analysis decisions reproduced here:
#   - GSE39582 is restricted to primary colorectal adenocarcinoma samples with
#     complete overall-survival information (expected N = 562).
#   - GSE17536 uses samples with complete overall-survival information
#     (expected N = 177).
#   - Optimal cutpoints are obtained with survminer::surv_cutpoint().
#   - The exact full-precision returned cutpoint is used; there is NO fallback
#     to the median and NO manually rounded cutpoint.
#   - score > cutoff = High; score <= cutoff = Low.
# ==============================================================================

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

# ------------------------------------------------------------------------------
# 1. Signature definitions
# ------------------------------------------------------------------------------

load(args$marker_list)  # expected object: marker_list
if (!exists("marker_list") || !is.list(marker_list)) {
  stop("marker_list.Rdata must contain a list object named marker_list")
}

# Manuscript-facing final name. Keep the original source signature unchanged as
# an alias so the provenance remains explicit.
if (!"dMMR_TLS_signature" %in% names(marker_list)) {
  if ("Stroma4_1" %in% names(marker_list)) {
    marker_list[["dMMR_TLS_signature"]] <- marker_list[["Stroma4_1"]]
  } else if ("dMMR_TLS_program" %in% names(marker_list)) {
    marker_list[["dMMR_TLS_signature"]] <- marker_list[["dMMR_TLS_program"]]
  } else {
    stop("Cannot construct dMMR_TLS_signature: neither Stroma4_1 nor dMMR_TLS_program was found in marker_list")
  }
}

required_signatures <- c("dMMR_TLS_signature", "Plasma_IgG")
missing_signatures <- setdiff(required_signatures, names(marker_list))
if (length(missing_signatures) > 0) {
  stop("Missing required signatures: ", paste(missing_signatures, collapse = ", "))
}

# ------------------------------------------------------------------------------
# 2. GEO download / input helpers
# ------------------------------------------------------------------------------

download_gse <- function(gse, data_dir) {
  expr_file <- file.path(data_dir, paste0(gse, "_exprSet.csv"))
  meta_file <- file.path(data_dir, paste0(gse, "_metadata.csv"))

  if (file.exists(expr_file) && file.exists(meta_file)) return(invisible(TRUE))

  message("Downloading ", gse, " from GEO...")
  eset_list <- GEOquery::getGEO(gse, destdir = data_dir, getGPL = FALSE)
  if (length(eset_list) < 1) stop("No ExpressionSet returned for ", gse)

  eset <- eset_list[[1]]
  write.csv(Biobase::exprs(eset), expr_file)
  write.csv(Biobase::pData(eset), meta_file)
  invisible(TRUE)
}

read_raw_metadata <- function(gse, data_dir) {
  path <- file.path(data_dir, paste0(gse, "_metadata.csv"))
  meta <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  if (ncol(meta) < 2) stop("Metadata file has too few columns: ", path)

  rownames(meta) <- as.character(meta[[1]])
  meta[[1]] <- NULL
  meta
}

read_survival_metadata <- function(gse, data_dir) {
  meta <- read_raw_metadata(gse, data_dir)

  if (gse == "GSE39582") {
    if (!"source_name_ch1" %in% colnames(meta)) {
      stop("GSE39582 metadata lacks source_name_ch1; cannot apply the primary-tumor restriction")
    }

    is_tumor <- grepl(
      "^Frozen tissue of primary colorectal Adenocarcinoma",
      meta$source_name_ch1,
      ignore.case = TRUE
    )

    if (ncol(meta) < 24) stop("GSE39582 metadata does not contain columns 23:24 used for OS")

    tmp <- meta[, c(23, 24), drop = FALSE]
    colnames(tmp) <- c("Event", "Time")

    clean <- data.frame(
      Event = sub("os.event: ", "", tmp$Event),
      Time = sub("os.delay \\(months\\): ", "", tmp$Time),
      is_tumor = is_tumor,
      row.names = rownames(tmp),
      stringsAsFactors = FALSE
    )

    os_complete <- clean$Event != "N/A" & clean$Time != "N/A" &
      !is.na(clean$Event) & !is.na(clean$Time)

    clean <- clean[os_complete & clean$is_tumor, , drop = FALSE]
    clean$Event <- suppressWarnings(as.numeric(clean$Event))
    clean$Time <- suppressWarnings(as.numeric(clean$Time))
    clean <- clean[is.finite(clean$Event) & is.finite(clean$Time), , drop = FALSE]

    if (nrow(clean) != 562) {
      warning("Expected GSE39582 tumor + OS-complete N=562; observed N=", nrow(clean))
    }

  } else if (gse == "GSE17536") {
    if (ncol(meta) < 18) stop("GSE17536 metadata does not contain columns 15 and 18 used for OS")

    tmp <- meta[, c(15, 18), drop = FALSE]
    colnames(tmp) <- c("Event", "Time")

    clean <- data.frame(
      Event = sub("overall_event \\(death from any cause\\): ", "", tmp$Event),
      Time = sub("overall survival follow-up time: ", "", tmp$Time),
      row.names = rownames(tmp),
      stringsAsFactors = FALSE
    )

    clean <- clean[clean$Event != "N/A" & clean$Time != "N/A" &
                     !is.na(clean$Event) & !is.na(clean$Time), , drop = FALSE]
    clean$Event <- ifelse(clean$Event == "no death", 0, 1)
    clean$Time <- suppressWarnings(as.numeric(clean$Time))
    clean$Event <- suppressWarnings(as.numeric(clean$Event))
    clean <- clean[is.finite(clean$Event) & is.finite(clean$Time), , drop = FALSE]

    if (nrow(clean) != 177) {
      warning("Expected GSE17536 OS-complete N=177; observed N=", nrow(clean))
    }

  } else {
    stop("Unsupported cohort: ", gse)
  }

  clean
}

map_probes_to_symbols <- function(expr_file, gpl_file) {
  expr <- data.table::fread(expr_file, data.table = FALSE)
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL

  gpl <- read.csv(gpl_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(gpl) < 11) stop("GPL570 annotation file must contain at least 11 columns")

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

run_ssgsea <- function(expr, gene_sets) {
  # Keep compatibility with both the older GSVA API used in the original
  # analysis and the newer parameter-object API.
  tryCatch(
    GSVA::gsva(expr, gene_sets, method = "ssgsea", parallel.sz = 1),
    error = function(e_old) {
      if (!exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) stop(e_old)
      param <- GSVA::ssgseaParam(expr, gene_sets)
      GSVA::gsva(param, verbose = FALSE)
    }
  )
}

xcell_to_sample_table <- function(xc, sample_ids) {
  xc <- as.matrix(xc)

  if (all(sample_ids %in% colnames(xc))) {
    out <- as.data.frame(t(xc[, sample_ids, drop = FALSE]), check.names = FALSE)
  } else if (all(sample_ids %in% rownames(xc))) {
    out <- as.data.frame(xc[sample_ids, , drop = FALSE], check.names = FALSE)
  } else {
    stop(
      "xCell output could not be aligned to samples. Row matches=",
      sum(sample_ids %in% rownames(xc)), ", column matches=",
      sum(sample_ids %in% colnames(xc)), ", expected=", length(sample_ids)
    )
  }

  out[sample_ids, , drop = FALSE]
}

# ------------------------------------------------------------------------------
# 3. Score one survival cohort
# ------------------------------------------------------------------------------

score_cohort <- function(gse) {
  if (args$download_if_missing) download_gse(gse, args$data_dir)

  clinical <- read_survival_metadata(gse, args$data_dir)
  expr <- map_probes_to_symbols(
    file.path(args$data_dir, paste0(gse, "_exprSet.csv")),
    args$gpl570
  )

  common <- intersect(rownames(clinical), colnames(expr))
  if (length(common) < 20) stop("Too few expression/survival matched samples for ", gse)

  clinical <- clinical[common, , drop = FALSE]
  expr <- expr[, common, drop = FALSE]

  message(gse, ": scoring ", length(common), " survival-eligible samples")

  gsva_scores <- run_ssgsea(expr, marker_list)
  gsva_df <- as.data.frame(t(gsva_scores[, common, drop = FALSE]), check.names = FALSE)

  if (!requireNamespace("xCell", quietly = TRUE)) {
    stop("Package xCell is required for the final B-cell/plasma-cell survival panels")
  }

  xcell_raw <- xCell::xCellAnalysis(
    expr,
    cell.types.use = c("B-cells", "Plasma cells", "Tregs")
  )
  xcell_df <- xcell_to_sample_table(xcell_raw, common)

  out <- cbind(
    clinical[, intersect(c("Time", "Event"), colnames(clinical)), drop = FALSE],
    gsva_df[common, , drop = FALSE],
    xcell_df[common, , drop = FALSE]
  )

  # Final manuscript/display aliases. Raw xCell labels are retained as source
  # columns; only downstream aliases use underscores / simplified names.
  if (!"B-cells" %in% colnames(out)) stop("xCell output lacks raw row/column 'B-cells'")
  if (!"Plasma cells" %in% colnames(out)) stop("xCell output lacks raw row/column 'Plasma cells'")

  out$B_cells <- as.numeric(out[["B-cells"]])
  out$Plasma <- as.numeric(out[["Plasma cells"]])

  if (!"dMMR_TLS_signature" %in% colnames(out)) {
    stop("ssGSEA output lacks dMMR_TLS_signature")
  }
  if (!"Plasma_IgG" %in% colnames(out)) {
    stop("ssGSEA output lacks Plasma_IgG")
  }

  # Compatibility alias for the existing Figure 5i Cox script. This does not
  # alter values and can be removed once that script is renamed consistently.
  out$dMMR_TLS_program <- out$dMMR_TLS_signature

  write.csv(
    out,
    file.path(args$output_dir, paste0(gse, "_clinical_signature_scores.csv")),
    row.names = TRUE,
    quote = FALSE
  )

  saveRDS(gsva_scores, file.path(args$output_dir, paste0(gse, "_signature_score_matrix.rds")))
  saveRDS(xcell_raw, file.path(args$output_dir, paste0(gse, "_xCell_score_matrix.rds")))

  out
}

# ------------------------------------------------------------------------------
# 4. Exact surv_cutpoint + KM helper
# ------------------------------------------------------------------------------

get_exact_cutpoint <- function(dat, variable) {
  cut <- survminer::surv_cutpoint(
    dat,
    time = "Time",
    event = "Event",
    variables = variable
  )

  cut_tbl <- summary(cut)
  if (!"cutpoint" %in% colnames(cut_tbl)) {
    stop("surv_cutpoint did not return a cutpoint for ", variable)
  }

  cutoff <- as.numeric(cut_tbl$cutpoint[1])
  if (!is.finite(cutoff)) stop("Non-finite cutpoint for ", variable)
  cutoff
}

run_km <- function(df, cohort, variable, panel, legend_title,
                   high_label = paste(variable, "High"),
                   low_label = paste(variable, "Low")) {
  if (!variable %in% colnames(df)) stop("Missing score column: ", variable, " in ", cohort)

  dat <- df[, c("Time", "Event", variable), drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat[[variable]] <- as.numeric(dat[[variable]])

  if (nrow(dat) < 20) stop("Too few complete samples for ", cohort, " / ", variable)
  if (length(unique(dat$Event)) < 2) stop("Survival event has <2 levels for ", cohort)

  cutoff <- get_exact_cutpoint(dat, variable)

  dat$Group <- ifelse(dat[[variable]] > cutoff, high_label, low_label)
  dat$Group <- factor(dat$Group, levels = c(high_label, low_label))

  fit <- survival::survfit(Surv(Time, Event) ~ Group, data = dat)
  lr <- survival::survdiff(Surv(Time, Event) ~ Group, data = dat)
  logrank_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  p <- survminer::ggsurvplot(
    fit,
    data = dat,
    title = cohort,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    surv.median.line = "hv",
    xlab = "Months",
    ylab = "Overall survival",
    legend.title = legend_title,
    legend.labs = c("High", "Low"),
    palette = c("#F94040", "#084583")
  )

  outfile <- file.path(args$output_dir, paste0(panel, "_", cohort, "_", variable, ".pdf"))
  pdf(outfile, width = 6, height = 6)
  print(p)
  dev.off()

  sample_table <- dat
  sample_table$sample_id <- rownames(dat)
  write.table(
    sample_table,
    file.path(args$output_dir, paste0(panel, "_sample_groups.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  stats <- data.frame(
    panel = panel,
    cohort = cohort,
    variable = variable,
    N = nrow(dat),
    cutoff = cutoff,
    High_N = sum(dat$Group == high_label),
    Low_N = sum(dat$Group == low_label),
    logrank_P = logrank_p,
    stringsAsFactors = FALSE
  )

  write.csv(
    stats,
    file.path(args$output_dir, paste0(panel, "_", cohort, "_", variable, "_stats.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  message(
    panel, " ", cohort, " ", variable,
    ": N=", nrow(dat),
    ", cutoff=", format(cutoff, digits = 17),
    ", High=", stats$High_N,
    ", Low=", stats$Low_N,
    ", log-rank P=", signif(logrank_p, 6)
  )

  invisible(list(fit = fit, cutoff = cutoff, stats = stats, data = dat))
}

# ------------------------------------------------------------------------------
# 5. Run FINAL Figure 5e-h mapping only
# ------------------------------------------------------------------------------

g39582 <- score_cohort("GSE39582")
g17536 <- score_cohort("GSE17536")

results <- list()

results[["Fig5e"]] <- run_km(
  df = g39582,
  cohort = "GSE39582",
  variable = "dMMR_TLS_signature",
  panel = "Fig5e",
  legend_title = "dMMR_TLS_signature",
  high_label = "dMMR_TLS_signature High",
  low_label = "dMMR_TLS_signature Low"
)

results[["Fig5f"]] <- run_km(
  df = g17536,
  cohort = "GSE17536",
  variable = "B_cells",
  panel = "Fig5f",
  legend_title = "B_cells",
  high_label = "B_cells High",
  low_label = "B_cells Low"
)

results[["Fig5g"]] <- run_km(
  df = g39582,
  cohort = "GSE39582",
  variable = "Plasma_IgG",
  panel = "Fig5g",
  legend_title = "Plasma_IgG",
  high_label = "Plasma_IgG High",
  low_label = "Plasma_IgG Low"
)

results[["Fig5h"]] <- run_km(
  df = g17536,
  cohort = "GSE17536",
  variable = "Plasma",
  panel = "Fig5h",
  legend_title = "Plasma",
  high_label = "Plasma High",
  low_label = "Plasma Low"
)

summary_table <- bind_rows(lapply(results, `[[`, "stats"))
write.csv(
  summary_table,
  file.path(args$output_dir, "Figure5e_h_KM_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(args$output_dir, "sessionInfo_Figure5e_h_survival.txt")
)

message("Done: final Figure 5e-h KM workflow written to ", args$output_dir)
