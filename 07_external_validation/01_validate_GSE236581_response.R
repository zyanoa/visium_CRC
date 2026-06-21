#!/usr/bin/env Rscript

# Project: CRC MMR spatial immune ecology
# Purpose: Validate dMMR TLS-program/plasma-cell-related signatures in the GSE236581 immunotherapy response cohort.
# Output: Boxplots and statistics for B-cells, Plasma cells, dMMR_TLS_program and Plasma_IgG scores.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(rstatix)
  library(ggpubr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    seurat_rds = NA_character_,
    marker_list = "signatures/marker_list.Rdata",
    output_dir = "results/07_external_validation/GSE236581_response",
    assay = "RNA",
    slot = "counts"
  )
  for (arg in args) {
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2 && kv[1] %in% names(out)) out[[kv[1]]] <- kv[2]
  }
  out
}

args <- parse_args()
if (is.na(args$seurat_rds)) stop("Please provide --seurat_rds=/path/to/GSE236581_NT.rds")
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

load(args$marker_list)  # expected object: marker_list
if (!exists("marker_list") || !is.list(marker_list)) stop("marker_list.Rdata must contain a list object named marker_list")
# Keep marker_list as the final signature source, but remove generic TLS entry
# and rename the internal Stroma4_1 score to a manuscript-facing name.
if ("TLS" %in% names(marker_list)) marker_list[["TLS"]] <- NULL
if ("Stroma4_1" %in% names(marker_list)) names(marker_list)[names(marker_list) == "Stroma4_1"] <- "dMMR_TLS_program"


obj <- readRDS(args$seurat_rds)

# Patient groups used in the original analysis.
cr_patients <- c("01", "04", "08", "09", "11", "12", "18", "20", "21", "22", "23", "26")
pr_patients <- c("03", "05", "14", "15", "17", "19", "24")
sd_patients <- c("02", "16", "25")
dmmr_patients <- c("01", "04", "08", "11", "12", "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26")
pmmr_patients <- c("02", "03", "05", "09", "14", "16")

assign_group <- function(ids, groups) {
  out <- rep(NA_character_, length(ids))
  for (nm in names(groups)) {
    pattern <- paste(groups[[nm]], collapse = "|")
    out[grepl(pattern, ids)] <- nm
  }
  out
}

cell_ids <- colnames(obj)
response <- assign_group(cell_ids, list(CR = cr_patients, PR = pr_patients, SD = sd_patients))
mmr_status <- assign_group(cell_ids, list(dMMR = dmmr_patients, pMMR = pmmr_patients))

obj$Response <- factor(response, levels = c("CR", "PR", "SD"), ordered = TRUE)
obj$MMR_Status <- factor(mmr_status, levels = c("dMMR", "pMMR"))

# xCell scores for B cells and plasma cells.
xcell_scores <- NULL
if (requireNamespace("xCell", quietly = TRUE)) {
  counts_mat <- Seurat::GetAssayData(obj, assay = args$assay, slot = args$slot)
  xcell_scores <- xCell::xCellAnalysis(counts_mat, cell.types.use = c("B-cells", "Plasma cells", "Tregs"))
  saveRDS(xcell_scores, file.path(args$output_dir, "GSE236581_xCell_scores.rds"))
} else {
  warning("Package xCell is not installed; B-cells and Plasma cells xCell scores will be skipped.")
}

# UCell scores for marker_list signatures.
if (!requireNamespace("UCell", quietly = TRUE)) stop("Package UCell is required for AddModuleScore_UCell")
obj <- UCell::AddModuleScore_UCell(obj, features = marker_list, name = NULL, assay = args$assay, slot = args$slot)
saveRDS(obj, file.path(args$output_dir, "GSE236581_with_UCell_scores.rds"))

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

score_df <- data.frame(cell_id = rownames(meta), Response = meta$Response, MMR_Status = meta$MMR_Status, check.names = FALSE)

# Add UCell scores. UCell versions may create either exact names or *_UCell columns.
for (sig in names(marker_list)) {
  candidates <- c(sig, paste0(sig, "_UCell"), paste0("UCell_", sig))
  col <- candidates[candidates %in% colnames(meta)][1]
  if (!is.na(col)) score_df[[sig]] <- meta[[col]]
}

if (!is.null(xcell_scores)) {
  xcell_t <- as.data.frame(t(xcell_scores), check.names = FALSE)
  xcell_t$cell_id <- rownames(xcell_t)
  score_df <- left_join(score_df, xcell_t[, intersect(c("cell_id", "B-cells", "Plasma cells", "Tregs"), colnames(xcell_t)), drop = FALSE], by = "cell_id")
}

write.csv(score_df, file.path(args$output_dir, "GSE236581_response_signature_scores.csv"), row.names = FALSE, quote = FALSE)

plot_one_score <- function(df, variable, y_label = "Score") {
  if (!variable %in% colnames(df)) return(NULL)
  plot_df <- df %>% filter(!is.na(Response), !is.na(.data[[variable]]))
  if (nrow(plot_df) == 0) return(NULL)
  stat <- plot_df %>%
    wilcox_test(as.formula(paste0("`", variable, "` ~ Response"))) %>%
    add_significance("p") %>%
    add_xy_position(x = "Response", dodge = 0.8)
  write.csv(stat, file.path(args$output_dir, paste0("stats_", make.names(variable), ".csv")), row.names = FALSE)
  p <- ggplot(plot_df, aes(x = Response, y = .data[[variable]], fill = Response)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    stat_pvalue_manual(stat, label = "p.signif", tip.length = 0.02, bracket.size = 0.5, step.increase = 0.08) +
    theme_classic(base_size = 14) +
    labs(title = "GSE236581", x = variable, y = y_label) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top")
  ggsave(file.path(args$output_dir, paste0("GSE236581_", make.names(variable), "_response_boxplot.pdf")), p, width = 4.5, height = 4.5)
  p
}

vars_to_plot <- intersect(c("dMMR_TLS_program", "Plasma_IgG", "B-cells", "Plasma cells"), colnames(score_df))
invisible(lapply(vars_to_plot, function(v) plot_one_score(score_df, v)))

message("Done: GSE236581 response validation results written to ", args$output_dir)
