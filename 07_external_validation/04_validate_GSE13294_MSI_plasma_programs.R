# ==============================================================================
# Script: 04_validate_GSE13294_MSI_plasma_programs.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   External MSI/MSS validation of plasma-cell enrichment in GSE13294.
#   This script reproduces the plasma-cell validation panels corresponding to:
#     1) xCell total plasma-cell score
#     2) ssGSEA scores for B0_Plasma_IgA and B1_Plasma_IgG programs
#
# Notes:
#   - This script is kept separate from Figure 5 response/survival validation.
#   - It requires a symbol-level expression matrix and sample metadata containing
#     an MSI/MSS group column.
#   - Figure 3D uses xCell built-in Plasma cells score.
#   - Figure 3E uses custom B0_Plasma_IgA and B1_Plasma_IgG signatures
#     derived from B-cell/plasma-cell subcluster marker genes.
#   - Plasma-subtype signatures should be provided as a CSV file with at least
#     columns: signature,gene
# ============================================================================== 

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(reshape2)
  library(ggpubr)
  library(GSVA)
})

option_list <- list(
  make_option("--expr_csv", type = "character", default = "data/external/GSE13294/GSE13294_exprSet_symbol.csv",
              help = "Gene-symbol expression matrix CSV. Rows are genes and columns are samples."),
  make_option("--metadata_csv", type = "character", default = "data/external/GSE13294/GSE13294_metadata.csv",
              help = "Metadata CSV containing sample groups."),
  make_option("--group_col", type = "character", default = "Group",
              help = "Column in metadata indicating MSS/MSI status."),
  make_option("--mss_label", type = "character", default = "MSS",
              help = "Label for MSS samples."),
  make_option("--msi_label", type = "character", default = "MSI",
              help = "Label for MSI samples."),
  make_option("--plasma_signature_csv", type = "character", default = "signatures/plasma_subtype_signatures.csv",
              help = "CSV with columns signature,gene; must contain B0_Plasma_IgA and B1_Plasma_IgG."),
  make_option("--output_dir", type = "character", default = "results/07_external_validation/GSE13294_MSI_plasma",
              help = "Output directory."),
  make_option("--run_xcell", action = "store_true", default = FALSE,
              help = "Run xCell plasma-cell scoring if xCell is installed."),
  make_option("--dataset_label", type = "character", default = "GSE13294",
              help = "Dataset label used in figure titles.")
)

args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

read_matrix_csv <- function(path) {
  if (!file.exists(path)) stop("Expression file not found: ", path)
  expr <- read.csv(path, row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  expr
}

read_metadata <- function(path, group_col, mss_label, msi_label) {
  if (!file.exists(path)) stop("Metadata file not found: ", path)
  meta <- read.csv(path, row.names = 1, check.names = FALSE)
  if (!group_col %in% colnames(meta)) {
    stop("Group column not found in metadata: ", group_col)
  }
  meta[[group_col]] <- factor(meta[[group_col]], levels = c(mss_label, msi_label))
  meta <- meta[!is.na(meta[[group_col]]), , drop = FALSE]
  meta
}

read_signature_csv <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Plasma-subtype signature file not found: ", path,
      "\nProvide a CSV with columns signature,gene containing B0_Plasma_IgA and B1_Plasma_IgG."
    )
  }
  sig <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("signature", "gene")
  if (!all(required_cols %in% colnames(sig))) {
    stop("Signature CSV must contain columns: signature,gene")
  }
  sig <- sig %>%
    filter(!is.na(signature), !is.na(gene), signature != "", gene != "")
  split(sig$gene, sig$signature) %>% lapply(unique)
}

plot_group_box <- function(df, feature_col, group_col, ylab, title, outfile,
                           ylim = NULL, x_text = TRUE) {
  plot_df <- df %>%
    filter(!is.na(.data[[feature_col]]), !is.na(.data[[group_col]]))

  plot_df$Feature <- feature_col

  # For a single feature, draw one box per group at the same x-position.
  p <- ggplot(plot_df, aes(x = Feature, y = .data[[feature_col]], fill = .data[[group_col]])) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, outlier.alpha = 0) +
    stat_compare_means(
      aes(group = .data[[group_col]]),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE
    ) +
    scale_fill_manual(values = c("#e6194b", "#3cb44b"), name = "Group") +
    labs(x = NULL, y = ylab, title = title) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      axis.text.x = element_text(size = 14, angle = ifelse(x_text, 0, 90), hjust = 0.5)
    )

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  ggsave(outfile, p, width = 4.5, height = 4.5)
  p
}

plot_multi_signature_box <- function(long_df, group_col, title, outfile) {
  p <- ggplot(long_df, aes(x = Signature, y = Score, fill = .data[[group_col]])) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, outlier.alpha = 0) +
    stat_compare_means(
      aes(group = .data[[group_col]]),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 5
    ) +
    scale_fill_manual(values = c("#e6194b", "#3cb44b"), name = "Group") +
    labs(x = NULL, y = "ssGSEA score", title = title) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(outfile, p, width = 5, height = 5)
  p
}

# ------------------------------------------------------------------------------
# 1. Load input data
# ------------------------------------------------------------------------------

expr <- read_matrix_csv(args$expr_csv)
metadata <- read_metadata(args$metadata_csv, args$group_col, args$mss_label, args$msi_label)

common_samples <- intersect(colnames(expr), rownames(metadata))
if (length(common_samples) == 0) {
  stop("No overlapping samples between expression matrix columns and metadata row names.")
}
expr <- expr[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

write.csv(metadata, file.path(args$output_dir, "GSE13294_metadata_used.csv"), quote = FALSE)

# ------------------------------------------------------------------------------
# 2. xCell total plasma-cell score
# ------------------------------------------------------------------------------

if (isTRUE(args$run_xcell)) {
  if (!requireNamespace("xCell", quietly = TRUE)) {
    warning("xCell is not installed; skipping xCell analysis.")
  } else {
    res_xcell <- xCell::xCellAnalysis(expr, cell.types.use = c("Plasma cells"))
    saveRDS(res_xcell, file.path(args$output_dir, "GSE13294_xCell_plasma_cells.rds"))

    xcell_df <- as.data.frame(t(res_xcell))
    xcell_df[[args$group_col]] <- metadata[[args$group_col]]
    xcell_df$Sample <- rownames(xcell_df)
    write.csv(xcell_df, file.path(args$output_dir, "GSE13294_xCell_plasma_cells_scores.csv"), row.names = FALSE, quote = FALSE)

    plot_df <- xcell_df %>% rename(Plasma_cells = `Plasma cells`)
    p <- plot_group_box(
      df = plot_df,
      feature_col = "Plasma_cells",
      group_col = args$group_col,
      ylab = "Cell composition",
      title = args$dataset_label,
      outfile = file.path(args$output_dir, "GSE13294_xCell_Plasma_cells_MSS_MSI.pdf")
    )
    ggsave(file.path(args$output_dir, "GSE13294_xCell_Plasma_cells_MSS_MSI.png"), p, width = 4.5, height = 4.5, dpi = 300)
  }
}

# ------------------------------------------------------------------------------
# 3. ssGSEA for IgA+ and IgG+ plasma-cell programs
# ------------------------------------------------------------------------------

plasma_signatures <- read_signature_csv(args$plasma_signature_csv)
needed <- c("B0_Plasma_IgA", "B1_Plasma_IgG")
missing <- setdiff(needed, names(plasma_signatures))
if (length(missing) > 0) {
  stop("Missing required plasma-subtype signatures: ", paste(missing, collapse = ", "))
}
plasma_signatures <- plasma_signatures[needed]

ssgsea_scores <- GSVA::gsva(expr, plasma_signatures, method = "ssgsea")
saveRDS(ssgsea_scores, file.path(args$output_dir, "GSE13294_plasma_subtype_ssGSEA_scores.rds"))

score_df <- as.data.frame(t(ssgsea_scores))
score_df[[args$group_col]] <- metadata[[args$group_col]]
score_df$Sample <- rownames(score_df)
write.csv(score_df, file.path(args$output_dir, "GSE13294_plasma_subtype_ssGSEA_scores.csv"), row.names = FALSE, quote = FALSE)

long_df <- score_df %>%
  pivot_longer(cols = all_of(needed), names_to = "Signature", values_to = "Score") %>%
  mutate(Signature = factor(Signature, levels = needed))

p <- plot_multi_signature_box(
  long_df = long_df,
  group_col = args$group_col,
  title = args$dataset_label,
  outfile = file.path(args$output_dir, "GSE13294_ssGSEA_Plasma_IgA_IgG_MSS_MSI.pdf")
)
ggsave(file.path(args$output_dir, "GSE13294_ssGSEA_Plasma_IgA_IgG_MSS_MSI.png"), p, width = 5, height = 5, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(args$output_dir, "sessionInfo_GSE13294_MSI_plasma_validation.txt"))
