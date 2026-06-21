#!/usr/bin/env Rscript
# Add SCENIC regulon AUC values to the Seurat object and plot spatial PAX5 activity.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
})

option_list <- list(
  make_option("--seurat_rds", type = "character", help = "Input Seurat object."),
  make_option("--auc_csv", type = "character", default = "results/06_regulon_analysis/SCENIC_regulon_AUC_matrix.csv"),
  make_option("--output_dir", type = "character", default = "results/06_regulon_analysis/PAX5_spatial"),
  make_option("--sections", type = "character", default = "PT18_1,PT55_2", help = "Comma-separated section IDs."),
  make_option("--regulon", type = "character", default = "PAX5(+)"),
  make_option("--max_cutoff", type = "double", default = 0.25)
)
args <- parse_args(OptionParser(option_list = option_list))

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(args$seurat_rds)
auc <- read.csv(args$auc_csv, row.names = 1, check.names = FALSE)
common_cells <- intersect(colnames(obj), rownames(auc))
if (length(common_cells) == 0) stop("No matched cells/spots between Seurat object and AUC table.")

obj <- obj[, common_cells]
auc <- auc[common_cells, , drop = FALSE]
obj@meta.data <- cbind(obj@meta.data, auc[rownames(obj@meta.data), , drop = FALSE])

if (!args$regulon %in% colnames(obj@meta.data)) {
  stop(paste0("Regulon not found in metadata: ", args$regulon))
}

sections <- unlist(strsplit(args$sections, ","))
for (section in sections) {
  p <- SpatialFeaturePlot(
    object = obj,
    images = section,
    features = args$regulon,
    pt.size.factor = 1.3,
    alpha = c(0.6, 1),
    max.cutoff = args$max_cutoff,
    slot = "data",
    combine = TRUE
  ) +
    scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      limits = c(0, args$max_cutoff)
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    ggtitle(paste0(section, "_", args$regulon))

  safe_regulon <- gsub("[()+]", "", args$regulon)
  ggsave(file.path(args$output_dir, paste0(section, "_", safe_regulon, ".pdf")), p, width = 6, height = 5)
}
