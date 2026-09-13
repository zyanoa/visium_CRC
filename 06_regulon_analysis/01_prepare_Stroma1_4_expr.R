#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(optparse)
})

option_list <- list(
  make_option(
    "--seurat_rds",
    type = "character",
    default = "data/processed/st_obj.rds",
    help = "Annotated Visium Seurat object containing metadata columns 'region' and 'type'."
  ),
  make_option(
    "--outdir",
    type = "character",
    default = "results/06_regulon_analysis/Stroma1_4/input",
    help = "Output directory for the SCENIC input matrix and metadata."
  ),
  make_option(
    "--assay",
    type = "character",
    default = "Spatial",
    help = "Assay containing raw counts. Default: Spatial."
  ),
  make_option(
    "--min_detected_fraction",
    type = "double",
    default = 0.01,
    help = "Minimum fraction of selected spots with non-zero expression. Default: 0.01."
  ),
  make_option(
    "--min_spots_floor",
    type = "integer",
    default = 20,
    help = "Minimum absolute number of selected spots with non-zero expression. Default: 20."
  ),
  make_option(
    "--min_counts",
    type = "double",
    default = 50,
    help = "Minimum total counts across selected spots. Default: 50."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$seurat_rds)) {
  stop("Seurat object not found: ", opt$seurat_rds)
}

message("Loading Seurat object: ", opt$seurat_rds)
obj <- readRDS(opt$seurat_rds)

required_meta <- c("region", "type")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

if (!opt$assay %in% names(obj@assays)) {
  stop("Assay not found: ", opt$assay)
}

target_regions <- c("Stroma1_IGLC1", "Stroma4_IGHM")
target_types <- c("dMMR", "pMMR")

type_chr <- as.character(obj@meta.data$type)
region_chr <- as.character(obj@meta.data$region)

cells_use <- rownames(obj@meta.data)[
  region_chr %in% target_regions & type_chr %in% target_types
]

if (length(cells_use) == 0) {
  stop("No spots found for Stroma1_IGLC1/Stroma4_IGHM in dMMR/pMMR samples.")
}

meta <- obj@meta.data[cells_use, , drop = FALSE]
meta$type_region <- paste0(as.character(meta$type), "_", as.character(meta$region))

group_levels <- c(
  "dMMR_Stroma1_IGLC1",
  "pMMR_Stroma1_IGLC1",
  "dMMR_Stroma4_IGHM",
  "pMMR_Stroma4_IGHM"
)
meta$type_region <- factor(meta$type_region, levels = group_levels)

if (anyNA(meta$type_region)) {
  bad <- unique(paste(meta$type, meta$region, sep = "_"))[is.na(unique(meta$type_region))]
  stop("Unexpected type/region combinations were selected. Check metadata values.")
}

message("Selected spot counts:")
print(table(meta$type_region, useNA = "ifany"))

# Use raw counts, matching the original Stroma1+Stroma4 pySCENIC workflow.
# The tryCatch keeps the script compatible with both Seurat v5 (layer) and
# older Seurat objects (slot).
expr_all <- tryCatch(
  Seurat::GetAssayData(obj, assay = opt$assay, layer = "counts"),
  error = function(e) Seurat::GetAssayData(obj, assay = opt$assay, slot = "counts")
)

missing_cells <- setdiff(cells_use, colnames(expr_all))
if (length(missing_cells) > 0) {
  stop("Selected metadata spots are missing from the count matrix. Example: ",
       paste(head(missing_cells, 10), collapse = ", "))
}

expr <- expr_all[, cells_use, drop = FALSE]

if (anyDuplicated(rownames(expr))) {
  dup <- unique(rownames(expr)[duplicated(rownames(expr))])
  stop("Duplicated gene names are not supported by this workflow. Examples: ",
       paste(head(dup, 10), collapse = ", "))
}

message("Raw expression matrix: ", nrow(expr), " genes x ", ncol(expr), " spots")

gene_detected <- Matrix::rowSums(expr > 0)
gene_counts <- Matrix::rowSums(expr)

min_spots <- max(
  as.integer(opt$min_spots_floor),
  ceiling(ncol(expr) * opt$min_detected_fraction)
)

keep_genes <- gene_detected >= min_spots & gene_counts >= opt$min_counts
expr_filt <- expr[keep_genes, , drop = FALSE]

message(
  "Gene filter: detected in >= ", min_spots,
  " spots and total counts >= ", opt$min_counts
)
message("Filtered expression matrix: ", nrow(expr_filt), " genes x ", ncol(expr_filt), " spots")

if (nrow(expr_filt) < 1000) {
  warning("Fewer than 1,000 genes remain after filtering; inspect the input and thresholds.")
}

expr_file <- file.path(opt$outdir, "Stroma1_4_exprMat.csv")
meta_file <- file.path(opt$outdir, "Stroma1_4_metadata.csv")
summary_file <- file.path(opt$outdir, "Stroma1_4_input_summary.tsv")

# The historical analysis used a gene-by-spot CSV before conversion to loom.
# This preserves that provenance for the public reproduction workflow.
expr_df <- as.data.frame(as.matrix(expr_filt), check.names = FALSE)
expr_df <- cbind(Gene = rownames(expr_df), expr_df)
write.csv(expr_df, expr_file, row.names = FALSE, quote = FALSE)
rm(expr_df)

gc()

meta_out <- meta
meta_out$cell_id <- rownames(meta_out)
write.csv(meta_out, meta_file, row.names = FALSE, quote = FALSE)

summary_df <- data.frame(
  metric = c(
    "selected_regions",
    "selected_types",
    "selected_spots",
    "raw_genes",
    "filtered_genes",
    "min_detected_fraction",
    "min_spots",
    "min_total_counts"
  ),
  value = c(
    paste(target_regions, collapse = ";"),
    paste(target_types, collapse = ";"),
    ncol(expr),
    nrow(expr),
    nrow(expr_filt),
    opt$min_detected_fraction,
    min_spots,
    opt$min_counts
  ),
  stringsAsFactors = FALSE
)
write.table(summary_df, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(opt$outdir, "sessionInfo_prepare_Stroma1_4.txt")
)

message("SCENIC input preparation completed.")
message("Expression: ", expr_file)
message("Metadata:   ", meta_file)
message("Summary:    ", summary_file)
