#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(optparse)
})

option_list <- list(
  make_option(
    "--seurat_rds",
    type = "character",
    default = "data/processed/st_obj.rds",
    help = "Annotated Visium Seurat object with spatial images and metadata columns region/type/SectionID."
  ),
  make_option(
    "--auc_csv",
    type = "character",
    default = "results/06_regulon_analysis/Stroma1_4/output/SCENIC_regulon_AUC_Stroma1_4.csv",
    help = "Spot-by-regulon AUCell table from step 04."
  ),
  make_option(
    "--outdir",
    type = "character",
    default = "results/06_regulon_analysis/Stroma1_4/output/spatial"
  ),
  make_option(
    "--sections",
    type = "character",
    default = "PT9_2,PT9_3,PT34_1,PT55_2",
    help = "Comma-separated representative sections. Defaults to the sections used in the historical SCENIC figure workflow."
  ),
  make_option(
    "--max_quantile",
    type = "double",
    default = 0.99,
    help = "Upper plotting quantile for each regulon. Default: 0.99."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$seurat_rds)) stop("Seurat object not found: ", opt$seurat_rds)
if (!file.exists(opt$auc_csv)) stop("AUCell table not found: ", opt$auc_csv)

obj <- readRDS(opt$seurat_rds)
required_meta <- c("region", "type", "SectionID")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

auc_df <- read.csv(opt$auc_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (!"cell_id" %in% colnames(auc_df)) stop("AUCell table must contain cell_id.")
if (anyDuplicated(auc_df$cell_id)) stop("Duplicated cell_id values in AUCell table.")
rownames(auc_df) <- auc_df$cell_id
auc_df$cell_id <- NULL

b_gc_tfs <- c("BCL6", "IRF8", "PAX5", "SPIB")
plasma_tfs <- c("IRF4", "XBP1", "PRDM1")

resolve_regulon <- function(tf, available) {
  preferred <- c(
    paste0(tf, "(+)"),
    tf,
    paste0(tf, "_extended(+)"),
    paste0(tf, "_extended")
  )
  hit <- preferred[preferred %in% available]
  if (length(hit) > 0) return(hit[[1]])
  pattern <- paste0("^", tf, "(_extended)?(\\(\\+\\))?$")
  hit <- grep(pattern, available, value = TRUE)
  if (length(hit) > 0) return(hit[[1]])
  NA_character_
}

target_tfs <- c(b_gc_tfs, plasma_tfs)
regulon_map <- setNames(
  vapply(target_tfs, resolve_regulon, character(1), available = colnames(auc_df)),
  target_tfs
)

missing_tfs <- names(regulon_map)[is.na(regulon_map)]
if (length(missing_tfs) > 0) {
  warning("Target regulons not found and will be skipped: ", paste(missing_tfs, collapse = ", "))
}
if (all(is.na(regulon_map))) stop("None of the target regulons were found in the AUCell table.")

common_cells <- intersect(colnames(obj), rownames(auc_df))
if (length(common_cells) == 0) stop("No overlapping spots between Seurat object and AUCell table.")

# Add AUCell scores as metadata columns. Scores are NA outside the Stroma1/Stroma4
# SCENIC input subset, so the original Spatial assay is never overwritten.
meta_col_for_tf <- character(0)
for (tf in names(regulon_map)[!is.na(regulon_map)]) {
  regulon <- regulon_map[[tf]]
  safe_col <- paste0("SCENIC_", tf)
  obj@meta.data[[safe_col]] <- NA_real_
  obj@meta.data[common_cells, safe_col] <- as.numeric(auc_df[common_cells, regulon])
  meta_col_for_tf[[tf]] <- safe_col
}

sections <- trimws(strsplit(opt$sections, ",", fixed = TRUE)[[1]])
sections <- sections[nzchar(sections)]
available_sections <- unique(as.character(obj$SectionID))
missing_sections <- setdiff(sections, available_sections)
if (length(missing_sections) > 0) {
  warning("Requested sections not found and will be skipped: ", paste(missing_sections, collapse = ", "))
}
sections <- intersect(sections, available_sections)
if (length(sections) == 0) stop("None of the requested sections are present in the Seurat object.")

plot_one_tf <- function(tf) {
  if (!tf %in% names(meta_col_for_tf)) return(NULL)
  feature <- meta_col_for_tf[[tf]]
  vals <- obj@meta.data[[feature]]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NULL)

  upper <- as.numeric(stats::quantile(vals, probs = opt$max_quantile, na.rm = TRUE, names = FALSE))
  if (!is.finite(upper) || upper <= 0) upper <- max(vals, na.rm = TRUE)

  plots <- SpatialFeaturePlot(
    object = obj,
    images = sections,
    features = feature,
    pt.size.factor = 1.3,
    alpha = c(0.6, 1),
    min.cutoff = 0,
    max.cutoff = upper,
    keep.scale = "all",
    combine = FALSE
  )

  plots <- lapply(seq_along(plots), function(i) {
    plots[[i]] +
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(9, "Spectral")),
        limits = c(0, upper),
        na.value = "grey90"
      ) +
      labs(title = paste0(sections[[i]], " - ", tf, " regulon")) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "right"
      )
  })

  combined <- patchwork::wrap_plots(plots, ncol = length(sections))
  outfile <- file.path(opt$outdir, paste0("SCENIC_AUCell_", tf, "_spatial.pdf"))
  ggsave(outfile, combined, width = 5 * length(sections), height = 5)
  combined
}

bgc_plots <- lapply(b_gc_tfs, plot_one_tf)
bgc_plots <- Filter(Negate(is.null), bgc_plots)
if (length(bgc_plots) > 0) {
  combined_bgc <- patchwork::wrap_plots(bgc_plots, ncol = 1)
  ggsave(
    file.path(opt$outdir, "Fig3g_SCENIC_AUCell_B_GC_regulons.pdf"),
    combined_bgc,
    width = 5 * length(sections),
    height = 5 * length(bgc_plots)
  )
}

plasma_plots <- lapply(plasma_tfs, plot_one_tf)
plasma_plots <- Filter(Negate(is.null), plasma_plots)
if (length(plasma_plots) > 0) {
  combined_plasma <- patchwork::wrap_plots(plasma_plots, ncol = 1)
  ggsave(
    file.path(opt$outdir, "Fig3i_SCENIC_AUCell_plasma_regulons.pdf"),
    combined_plasma,
    width = 5 * length(sections),
    height = 5 * length(plasma_plots)
  )
}

write.table(
  data.frame(
    requested_tf = names(regulon_map),
    regulon = unname(regulon_map),
    metadata_column = ifelse(names(regulon_map) %in% names(meta_col_for_tf),
                             unname(meta_col_for_tf[names(regulon_map)]), NA_character_),
    stringsAsFactors = FALSE
  ),
  file.path(opt$outdir, "SCENIC_spatial_regulon_mapping.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(opt$outdir, "sessionInfo_SCENIC_spatial_Stroma1_4.txt")
)

message("SCENIC spatial plotting completed: ", opt$outdir)
