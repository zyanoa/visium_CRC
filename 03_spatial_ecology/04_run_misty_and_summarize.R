# ==============================================================================
# Script: 04_misty_run_and_summary.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Merge cell2location abundance matrices
#   2. Add cell2location abundance as a new Seurat assay
#   3. Run MISTy for each Visium section
#   4. Summarize MISTy improvements and interactions across sections
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(mistyR)
})

source("00_utils/misty_utilities.R")

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

# ------------------------------------------------------------------------------
# 1. Input and output
# ------------------------------------------------------------------------------

input_rds <- file.path(PROCESSED_DATA_DIR, "st_obj.rds")
cell2location_dir <- file.path(RESULTS_DIR, "02_celltype_annotation", "cell2location")
outdir <- file.path(RESULTS_DIR, "03_spatial_ecology", "misty")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

st_obj <- readRDS(input_rds)

# ------------------------------------------------------------------------------
# 2. Helper function to run colocalization
# ------------------------------------------------------------------------------

run_colocalization <- function(slide_obj, image_id, assay, useful_features, out_label, misty_out_alias) {
  view_assays <- list(
    main = assay,
    juxta = assay,
    para = assay
  )

  view_features <- list(
    main = useful_features,
    juxta = useful_features,
    para = useful_features
  )

  view_types <- list(
    main = "intra",
    juxta = "juxta",
    para = "para"
  )

  view_params <- list(
    main = NULL,
    juxta = 2,
    para = 5
  )

  misty_out <- file.path(misty_out_alias, paste0(out_label, "_", assay))

  run_misty_seurat(
    visium_slide = slide_obj,
    image_id = image_id,
    view_assays = view_assays,
    view_features = view_features,
    view_types = view_types,
    view_params = view_params,
    spot_ids = NULL,
    out_alias = misty_out
  )

  return(misty_out)
}

# ------------------------------------------------------------------------------
# 3. Read and merge cell2location abundance matrices
# ------------------------------------------------------------------------------

sample_dirs <- list.dirs(cell2location_dir, full.names = TRUE, recursive = FALSE)
abundance_list <- list()

for (sub_dir in sample_dirs) {
  sample_name <- basename(sub_dir)
  abundance_file <- file.path(sub_dir, paste0(sample_name, ".spatial.deconvolution.csv"))

  if (file.exists(abundance_file)) {
    tmp_df <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
    abundance_list[[sample_name]] <- tmp_df
  }
}

if (length(abundance_list) == 0) {
  stop("No cell2location abundance files found.")
}

common_columns <- Reduce(intersect, lapply(abundance_list, colnames))
merged_abundance <- do.call(rbind, lapply(abundance_list, `[`, common_columns))
rownames(merged_abundance) <- gsub("\\.", "_", rownames(merged_abundance))
merged_abundance <- merged_abundance[colnames(st_obj), common_columns, drop = FALSE]

# Preserve original cell type label if needed
if ("Tgd" %in% colnames(merged_abundance)) {
  colnames(merged_abundance)[colnames(merged_abundance) == "Tgd"] <- "γδ_T"
}

st_obj[["cell2location"]] <- CreateAssayObject(counts = t(merged_abundance))

# ------------------------------------------------------------------------------
# 4. Run MISTy for each section
# ------------------------------------------------------------------------------

section_ids <- c(
  "PT18_1", "PT18_2", "PT18_3",
  "PT32_1", "PT32_2", "PT32_3",
  "PT34_1", "PT34_2",
  "PT55_1", "PT55_2", "PT55_3",
  "PT57_1", "PT57_2", "PT57_3", "PT57_4",
  "PT9_1", "PT9_2", "PT9_3"
)

for (slide_id in section_ids) {
  slide_obj <- subset(st_obj, subset = SectionID == slide_id)
  DefaultAssay(slide_obj) <- "cell2location"

  sample_outdir <- file.path(outdir, slide_id)
  if (!dir.exists(sample_outdir)) dir.create(sample_outdir, recursive = TRUE)

  useful_features <- rownames(slide_obj[["cell2location"]])

  misty_out <- run_colocalization(
    slide_obj = slide_obj,
    image_id = slide_id,
    assay = "cell2location",
    useful_features = useful_features,
    out_label = slide_id,
    misty_out_alias = sample_outdir
  )

  misty_res_slide <- collect_results(misty_out)

  plot_folder <- paste0(misty_out, "_plots")
  if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)

  views_and_cutoffs <- list(
    list(view = "intra", cutoff = 0),
    list(view = "intra", cutoff = 0.5),
    list(view = "juxta_2", cutoff = 0),
    list(view = "juxta_2", cutoff = 0.5),
    list(view = "para_5", cutoff = 0),
    list(view = "para_5", cutoff = 0.5)
  )

  for (vc in views_and_cutoffs) {
    pdf(file = file.path(plot_folder, paste0(slide_id, "_", vc$view, "_improvement_stats.pdf")))
    mistyR::plot_improvement_stats(misty_res_slide)
    dev.off()

    pdf(file = file.path(plot_folder, paste0(slide_id, "_", vc$view, "_view_contributions.pdf")))
    mistyR::plot_view_contributions(misty_res_slide)
    dev.off()

    pdf(file = file.path(plot_folder, paste0(slide_id, "_", vc$view, "_interaction_heatmap_cutoff", vc$cutoff, ".pdf")))
    mistyR::plot_interaction_heatmap(misty_res_slide, vc$view, cutoff = vc$cutoff)
    dev.off()

    pdf(file = file.path(plot_folder, paste0(slide_id, "_", vc$view, "_interaction_communities_cutoff", vc$cutoff, ".pdf")))
    mistyR::plot_interaction_communities(misty_res_slide, vc$view, cutoff = vc$cutoff)
    dev.off()
  }
}

# ------------------------------------------------------------------------------
# 5. Aggregate MISTy results across sections
# ------------------------------------------------------------------------------

misty_runs <- list.files(outdir, full.names = TRUE)
misty_runs <- misty_runs[file.info(misty_runs)$isdir]

misty_result_dirs <- unlist(lapply(misty_runs, function(x) {
  list.files(x, full.names = TRUE)
}))
misty_result_dirs <- misty_result_dirs[file.info(misty_result_dirs)$isdir]

misty_res <- collect_results(misty_result_dirs)

saveRDS(
  misty_res,
  file = file.path(outdir, "misty_results_aggregated.rds")
)

if (!is.null(misty_res$improvements)) {
  write.csv(
    misty_res$improvements,
    file = file.path(outdir, "misty_improvements_all_sections.csv"),
    row.names = FALSE
  )
}

if (!is.null(misty_res$contributions)) {
  write.csv(
    misty_res$contributions,
    file = file.path(outdir, "misty_contributions_all_sections.csv"),
    row.names = FALSE
  )
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo_04_misty_run_and_summary.txt")
)