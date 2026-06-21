# ==============================================================================
# Script: 02_MIA_run.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Load precomputed spatial-region and scRNA cell-type marker tables
#   2. Run multimodal intersection analysis (MIA)
#   3. Export enrichment results and heatmaps
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

# ------------------------------------------------------------------------------
# 1. Input and output
# ------------------------------------------------------------------------------

input_region_markers <- file.path(RESULTS_DIR, "02_celltype_annotation", "MIA", "region_specific_for_MIA.csv")
input_celltype_markers <- file.path(RESULTS_DIR, "02_celltype_annotation", "MIA", "celltype_specific_for_MIA.csv")
mia_function_file <- file.path("scripts", "00_utils", "MIA_function.R")

outdir <- file.path(RESULTS_DIR, "02_celltype_annotation", "MIA")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

sample_name <- "dMMR_pMMR_all_celltypes"

# ------------------------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------------------------

region_specific <- read.csv(input_region_markers, stringsAsFactors = FALSE)
celltype_specific <- read.csv(input_celltype_markers, stringsAsFactors = FALSE)

stopifnot(all(c("region", "gene") %in% colnames(region_specific)))
stopifnot(all(c("celltype", "gene") %in% colnames(celltype_specific)))

# ------------------------------------------------------------------------------
# 3. Load MIA function
# ------------------------------------------------------------------------------

source(mia_function_file)

# ------------------------------------------------------------------------------
# 4. Run MIA
# ------------------------------------------------------------------------------

mia_result <- zhao_MIA(
  sp.diff = region_specific,
  sc.diff = celltype_specific,
  sample = sample_name,
  outdir = outdir
)

# ------------------------------------------------------------------------------
# 5. Save serialized result object if returned
# ------------------------------------------------------------------------------

saveRDS(
  mia_result,
  file = file.path(outdir, paste0(sample_name, "_MIA_result.rds"))
)

# ------------------------------------------------------------------------------
# 6. Save session information
# ------------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo_02_MIA_run.txt")
)

message("MIA analysis completed successfully.")
message("Saved outputs to: ", outdir)