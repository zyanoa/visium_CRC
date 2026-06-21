# ==============================================================================
# Project-level configuration
# ==============================================================================
# Run scripts from the repository root whenever possible.
# You can override the default locations without editing scripts:
#   export CRC_MMR_ST_PROJECT=/path/to/CRC_MMR_spatial_TLS
#   export CRC_MMR_ST_DATA=/path/to/local/data
#   export CRC_MMR_ST_RAW=/path/to/raw/SpaceRanger_outputs
# ==============================================================================

PROJECT_DIR <- Sys.getenv("CRC_MMR_ST_PROJECT", unset = normalizePath(".", mustWork = FALSE))
DATA_DIR <- Sys.getenv("CRC_MMR_ST_DATA", unset = file.path(PROJECT_DIR, "data"))
RAW_DATA_DIR <- Sys.getenv("CRC_MMR_ST_RAW", unset = file.path(DATA_DIR, "raw"))
PROCESSED_DATA_DIR <- file.path(DATA_DIR, "processed")
EXTERNAL_DATA_DIR <- file.path(DATA_DIR, "external")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIGURES_DIR <- file.path(PROJECT_DIR, "figures")
METADATA_DIR <- file.path(PROJECT_DIR, "metadata")
SIGNATURE_DIR <- file.path(METADATA_DIR, "signatures")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

load_project_config_message <- function() {
  message("[config] PROJECT_DIR: ", PROJECT_DIR)
  message("[config] DATA_DIR: ", DATA_DIR)
  message("[config] RAW_DATA_DIR: ", RAW_DATA_DIR)
  message("[config] RESULTS_DIR: ", RESULTS_DIR)
}
