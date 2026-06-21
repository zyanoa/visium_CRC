# ==============================================================================
# Script: misty_utilities.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   Utility functions for running MISTy on Visium Seurat objects
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(mistyR)
  library(furrr)
  library(R.utils)
  library(readr)
  library(stringr)
  library(tidyr)
  library(dplyr)
})

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

# ------------------------------------------------------------------------------
# 1. Run MISTy from a Seurat object
# ------------------------------------------------------------------------------

run_misty_seurat <- function(
    visium_slide,
    image_id,
    view_assays,
    view_features = NULL,
    view_types,
    view_params,
    spot_ids = NULL,
    out_alias = "results"
) {
  geometry <- GetTissueCoordinates(
    visium_slide,
    cols = c("row", "col"),
    scale = NULL,
    image = image_id
  )

  view_data <- purrr::map(
    view_assays,
    extract_seurat_data,
    geometry = geometry,
    visium_slide = visium_slide
  )

  build_misty_pipeline(
    view_data = view_data,
    view_features = view_features,
    view_types = view_types,
    view_params = view_params,
    geometry = geometry,
    spot_ids = spot_ids,
    out_alias = out_alias
  )
}

# ------------------------------------------------------------------------------
# 2. Extract assay data aligned to geometry
# ------------------------------------------------------------------------------

extract_seurat_data <- function(visium_slide, assay, geometry) {
  data <- GetAssayData(visium_slide, assay = assay) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()

  data <- data[match(rownames(geometry), rownames(data)), , drop = FALSE]
  return(data)
}

# ------------------------------------------------------------------------------
# 3. Filter features of interest
# ------------------------------------------------------------------------------

filter_data_features <- function(data, features = NULL) {
  if (is.null(features)) {
    features <- colnames(data)
  }

  data %>%
    rownames_to_column("rowname") %>%
    dplyr::select(rowname, all_of(features)) %>%
    rename_with(make.names) %>%
    column_to_rownames("rowname")
}

# ------------------------------------------------------------------------------
# 4. Build default MISTy views
# ------------------------------------------------------------------------------

create_default_views <- function(
    data,
    view_type,
    view_param,
    view_name,
    spot_ids,
    geometry
) {
  mistyR::clear_cache()

  view_data_init <- create_initial_view(data)

  if (!(view_type %in% c("intra", "para", "juxta"))) {
    view_type <- "intra"
  }

  if (view_type == "intra") {
    data_red <- view_data_init[["intraview"]]$data %>%
      rownames_to_column("rowname") %>%
      filter(rowname %in% spot_ids) %>%
      dplyr::select(-rowname)
  }

  if (view_type == "para") {
    view_data_tmp <- view_data_init %>%
      add_paraview(geometry, l = view_param)

    data_ix <- paste0("paraview.", view_param)
    data_red <- view_data_tmp[[data_ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot_ids) %>%
      dplyr::select(-rowname)
  }

  if (view_type == "juxta") {
    view_data_tmp <- view_data_init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view_param
      )

    data_ix <- paste0("juxtaview.", view_param)
    data_red <- view_data_tmp[[data_ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot_ids) %>%
      dplyr::select(-rowname)
  }

  if (is.null(view_param)) {
    misty_view <- create_view(view_name, data_red)
  } else {
    misty_view <- create_view(paste0(view_name, "_", view_param), data_red)
  }

  return(misty_view)
}

# ------------------------------------------------------------------------------
# 5. Build and run a MISTy pipeline
# ------------------------------------------------------------------------------

build_misty_pipeline <- function(
    view_data,
    view_features,
    view_types,
    view_params,
    geometry,
    spot_ids = NULL,
    out_alias = "default"
) {
  if (is.null(spot_ids)) {
    spot_ids <- rownames(view_data[[1]])
  }

  view_data_filt <- purrr::map2(view_data, view_features, filter_data_features)

  views_main <- create_initial_view(
    view_data_filt[[1]] %>%
      rownames_to_column("rowname") %>%
      filter(rowname %in% spot_ids) %>%
      dplyr::select(-rowname)
  )

  view_names <- names(view_data_filt)

  all_views <- purrr::pmap(
    list(
      view_data_filt[-1],
      view_types[-1],
      view_params[-1],
      view_names[-1]
    ),
    create_default_views,
    spot_ids = spot_ids,
    geometry = geometry
  )

  pipeline_views <- add_views(
    views_main,
    unlist(all_views, recursive = FALSE)
  )

  run_misty(pipeline_views, out_alias, cached = FALSE)
}

# ------------------------------------------------------------------------------
# 6. Collect MISTy results
# ------------------------------------------------------------------------------

collect_results_v2 <- function(folders) {
  samples <- R.utils::getAbsolutePath(folders)

  message("\nCollecting improvements")
  improvements <- samples %>%
    furrr::future_map_dfr(function(sample) {
      performance <- readr::read_table2(
        paste0(sample, .Platform$file.sep, "performance.txt"),
        na = c("", "NA", "NaN"),
        col_types = readr::cols()
      ) %>%
        dplyr::distinct()

      performance %>%
        dplyr::mutate(
          sample = sample,
          gain.RMSE = 100 * (.data$intra.RMSE - .data$multi.RMSE) / .data$intra.RMSE,
          gain.R2 = 100 * (.data$multi.R2 - .data$intra.R2)
        )
    }, .progress = TRUE) %>%
    tidyr::pivot_longer(
      -c(.data$sample, .data$target),
      names_to = "measure"
    )

  message("\nCollecting contributions")
  contributions <- samples %>%
    furrr::future_map_dfr(function(sample) {
      coefficients <- readr::read_table2(
        paste0(sample, .Platform$file.sep, "coefficients.txt"),
        na = c("", "NA", "NaN"),
        col_types = readr::cols()
      ) %>%
        dplyr::distinct()

      coefficients %>%
        dplyr::mutate(sample = sample, .after = "target") %>%
        tidyr::pivot_longer(
          cols = -c(.data$sample, .data$target),
          names_to = "view"
        )
    }, .progress = TRUE)

  improvements.stats <- improvements %>%
    dplyr::filter(!stringr::str_starts(.data$measure, "p\\.")) %>%
    dplyr::group_by(.data$target, .data$measure) %>%
    dplyr::summarise(
      mean = mean(.data$value),
      sd = stats::sd(.data$value),
      cv = .data$sd / .data$mean,
      .groups = "drop"
    )

  list(
    improvements = improvements,
    contributions = contributions,
    improvements.stats = improvements.stats
  )
}