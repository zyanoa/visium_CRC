#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
})

option_list <- list(
  make_option("--loom", type = "character", help = "pySCENIC AUCell loom: sce_SCENIC.loom"),
  make_option("--metadata", type = "character", help = "Stroma1_4_metadata.csv from step 01."),
  make_option("--outdir", type = "character", default = "results/06_regulon_analysis/Stroma1_4/output")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$loom) || !file.exists(opt$loom)) stop("Missing --loom input.")
if (is.null(opt$metadata) || !file.exists(opt$metadata)) stop("Missing --metadata input.")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

required_groups <- c(
  "dMMR_Stroma1_IGLC1",
  "pMMR_Stroma1_IGLC1",
  "dMMR_Stroma4_IGHM",
  "pMMR_Stroma4_IGHM"
)

b_gc_tfs <- c("BCL6", "IRF8", "PAX5", "SPIB")
plasma_tfs <- c("IRF4", "XBP1", "PRDM1")

meta <- read.csv(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("cell_id", "type_region") %in% colnames(meta))) {
  stop("Metadata must contain columns: cell_id, type_region")
}
if (anyDuplicated(meta$cell_id)) stop("Duplicated cell_id values in metadata.")

meta <- meta[meta$type_region %in% required_groups, , drop = FALSE]
missing_groups <- setdiff(required_groups, unique(meta$type_region))
if (length(missing_groups) > 0) {
  stop("Missing required type_region groups: ", paste(missing_groups, collapse = ", "))
}

sce <- open_loom(opt$loom)
on.exit(close_loom(sce), add = TRUE)

regulon_auc <- get_regulons_AUC(sce, column.attr.name = "RegulonsAUC")
auc <- getAUC(regulon_auc)

common_cells <- intersect(colnames(auc), meta$cell_id)
if (length(common_cells) == 0) {
  stop("No overlapping spots between SCENIC AUCell output and metadata.")
}

auc <- auc[, common_cells, drop = FALSE]
meta <- meta[match(common_cells, meta$cell_id), , drop = FALSE]
stopifnot(identical(meta$cell_id, colnames(auc)))

cell_annotation <- factor(meta$type_region, levels = required_groups)
rss <- calcRSS(AUC = auc, cellAnnotation = cell_annotation)
rss <- as.data.frame(rss, check.names = FALSE)
rss <- rss[, required_groups, drop = FALSE]
rss <- rss[complete.cases(rss), , drop = FALSE]

write.csv(
  rss,
  file.path(opt$outdir, "SCENIC_RSS_all_regulons_Stroma1_4.csv"),
  quote = FALSE
)

auc_out <- as.data.frame(t(auc), check.names = FALSE)
auc_out <- cbind(cell_id = rownames(auc_out), auc_out)
write.csv(
  auc_out,
  file.path(opt$outdir, "SCENIC_regulon_AUC_Stroma1_4.csv"),
  row.names = FALSE,
  quote = FALSE
)

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

requested_tfs <- c(b_gc_tfs, plasma_tfs)
resolution <- data.frame(
  requested_tf = requested_tfs,
  regulon = vapply(requested_tfs, resolve_regulon, character(1), available = rownames(rss)),
  stringsAsFactors = FALSE
)
resolution$status <- ifelse(is.na(resolution$regulon), "MISSING", "FOUND")
write.table(
  resolution,
  file.path(opt$outdir, "SCENIC_target_regulon_resolution.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (all(is.na(resolution$regulon))) {
  stop("None of the manuscript target regulons were found in the pySCENIC output.")
}
if (any(is.na(resolution$regulon))) {
  warning(
    "Some manuscript target regulons were not found: ",
    paste(resolution$requested_tf[is.na(resolution$regulon)], collapse = ", ")
  )
}

make_target_table <- function(tfs, groups, label) {
  map <- resolution[match(tfs, resolution$requested_tf), , drop = FALSE]
  map <- map[!is.na(map$regulon), , drop = FALSE]
  if (nrow(map) == 0) return(NULL)

  out <- rss[map$regulon, groups, drop = FALSE]
  out <- cbind(
    requested_tf = map$requested_tf,
    regulon = map$regulon,
    as.data.frame(out, check.names = FALSE)
  )
  write.csv(
    out,
    file.path(opt$outdir, paste0(label, "_target_RSS.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  out
}

make_target_table(
  b_gc_tfs,
  c("dMMR_Stroma4_IGHM", "pMMR_Stroma4_IGHM"),
  "Fig3f_Stroma4_B_GC"
)
make_target_table(
  plasma_tfs,
  c("dMMR_Stroma1_IGLC1", "pMMR_Stroma1_IGLC1"),
  "Fig3h_Stroma1_plasma"
)

make_ranking_plot <- function(rss_mat, groups, target_tfs, title, outfile) {
  target_map <- resolution[match(target_tfs, resolution$requested_tf), , drop = FALSE]
  target_regulons <- na.omit(target_map$regulon)

  df <- rss_mat[, groups, drop = FALSE] %>%
    tibble::rownames_to_column("regulon") %>%
    pivot_longer(-regulon, names_to = "group", values_to = "RSS") %>%
    group_by(group) %>%
    arrange(desc(RSS), .by_group = TRUE) %>%
    mutate(rank = row_number(), target = regulon %in% target_regulons) %>%
    ungroup()

  label_map <- setNames(target_map$requested_tf, target_map$regulon)
  highlight <- df %>%
    filter(target) %>%
    mutate(label = unname(label_map[regulon]))

  p <- ggplot(df, aes(x = rank, y = RSS)) +
    geom_line(linewidth = 0.55) +
    geom_point(data = highlight, size = 2.4) +
    ggrepel::geom_text_repel(
      data = highlight,
      aes(label = label),
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    facet_wrap(~group, nrow = 1) +
    theme_classic(base_size = 12) +
    labs(title = title, x = "Regulon rank", y = "Regulon specificity score (RSS)") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold")
    )

  ggsave(outfile, p, width = 8, height = 4.5)
  p
}

make_ranking_plot(
  rss,
  c("dMMR_Stroma4_IGHM", "pMMR_Stroma4_IGHM"),
  b_gc_tfs,
  "B-cell / germinal-center regulons in Stroma4_IGHM",
  file.path(opt$outdir, "Fig3f_SCENIC_RSS_Stroma4_B_GC.pdf")
)

make_ranking_plot(
  rss,
  c("dMMR_Stroma1_IGLC1", "pMMR_Stroma1_IGLC1"),
  plasma_tfs,
  "Plasma-cell regulons in Stroma1_IGLC1",
  file.path(opt$outdir, "Fig3h_SCENIC_RSS_Stroma1_plasma.pdf")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(opt$outdir, "sessionInfo_SCENIC_RSS_Stroma1_4.txt")
)

message("SCENIC RSS analysis completed: ", opt$outdir)
