#!/usr/bin/env Rscript

source("R/encoding.R")
source("R/preprocess.R")
source("R/tree_utils.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

alignment <- read_alignment_dnastringset(config$paths$alignment)
tree <- read_tree_any(config$paths$tree)
matched <- match_prune_reorder_alignment_tree(alignment, tree)

aln_int <- encode_alignment_int(matched$alignment)
clades <- generate_clade_masks(
  matched$tree,
  min_clade_size = config$analysis$min_clade_size,
  max_clade_frac = config$analysis$max_clade_frac
)
node_metadata <- calculate_node_metadata(matched$tree, clades)
polymorphic_sites <- find_polymorphic_sites(aln_int)
site_map <- make_site_map(aln_int)

precomputed <- list(
  params = config$analysis,
  paths = config$paths,
  tree = matched$tree,
  alignment_names = rownames(aln_int),
  aln_int = aln_int,
  site_map = site_map,
  clades = clades,
  target_mask = clades$target_mask,
  node_metadata = node_metadata,
  polymorphic_sites = polymorphic_sites,
  dropped_alignment = matched$dropped_alignment,
  dropped_tree = matched$dropped_tree
)

saveRDS(precomputed, file.path(processed_dir, "precomputed.rds"))

message("Preprocessing complete")
message("  tips retained: ", nrow(aln_int))
message("  sites: ", ncol(aln_int))
message("  eligible internal nodes: ", nrow(clades$target_mask))
message("  polymorphic sites: ", length(polymorphic_sites))
message("  dropped alignment sequences: ", length(matched$dropped_alignment))
message("  dropped tree tips: ", length(matched$dropped_tree))
message("  wrote: ", file.path(processed_dir, "precomputed.rds"))
