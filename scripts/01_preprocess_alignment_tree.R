#!/usr/bin/env Rscript

# Stage 01: preprocess the configured whole-genome alignment and MCC tree.
# Inputs: config/config.yml paths.alignment, paths.tree, and analysis clade filters.
# Outputs: data/processed/precomputed.rds with encoded alignment, pruned tree,
# broad diagnostic clades, target clades, node metadata, site map, and
# polymorphic-site coordinates.
# Run directly from the repository root before all downstream stages.

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
all_clades <- generate_clade_masks(
  matched$tree,
  min_clade_size = config$analysis$min_clade_size,
  max_clade_frac = config$analysis$max_clade_frac
)
all_node_metadata <- calculate_node_metadata(matched$tree, all_clades)
target_clades <- resolve_target_clades(
  matched$tree,
  all_clades,
  config$analysis$target_clades,
  min_posterior_support = config$analysis$min_target_posterior_support,
  min_clade_size = config$analysis$min_clade_size,
  all_node_metadata = all_node_metadata
)
target_node_metadata <- calculate_node_metadata(matched$tree, target_clades)
polymorphic_sites <- find_polymorphic_sites(aln_int)
site_map <- make_site_map(aln_int)

precomputed <- list(
  params = config$analysis,
  paths = config$paths,
  tree = matched$tree,
  alignment_names = rownames(aln_int),
  aln_int = aln_int,
  site_map = site_map,
  all_clades = all_clades,
  all_clade_mask = all_clades$all_clade_mask,
  all_node_metadata = all_node_metadata,
  target_clades = target_clades,
  target_clade_mask = target_clades$target_clade_mask,
  target_node_metadata = target_node_metadata,
  # Compatibility aliases for downstream stages that have not yet been
  # migrated to the new terminology.
  clades = all_clades,
  target_mask = all_clades$all_clade_mask,
  node_metadata = all_node_metadata,
  polymorphic_sites = polymorphic_sites,
  dropped_alignment = matched$dropped_alignment,
  dropped_tree = matched$dropped_tree
)

saveRDS(precomputed, file.path(processed_dir, "precomputed.rds"))

message("Preprocessing complete")
message("  tips retained: ", nrow(aln_int))
message("  sites: ", ncol(aln_int))
message("  diagnostic clades: ", nrow(all_clades$all_clade_mask))
message("  target clades: ", nrow(target_clades$target_clade_mask))
message("  target node IDs: ", paste(target_node_metadata$node_id, collapse = ", "))
message("  target clade sizes: ", paste(target_node_metadata$clade_size, collapse = ", "))
message(
  "  target posterior supports: ",
  paste(sprintf("%.6f", target_node_metadata$posterior_support), collapse = ", ")
)
message("  polymorphic sites: ", length(polymorphic_sites))
message("  dropped alignment sequences: ", length(matched$dropped_alignment))
message("  dropped tree tips: ", length(matched$dropped_tree))
message("  wrote: ", file.path(processed_dir, "precomputed.rds"))
