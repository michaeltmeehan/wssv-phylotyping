#!/usr/bin/env Rscript

# Stage 02: catalogue allele frequencies for every polymorphic site
# across every retained clade/node in the reference tree.
#
# Input:
#   data/processed/precomputed.rds
#
# Output:
#   data/processed/site_node_alleles.csv
#
# For each node x polymorphic-site combination, records:
#   - clade size
#   - number of observed/called taxa
#   - A/C/G/T counts
#   - A/C/G/T frequencies among called taxa
#   - call rate within the clade
#
# Missing/ambiguous alignment states are encoded as 0 and contribute to
# n_taxa but not n_called or allele frequencies.
#
# This stage is deliberately target-agnostic. It describes the allele
# distribution across the full retained tree and performs no SNP scoring,
# feature selection, or classification.

source("R/encoding.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

precomputed_path <- file.path(processed_dir, "precomputed.rds")

if (!file.exists(precomputed_path)) {
  stop(
    "Precomputed data not found: ", precomputed_path,
    "\nRun scripts/01_preprocess_alignment_tree.R first.",
    call. = FALSE
  )
}

precomputed <- readRDS(precomputed_path)

aln_int <- precomputed$aln_int
clade_mask <- precomputed$all_clade_mask
node_metadata <- precomputed$all_node_metadata
polymorphic_sites <- precomputed$polymorphic_sites
site_map <- precomputed$site_map

if (!identical(colnames(clade_mask), rownames(aln_int))) {
  stop(
    "Clade-mask taxa do not match the alignment ordering.",
    call. = FALSE
  )
}

if (length(polymorphic_sites) == 0L) {
  stop("No polymorphic sites were identified in stage 01.", call. = FALSE)
}

node_ids <- as.integer(rownames(clade_mask))

if (anyNA(node_ids)) {
  stop("Could not recover node IDs from clade-mask row names.", call. = FALSE)
}

n_nodes <- nrow(clade_mask)
n_sites <- length(polymorphic_sites)

results <- vector("list", n_nodes)

for (i in seq_len(n_nodes)) {
  
  in_clade <- clade_mask[i, ]
  clade_alignment <- aln_int[in_clade, polymorphic_sites, drop = FALSE]
  
  n_taxa <- sum(in_clade)
  
  n_called <- colSums(clade_alignment != 0L)
  
  n_A <- colSums(clade_alignment == 1L)
  n_C <- colSums(clade_alignment == 2L)
  n_G <- colSums(clade_alignment == 3L)
  n_T <- colSums(clade_alignment == 4L)
  
  p_A <- ifelse(n_called > 0L, n_A / n_called, NA_real_)
  p_C <- ifelse(n_called > 0L, n_C / n_called, NA_real_)
  p_G <- ifelse(n_called > 0L, n_G / n_called, NA_real_)
  p_T <- ifelse(n_called > 0L, n_T / n_called, NA_real_)
  
  results[[i]] <- data.frame(
    node_id = node_ids[[i]],
    site = polymorphic_sites,
    alignment_position = site_map$alignment_position[polymorphic_sites],
    
    n_taxa = n_taxa,
    n_called = n_called,
    n_missing = n_taxa - n_called,
    call_rate = n_called / n_taxa,
    
    n_A = n_A,
    n_C = n_C,
    n_G = n_G,
    n_T = n_T,
    
    p_A = p_A,
    p_C = p_C,
    p_G = p_G,
    p_T = p_T,
    
    stringsAsFactors = FALSE
  )
}

site_node_alleles <- do.call(rbind, results)

# Add useful tree metadata while keeping this table target-agnostic.
metadata_columns <- intersect(
  c(
    "node_id",
    "clade_size",
    "posterior_support",
    "depth",
    "normalized_depth",
    "balance"
  ),
  names(node_metadata)
)

if ("node_id" %in% metadata_columns) {
  
  metadata <- node_metadata[, metadata_columns, drop = FALSE]
  
  site_node_alleles <- merge(
    site_node_alleles,
    metadata,
    by = "node_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  site_node_alleles <- site_node_alleles[
    order(
      match(site_node_alleles$node_id, node_ids),
      match(site_node_alleles$site, polymorphic_sites)
    ),
  ]
}

rownames(site_node_alleles) <- NULL

output_path <- file.path(processed_dir, "site_node_alleles.csv")

write.csv(
  site_node_alleles,
  output_path,
  row.names = FALSE
)

message("Site-by-node allele catalogue complete")
message("  nodes: ", n_nodes)
message("  polymorphic sites: ", n_sites)
message("  node x site records: ", nrow(site_node_alleles))
message("  wrote: ", output_path)


rds_output_path <- file.path(processed_dir, "site_node_alleles.rds")

saveRDS(
  site_node_alleles,
  rds_output_path
)

message("  wrote: ", rds_output_path)
