#!/usr/bin/env Rscript

# Stage 02: score polymorphic SNPs against eligible MCC-tree clades.
# Inputs: data/processed/precomputed.rds and analysis SNP-scoring thresholds.
# Outputs: outputs/tables/site_node_scores.rds, one row per informative
# site-node rule that passes the configured observation and gain filters.
# Run directly after scripts/01_preprocess_alignment_tree.R.

source("R/encoding.R")
source("R/snp_scoring.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

precomputed_path <- file.path(processed_dir, "precomputed.rds")
if (!file.exists(precomputed_path)) {
  stop("Missing ", precomputed_path, ". Run scripts/01_preprocess_alignment_tree.R first.", call. = FALSE)
}

pre <- readRDS(precomputed_path)
sites <- pre$polymorphic_sites

message("Scoring SNPs")
message("  tips: ", nrow(pre$aln_int))
message("  sites: ", ncol(pre$aln_int))
message("  eligible internal nodes: ", nrow(pre$target_mask))
message("  polymorphic sites: ", length(sites))

site_node_scores <- score_snp_sites(
  aln_int = pre$aln_int,
  target_mask = pre$target_mask,
  node_metadata = pre$node_metadata,
  sites = sites,
  min_total_obs = config$analysis$min_total_obs,
  min_side_obs = config$analysis$min_side_obs,
  min_site_maf = config$analysis$min_site_maf,
  min_gain_norm = config$analysis$min_gain_norm,
  chunk_size = config$analysis$chunk_size,
  progress = TRUE
)

score_path <- file.path(table_dir, "site_node_scores.rds")
saveRDS(site_node_scores, score_path)

message("SNP scoring complete")
message("  scored node-site pairs: ", nrow(site_node_scores))
message("  wrote: ", score_path)
