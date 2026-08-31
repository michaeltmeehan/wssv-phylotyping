#!/usr/bin/env Rscript

# Stage 02: score polymorphic SNPs against eligible MCC-tree clades.
# Inputs: data/processed/precomputed.rds and analysis SNP-scoring thresholds.
# Outputs:
# - outputs/tables/site_node_scores_all.rds for the full diagnostic clade set
# - outputs/tables/site_node_scores_targets.rds for the restricted target set
# - outputs/tables/site_node_scores.rds as a transitional alias for the full set
# Each file contains one row per informative site-node rule that passes the
# configured observation and gain filters.
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

score_clade_set <- function(label, clade_mask, node_metadata, output_path) {
  message("  scoring ", label, " clades: ", nrow(clade_mask))
  scores <- score_snp_sites(
    aln_int = pre$aln_int,
    target_mask = clade_mask,
    node_metadata = node_metadata,
    sites = sites,
    min_total_obs = config$analysis$min_total_obs,
    min_side_obs = config$analysis$min_side_obs,
    min_site_maf = config$analysis$min_site_maf,
    min_gain_norm = config$analysis$min_gain_norm,
    chunk_size = config$analysis$chunk_size,
    progress = interactive()
  )
  saveRDS(scores, output_path)
  scores
}

message("Scoring SNPs")
message("  tips: ", nrow(pre$aln_int))
message("  sites: ", ncol(pre$aln_int))
message("  polymorphic sites: ", length(sites))
message("  diagnostic nodes scored: ", nrow(pre$all_clade_mask))
message("  target nodes scored: ", nrow(pre$target_clade_mask))

all_score_path <- file.path(table_dir, "site_node_scores_all.rds")
target_score_path <- file.path(table_dir, "site_node_scores_targets.rds")
legacy_score_path <- file.path(table_dir, "site_node_scores.rds")

site_node_scores_all <- score_clade_set(
  label = "diagnostic",
  clade_mask = pre$all_clade_mask,
  node_metadata = pre$all_node_metadata,
  output_path = all_score_path
)
site_node_scores_targets <- score_clade_set(
  label = "target",
  clade_mask = pre$target_clade_mask,
  node_metadata = pre$target_node_metadata,
  output_path = target_score_path
)

saveRDS(site_node_scores_all, legacy_score_path)

message("SNP scoring complete")
message("  diagnostic site-node combinations retained: ", nrow(site_node_scores_all))
message("  target site-node combinations retained: ", nrow(site_node_scores_targets))
message("  wrote: ", all_score_path)
message("  wrote: ", target_score_path)
message("  legacy alias written: ", legacy_score_path)
