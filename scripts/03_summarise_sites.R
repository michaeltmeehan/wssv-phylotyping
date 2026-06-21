#!/usr/bin/env Rscript

# Stage 03: summarise informative SNP and node-level scoring results.
# Inputs: data/processed/precomputed.rds and outputs/tables/site_node_scores.rds.
# Outputs: outputs/tables/site_summary.csv and outputs/tables/node_summary.csv.
# Uses the same config paths as upstream scripts; no new analysis thresholds.
# Run directly after scripts/02_score_snps.R.

source("R/snp_scoring.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")

precomputed_path <- file.path(processed_dir, "precomputed.rds")
score_path <- file.path(table_dir, "site_node_scores.rds")
if (!file.exists(precomputed_path)) {
  stop("Missing ", precomputed_path, ". Run scripts/01_preprocess_alignment_tree.R first.", call. = FALSE)
}
if (!file.exists(score_path)) {
  stop("Missing ", score_path, ". Run scripts/02_score_snps.R first.", call. = FALSE)
}

pre <- readRDS(precomputed_path)
site_node_scores <- readRDS(score_path)

site_summary <- make_site_summary(site_node_scores)
node_summary <- make_node_summary(site_node_scores, pre$node_metadata)

write.csv(site_summary, file.path(table_dir, "site_summary.csv"), row.names = FALSE)
write.csv(node_summary, file.path(table_dir, "node_summary.csv"), row.names = FALSE)

message("Summaries complete")
message("  tips: ", nrow(pre$aln_int))
message("  sites: ", ncol(pre$aln_int))
message("  eligible internal nodes: ", nrow(pre$target_mask))
message("  polymorphic sites: ", length(pre$polymorphic_sites))
message("  scored node-site pairs: ", nrow(site_node_scores))
message("  site summary rows: ", nrow(site_summary))
message("  node summary rows: ", nrow(node_summary))
