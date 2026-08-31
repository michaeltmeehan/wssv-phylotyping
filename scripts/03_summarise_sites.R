#!/usr/bin/env Rscript

# Stage 03: summarise informative SNP and node-level scoring results.
# Inputs: data/processed/precomputed.rds plus the stage-02 all/target score
# tables.
# Outputs:
# - outputs/tables/site_summary_all.csv
# - outputs/tables/node_summary_all.csv
# - outputs/tables/site_summary_targets.csv
# - outputs/tables/node_summary_targets.csv
# - outputs/tables/target_clade_diagnostics.csv
# - outputs/tables/target_clade_strongest_snps.csv
# Legacy compatibility aliases are also written to outputs/tables/site_summary.csv
# and outputs/tables/node_summary.csv for the full diagnostic set.
# Use tree posterior support as descriptive metadata only; site ranking is based
# on SNP discriminatory gain, not phylogenetic support or later classifier
# evidence.
# Run directly after scripts/02_score_snps.R.

source("R/snp_scoring.R")
source("R/target_diagnostics.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")

precomputed_path <- file.path(processed_dir, "precomputed.rds")
all_score_path <- file.path(table_dir, "site_node_scores_all.rds")
target_score_path <- file.path(table_dir, "site_node_scores_targets.rds")
if (!file.exists(precomputed_path)) {
  stop("Missing ", precomputed_path, ". Run scripts/01_preprocess_alignment_tree.R first.", call. = FALSE)
}
if (!file.exists(all_score_path)) {
  stop("Missing ", all_score_path, ". Run scripts/02_score_snps.R first.", call. = FALSE)
}
if (!file.exists(target_score_path)) {
  stop("Missing ", target_score_path, ". Run scripts/02_score_snps.R first.", call. = FALSE)
}

pre <- readRDS(precomputed_path)
site_node_scores_all <- readRDS(all_score_path)
site_node_scores_targets <- readRDS(target_score_path)

site_summary_all <- make_site_summary(
  site_node_scores_all,
  site_map = pre$site_map,
  n_tip = nrow(pre$aln_int),
  helped_count_name = "nodes_helped"
)
node_summary_all <- make_node_summary(site_node_scores_all, pre$all_node_metadata)
site_summary_targets <- make_site_summary(
  site_node_scores_targets,
  site_map = pre$site_map,
  n_tip = nrow(pre$aln_int),
  helped_count_name = "target_clades_helped"
)
node_summary_targets <- make_node_summary(site_node_scores_targets, pre$target_node_metadata)
target_checkpoint <- make_target_clade_checkpoint(pre, site_node_scores_targets, top_n = 3L)

write.csv(site_summary_all, file.path(table_dir, "site_summary_all.csv"), row.names = FALSE)
write.csv(node_summary_all, file.path(table_dir, "node_summary_all.csv"), row.names = FALSE)
write.csv(site_summary_targets, file.path(table_dir, "site_summary_targets.csv"), row.names = FALSE)
write.csv(node_summary_targets, file.path(table_dir, "node_summary_targets.csv"), row.names = FALSE)
write.csv(target_checkpoint$summary, file.path(table_dir, "target_clade_diagnostics.csv"), row.names = FALSE)
write.csv(target_checkpoint$strongest_snps, file.path(table_dir, "target_clade_strongest_snps.csv"), row.names = FALSE)

# Backward-compatible aliases for downstream scripts and report generation that
# still expect the legacy full-diagnostic filenames.
write.csv(site_summary_all, file.path(table_dir, "site_summary.csv"), row.names = FALSE)
write.csv(node_summary_all, file.path(table_dir, "node_summary.csv"), row.names = FALSE)

message("Summaries complete")
message("  tips: ", nrow(pre$aln_int))
message("  sites: ", ncol(pre$aln_int))
message("  eligible internal nodes: ", nrow(pre$target_mask))
message("  polymorphic sites: ", length(pre$polymorphic_sites))
message("  scored node-site pairs (all): ", nrow(site_node_scores_all))
message("  scored node-site pairs (targets): ", nrow(site_node_scores_targets))
message("  site summary rows (all): ", nrow(site_summary_all))
message("  node summary rows (all): ", nrow(node_summary_all))
message("  site summary rows (targets): ", nrow(site_summary_targets))
message("  node summary rows (targets): ", nrow(node_summary_targets))
message("  target-clade diagnostics rows: ", nrow(target_checkpoint$summary))
message("  target-clade strongest SNP rows: ", nrow(target_checkpoint$strongest_snps))
message("  wrote: ", file.path(table_dir, "target_clade_diagnostics.csv"))
message("  wrote: ", file.path(table_dir, "target_clade_strongest_snps.csv"))
