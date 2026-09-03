#!/usr/bin/env Rscript

# Stage 04: identify and characterise candidate amplicons from informative
# target-clade SNP geography.
#
# Inputs:
# - data/processed/precomputed.rds
# - outputs/tables/site_node_scores_targets.rds
# - analysis.amplicons settings from config/config.yml
#
# Outputs:
# - outputs/tables/candidate_amplicons.csv
# - outputs/tables/candidate_amplicon_target_scores.csv
# - outputs/tables/candidate_amplicon_snps.csv
# - matching .rds files for the same tables
# - outputs/figures/candidate_amplicons_spatial.png
#
# Stage 04 now answers: where are assay-feasible genomic regions containing
# discriminatory signal for the primary clades? It does not select a final
# panel or design primers.

source("R/window_scoring.R")
source("R/script_utils.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")
figure_dir <- file.path(config$paths$outputs, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

precomputed_path <- file.path(processed_dir, "precomputed.rds")
if (!file.exists(precomputed_path)) {
  stop("Missing ", precomputed_path, ". Run scripts/01_preprocess_alignment_tree.R first.", call. = FALSE)
}

score_path <- file.path(table_dir, "site_node_scores_targets.rds")
if (!file.exists(score_path)) {
  stop("Missing ", score_path, ". Run scripts/02_score_snps.R first.", call. = FALSE)
}

pre <- readRDS(precomputed_path)
site_node_scores_targets <- readRDS(score_path)
primary_target_nodes <- resolve_primary_target_nodes(config)
alignment_length <- ncol(pre$aln_int)

message("Scoring candidate amplicons")
message("  alignment length: ", alignment_length)
message("  target SNP rows read: ", nrow(site_node_scores_targets))
message("  primary target nodes: ", paste(primary_target_nodes, collapse = ", "))

target_rows <- site_node_scores_targets[site_node_scores_targets$node_id %in% primary_target_nodes, , drop = FALSE]
ignored_rows <- nrow(site_node_scores_targets) - nrow(target_rows)
if (ignored_rows > 0L) {
  message("  ignored non-primary diagnostic rows: ", ignored_rows)
}

scored <- score_candidate_amplicons(
  site_node_scores = target_rows,
  polymorphic_sites = pre$polymorphic_sites,
  aln_int = pre$aln_int,
  site_map = pre$site_map,
  config = config
)

candidate_amplicons_path <- file.path(table_dir, "candidate_amplicons.csv")
candidate_target_scores_path <- file.path(table_dir, "candidate_amplicon_target_scores.csv")
candidate_snps_path <- file.path(table_dir, "candidate_amplicon_snps.csv")

write.csv(scored$candidate_amplicons, candidate_amplicons_path, row.names = FALSE)
write.csv(scored$candidate_amplicon_target_scores, candidate_target_scores_path, row.names = FALSE)
write.csv(scored$candidate_amplicon_snps, candidate_snps_path, row.names = FALSE)

saveRDS(scored$candidate_amplicons, file.path(table_dir, "candidate_amplicons.rds"))
saveRDS(scored$candidate_amplicon_target_scores, file.path(table_dir, "candidate_amplicon_target_scores.rds"))
saveRDS(scored$candidate_amplicon_snps, file.path(table_dir, "candidate_amplicon_snps.rds"))

if (nrow(scored$candidate_amplicons) > 0L) {
  plot_candidate_amplicons_spatial(
    candidate_amplicons = scored$candidate_amplicons,
    candidate_amplicon_snps = scored$candidate_amplicon_snps,
    primary_target_nodes = primary_target_nodes,
    output_path = file.path(figure_dir, "candidate_amplicons_spatial.png"),
    alignment_length = alignment_length
  )
}

message("Candidate amplicon scoring complete")
message("  candidate amplicons retained: ", nrow(scored$candidate_amplicons))
message("  candidate-target rows produced: ", nrow(scored$candidate_amplicon_target_scores))
message("  candidate-SNP rows produced: ", nrow(scored$candidate_amplicon_snps))
message("  wrote: ", candidate_amplicons_path)
message("  wrote: ", candidate_target_scores_path)
message("  wrote: ", candidate_snps_path)
if (nrow(scored$candidate_amplicons) > 0L) {
  message("  wrote: ", file.path(figure_dir, "candidate_amplicons_spatial.png"))
}
