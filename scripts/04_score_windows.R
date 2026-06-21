#!/usr/bin/env Rscript

# Stage 04: build fixed and optional SNP-centred genomic windows, then
# aggregate informative SNP scores into candidate marker-window summaries.
# Inputs: data/processed/precomputed.rds, outputs/tables/site_node_scores.rds,
# and analysis.windows settings from config/config.yml.
# Outputs: candidate_windows, window_summary, and window_node_summary as CSV/RDS
# files under outputs/tables/.
# Run directly after scripts/02_score_snps.R.

source("R/window_scoring.R")

parse_window_steps <- function(steps) {
  if (is.null(steps)) {
    return(NULL)
  }
  if (is.list(steps)) {
    out <- as.integer(unlist(steps))
    names(out) <- names(steps)
    return(out)
  }
  as.integer(steps)
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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
alignment_length <- ncol(pre$aln_int)
window_cfg <- config$analysis$windows
if (is.null(window_cfg)) {
  window_cfg <- list(widths = config$analysis$candidate_window_sizes)
}

widths <- as.integer(window_cfg$widths)
if (length(widths) == 0L || any(is.na(widths))) {
  stop("Configure analysis.windows.widths with one or more positive integers.", call. = FALSE)
}

fixed_windows <- generate_fixed_windows(
  alignment_length = alignment_length,
  widths = widths,
  step = parse_window_steps(window_cfg$steps),
  overlap = window_cfg$overlap
)

snp_centered_windows <- empty_windows()
if (isTRUE(window_cfg$snp_centered$enabled)) {
  snp_centered_widths <- window_cfg$snp_centered$widths
  if (is.null(snp_centered_widths)) {
    snp_centered_widths <- widths
  }
  snp_centered_windows <- generate_snp_centered_windows(
    sites = unique(site_node_scores$site),
    alignment_length = alignment_length,
    widths = as.integer(snp_centered_widths)
  )
}

all_windows <- unique(rbind(
  fixed_windows[-match("window_id", names(fixed_windows))],
  snp_centered_windows[-match("window_id", names(snp_centered_windows))]
))
all_windows <- add_window_ids(all_windows)

min_informative_snps <- window_cfg$min_informative_snps
if (is.null(min_informative_snps)) {
  min_informative_snps <- 1L
}

message("Scoring candidate windows")
message("  alignment length: ", alignment_length)
message("  polymorphic SNPs: ", length(pre$polymorphic_sites))
message("  informative SNPs: ", length(unique(site_node_scores$site)))
message("  generated windows: ", nrow(all_windows))

scored <- score_candidate_windows(
  windows = all_windows,
  polymorphic_sites = pre$polymorphic_sites,
  site_node_scores = site_node_scores,
  aln_int = pre$aln_int,
  min_informative_snps = as.integer(min_informative_snps)
)

candidate_windows_path <- file.path(table_dir, "candidate_windows.csv")
window_summary_path <- file.path(table_dir, "window_summary.csv")
window_node_summary_path <- file.path(table_dir, "window_node_summary.csv")
write.csv(scored$candidate_windows, candidate_windows_path, row.names = FALSE)
write.csv(scored$window_summary, window_summary_path, row.names = FALSE)
write.csv(scored$window_node_summary, window_node_summary_path, row.names = FALSE)

saveRDS(scored$candidate_windows, file.path(table_dir, "candidate_windows.rds"))
saveRDS(scored$window_summary, file.path(table_dir, "window_summary.rds"))
saveRDS(scored$window_node_summary, file.path(table_dir, "window_node_summary.rds"))

assigned <- assign_snps_to_windows(unique(site_node_scores$site), scored$candidate_windows)
message("Window scoring complete")
message("  candidate windows retained: ", nrow(scored$candidate_windows))
message("  candidate windows filtered out: ", nrow(all_windows) - nrow(scored$candidate_windows))
message("  informative SNP assignments retained: ", nrow(assigned))
message("  window-node rows produced: ", nrow(scored$window_node_summary))
message("  wrote: ", candidate_windows_path)
message("  wrote: ", window_summary_path)
message("  wrote: ", window_node_summary_path)
