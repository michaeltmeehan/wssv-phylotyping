#!/usr/bin/env Rscript

# Stage 06: validate how much complete-genome signal is retained by selected
# panel sizes and by random/baseline fragment sets.
# Inputs: precomputed data, site_node_scores_all with fallback to the legacy
# site_node_scores alias, window_summary, selected_panel, and analysis.validation
# settings from config/config.yml.
# Outputs: validation_panel_summary, validation_tip_summary,
# validation_fragment_summary, and validation_baseline_summary as CSV/RDS files.
# This remains in the legacy assay-oriented workflow pending the later stage-04+
# redesign.
# Run directly after scripts/05_select_marker_panel.R.

source("R/validation.R")
source("R/script_utils.R")

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
site_node_scores_all <- read_table_output(table_dir, "site_node_scores_all", "scripts/02_score_snps.R", fallback_stems = "site_node_scores")
window_summary <- read_table_output(table_dir, "window_summary", "scripts/04_score_windows.R")
selected_panel <- read_table_output(table_dir, "selected_panel", "scripts/05_select_marker_panel.R")

validation_cfg <- config$analysis$validation
if (is.null(validation_cfg)) {
  validation_cfg <- list()
}

panel_sizes <- as.integer(validation_cfg$panel_sizes)
if (length(panel_sizes) == 0L || all(is.na(panel_sizes))) {
  panel_sizes <- as.integer(config$analysis$panel_selection$panel_sizes)
}
panel_sizes <- sort(unique(panel_sizes[!is.na(panel_sizes) & panel_sizes > 0L]))
requested_panel_sizes <- panel_sizes
panel_sizes <- panel_sizes[panel_sizes <= nrow(selected_panel)]

fragment_lengths <- as.integer(validation_cfg$fragment_lengths)
fragment_lengths <- fragment_lengths[!is.na(fragment_lengths) & fragment_lengths > 0L]
baseline_replicates <- validation_cfg$random_baseline_replicates
if (is.null(baseline_replicates)) {
  baseline_replicates <- 0L
}
min_informative_sites <- validation_cfg$min_informative_sites
if (is.null(min_informative_sites)) {
  min_informative_sites <- 1L
}
min_informative_nodes <- validation_cfg$min_informative_nodes
if (is.null(min_informative_nodes)) {
  min_informative_nodes <- 1L
}
seed <- validation_cfg$random_seed
if (is.null(seed)) {
  seed <- 1L
}
set.seed(as.integer(seed))

target_clade_mask <- pre$target_clade_mask
if (!is.null(target_clade_mask) && is.null(rownames(target_clade_mask)) && "node_id" %in% names(pre$target_node_metadata)) {
  rownames(target_clade_mask) <- pre$target_node_metadata$node_id
}

tip_rows <- lapply(panel_sizes, function(size) {
  evaluate_panel_signal(
    aln_int = pre$aln_int,
    windows = selected_panel,
    site_node_scores = site_node_scores_all,
    target_mask = target_clade_mask,
    panel_size = size,
    method = "selected_panel",
    min_informative_sites = as.integer(min_informative_sites),
    min_informative_nodes = as.integer(min_informative_nodes)
  )
})
validation_tip_summary <- if (length(tip_rows) == 0L) data.frame() else do.call(rbind, tip_rows)
validation_panel_summary <- summarise_validation_by_panel(validation_tip_summary)

fragments <- make_random_fragments(
  alignment_length = ncol(pre$aln_int),
  fragment_lengths = fragment_lengths,
  n_per_length = 1L,
  seed = as.integer(seed),
  label = "random_fixed"
)
validation_fragment_summary <- summarise_fragment_signal(
  aln_int = pre$aln_int,
  fragments = fragments,
  site_node_scores = site_node_scores_all,
  min_informative_sites = as.integer(min_informative_sites),
  min_informative_nodes = as.integer(min_informative_nodes)
)

baseline_rows <- list()
if (baseline_replicates > 0L && length(panel_sizes) > 0L) {
  for (size in panel_sizes) {
    panel_subset <- selected_panel[seq_len(size), , drop = FALSE]
    random_windows <- make_matched_random_windows(
      panel_windows = panel_subset,
      alignment_length = ncol(pre$aln_int),
      n_replicates = as.integer(baseline_replicates),
      seed = as.integer(seed) + size,
      label = "random_matched"
    )
    for (replicate in seq_len(as.integer(baseline_replicates))) {
      rep_windows <- random_windows[random_windows$replicate == replicate, , drop = FALSE]
      baseline_rows[[length(baseline_rows) + 1L]] <- evaluate_panel_signal(
        aln_int = pre$aln_int,
        windows = rep_windows,
        site_node_scores = site_node_scores_all,
        target_mask = target_clade_mask,
        panel_size = nrow(rep_windows),
        method = "random_matched",
        replicate = replicate,
        min_informative_sites = as.integer(min_informative_sites),
        min_informative_nodes = as.integer(min_informative_nodes)
      )
    }
  }
}
baseline_tip_summary <- if (length(baseline_rows) == 0L) data.frame() else do.call(rbind, baseline_rows)
validation_baseline_summary <- summarise_validation_by_panel(baseline_tip_summary)

top_rows <- lapply(panel_sizes, function(size) {
  top_windows <- window_summary[order(-window_summary$total_weighted_gain, -window_summary$informative_snps, window_summary$start), , drop = FALSE]
  evaluate_panel_signal(
    aln_int = pre$aln_int,
    windows = top_windows,
    site_node_scores = site_node_scores_all,
    target_mask = target_clade_mask,
    panel_size = size,
    method = "top_weighted_gain",
    min_informative_sites = as.integer(min_informative_sites),
    min_informative_nodes = as.integer(min_informative_nodes)
  )
})
if (length(top_rows) > 0L) {
  validation_baseline_summary <- rbind(validation_baseline_summary, summarise_validation_by_panel(do.call(rbind, top_rows)))
}

write.csv(validation_panel_summary, file.path(table_dir, "validation_panel_summary.csv"), row.names = FALSE)
write.csv(validation_tip_summary, file.path(table_dir, "validation_tip_summary.csv"), row.names = FALSE)
write.csv(validation_fragment_summary, file.path(table_dir, "validation_fragment_summary.csv"), row.names = FALSE)
write.csv(validation_baseline_summary, file.path(table_dir, "validation_baseline_summary.csv"), row.names = FALSE)

saveRDS(validation_panel_summary, file.path(table_dir, "validation_panel_summary.rds"))
saveRDS(validation_tip_summary, file.path(table_dir, "validation_tip_summary.rds"))
saveRDS(validation_fragment_summary, file.path(table_dir, "validation_fragment_summary.rds"))
saveRDS(validation_baseline_summary, file.path(table_dir, "validation_baseline_summary.rds"))

message("Panel validation scaffold complete")
message("  requested panel sizes: ", format_count_list(requested_panel_sizes))
message("  panel sizes assessed: ", format_count_list(panel_sizes))
message("  panel sizes skipped because unavailable: ", format_count_list(setdiff(requested_panel_sizes, panel_sizes)))
message("  tips assessed: ", nrow(pre$aln_int))
message("  artificial fragments generated: ", nrow(fragments))
message("  random baseline comparisons completed: ", nrow(validation_baseline_summary[validation_baseline_summary$method == "random_matched", , drop = FALSE]))
message("  wrote validation tables under: ", table_dir)
