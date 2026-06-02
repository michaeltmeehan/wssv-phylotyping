#!/usr/bin/env Rscript

source("R/reporting.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
output_dir <- config$paths$outputs
table_dir <- file.path(output_dir, "tables")
report_dir <- file.path(output_dir, "reports")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

pre_path <- file.path(processed_dir, "precomputed.rds")
pre <- if (file.exists(pre_path)) readRDS(pre_path) else NULL

tables <- list(
  site_node_scores = read_optional_output(table_dir, "site_node_scores"),
  site_summary = read_optional_output(table_dir, "site_summary"),
  node_summary = read_optional_output(table_dir, "node_summary"),
  candidate_windows = read_optional_output(table_dir, "candidate_windows"),
  window_summary = read_optional_output(table_dir, "window_summary"),
  window_node_summary = read_optional_output(table_dir, "window_node_summary"),
  selected_panel = read_optional_output(table_dir, "selected_panel"),
  panel_summary = read_optional_output(table_dir, "panel_summary"),
  panel_node_coverage = read_optional_output(table_dir, "panel_node_coverage"),
  validation_panel_summary = read_optional_output(table_dir, "validation_panel_summary"),
  validation_tip_summary = read_optional_output(table_dir, "validation_tip_summary"),
  validation_fragment_summary = read_optional_output(table_dir, "validation_fragment_summary"),
  validation_baseline_summary = read_optional_output(table_dir, "validation_baseline_summary"),
  classifier_training_summary = read_optional_output(table_dir, "classifier_training_summary"),
  classifier_tip_predictions = read_optional_output(table_dir, "classifier_tip_predictions"),
  classifier_node_evidence = read_optional_output(table_dir, "classifier_node_evidence"),
  partial_classifications = read_optional_output(table_dir, "partial_classifications")
)

data <- lapply(tables, function(x) x$data)
missing_note <- function(stem) paste0("Missing output: `", file.path(table_dir, paste0(stem, ".csv")), "` or matching `.rds`.")
safe_nrow <- function(x) if (is.null(x)) NA_integer_ else nrow(x)
safe_unique <- function(x, col) if (is.null(x) || !col %in% names(x)) NA_integer_ else length(unique(x[[col]]))
status_counts <- function(x) {
  if (is.null(x) || !"status" %in% names(x)) data.frame() else count_values(x$status)
}

partial_dir <- config$paths$partials
if (is.null(partial_dir)) {
  partial_dir <- config$analysis$partials$input_dir
}
partial_ext <- config$analysis$partials$fasta_extensions
if (is.null(partial_ext)) {
  partial_ext <- c("fa", "fas", "fasta", "fna")
}
partial_files <- if (!is.null(partial_dir) && dir.exists(partial_dir)) {
  list.files(partial_dir, pattern = paste0("\\.(", paste(partial_ext, collapse = "|"), ")$"), ignore.case = TRUE)
} else {
  character()
}

if (!is.null(data$site_summary)) {
  top_sites <- data$site_summary[order(-data$site_summary$weighted_gain_sum, -data$site_summary$nodes_helped), , drop = FALSE]
  invisible(write_report_extract(top_sites, file.path(table_dir, "report_top_sites.csv"), n = 50L))
}
if (!is.null(data$window_summary)) {
  top_windows <- data$window_summary[order(-data$window_summary$total_weighted_gain, -data$window_summary$nodes_covered), , drop = FALSE]
  invisible(write_report_extract(top_windows, file.path(table_dir, "report_top_windows.csv"), n = 50L))
}
invisible(write_report_extract(data$selected_panel, file.path(table_dir, "report_selected_panel.csv")))
invisible(write_report_extract(data$validation_panel_summary, file.path(table_dir, "report_validation_summary.csv")))
invisible(write_report_extract(data$classifier_training_summary, file.path(table_dir, "report_classifier_summary.csv")))

final_covered_nodes <- if (is.null(data$selected_panel) || nrow(data$selected_panel) == 0L) character() else {
  ids <- unlist(strsplit(paste(data$selected_panel$new_node_ids, collapse = ";"), ";", fixed = TRUE))
  unique(ids[nzchar(ids)])
}
eligible_node_ids <- if (!is.null(pre) && !is.null(pre$node_metadata) && "node_id" %in% names(pre$node_metadata)) {
  as.character(pre$node_metadata$node_id)
} else {
  character()
}
uncovered_nodes <- setdiff(eligible_node_ids, final_covered_nodes)

lines <- c(
  "# WSSV Phylotyping Exploratory Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This report summarizes existing pipeline outputs. It does not rerun SNP scoring, window scoring, marker-panel selection, validation, classifier training, or partial-genome classification.",
  "",
  "## Pipeline Summary",
  "",
  paste0("- Input alignment: `", config$paths$alignment, "`"),
  paste0("- Input tree: `", config$paths$tree, "`"),
  paste0("- Retained tips: ", if (is.null(pre)) "missing precomputed output" else nrow(pre$aln_int)),
  paste0("- Alignment length: ", if (is.null(pre)) "missing precomputed output" else ncol(pre$aln_int)),
  paste0("- Eligible internal nodes: ", if (is.null(pre)) "missing precomputed output" else nrow(pre$target_mask)),
  paste0("- Polymorphic sites: ", if (is.null(pre)) "missing precomputed output" else length(pre$polymorphic_sites)),
  paste0("- Informative/scored node-site pairs: ", safe_nrow(data$site_node_scores)),
  paste0("- Candidate windows: ", safe_nrow(data$candidate_windows)),
  paste0("- Selected panel windows: ", safe_nrow(data$selected_panel)),
  paste0("- Partial FASTA files present: ", if (length(partial_files) > 0L) "yes" else "no"),
  "",
  "## SNP Scoring Summary",
  "",
  "Top SNPs by weighted gain:",
  "",
  table_or_note(data$site_summary, c("site", "nodes_helped", "weighted_gain_sum", "gain_max", "normalized_gain_max"), missing_note = missing_note("site_summary")),
  "",
  "Top nodes by number of helpful SNPs:",
  "",
  table_or_note(if (is.null(data$node_summary)) NULL else data$node_summary[order(-data$node_summary$sites_helpful), , drop = FALSE], c("node_id", "sites_helpful", "best_site", "normalized_gain_max", "clade_size", "depth"), missing_note = missing_note("node_summary")),
  "",
  interpret_signal_concentration(data$site_summary),
  "",
  "Missingness is available mainly in window-level outputs; site-level missingness is not currently written in `site_summary.csv`.",
  "",
  "## Window Scoring Summary",
  "",
  "Top windows by weighted score:",
  "",
  table_or_note(if (is.null(data$window_summary)) NULL else data$window_summary[order(-data$window_summary$total_weighted_gain), , drop = FALSE], c("window_id", "start", "end", "width", "informative_snps", "nodes_covered", "total_weighted_gain", "missing_fraction"), missing_note = missing_note("window_summary")),
  "",
  "Top windows by node coverage:",
  "",
  table_or_note(if (is.null(data$window_summary)) NULL else data$window_summary[order(-data$window_summary$nodes_covered, -data$window_summary$total_weighted_gain), , drop = FALSE], c("window_id", "start", "end", "width", "informative_snps", "nodes_covered", "total_weighted_gain"), missing_note = missing_note("window_summary")),
  "",
  "Informative SNPs per window:",
  "",
  table_or_note(if (is.null(data$window_summary)) NULL else numeric_summary(data$window_summary$informative_snps), missing_note = missing_note("window_summary")),
  "",
  "Nodes covered per window:",
  "",
  table_or_note(if (is.null(data$window_summary)) NULL else numeric_summary(data$window_summary$nodes_covered), missing_note = missing_note("window_summary")),
  "",
  "The selected-panel section below is the clearest redundancy check: windows with positive marginal gain add complementary node coverage, while high total scores alone may still overlap in phylogenetic signal.",
  "",
  "## Selected Marker Panel Summary",
  "",
  "Selected candidate marker regions:",
  "",
  table_or_note(data$selected_panel, c("selection_step", "window_id", "start", "end", "width", "marginal_gain", "cumulative_gain", "newly_covered_nodes", "total_covered_nodes"), missing_note = missing_note("selected_panel")),
  "",
  "Panel-level summaries for requested sizes:",
  "",
  table_or_note(panel_size_subset(data$panel_summary), c("method", "panel_size", "windows_selected", "cumulative_gain", "total_window_score", "nodes_covered", "newly_covered_nodes_last_step"), missing_note = missing_note("panel_summary")),
  "",
  paste0("Covered nodes in the selected panel: ", if (length(final_covered_nodes) == 0L) "not available" else paste(final_covered_nodes, collapse = ", ")),
  paste0("Eligible nodes not covered by the selected panel: ", if (length(uncovered_nodes) == 0L) "none identified" else paste(uncovered_nodes, collapse = ", ")),
  "",
  "These are candidate marker regions from tree-informed signal summaries. They are not validated primer assays, and primer-feasible flanks have not been assessed by this reporting script.",
  "",
  "## Validation Summary",
  "",
  "Selected-panel validation summarizes signal availability and true-path support in complete genomes:",
  "",
  table_or_note(data$validation_panel_summary, c("method", "panel_size", "panel_width", "tips_assessed", "mean_observed_informative_sites", "mean_observed_informative_nodes", "mean_supported_true_nodes", "resolved_fraction"), missing_note = missing_note("validation_panel_summary")),
  "",
  "Tip-level signal availability summary:",
  "",
  table_or_note(if (is.null(data$validation_tip_summary)) NULL else numeric_summary(data$validation_tip_summary$observed_informative_sites), missing_note = missing_note("validation_tip_summary")),
  "",
  "Artificial-fragment signal summary:",
  "",
  table_or_note(data$validation_fragment_summary, c("fragment_id", "width", "start", "end", "mean_observed_informative_sites", "mean_observed_informative_nodes", "resolved_fraction"), missing_note = missing_note("validation_fragment_summary")),
  "",
  "Selected panel versus top-window and random-window baselines:",
  "",
  table_or_note(data$validation_baseline_summary, c("method", "replicate", "panel_size", "panel_width", "mean_observed_informative_sites", "mean_observed_informative_nodes", "mean_supported_true_nodes", "resolved_fraction"), n = 20L, missing_note = missing_note("validation_baseline_summary")),
  "",
  "Interpretation: signal availability counts observed informative sites/nodes; true-path support counts informative nodes compatible with each known training tip path. These summaries are not independent classifier accuracy unless a validation table explicitly implements that design.",
  "",
  "## Classifier Sanity-Check Summary",
  "",
  "Training-set/internal summary:",
  "",
  table_or_note(data$classifier_training_summary, missing_note = missing_note("classifier_training_summary")),
  "",
  "Prediction status counts:",
  "",
  table_or_note(status_counts(data$classifier_tip_predictions), missing_note = missing_note("classifier_tip_predictions")),
  "",
  "Informative sites used per training prediction:",
  "",
  table_or_note(if (is.null(data$classifier_tip_predictions)) NULL else numeric_summary(data$classifier_tip_predictions$observed_informative_sites), missing_note = missing_note("classifier_tip_predictions")),
  "",
  "Assigned-depth summary:",
  "",
  table_or_note(if (is.null(data$classifier_tip_predictions) || !"assigned_depth" %in% names(data$classifier_tip_predictions)) NULL else numeric_summary(data$classifier_tip_predictions$assigned_depth), missing_note = missing_note("classifier_tip_predictions")),
  "",
  "## Classifier Conflict Diagnostics",
  "",
  "Classifier statuses now separate absent evidence, weak evidence, nested true-path-compatible support, and strong incompatible off-path support. Nested support along a single tree path is treated as coherent evidence rather than conflict; strong support for incompatible clades remains conservative and is reported as conflicting.",
  "",
  "Conflict reason counts:",
  "",
  table_or_note(if (is.null(data$classifier_tip_predictions) || !"conflict_reason" %in% names(data$classifier_tip_predictions)) NULL else count_values(data$classifier_tip_predictions$conflict_reason), missing_note = missing_note("classifier_tip_predictions")),
  "",
  "Truth-path assignment counts:",
  "",
  table_or_note(if (is.null(data$classifier_tip_predictions) || !"assigned_on_true_path" %in% names(data$classifier_tip_predictions)) NULL else count_values(data$classifier_tip_predictions$assigned_on_true_path[!is.na(data$classifier_tip_predictions$assigned_node)]), missing_note = missing_note("classifier_tip_predictions")),
  "",
  "Training tips with strongest off-path support:",
  "",
  table_or_note(
    if (is.null(data$classifier_tip_predictions) || !"strongest_off_path_support" %in% names(data$classifier_tip_predictions)) NULL else data$classifier_tip_predictions[order(-data$classifier_tip_predictions$strongest_off_path_support), , drop = FALSE],
    c("query_id", "status", "assigned_node", "assigned_on_true_path", "deepest_supported_true_path_node", "strongest_on_path_node", "strongest_on_path_support", "strongest_off_path_node", "strongest_off_path_support", "on_path_nodes_supported", "off_path_nodes_supported", "conflict_reason"),
    n = 15L,
    missing_note = missing_note("classifier_tip_predictions")
  ),
  "",
  "Interpretation: `on_path_nodes_supported` counts strong supported nodes on the known training-tip MCC path; `off_path_nodes_supported` counts strong supported nodes outside that path. Weak off-path evidence is shown by strongest support columns but does not automatically make a prediction conflicting.",
  "",
  "These predictions are internal sanity checks on training genomes, not independent validation.",
  "",
  "## Partial-Genome Classification Summary",
  "",
  if (length(partial_files) == 0L) "The workflow ran successfully, but no external partial FASTA files were available in the configured partial input directory." else paste0("Partial FASTA files found: ", paste(partial_files, collapse = ", ")),
  "",
  "Partial classification rows:",
  "",
  table_or_note(data$partial_classifications, c("sequence_id", "mapped", "mapped_start", "mapped_end", "observed_non_missing_sites", "observed_informative_sites", "assigned_node", "support_score", "status", "note"), missing_note = missing_note("partial_classifications")),
  "",
  "Partial status counts:",
  "",
  table_or_note(status_counts(data$partial_classifications), missing_note = missing_note("partial_classifications")),
  "",
  "Mapped/unmapped counts:",
  "",
  table_or_note(if (is.null(data$partial_classifications) || !"mapped" %in% names(data$partial_classifications)) NULL else count_values(data$partial_classifications$mapped), missing_note = missing_note("partial_classifications")),
  "",
  "Assigned nodes among partial classifications:",
  "",
  table_or_note(if (is.null(data$partial_classifications) || !"assigned_node" %in% names(data$partial_classifications)) NULL else count_values(data$partial_classifications$assigned_node), missing_note = missing_note("partial_classifications")),
  "",
  "## Warnings And Next Actions",
  "",
  "- Inspect selected panel coordinates biologically.",
  "- Check primer-feasible flanks around candidate windows.",
  "- Compare selected regions with known WSSV marker and deletion regions.",
  "- Add external partial sequences when available.",
  "- Run partial classification only after coordinate mapping is confirmed.",
  "- Consider leave-one-out classifier validation if independent classifier accuracy is required."
)

report_path <- file.path(report_dir, "wssv_phylotyping_report.md")
writeLines(lines, report_path)

message("Exploratory report complete")
message("  wrote: ", report_path)
message("  optional extracts written under: ", table_dir)
