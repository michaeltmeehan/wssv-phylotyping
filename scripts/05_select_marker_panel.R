#!/usr/bin/env Rscript

# Stage 05: greedily select complementary marker windows for requested panel
# sizes and compare simple ranked-window alternatives.
# Inputs: outputs/tables/window_summary.* and window_node_summary.*, plus
# analysis.panel_selection settings from config/config.yml.
# Outputs: selected_panel, selected_panel_steps, panel_node_coverage, and
# panel_summary as CSV/RDS files under outputs/tables/.
# This remains part of the legacy assay-oriented panel-selection path that will
# be revisited after the stage 01-03 naming cleanup.
# Run directly after scripts/04_score_windows.R.

source("R/panel_selection.R")
source("R/script_utils.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
table_dir <- file.path(config$paths$outputs, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

window_summary <- read_table_output(table_dir, "window_summary", "scripts/04_score_windows.R")
window_node_summary <- read_table_output(table_dir, "window_node_summary", "scripts/04_score_windows.R")

panel_cfg <- config$analysis$panel_selection
if (is.null(panel_cfg)) {
  panel_cfg <- list()
}

panel_sizes <- as.integer(panel_cfg$panel_sizes)
if (length(panel_sizes) == 0L || all(is.na(panel_sizes))) {
  panel_sizes <- c(1L, 2L, 3L, 5L, 10L)
}
panel_sizes <- sort(unique(panel_sizes[!is.na(panel_sizes) & panel_sizes > 0L]))
max_panel_size <- max(panel_sizes)

min_informative_snps <- panel_cfg$min_informative_snps
if (is.null(min_informative_snps)) {
  min_informative_snps <- 1L
}
min_total_weighted_gain <- panel_cfg$min_total_weighted_gain
if (is.null(min_total_weighted_gain)) {
  min_total_weighted_gain <- 0
}
min_marginal_gain <- panel_cfg$min_marginal_gain
if (is.null(min_marginal_gain)) {
  min_marginal_gain <- 0
}
max_allowed_overlap <- panel_cfg$max_allowed_overlap
if (is.null(max_allowed_overlap)) {
  max_allowed_overlap <- 1
}
repeated_node_credit <- panel_cfg$repeated_node_credit
if (is.null(repeated_node_credit)) {
  repeated_node_credit <- 0
}
permitted_widths <- panel_cfg$permitted_widths
if (!is.null(permitted_widths)) {
  permitted_widths <- as.integer(unlist(permitted_widths))
}
cap_repeated_coverage <- panel_cfg$cap_repeated_coverage
if (is.null(cap_repeated_coverage)) {
  cap_repeated_coverage <- TRUE
}

filtered <- filter_candidate_windows(
  window_summary = window_summary,
  min_informative_snps = as.integer(min_informative_snps),
  min_total_weighted_gain = min_total_weighted_gain,
  max_missing_fraction = panel_cfg$max_missing_fraction,
  permitted_widths = permitted_widths,
  exclude_near_duplicates = isTRUE(panel_cfg$exclude_near_duplicates),
  near_duplicate_max_overlap = max_allowed_overlap,
  near_duplicate_min_distance = panel_cfg$min_distance
)

coverage <- calculate_window_node_coverage(
  window_node_summary = window_node_summary,
  retained_window_ids = filtered$window_id
)

selection <- greedy_select_panel(
  window_summary = filtered,
  node_coverage = coverage,
  max_panel_size = max_panel_size,
  min_marginal_gain = min_marginal_gain,
  max_allowed_overlap = max_allowed_overlap,
  min_distance = panel_cfg$min_distance,
  repeated_node_credit = repeated_node_credit,
  cap_repeated_coverage = isTRUE(cap_repeated_coverage)
)

greedy_summary <- summarise_selected_panels(
  selected_panel = selection$selected_panel,
  node_coverage = coverage,
  panel_sizes = panel_sizes,
  label = "greedy"
)
alternative_summary <- compare_panel_alternatives(
  window_summary = filtered,
  node_coverage = coverage,
  panel_sizes = panel_sizes,
  max_allowed_overlap = max_allowed_overlap,
  min_distance = panel_cfg$min_distance
)
panel_summary <- rbind(greedy_summary, alternative_summary)
panel_node_coverage <- summarise_panel_node_coverage_by_size(
  selected_panel = selection$selected_panel,
  node_coverage = coverage,
  panel_sizes = panel_sizes,
  label = "greedy"
)

selected_panel_path <- file.path(table_dir, "selected_panel.csv")
selected_panel_steps_path <- file.path(table_dir, "selected_panel_steps.csv")
panel_node_coverage_path <- file.path(table_dir, "panel_node_coverage.csv")
panel_summary_path <- file.path(table_dir, "panel_summary.csv")

write.csv(selection$selected_panel, selected_panel_path, row.names = FALSE)
write.csv(selection$steps, selected_panel_steps_path, row.names = FALSE)
write.csv(panel_node_coverage, panel_node_coverage_path, row.names = FALSE)
write.csv(panel_summary, panel_summary_path, row.names = FALSE)

saveRDS(selection$selected_panel, file.path(table_dir, "selected_panel.rds"))
saveRDS(selection$steps, file.path(table_dir, "selected_panel_steps.rds"))
saveRDS(panel_node_coverage, file.path(table_dir, "panel_node_coverage.rds"))
saveRDS(panel_summary, file.path(table_dir, "panel_summary.rds"))

message("Marker panel selection complete")
message("  candidate windows considered: ", nrow(window_summary))
message("  retained after filtering: ", nrow(filtered))
message("  filtered out: ", nrow(window_summary) - nrow(filtered))
message("  selected windows: ", nrow(selection$selected_panel))
message("  requested panel sizes: ", paste(panel_sizes, collapse = ", "))
message("  selected panel sizes available: ", paste(panel_sizes[panel_sizes <= nrow(selection$selected_panel)], collapse = ", "))
message("  cumulative greedy score: ", sum(selection$selected_panel$marginal_gain, na.rm = TRUE))
message("  nodes covered: ", nrow(selection$node_coverage))
message("  wrote: ", selected_panel_path)
message("  wrote: ", selected_panel_steps_path)
message("  wrote: ", panel_node_coverage_path)
message("  wrote: ", panel_summary_path)
