#!/usr/bin/env Rscript

# Stage 05: select and compare small candidate amplicon panels from the
# stage-04 amplicon outputs.
#
# Inputs:
# - outputs/tables/candidate_amplicons.*
# - outputs/tables/candidate_amplicon_target_scores.*
# - optional outputs/tables/candidate_amplicon_snps.*
# - analysis.panel_selection settings from config/config.yml
#
# Outputs:
# - outputs/tables/panel_candidates.csv
# - outputs/tables/panel_target_scores.csv
# - outputs/tables/pareto_panels_by_size.csv
# - outputs/tables/cross_size_pareto_panels.csv
# - outputs/tables/recommended_provisional_panel.csv
# - outputs/tables/recommended_panel_alternatives.csv
# - outputs/tables/panel_size_comparison.csv
# - outputs/tables/panel_summary.csv (compatibility alias)
# - outputs/tables/panel_search_summary.csv
# - outputs/figures/recommended_provisional_panel_panels.png
# - outputs/figures/selected_panel_panels.png (compatibility alias)

source("R/panel_selection.R")
source("R/script_utils.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

read_optional_table <- function(table_dir, stem) {
  rds_path <- file.path(table_dir, paste0(stem, ".rds"))
  csv_path <- file.path(table_dir, paste0(stem, ".csv"))
  if (file.exists(rds_path) || file.exists(csv_path)) {
    read_table_output(table_dir, stem, "scripts/04_score_windows.R")
  } else {
    NULL
  }
}

plot_recommended_panel_overview <- function(recommended_panel, panel_size_comparison, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(output_path, width = 1800, height = 1400, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (is.null(recommended_panel) || nrow(recommended_panel) == 0L) {
    graphics::plot.new()
    graphics::title(main = "Recommended provisional panel", sub = "No panel available")
    return(invisible(NULL))
  }

  targets <- sort(as.integer(gsub("^gain_", "", grep("^gain_[0-9]+$", names(recommended_panel), value = TRUE))))
  gain_cols <- paste0("gain_", targets)
  rel_cols <- paste0("relative_gain_", targets)
  panel <- recommended_panel[order(recommended_panel$start, recommended_panel$end, recommended_panel$candidate_id), , drop = FALSE]
  panel_index <- seq_len(nrow(panel))
  best_target <- targets[apply(panel[gain_cols], 1L, which.max)]
  colors <- c("grey60", "#1b9e77", "#d95f02", "#7570b3", "#e7298a")
  target_colors <- stats::setNames(colors[seq_along(targets) + 1L], as.character(targets))

  graphics::par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  graphics::plot(
    NA,
    xlim = c(min(panel$start), max(panel$end)),
    ylim = c(0.5, nrow(panel) + 0.5),
    xlab = "Genome position",
    ylab = "Panel amplicon",
    yaxt = "n",
    main = "Recommended provisional amplicon panel"
  )
  graphics::axis(2, at = panel_index, labels = panel$candidate_id, las = 2, cex.axis = 0.7)
  for (i in seq_len(nrow(panel))) {
    col <- target_colors[[as.character(best_target[[i]])]]
    if (is.null(col)) col <- "grey70"
    graphics::rect(
      panel$start[[i]], i - 0.35, panel$end[[i]], i + 0.35,
      col = grDevices::adjustcolor(col, alpha.f = 0.35),
      border = col,
      lwd = 2
    )
    graphics::text(
      x = panel$start[[i]] + (panel$end[[i]] - panel$start[[i]]) / 2,
      y = i + 0.48,
      labels = paste0("len=", panel$length[[i]]),
      cex = 0.65
    )
  }

  rel_matrix <- as.matrix(panel[rel_cols])
  graphics::plot(
    NA,
    xlim = c(0.5, length(targets) + 0.5),
    ylim = c(0.5, nrow(panel) + 0.5),
    xlab = "Primary target",
    ylab = "Panel amplicon",
    yaxt = "n",
    xaxt = "n",
    main = "Target-relative gain heatmap"
  )
  graphics::axis(1, at = seq_along(targets), labels = targets)
  graphics::axis(2, at = seq_len(nrow(panel)), labels = panel$candidate_id, las = 2, cex.axis = 0.7)
  ramp <- grDevices::colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(50)
  for (i in seq_len(nrow(panel))) {
    for (j in seq_along(targets)) {
      value <- rel_matrix[[i, j]]
      if (is.na(value)) value <- 0
      fill <- ramp[max(1L, min(length(ramp), as.integer(round(value * (length(ramp) - 1L))) + 1L))]
      graphics::rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = fill, border = "white")
      graphics::text(j, i, labels = formatC(value, digits = 2, format = "f"), cex = 0.7)
    }
  }
  graphics::box()
  graphics::mtext("Target-relative gain", side = 3, line = 0.5)

  if (!is.null(panel_size_comparison) && nrow(panel_size_comparison) > 0L) {
    valid <- panel_size_comparison[!is.na(panel_size_comparison$best_panel_id), , drop = FALSE]
    graphics::par(mfrow = c(1, 1), mar = c(5, 4, 3, 4))
    if (nrow(valid) > 0L) {
      graphics::plot(
        valid$panel_size,
        valid$min_absolute_target_gain,
        type = "b",
        pch = 19,
        col = "#1b9e77",
        ylim = range(c(valid$min_absolute_target_gain, valid$min_target_relative_gain, valid$capped_target_redundancy), na.rm = TRUE),
        xlab = "Panel size",
        ylab = "Best-by-size summary",
        main = "Best panel by size"
      )
      graphics::lines(valid$panel_size, valid$min_target_relative_gain, type = "b", pch = 17, col = "#d95f02")
      graphics::lines(valid$panel_size, valid$capped_target_redundancy, type = "b", pch = 15, col = "#7570b3")
      graphics::legend(
        "bottomright",
        legend = c("Min absolute gain", "Min relative gain", "Capped redundancy"),
        col = c("#1b9e77", "#d95f02", "#7570b3"),
        pch = c(19, 17, 15),
        lty = 1,
        bty = "n",
        cex = 0.8
      )
      if (any(valid$recommended_provisional_size, na.rm = TRUE)) {
        recommended <- valid[valid$recommended_provisional_size, , drop = FALSE][1L, , drop = FALSE]
        graphics::abline(v = recommended$panel_size, lty = 2, col = "grey40")
        graphics::mtext(paste0("Recommended size: ", recommended$panel_size), side = 3, line = 0.2, adj = 0)
      }
      graphics::text(valid$panel_size, valid$min_absolute_target_gain, labels = valid$best_panel_id, pos = 3, cex = 0.65)
    } else {
      graphics::plot.new()
      graphics::title(main = "Best panel by size", sub = "No feasible panel sizes available")
    }
  }
}

config <- yaml::read_yaml("config/config.yml")
table_dir <- file.path(config$paths$outputs, "tables")
figure_dir <- file.path(config$paths$outputs, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

candidate_amplicons <- read_table_output(table_dir, "candidate_amplicons", "scripts/04_score_windows.R")
candidate_amplicon_target_scores <- read_table_output(table_dir, "candidate_amplicon_target_scores", "scripts/04_score_windows.R")
candidate_amplicon_snps <- read_optional_table(table_dir, "candidate_amplicon_snps")

selection <- select_amplicon_panels(
  candidate_amplicons = candidate_amplicons,
  candidate_amplicon_target_scores = candidate_amplicon_target_scores,
  candidate_amplicon_snps = candidate_amplicon_snps,
  config = config
)

panel_candidates_path <- file.path(table_dir, "panel_candidates.csv")
panel_target_scores_path <- file.path(table_dir, "panel_target_scores.csv")
pareto_panels_by_size_path <- file.path(table_dir, "pareto_panels_by_size.csv")
cross_size_pareto_panels_path <- file.path(table_dir, "cross_size_pareto_panels.csv")
recommended_panel_path <- file.path(table_dir, "recommended_provisional_panel.csv")
recommended_panel_alternatives_path <- file.path(table_dir, "recommended_panel_alternatives.csv")
selected_panel_path <- file.path(table_dir, "selected_panel.csv")
pareto_panels_path <- file.path(table_dir, "pareto_panels.csv")
panel_size_comparison_path <- file.path(table_dir, "panel_size_comparison.csv")
panel_summary_path <- file.path(table_dir, "panel_summary.csv")
panel_search_summary_path <- file.path(table_dir, "panel_search_summary.csv")

write.csv(selection$panel_candidates, panel_candidates_path, row.names = FALSE)
write.csv(selection$panel_target_scores, panel_target_scores_path, row.names = FALSE)
write.csv(selection$pareto_panels_by_size, pareto_panels_by_size_path, row.names = FALSE)
write.csv(selection$cross_size_pareto_panels, cross_size_pareto_panels_path, row.names = FALSE)
write.csv(selection$recommended_provisional_panel, recommended_panel_path, row.names = FALSE)
write.csv(selection$recommended_panel_alternatives, recommended_panel_alternatives_path, row.names = FALSE)
write.csv(selection$selected_panel, selected_panel_path, row.names = FALSE)
write.csv(selection$pareto_panels, pareto_panels_path, row.names = FALSE)
write.csv(selection$panel_size_comparison, panel_size_comparison_path, row.names = FALSE)
write.csv(selection$panel_size_comparison, panel_summary_path, row.names = FALSE)
write.csv(panel_search_summary(
  selection$candidate_count_before,
  selection$candidate_count_after_filter,
  selection$candidate_count_after_reduction,
  selection$panel_size_comparison,
  selection$config
), panel_search_summary_path, row.names = FALSE)

saveRDS(selection$panel_candidates, file.path(table_dir, "panel_candidates.rds"))
saveRDS(selection$panel_target_scores, file.path(table_dir, "panel_target_scores.rds"))
saveRDS(selection$pareto_panels_by_size, file.path(table_dir, "pareto_panels_by_size.rds"))
saveRDS(selection$cross_size_pareto_panels, file.path(table_dir, "cross_size_pareto_panels.rds"))
saveRDS(selection$recommended_provisional_panel, file.path(table_dir, "recommended_provisional_panel.rds"))
saveRDS(selection$recommended_panel_alternatives, file.path(table_dir, "recommended_panel_alternatives.rds"))
saveRDS(selection$selected_panel, file.path(table_dir, "selected_panel.rds"))
saveRDS(selection$pareto_panels, file.path(table_dir, "pareto_panels.rds"))
saveRDS(selection$panel_size_comparison, file.path(table_dir, "panel_size_comparison.rds"))
saveRDS(selection$panel_size_comparison, file.path(table_dir, "panel_summary.rds"))
saveRDS(panel_search_summary(
  selection$candidate_count_before,
  selection$candidate_count_after_filter,
  selection$candidate_count_after_reduction,
  selection$panel_size_comparison,
  selection$config
), file.path(table_dir, "panel_search_summary.rds"))

plot_recommended_panel_overview(
  selection$recommended_provisional_panel,
  selection$panel_size_comparison,
  file.path(figure_dir, "recommended_provisional_panel_panels.png")
)
file.copy(
  file.path(figure_dir, "recommended_provisional_panel_panels.png"),
  file.path(figure_dir, "selected_panel_panels.png"),
  overwrite = TRUE
)

size_table <- selection$panel_size_comparison
if (nrow(size_table) > 0L) {
  display_cols <- c(
    "panel_size", "search_method", "panels_evaluated", "within_size_pareto_count",
    "min_absolute_target_gain", "min_target_relative_gain", "capped_target_redundancy",
    "total_length"
  )
  display_cols <- display_cols[display_cols %in% names(size_table)]
  message("Stage 05 per-size summary")
  print(size_table[display_cols], row.names = FALSE)
}

recommended_size <- selection$recommended_panel_size
recommended_panel_ids <- if (nrow(selection$recommended_provisional_panel) == 0L) {
  "n/a"
} else {
  paste(selection$recommended_provisional_panel$candidate_id, collapse = ";")
}

message("Stage 05 panel selection complete")
message("  candidates before filtering: ", selection$candidate_count_before)
message("  candidates after filtering: ", selection$candidate_count_after_filter)
message("  candidates after reduction: ", selection$candidate_count_after_reduction)
message("  total panels evaluated: ", selection$panels_evaluated)
message("  recommended provisional size: ", if (is.na(recommended_size)) "n/a" else recommended_size)
message("  recommended provisional panel ID: ", if (is.na(recommended_size) || nrow(selection$panel_size_comparison) == 0L) "n/a" else selection$panel_size_comparison$best_panel_id[match(recommended_size, selection$panel_size_comparison$panel_size)])
message("  recommended provisional candidate IDs: ", recommended_panel_ids)
message("  recommendation reason: ", selection$recommended_reason)
message("  wrote: ", panel_candidates_path)
message("  wrote: ", panel_target_scores_path)
message("  wrote: ", pareto_panels_by_size_path)
message("  wrote: ", cross_size_pareto_panels_path)
message("  wrote: ", recommended_panel_path)
message("  wrote: ", recommended_panel_alternatives_path)
message("  wrote: ", selected_panel_path)
message("  wrote: ", pareto_panels_path)
message("  wrote: ", panel_size_comparison_path)
message("  wrote: ", panel_search_summary_path)
message("  wrote: ", file.path(figure_dir, "recommended_provisional_panel_panels.png"))
message("  wrote: ", file.path(figure_dir, "selected_panel_panels.png"))
