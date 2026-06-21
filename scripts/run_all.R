#!/usr/bin/env Rscript

# Convenience wrapper: run the numbered workflow scripts in dependency order.
# Inputs: config/config.yml and all raw files referenced by that config. The
# config path is fixed because the numbered scripts read config/config.yml.
# Outputs: all standard processed objects, output tables, model files, partial
# classifications, and the final Markdown report produced by stages 01-09.
# Run directly from the repository root for a baseline end-to-end analysis. This
# wrapper overwrites/regenerates the normal pipeline outputs under configured
# paths; copy outputs first if you need to preserve a previous run.

config_path <- "config/config.yml"
if (!file.exists(config_path)) {
  stop("Missing required config file: ", config_path, call. = FALSE)
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read ", config_path, ".", call. = FALSE)
}
config <- yaml::read_yaml(config_path)
required_inputs <- c(config$paths$alignment, config$paths$tree)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input file(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}
message("Using config: ", config_path)
message("Normal outputs under ", config$paths$outputs, " and ", config$paths$processed, " may be overwritten/regenerated.")

stages <- file.path(
  "scripts",
  c(
    "01_preprocess_alignment_tree.R",
    "02_score_snps.R",
    "03_summarise_sites.R",
    "04_score_windows.R",
    "05_select_marker_panel.R",
    "06_validate_panel.R",
    "07_train_classifier.R",
    "08_classify_partials.R",
    "09_generate_report.R"
  )
)

for (stage in stages) {
  if (!file.exists(stage)) {
    stop("Missing workflow stage: ", stage, call. = FALSE)
  }
  message("\n== Running ", stage, " ==")
  status <- system2(file.path(R.home("bin"), "Rscript"), stage)
  if (!identical(status, 0L)) {
    stop("Stage failed: ", stage, call. = FALSE)
  }
}

message("\nWorkflow complete.")
