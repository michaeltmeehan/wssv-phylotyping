#!/usr/bin/env Rscript

# Stage 08: map external partial FASTA records to alignment coordinates and
# classify them with the trained panel classifier and/or opportunistic SNP
# evidence, depending on config.
# Inputs: outputs/models/wssv_classifier.rds, data/processed/precomputed.rds,
# optional site_node_scores for opportunistic/auto mode, raw partial FASTA files,
# and analysis.partials settings, including reference_id, from config/config.yml.
# Outputs: partial classification, evidence, site-evidence, and mapping
# diagnostics tables under outputs/tables/.
# Run directly after scripts/07_train_classifier.R when partial genomes are
# available. It also writes empty well-formed tables if none are present.

source("R/encoding.R")
source("R/preprocess.R")
source("R/classifier.R")
source("R/partials.R")
source("R/script_utils.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")
model_dir <- file.path(config$paths$outputs, "models")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

safe_write_csv <- function(x, path) {
  ok <- tryCatch({
    write.csv(x, path, row.names = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (ok) {
    return(invisible(path))
  }
  fallback <- file.path(
    dirname(path),
    paste0(tools::file_path_sans_ext(basename(path)), "_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv")
  )
  write.csv(x, fallback, row.names = FALSE)
  warning("Could not write locked output ", path, "; wrote ", fallback, " instead.", call. = FALSE)
  invisible(fallback)
}

safe_save_rds <- function(x, path) {
  ok <- tryCatch({
    saveRDS(x, path)
    TRUE
  }, error = function(e) FALSE)
  if (ok) {
    return(invisible(path))
  }
  fallback <- file.path(
    dirname(path),
    paste0(tools::file_path_sans_ext(basename(path)), "_", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")
  )
  saveRDS(x, fallback)
  warning("Could not write locked output ", path, "; wrote ", fallback, " instead.", call. = FALSE)
  invisible(fallback)
}

precomputed_path <- file.path(processed_dir, "precomputed.rds")
classifier_path <- file.path(model_dir, "wssv_classifier.rds")
site_node_scores_path <- file.path(table_dir, "site_node_scores.rds")
if (!file.exists(precomputed_path)) {
  stop("Missing ", precomputed_path, ". Run scripts/01_preprocess_alignment_tree.R first.", call. = FALSE)
}
if (!file.exists(classifier_path)) {
  stop("Missing ", classifier_path, ". Run scripts/07_train_classifier.R first.", call. = FALSE)
}

pre <- readRDS(precomputed_path)
classifier <- load_classifier(classifier_path)

partial_cfg <- value_or(config$analysis$partials, list())
partial_dir <- value_or(partial_cfg$input_dir, value_or(config$paths$partials, "data/raw/partials/"))
partial_reference_id <- value_or(partial_cfg$reference_id, "CN01_1994")
extensions <- value_or(partial_cfg$fasta_extensions, partial_default_extensions)
mapping_mode <- value_or(partial_cfg$mapping_mode, "auto")
classification_mode <- value_or(partial_cfg$classification_mode, "panel")
pairwise_local_thresholds <- value_or(partial_cfg$pairwise_local, list())
min_observed_informative_sites <- as.integer(value_or(partial_cfg$min_observed_informative_sites, classifier$settings$min_total_informative_sites))
opportunistic_settings <- list(
  min_observed_informative_sites = min_observed_informative_sites,
  min_support = as.numeric(value_or(partial_cfg$min_support, classifier$settings$min_support)),
  min_support_margin = as.numeric(value_or(partial_cfg$min_support_margin, classifier$settings$support_margin)),
  conflict_support_threshold = as.numeric(value_or(partial_cfg$conflict_support_threshold, classifier$settings$max_conflict)),
  use_weighted_support = isTRUE(value_or(partial_cfg$use_weighted_support, TRUE))
)
write_unresolved_evidence <- isTRUE(value_or(partial_cfg$write_unresolved_evidence, FALSE))
site_node_scores <- NULL
if (classification_mode %in% c("opportunistic", "auto")) {
  if (!file.exists(site_node_scores_path)) {
    stop("Missing ", site_node_scores_path, ". Run scripts/02_score_snps.R first.", call. = FALSE)
  }
  site_node_scores <- readRDS(site_node_scores_path)
}

partials <- read_partial_fastas(partial_dir, extensions)
classification_path <- file.path(table_dir, "partial_classifications.csv")
evidence_path <- file.path(table_dir, "partial_node_evidence.csv")
panel_classification_path <- file.path(table_dir, "partial_classifications_panel.csv")
opportunistic_classification_path <- file.path(table_dir, "partial_classifications_opportunistic.csv")
opportunistic_evidence_path <- file.path(table_dir, "partial_opportunistic_node_evidence.csv")
opportunistic_site_evidence_path <- file.path(table_dir, "partial_opportunistic_site_evidence.csv")
opportunistic_site_summary_path <- file.path(table_dir, "partial_opportunistic_site_summary.csv")
region_diagnostics_path <- file.path(table_dir, "partial_region_diagnostics.csv")

if (nrow(partials) == 0L) {
  empty <- empty_partial_classification_table()
  classification_path <- safe_write_csv(empty, classification_path)
  panel_classification_path <- safe_write_csv(empty, panel_classification_path)
  opportunistic_classification_path <- safe_write_csv(empty, opportunistic_classification_path)
  opportunistic_evidence_path <- safe_write_csv(empty_partial_opportunistic_node_evidence_table(), opportunistic_evidence_path)
  opportunistic_site_evidence_path <- safe_write_csv(empty_partial_opportunistic_site_evidence_table(), opportunistic_site_evidence_path)
  opportunistic_site_summary_path <- safe_write_csv(empty_partial_opportunistic_site_summary_table(), opportunistic_site_summary_path)
  region_diagnostics_path <- safe_write_csv(empty_partial_region_diagnostics_table(), region_diagnostics_path)
  safe_save_rds(empty, file.path(table_dir, "partial_classifications.rds"))
  safe_save_rds(empty_partial_opportunistic_site_evidence_table(), file.path(table_dir, "partial_opportunistic_site_evidence.rds"))
  safe_save_rds(empty_partial_opportunistic_site_summary_table(), file.path(table_dir, "partial_opportunistic_site_summary.rds"))
  message("No partial FASTA files found in ", partial_dir)
  message("  partial records read: 0")
  message("  mapped records: 0")
  message("  unmapped records: 0")
  message("  resolved records: 0")
  message("  unresolved or conservative non-resolved records: 0")
  message("  wrote empty classification table: ", classification_path)
  quit(save = "no", status = 0L)
}

reference_sequence <- reference_sequence_from_alignment(pre$aln_int, reference_id = partial_reference_id)
mapped <- map_partial_sequences(
  partials,
  alignment_length = ncol(pre$aln_int),
  reference_sequence = reference_sequence,
  mapping_mode = mapping_mode,
  pairwise_local_thresholds = pairwise_local_thresholds
)

classified <- classify_mapped_partials(
  mapped,
  classifier,
  min_total_informative_sites = min_observed_informative_sites,
  write_unresolved_evidence = write_unresolved_evidence,
  classification_mode = classification_mode,
  site_node_scores = site_node_scores,
  opportunistic_settings = opportunistic_settings
)

classification_path <- safe_write_csv(classified$classifications, classification_path)
safe_save_rds(classified$classifications, file.path(table_dir, "partial_classifications.rds"))
panel_classification_path <- safe_write_csv(classified$panel_classifications, panel_classification_path)
safe_save_rds(classified$panel_classifications, file.path(table_dir, "partial_classifications_panel.rds"))
opportunistic_classification_path <- safe_write_csv(classified$opportunistic_classifications, opportunistic_classification_path)
safe_save_rds(classified$opportunistic_classifications, file.path(table_dir, "partial_classifications_opportunistic.rds"))
opportunistic_evidence_path <- safe_write_csv(classified$opportunistic_node_evidence, opportunistic_evidence_path)
safe_save_rds(classified$opportunistic_node_evidence, file.path(table_dir, "partial_opportunistic_node_evidence.rds"))
opportunistic_site_evidence_path <- safe_write_csv(classified$opportunistic_site_evidence, opportunistic_site_evidence_path)
safe_save_rds(classified$opportunistic_site_evidence, file.path(table_dir, "partial_opportunistic_site_evidence.rds"))
opportunistic_site_summary_path <- safe_write_csv(classified$opportunistic_site_summary, opportunistic_site_summary_path)
safe_save_rds(classified$opportunistic_site_summary, file.path(table_dir, "partial_opportunistic_site_summary.rds"))
region_diagnostics_path <- safe_write_csv(classified$region_diagnostics, region_diagnostics_path)
safe_save_rds(classified$region_diagnostics, file.path(table_dir, "partial_region_diagnostics.rds"))

if (nrow(classified$node_evidence) > 0L) {
  write.csv(classified$node_evidence, evidence_path, row.names = FALSE)
  saveRDS(classified$node_evidence, file.path(table_dir, "partial_node_evidence.rds"))
}

message("Partial classification complete")
message("  partial records read: ", nrow(partials))
message("  classification mode: ", classification_mode)
message("  mapped records: ", sum(mapped$mapped))
message("  panel-informative records: ", sum(classified$panel_classifications$mapped & classified$panel_classifications$observed_informative_sites > 0L, na.rm = TRUE))
message("  opportunistically informative records: ", sum(classified$opportunistic_classifications$mapped & classified$opportunistic_classifications$observed_informative_sites > 0L, na.rm = TRUE))
message("  resolved panel records: ", sum(classified$panel_classifications$status == "resolved"))
message("  resolved opportunistic records: ", sum(classified$opportunistic_classifications$status == "resolved_opportunistic"))
message("  final resolved records: ", sum(classified$classifications$status %in% c("resolved", "resolved_opportunistic")))
message("  weak records: ", sum(classified$classifications$status == "weak_support"))
message("  conflicting records: ", sum(classified$classifications$status == "conflicting"))
message("  no_informative_sites records: ", sum(classified$classifications$status == "no_informative_sites"))
message("  no_observed_informative_sites records: ", sum(classified$classifications$status == "no_observed_informative_sites"))
message("  unmapped records: ", sum(classified$classifications$status == "unmapped"))
message("  node evidence rows written: ", nrow(classified$node_evidence))
message("  opportunistic site evidence rows written: ", nrow(classified$opportunistic_site_evidence))
message("  classification table: ", classification_path)
message("  panel classification table: ", panel_classification_path)
message("  opportunistic classification table: ", opportunistic_classification_path)
message("  opportunistic node evidence table: ", opportunistic_evidence_path)
message("  opportunistic site evidence table: ", opportunistic_site_evidence_path)
message("  opportunistic site summary table: ", opportunistic_site_summary_path)
message("  region diagnostics table: ", region_diagnostics_path)
if (nrow(classified$node_evidence) > 0L) {
  message("  node evidence table: ", evidence_path)
}
