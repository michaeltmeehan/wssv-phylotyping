#!/usr/bin/env Rscript

source("R/encoding.R")
source("R/preprocess.R")
source("R/classifier.R")
source("R/partials.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

value_or <- function(x, default) {
  if (is.null(x)) default else x
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")
model_dir <- file.path(config$paths$outputs, "models")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

precomputed_path <- file.path(processed_dir, "precomputed.rds")
classifier_path <- file.path(model_dir, "wssv_classifier.rds")
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
extensions <- value_or(partial_cfg$fasta_extensions, partial_default_extensions)
mapping_mode <- value_or(partial_cfg$mapping_mode, "auto")
min_observed_informative_sites <- as.integer(value_or(partial_cfg$min_observed_informative_sites, classifier$settings$min_total_informative_sites))
write_unresolved_evidence <- isTRUE(value_or(partial_cfg$write_unresolved_evidence, FALSE))

partials <- read_partial_fastas(partial_dir, extensions)
classification_path <- file.path(table_dir, "partial_classifications.csv")
evidence_path <- file.path(table_dir, "partial_node_evidence.csv")

if (nrow(partials) == 0L) {
  empty <- empty_partial_classification_table()
  write.csv(empty, classification_path, row.names = FALSE)
  saveRDS(empty, file.path(table_dir, "partial_classifications.rds"))
  message("No partial FASTA files found in ", partial_dir)
  message("  wrote empty classification table: ", classification_path)
  quit(save = "no", status = 0L)
}

reference_sequence <- reference_sequence_from_alignment(pre$aln_int)
mapped <- map_partial_sequences(
  partials,
  alignment_length = ncol(pre$aln_int),
  reference_sequence = reference_sequence,
  mapping_mode = mapping_mode
)

classified <- classify_mapped_partials(
  mapped,
  classifier,
  min_total_informative_sites = min_observed_informative_sites,
  write_unresolved_evidence = write_unresolved_evidence
)

write.csv(classified$classifications, classification_path, row.names = FALSE)
saveRDS(classified$classifications, file.path(table_dir, "partial_classifications.rds"))

if (nrow(classified$node_evidence) > 0L) {
  write.csv(classified$node_evidence, evidence_path, row.names = FALSE)
  saveRDS(classified$node_evidence, file.path(table_dir, "partial_node_evidence.rds"))
}

message("Partial classification complete")
message("  partial records read: ", nrow(partials))
message("  mapped records: ", sum(mapped$mapped))
message("  classification table: ", classification_path)
if (nrow(classified$node_evidence) > 0L) {
  message("  node evidence table: ", evidence_path)
}
