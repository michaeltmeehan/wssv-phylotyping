#!/usr/bin/env Rscript

source("R/encoding.R")
source("R/classifier.R")
source("R/script_utils.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed
table_dir <- file.path(config$paths$outputs, "tables")
model_dir <- file.path(config$paths$outputs, "models")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

precomputed_path <- file.path(processed_dir, "precomputed.rds")
if (!file.exists(precomputed_path)) {
  stop("Missing ", precomputed_path, ". Run scripts/01_preprocess_alignment_tree.R first.", call. = FALSE)
}

pre <- readRDS(precomputed_path)
site_node_scores <- read_table_output(table_dir, "site_node_scores", "scripts/02_score_snps.R")
selected_panel <- read_table_output(table_dir, "selected_panel", "scripts/05_select_marker_panel.R")

classifier_cfg <- config$analysis$classifier
if (is.null(classifier_cfg)) {
  classifier_cfg <- list()
}

use_selected_panel <- isTRUE(value_or(classifier_cfg$use_selected_panel, TRUE))
panel_size <- classifier_cfg$panel_size
if (is.null(panel_size) || is.na(panel_size)) {
  panel_size <- nrow(selected_panel)
}

message("Training tree-path classifier")
message("  training tips: ", nrow(pre$aln_int))
message("  use selected panel: ", use_selected_panel)
message("  panel windows used: ", if (use_selected_panel) min(as.integer(panel_size), nrow(selected_panel)) else 0L)

classifier <- train_classifier(
  aln_int = pre$aln_int,
  target_mask = pre$target_mask,
  node_metadata = pre$node_metadata,
  site_node_scores = site_node_scores,
  selected_panel = selected_panel,
  use_selected_panel = use_selected_panel,
  panel_size = as.integer(panel_size),
  min_sites_per_node = as.integer(value_or(classifier_cfg$min_sites_per_node, 1L)),
  min_total_informative_sites = as.integer(value_or(classifier_cfg$min_total_informative_sites, 1L)),
  min_support = as.numeric(value_or(classifier_cfg$min_support, 0.8)),
  max_conflict = as.numeric(value_or(classifier_cfg$max_conflict, 0.2)),
  support_margin = as.numeric(value_or(classifier_cfg$support_margin, 0.05))
)

predictions <- vector("list", nrow(pre$aln_int))
evidence_rows <- vector("list", nrow(pre$aln_int))
for (i in seq_len(nrow(pre$aln_int))) {
  tip <- rownames(pre$aln_int)[[i]]
  pred <- classify_tree_path(pre$aln_int[i, ], classifier)
  predictions[[i]] <- classifier_prediction_row(pred, tip)
  ev <- pred$evidence
  ev$query_id <- tip
  evidence_rows[[i]] <- ev
}
tip_predictions <- do.call(rbind, predictions)
node_evidence <- do.call(rbind, evidence_rows)
node_evidence <- node_evidence[, c("query_id", setdiff(names(node_evidence), "query_id")), drop = FALSE]

training_summary <- data.frame(
  check_type = "training_set_internal_sanity_check",
  tips_checked = nrow(tip_predictions),
  resolved = sum(tip_predictions$status == "resolved"),
  weak_support = sum(tip_predictions$status == "weak_support"),
  no_informative_sites = sum(tip_predictions$status == "no_informative_sites"),
  conflicting = sum(tip_predictions$status == "conflicting"),
  unresolved = sum(tip_predictions$status == "unresolved"),
  mean_observed_informative_sites = mean(tip_predictions$observed_informative_sites),
  median_observed_informative_sites = stats::median(tip_predictions$observed_informative_sites),
  mean_informative_nodes_evaluated = mean(tip_predictions$informative_nodes_evaluated),
  classifier_nodes = nrow(classifier$node_table),
  classifier_rules = nrow(classifier$rules),
  classifier_informative_sites = length(classifier$informative_sites),
  independent_validation = FALSE
)

classifier_path <- file.path(model_dir, "wssv_classifier.rds")
save_classifier(classifier, classifier_path)
write.csv(training_summary, file.path(table_dir, "classifier_training_summary.csv"), row.names = FALSE)
write.csv(tip_predictions, file.path(table_dir, "classifier_tip_predictions.csv"), row.names = FALSE)
if (nrow(node_evidence) <= 1000000L) {
  write.csv(node_evidence, file.path(table_dir, "classifier_node_evidence.csv"), row.names = FALSE)
}

message("Classifier training complete")
message("  wrote model: ", classifier_path)
message("  candidate node-site score rows read: ", nrow(site_node_scores))
message("  rules: ", nrow(classifier$rules))
message("  classifier nodes retained: ", nrow(classifier$node_table))
message("  informative sites: ", length(classifier$informative_sites))
message("  resolved training-set sanity checks: ", training_summary$resolved, " of ", training_summary$tips_checked)
message("  note: predictions are internal training-set checks, not independent validation")
