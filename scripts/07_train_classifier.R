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

training_tip_diagnostics <- function(tip, prediction, target_mask, node_metadata, min_support, max_conflict) {
  tip_col <- match(tip, colnames(target_mask))
  if (is.na(tip_col)) {
    stop("Tip not found in target mask: ", tip, call. = FALSE)
  }
  metadata <- node_metadata
  metadata$node_id <- as.character(metadata$node_id)
  true_ids <- as.character(rownames(target_mask)[target_mask[, tip_col]])
  true_meta <- metadata[metadata$node_id %in% true_ids, , drop = FALSE]
  true_meta <- true_meta[order(true_meta$depth, true_meta$clade_size, true_meta$node_id), , drop = FALSE]
  true_path <- true_meta$node_id

  evidence <- prediction$evidence
  evidence$node_id <- as.character(evidence$node_id)
  supported <- supported_evidence(evidence, min_support, max_conflict)
  on_path <- supported[supported$node_id %in% true_path, , drop = FALSE]
  off_path <- supported[!supported$node_id %in% true_path, , drop = FALSE]
  deepest <- if (nrow(on_path) == 0L) NA_character_ else {
    on_path[order(-on_path$depth, -on_path$support, on_path$conflict), "node_id"][[1L]]
  }
  assigned <- prediction$assigned_node
  assigned_on_path <- !is.na(assigned) && assigned %in% true_path
  strongest_on <- evidence[evidence$node_id %in% true_path & !is.na(evidence$support), , drop = FALSE]
  strongest_off <- evidence[!evidence$node_id %in% true_path & !is.na(evidence$support), , drop = FALSE]
  strongest_on <- if (nrow(strongest_on) == 0L) strongest_on else strongest_on[order(-strongest_on$support, strongest_on$conflict, -strongest_on$depth), , drop = FALSE]
  strongest_off <- if (nrow(strongest_off) == 0L) strongest_off else strongest_off[order(-strongest_off$support, strongest_off$conflict, -strongest_off$depth), , drop = FALSE]
  off_reason <- if (nrow(off_path) == 0L) {
    NA_character_
  } else if (nrow(on_path) == 0L) {
    "strong_off_path_support_without_true_path_support"
  } else {
    "strong_off_path_support_competes_with_true_path"
  }
  conflict_reason <- prediction$conflict_reason
  if (prediction$status == "conflicting" && !is.na(off_reason)) {
    conflict_reason <- off_reason
  } else if (!is.na(assigned) && !assigned_on_path) {
    conflict_reason <- "assigned_off_true_path_nested_offshoot_support"
  }

  data.frame(
    query_id = tip,
    true_mcc_path = paste(true_path, collapse = ";"),
    true_mcc_path_depth = length(true_path),
    deepest_supported_true_path_node = deepest,
    assigned_node = assigned,
    assigned_on_true_path = assigned_on_path,
    strongest_on_path_node = if (nrow(strongest_on) > 0L) strongest_on$node_id[[1L]] else NA_character_,
    strongest_on_path_support = if (nrow(strongest_on) > 0L) strongest_on$support[[1L]] else NA_real_,
    strongest_off_path_node = if (nrow(strongest_off) > 0L) strongest_off$node_id[[1L]] else NA_character_,
    strongest_off_path_support = if (nrow(strongest_off) > 0L) strongest_off$support[[1L]] else NA_real_,
    on_path_nodes_supported = nrow(on_path),
    off_path_nodes_supported = nrow(off_path),
    conflict_reason = conflict_reason
  )
}

predictions <- vector("list", nrow(pre$aln_int))
evidence_rows <- vector("list", nrow(pre$aln_int))
diagnostic_rows <- vector("list", nrow(pre$aln_int))
for (i in seq_len(nrow(pre$aln_int))) {
  tip <- rownames(pre$aln_int)[[i]]
  pred <- classify_tree_path(pre$aln_int[i, ], classifier)
  predictions[[i]] <- classifier_prediction_row(pred, tip)
  diagnostic_rows[[i]] <- training_tip_diagnostics(
    tip,
    pred,
    pre$target_mask,
    pre$node_metadata,
    classifier$settings$min_support,
    classifier$settings$max_conflict
  )
  ev <- pred$evidence
  ev$query_id <- tip
  evidence_rows[[i]] <- ev
}
tip_predictions <- do.call(rbind, predictions)
tip_diagnostics <- do.call(rbind, diagnostic_rows)
tip_predictions <- merge(
  tip_predictions,
  tip_diagnostics[, setdiff(names(tip_diagnostics), c("assigned_node", "conflict_reason")), drop = FALSE],
  by = "query_id",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_diagnostic")
)
tip_predictions$conflict_reason <- tip_diagnostics$conflict_reason[match(tip_predictions$query_id, tip_diagnostics$query_id)]
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
  assigned_on_true_path = sum(tip_predictions$assigned_on_true_path, na.rm = TRUE),
  assigned_off_true_path = sum(!tip_predictions$assigned_on_true_path & !is.na(tip_predictions$assigned_node), na.rm = TRUE),
  mean_on_path_nodes_supported = mean(tip_predictions$on_path_nodes_supported),
  mean_off_path_nodes_supported = mean(tip_predictions$off_path_nodes_supported),
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
