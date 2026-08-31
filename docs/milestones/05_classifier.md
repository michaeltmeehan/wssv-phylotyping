# Milestone 5: Classifier

## Goal

Train a conservative tree-path classifier for assigning complete or partial
encoded WSSV sequences to MCC-tree clades using only observed informative sites.

## Expected Inputs

- Preprocessed alignment/tree object.
- Selected marker panel.
- SNP scoring outputs.
- Validation design from Milestone 4.

## Expected Outputs

- `outputs/models/wssv_classifier.rds`.
- `outputs/tables/classifier_training_summary.csv`.
- `outputs/tables/classifier_tip_predictions.csv`.
- `outputs/tables/classifier_node_evidence.csv` when the table is compact enough.
- Prediction diagnostics including assigned node, depth, observed informative
  sites, nodes evaluated, classifier support, competitor support, status, and
  per-node evidence.

## Success Criteria

- Classifiers use only observed informative sites in a query sequence.
- The tree-path classifier stops when classifier support is insufficient.
- Predictions include assigned clade, classifier support, alternative clades, sites used, and status.
- Training-set predictions are labelled as internal sanity checks, not
  independent validation.
