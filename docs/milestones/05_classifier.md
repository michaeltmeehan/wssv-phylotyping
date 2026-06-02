# Milestone 5: Classifier

## Goal

Train conservative probabilistic and tree-path classifiers for assigning complete or partial WSSV sequences to MCC-tree clades.

## Expected Inputs

- Preprocessed alignment/tree object.
- Selected marker panel.
- SNP scoring outputs.
- Validation design from Milestone 4.

## Expected Outputs

- Naive probabilistic classifier object.
- Tree-path classifier object.
- Classifier performance summary.
- Prediction output schema.

## Success Criteria

- Classifiers use only observed informative sites in a query sequence.
- The tree-path classifier stops when support is insufficient.
- Predictions include assigned clade, support, alternative clades, sites used, and status.

