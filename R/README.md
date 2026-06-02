# R Function Files

The `R/` directory will hold reusable functions shared by the numbered workflow scripts.

- `preprocess.R`: alignment/tree loading, matching, pruning, and coordinate preparation.
- `tree_utils.R`: MCC-tree node, clade, depth, balance, weight, and target-mask helpers.
- `encoding.R`: alignment and allele encoding helpers.
- `snp_scoring.R`: per-site/per-node SNP gain scoring and summary helpers.
- `window_scoring.R`: candidate window generation and scoring helpers.
- `panel_selection.R`: complementary marker panel selection routines.
- `validation.R`: selected-window extraction, complete-genome masking,
  synthetic fragment generation, observed signal counting, and baseline
  panel-summary helpers.
- `classifier.R`: conservative tree-path classifier training, evidence, prediction,
  and save/load helpers.
- `partials.R`: partial-sequence mapping and extraction helpers.
- `plotting.R`: shared plotting helpers.

Milestones 1-5 are implemented through conservative classifier training and
training-set sanity checks. Partial-genome ingestion/classification remains a
placeholder until its milestone begins.
