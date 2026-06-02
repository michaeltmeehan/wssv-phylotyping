# R Function Files

The `R/` directory will hold reusable functions shared by the numbered workflow scripts.

- `preprocess.R`: alignment/tree loading, matching, pruning, and coordinate preparation.
- `tree_utils.R`: MCC-tree node, clade, depth, balance, weight, and target-mask helpers.
- `encoding.R`: alignment and allele encoding helpers.
- `snp_scoring.R`: per-site/per-node SNP gain scoring and summary helpers.
- `window_scoring.R`: candidate window generation and scoring helpers.
- `panel_selection.R`: complementary marker panel selection routines.
- `validation.R`: leave-one-out, synthetic partial, and tree-recovery validation helpers.
- `classifier.R`: classifier training and prediction helpers.
- `partials.R`: partial-sequence mapping and extraction helpers.
- `plotting.R`: shared plotting helpers.

Milestone 1 logic is implemented for preprocessing and SNP scoring only.
Later-stage files remain placeholders until their milestones begin.
