# R Function Files

The `R/` directory will hold reusable functions shared by the numbered workflow scripts.

- `preprocess.R`: alignment/tree loading, matching, pruning, and coordinate preparation.
- `tree_utils.R`: MCC-tree node, posterior support, broad diagnostic clade, target clade, depth, balance, legacy weight, and mask helpers.
- `encoding.R`: alignment and allele encoding helpers.
- `snp_scoring.R`: per-site/per-node SNP gain scoring and summary helpers.
- `window_scoring.R`: candidate window generation and scoring helpers.
- `panel_selection.R`: complementary marker panel selection routines.
- `validation.R`: selected-window extraction, complete-genome masking,
  synthetic fragment generation, observed signal counting, and baseline
  panel-summary helpers.
- `classifier.R`: conservative tree-path classifier training, evidence, prediction,
  and save/load helpers.
- `partials.R`: partial FASTA ingestion, conservative coordinate mapping,
  optional local-alignment mapping, partial encoding, panel and opportunistic
  partial-genome classification, interval diagnostics, and per-node evidence
  helpers.
- `script_utils.R`: small shared helpers used by numbered workflow scripts.
- `plotting.R`: shared plotting helpers.

Milestones 0-6 are implemented through conservative partial-genome
classification with explicit unmapped, weak, conflicting, unresolved,
no-observed-signal, and no-informative-site outcomes. Partial classification
supports panel, opportunistic, and auto modes; opportunistic classifications use
all scored informative SNPs within each mapped interval and are reported
separately from designed-panel classifier calls. The `weight` column in node
metadata is kept as a legacy compatibility field for now; use
`posterior_support` for tree credibility and SNP `gain` / `normalized_gain` for
stage-02 marker scoring.
