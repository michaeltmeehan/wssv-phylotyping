# Workflow Scripts

The `scripts/` directory contains numbered entry-point scripts for the staged workflow:

1. `01_preprocess_alignment_tree.R`: prepare alignment, tree, coordinate, and clade objects.
2. `02_score_snps.R`: score SNPs against MCC-tree clades.
3. `03_summarise_sites.R`: summarize and rank informative sites.
4. `04_score_windows.R`: aggregate SNP signal into genomic windows.
5. `05_select_marker_panel.R`: choose complementary marker windows.
6. `06_validate_panel.R`: validate reduced-region placement.
7. `07_train_classifier.R`: train probabilistic and tree-path classifiers.
8. `08_classify_partials.R`: classify partial genomes where signal is sufficient.

Milestone 1 scripts are runnable from the repository root:

- `01_preprocess_alignment_tree.R` reads configured raw inputs, matches tree
  tips to alignment names, encodes bases, derives eligible node target masks,
  and writes `data/processed/precomputed.rds`.
- `02_score_snps.R` scores polymorphic sites against eligible clades and writes
  `outputs/tables/site_node_scores.rds`.
- `03_summarise_sites.R` writes `outputs/tables/site_summary.csv` and
  `outputs/tables/node_summary.csv`.

Later-stage scripts remain placeholders until their milestones begin.
