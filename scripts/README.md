# Workflow Scripts

The `scripts/` directory contains numbered entry-point scripts for the staged workflow:

1. `01_preprocess_alignment_tree.R`: prepare alignment, tree, coordinate, and clade objects.
2. `02_score_snps.R`: score SNPs against MCC-tree clades.
3. `03_summarise_sites.R`: summarize and rank informative sites.
4. `04_score_windows.R`: aggregate SNP signal into genomic windows.
5. `05_select_marker_panel.R`: choose complementary marker windows.
6. `06_validate_panel.R`: validate reduced-region placement.
7. `07_train_classifier.R`: train the conservative tree-path classifier.
8. `08_classify_partials.R`: classify partial genomes where signal is sufficient.

Milestone scripts are runnable from the repository root:

- `01_preprocess_alignment_tree.R` reads configured raw inputs, matches tree
  tips to alignment names, encodes bases, derives eligible node target masks,
  and writes `data/processed/precomputed.rds`.
- `02_score_snps.R` scores polymorphic sites against eligible clades and writes
  `outputs/tables/site_node_scores.rds`.
- `03_summarise_sites.R` writes `outputs/tables/site_summary.csv` and
  `outputs/tables/node_summary.csv`.
- `04_score_windows.R` generates fixed and optional SNP-centred windows from
  configured widths, aggregates site-node scores, and writes
  `outputs/tables/candidate_windows.csv`,
  `outputs/tables/window_summary.csv`, and
  `outputs/tables/window_node_summary.csv`.
- `05_select_marker_panel.R` filters scored windows, greedily selects
  complementary non-redundant windows, compares simple weighted-gain baselines,
  and writes `outputs/tables/selected_panel.csv`,
  `outputs/tables/selected_panel_steps.csv`,
  `outputs/tables/panel_node_coverage.csv`, and
  `outputs/tables/panel_summary.csv`.
- `06_validate_panel.R` evaluates whether selected panel regions retain
  informative complete-genome signal, compares top-weighted and random matched
  window baselines, summarizes random fixed-length artificial fragments, and
  writes `outputs/tables/validation_panel_summary.csv`,
  `outputs/tables/validation_tip_summary.csv`,
  `outputs/tables/validation_fragment_summary.csv`, and
  `outputs/tables/validation_baseline_summary.csv`.
- `07_train_classifier.R` trains an interpretable tree-path classifier from the
  selected marker panel or full informative-site catalogue, performs
  training-genome internal sanity checks, and writes
  `outputs/models/wssv_classifier.rds`,
  `outputs/tables/classifier_training_summary.csv`,
  `outputs/tables/classifier_tip_predictions.csv`, and, when compact enough,
  `outputs/tables/classifier_node_evidence.csv`. These predictions are
  training-set checks, not independent validation.
- `08_classify_partials.R` reads external partial FASTA files from the configured
  raw partial input directory, conservatively maps aligned full-length records
  or exact reference substrings to alignment coordinates, classifies mapped
  records with the trained classifier, and writes
  `outputs/tables/partial_classifications.csv` plus
  `outputs/tables/partial_node_evidence.csv` when node evidence is available.
  If no partial FASTA files are present, it writes an empty well-formed
  classification table.
