# Workflow Scripts

The `scripts/` directory contains numbered entry-point scripts for the staged workflow:

0. `run_all.R`: convenience wrapper that runs the numbered workflow from `01` through `09`.
1. `01_preprocess_alignment_tree.R`: prepare alignment, tree, coordinate, and clade objects.
2. `02_score_snps.R`: score SNPs against MCC-tree clades.
3. `03_summarise_sites.R`: summarize and rank informative sites.
4. `04_score_windows.R`: aggregate SNP signal into genomic windows.
5. `05_select_panel.R`: compare small candidate amplicon panels from the stage-04 outputs.
6. `06_validate_panel.R`: validate reduced-region placement.
7. `07_train_classifier.R`: train the conservative tree-path classifier.
8. `08_classify_partials.R`: classify partial genomes where signal is sufficient.
9. `09_generate_report.R`: generate a lightweight exploratory Markdown report from existing outputs.

Milestone scripts are runnable from the repository root:

- `run_all.R` executes the full baseline workflow in dependency order. Use this
  when you want a single command for a configured run:
  `Rscript scripts/run_all.R`.
- `01_preprocess_alignment_tree.R` reads configured raw inputs, matches tree
  tips to alignment names, encodes bases, derives the full diagnostic clade set
  plus the explicitly configured target clades, and writes
  `data/processed/precomputed.rds`.
- `02_score_snps.R` scores polymorphic sites twice: once against the full
  diagnostic clade set and once against the restricted target clade set. It
  writes `outputs/tables/site_node_scores_all.rds`,
  `outputs/tables/site_node_scores_targets.rds`, and the transitional alias
  `outputs/tables/site_node_scores.rds` for the full diagnostic set. Future
  primer-panel design should use the target-specific output.
- `03_summarise_sites.R` writes
  `outputs/tables/site_summary_all.csv`,
  `outputs/tables/node_summary_all.csv`,
  `outputs/tables/site_summary_targets.csv`,
  `outputs/tables/node_summary_targets.csv`, the target-clade checkpoint tables
  `outputs/tables/target_clade_diagnostics.csv` and
  `outputs/tables/target_clade_strongest_snps.csv`, plus the legacy
  compatibility aliases `outputs/tables/site_summary.csv` and
  `outputs/tables/node_summary.csv` for the full diagnostic summaries.
- `04_score_windows.R` generates candidate amplicons from informative target
  SNP geography, summarises per-target evidence, and writes
  `outputs/tables/candidate_amplicons.csv`,
  `outputs/tables/candidate_amplicon_target_scores.csv`,
  `outputs/tables/candidate_amplicon_snps.csv`, and the spatial diagnostic
  figure `outputs/figures/candidate_amplicons_spatial.png`.
- `05_select_panel.R` compares candidate amplicon panels using transparent
  multi-objective criteria, evaluates each allowed panel size independently,
  retains within-size and cross-size Pareto sets, and writes
  `outputs/tables/panel_candidates.csv`,
  `outputs/tables/panel_target_scores.csv`,
  `outputs/tables/pareto_panels_by_size.csv`,
  `outputs/tables/cross_size_pareto_panels.csv`,
  `outputs/tables/recommended_provisional_panel.csv`,
  `outputs/tables/recommended_panel_alternatives.csv`,
  `outputs/tables/panel_size_comparison.csv`,
  `outputs/tables/panel_summary.csv`, and
  compatibility aliases `outputs/tables/selected_panel.csv` and
  `outputs/tables/pareto_panels.csv`, plus
  `outputs/figures/recommended_provisional_panel_panels.png`.
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
  raw partial input directory, conservatively maps aligned full-length records,
  exact reference substrings, or accepted local alignments to alignment
  coordinates, and classifies mapped records in configured `panel`,
  `opportunistic`, or `auto` mode. Panel mode preserves the trained selected
  panel classifier; opportunistic mode uses all scored informative SNPs in each
  mapped interval; auto mode falls back to opportunistic mode only when panel
  sites are not observed. It writes `outputs/tables/partial_classifications.csv`,
  separate panel/opportunistic classification tables, opportunistic node
  evidence, and interval diagnostics. If no partial FASTA files are present, it
  writes empty well-formed classification tables.
- `09_generate_report.R` reads existing processed objects and output tables,
  writes `outputs/reports/wssv_phylotyping_report.md`, and writes optional
  compact report extract tables under `outputs/tables/`. Missing upstream
  outputs are reported in the Markdown instead of causing avoidable failures.
