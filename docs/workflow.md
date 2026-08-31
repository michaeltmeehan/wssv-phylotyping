# Workflow

This pipeline turns a WSSV whole-genome alignment and MCC tree into scored SNPs, candidate marker windows, selected marker panels, validation summaries, a trained classifier, optional partial-genome classifications, and a report.

## Pipeline Diagram

```text
data/raw/alignment/*.fasta + data/raw/tree/*
        |
        v
01_preprocess_alignment_tree.R
        |
        v
data/processed/precomputed.rds
        |
        v
02_score_snps.R -> site_node_scores.rds
        |
        v
03_summarise_sites.R -> site_summary.csv, node_summary.csv
        |
        v
04_score_windows.R -> candidate_windows.csv, window_summary.csv, window_node_summary.csv
        |
        v
05_select_marker_panel.R -> selected_panel.csv, panel_summary.csv
        |
        +--> 06_validate_panel.R -> validation_*.csv
        |
        v
07_train_classifier.R -> outputs/models/wssv_classifier.rds, classifier_*.csv
        |
        v
08_classify_partials.R + data/raw/partials/*.fasta -> partial_*.csv
        |
        v
09_generate_report.R -> outputs/reports/wssv_phylotyping_report.md
```

## Stages

| Stage | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| `scripts/01_preprocess_alignment_tree.R` | Reads the alignment and tree, keeps matching tips, encodes bases, finds polymorphic sites, builds the broad diagnostic clade masks, and resolves the explicitly configured target clades. | `config/config.yml`, `paths.alignment`, `paths.tree` | `data/processed/precomputed.rds` with `all_clade_mask`, `target_clade_mask`, `all_node_metadata`, and `target_node_metadata` |
| `scripts/02_score_snps.R` | Scores each polymorphic SNP against each eligible clade. | `precomputed.rds`, `analysis.min_*`, `analysis.chunk_size` | `outputs/tables/site_node_scores.rds` |
| `scripts/03_summarise_sites.R` | Makes student-readable SNP and node summaries. | `precomputed.rds`, `site_node_scores.rds` | `site_summary.csv`, `node_summary.csv` |
| `scripts/04_score_windows.R` | Creates fixed and SNP-centred windows and aggregates SNP scores into window scores. | `site_node_scores.rds`, `analysis.windows` | `candidate_windows.csv`, `window_summary.csv`, `window_node_summary.csv` |
| `scripts/05_select_marker_panel.R` | Filters windows and greedily selects complementary marker windows. | `window_summary.*`, `window_node_summary.*`, `analysis.panel_selection` | `selected_panel.csv`, `selected_panel_steps.csv`, `panel_node_coverage.csv`, `panel_summary.csv` |
| `scripts/06_validate_panel.R` | Checks complete-genome signal retained by selected panels and random/baseline fragments. | `selected_panel.*`, `site_node_scores.*`, `analysis.validation` | `validation_panel_summary.csv`, `validation_tip_summary.csv`, `validation_fragment_summary.csv`, `validation_baseline_summary.csv` |
| `scripts/07_train_classifier.R` | Builds an interpretable tree-path classifier and runs training-tip sanity checks. | `selected_panel.*`, `site_node_scores.*`, `analysis.classifier` | `outputs/models/wssv_classifier.rds`, `classifier_training_summary.csv`, `classifier_tip_predictions.csv`, optional `classifier_node_evidence.csv` |
| `scripts/08_classify_partials.R` | Reads partial FASTA files, maps them to alignment coordinates, and classifies mapped records. | `wssv_classifier.rds`, `data/raw/partials/`, `analysis.partials` | `partial_classifications.csv`, panel/opportunistic classification tables, evidence tables, region diagnostics |
| `scripts/09_generate_report.R` | Summarises existing outputs into a Markdown report. It does not rerun analysis. | Existing output tables and `config/config.yml` | `outputs/reports/wssv_phylotyping_report.md`, `report_*.csv` extracts |

## Running The Workflow

From the repository root:

```sh
Rscript scripts/run_all.R
```

For careful debugging or parameter tuning, run scripts one at a time. If a downstream script reports a missing file, run the upstream stage named in the error message.

## Important Rules

- `data/raw/` is input storage. Do not generate or modify raw data there.
- `data/processed/` and `outputs/` are generated and can be regenerated from the raw inputs and config.
- The report in `outputs/reports/` is a summary of existing outputs, not a substitute for inspecting tables.
