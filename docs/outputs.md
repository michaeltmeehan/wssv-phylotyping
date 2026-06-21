# Outputs Catalogue

Generated files are written under `data/processed/` and `outputs/`. They are safe to delete and regenerate if the raw inputs and config are preserved. Do not delete raw files under `data/raw/`.

| Output | Created by | Contains | Intermediate or final | When to inspect | Safe to delete/regenerate |
|---|---|---|---|---|---|
| `data/processed/precomputed.rds` | Stage 01 | Encoded alignment, pruned tree, clade masks, node metadata, polymorphic sites. | Intermediate | After changing raw alignment/tree or clade filters. | Yes |
| `outputs/tables/site_node_scores.rds` | Stage 02 | Informative SNP rules by site and node. | Intermediate | After changing SNP scoring thresholds. | Yes |
| `outputs/tables/site_summary.csv` | Stage 03 | Per-site score summaries. | Diagnostic/final summary | To find strongest SNPs and signal concentration. | Yes |
| `outputs/tables/node_summary.csv` | Stage 03 | Per-node score summaries. | Diagnostic/final summary | To see which clades have informative SNPs. | Yes |
| `outputs/tables/candidate_windows.csv` | Stage 04 | All generated candidate windows. | Intermediate | To check window count, widths, and coordinates. | Yes |
| `outputs/tables/window_summary.csv` | Stage 04 | Window-level signal, missingness, and node coverage summaries. | Diagnostic/intermediate | Before panel selection and when tuning windows. | Yes |
| `outputs/tables/window_node_summary.csv` | Stage 04 | Window-by-node score summaries. | Intermediate | When diagnosing node coverage. | Yes |
| `outputs/tables/selected_panel.csv` | Stage 05 | Greedy selected marker windows in order. | Final candidate panel | Always inspect after panel tuning. | Yes |
| `outputs/tables/selected_panel_steps.csv` | Stage 05 | Marginal gain and cumulative gain for each selected step. | Diagnostic | To check whether later windows add real signal. | Yes |
| `outputs/tables/panel_node_coverage.csv` | Stage 05 | Node coverage for selected panel sizes. | Diagnostic | To identify covered and uncovered nodes. | Yes |
| `outputs/tables/panel_summary.csv` | Stage 05 | Panel-size summaries and alternative ranking comparisons. | Final summary | To compare panel sizes and selection settings. | Yes |
| `outputs/tables/validation_panel_summary.csv` | Stage 06 | Complete-genome validation summaries by panel size. | Final diagnostic | To compare parameter sets. | Yes |
| `outputs/tables/validation_tip_summary.csv` | Stage 06 | Per-tip retained signal summaries. | Diagnostic | To find genomes with poor retained signal. | Yes |
| `outputs/tables/validation_fragment_summary.csv` | Stage 06 | Random fixed-fragment signal summaries. | Diagnostic baseline | To compare selected panels to generic fragments. | Yes |
| `outputs/tables/validation_baseline_summary.csv` | Stage 06 | Random matched-window baseline summaries. | Diagnostic baseline | To check selected panel benefit over random windows. | Yes |
| `outputs/models/wssv_classifier.rds` | Stage 07 | Saved classifier object. | Final model for current config | Before partial classification; archive with experiment outputs if needed. | Yes, but rerun stage 07 before classifying. |
| `outputs/tables/classifier_training_summary.csv` | Stage 07 | Training-tip status counts and evidence means. | Diagnostic | To check classifier strictness. | Yes |
| `outputs/tables/classifier_tip_predictions.csv` | Stage 07 | Training-tip predictions and path diagnostics. | Diagnostic | To inspect resolved, weak, and conflicting complete-genome calls. | Yes |
| `outputs/tables/classifier_node_evidence.csv` | Stage 07 | Node evidence for training tips when compact enough. | Diagnostic | To investigate specific calls. | Yes |
| `outputs/tables/partial_classifications.csv` | Stage 08 | Final partial-genome classification table. | Final output for partials | Always inspect when partials are used. | Yes |
| `outputs/tables/partial_classifications_panel.csv` | Stage 08 | Panel-only partial classifications. | Diagnostic/final | To separate designed-panel results from fallback results. | Yes |
| `outputs/tables/partial_classifications_opportunistic.csv` | Stage 08 | Opportunistic partial classifications. | Diagnostic/final | To interpret auto/opportunistic calls cautiously. | Yes |
| `outputs/tables/partial_node_evidence.csv` | Stage 08 | Panel node evidence for partial records. | Diagnostic | To inspect weak/conflicting/resolved partial calls. | Yes |
| `outputs/tables/partial_opportunistic_node_evidence.csv` | Stage 08 | Opportunistic node evidence. | Diagnostic | To audit opportunistic assignments. | Yes |
| `outputs/tables/partial_opportunistic_site_evidence.csv` | Stage 08 | Site-by-node opportunistic evidence. | Diagnostic | To see which SNPs drove opportunistic calls. | Yes |
| `outputs/tables/partial_opportunistic_site_summary.csv` | Stage 08 | Per-record opportunistic site summary. | Diagnostic | To detect off-path support or sparse evidence. | Yes |
| `outputs/tables/partial_region_diagnostics.csv` | Stage 08 | Partial mapping interval and observed-site diagnostics. | Diagnostic | First place to check for partial problems. | Yes |
| `outputs/reports/wssv_phylotyping_report.md` | Stage 09 | Human-readable summary of existing outputs. | Final report summary | After a baseline or experiment run. | Yes |
| `outputs/tables/report_*.csv` | Stage 09 | Compact report extracts such as top sites/windows. | Convenience outputs | When reviewing report inputs. | Yes |

Most tables are also written as `.rds` files where downstream scripts need exact R types. Use CSV files for manual review; use RDS files for pipeline reruns.
