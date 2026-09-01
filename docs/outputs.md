# Outputs Catalogue

Generated files are written under `data/processed/` and `outputs/`. They are safe to delete and regenerate if the raw inputs and config are preserved. Do not delete raw files under `data/raw/`.

| Output | Created by | Contains | Intermediate or final | When to inspect | Safe to delete/regenerate |
|---|---|---|---|---|---|
| `data/processed/precomputed.rds` | Stage 01 | Encoded alignment, pruned tree, clade masks, node metadata, polymorphic sites. | Intermediate | After changing raw alignment/tree or clade filters. | Yes |
| `outputs/tables/site_node_scores_all.rds` | Stage 02 | Informative SNP rules by site and node for the full diagnostic clade set. | Intermediate | After changing SNP scoring thresholds. | Yes |
| `outputs/tables/site_node_scores_targets.rds` | Stage 02 | Informative SNP rules by site and node for the restricted target clade set. | Intermediate | After changing SNP scoring thresholds or target clades. | Yes |
| `outputs/tables/site_node_scores.rds` | Stage 02 | Transitional alias of `site_node_scores_all.rds` kept for downstream compatibility. | Intermediate | During the migration period before downstream stages are repointed. | Yes |
| `outputs/tables/site_summary_all.csv` | Stage 03 | Per-site score summaries for the full diagnostic clade set. | Diagnostic/final summary | To find strongest SNPs and signal concentration across all eligible clades. | Yes |
| `outputs/tables/node_summary_all.csv` | Stage 03 | Per-node score summaries for the full diagnostic clade set. | Diagnostic/final summary | To see which broad clades have informative SNPs. | Yes |
| `outputs/tables/site_summary_targets.csv` | Stage 03 | Per-site score summaries for the restricted target clade set. | Diagnostic/final summary | To inspect target-clade SNP candidates and their observed support. | Yes |
| `outputs/tables/node_summary_targets.csv` | Stage 03 | Per-node score summaries for the restricted target clade set. | Diagnostic/final summary | To inspect target clades individually. | Yes |
| `outputs/tables/target_clade_diagnostics.csv` | Stage 03 | Compact one-row-per-target checkpoint for deciding whether the configured target clades are suitable for primer-panel work. | Diagnostic checkpoint | Read before window or primer design to check posterior support, descendant size, informative-site counts, strongest gain, and best-site missingness. | Yes |
| `outputs/tables/target_clade_strongest_snps.csv` | Stage 03 | Strongest informative SNPs for each configured target clade. | Diagnostic checkpoint | Read alongside `target_clade_diagnostics.csv` to inspect the exact SNPs driving each target-clade signal. | Yes |
| `outputs/tables/candidate_amplicons.csv` | Stage 04 | Candidate amplicon intervals with per-candidate completeness, flank, and target-coverage summaries. | Diagnostic/intermediate | To inspect which genomic regions have assay-feasible discriminatory signal. | Yes |
| `outputs/tables/candidate_amplicon_target_scores.csv` | Stage 04 | One row per candidate amplicon and primary target clade with target-specific gain summaries. | Diagnostic/intermediate | To inspect which clades contribute signal to each candidate region. | Yes |
| `outputs/tables/candidate_amplicon_snps.csv` | Stage 04 | One row per candidate amplicon, target clade, and informative SNP association. | Diagnostic/intermediate | To audit the exact SNPs supporting each region. | Yes |
| `outputs/figures/candidate_amplicons_spatial.png` | Stage 04 | Simple spatial diagnostic of informative SNPs and candidate amplicons. | Diagnostic figure | To see SNP clustering along the genome at a glance. | Yes |
| `outputs/tables/panel_candidates.csv` | Stage 05 | One row per evaluated candidate panel, including explicit per-target evidence vectors, redundancy metrics, feasibility summaries, within-size Pareto status, cross-size Pareto status, and deterministic within-size rank. | Final diagnostic summary | To compare candidate panels transparently. | Yes |
| `outputs/tables/panel_target_scores.csv` | Stage 05 | One row per evaluated panel and primary target, with the best supporting amplicon, best and second-best normalized gains, relative gain, and support counts. | Final diagnostic summary | To inspect how each panel supports each target. | Yes |
| `outputs/tables/pareto_panels_by_size.csv` | Stage 05 | Non-dominated subset of evaluated panels within each fixed panel size. | Final diagnostic summary | To review the within-size Pareto frontiers. | Yes |
| `outputs/tables/cross_size_pareto_panels.csv` | Stage 05 | Non-dominated subset of evaluated panels when panel size is treated as a separate minimisation objective. | Secondary diagnostic summary | To review trade-offs across sizes. | Yes |
| `outputs/tables/pareto_panels.csv` | Stage 05 | Compatibility alias of `cross_size_pareto_panels.csv`. | Compatibility summary | For legacy downstream scripts. | Yes |
| `outputs/tables/recommended_provisional_panel.csv` | Stage 05 | Amplicons from the recommended provisional panel, with per-target evidence, feasibility metrics, and recommendation metadata. | Final recommended panel | Always inspect after panel tuning. | Yes |
| `outputs/tables/selected_panel.csv` | Stage 05 | Compatibility alias of `recommended_provisional_panel.csv`. | Compatibility summary | For legacy downstream scripts. | Yes |
| `outputs/tables/recommended_panel_alternatives.csv` | Stage 05 | Top alternative panels at the recommended size. | Final diagnostic summary | To keep backup choices if primer design later fails. | Yes |
| `outputs/tables/panel_size_comparison.csv` | Stage 05 | Best panel of each size together with search method, evaluation counts, Pareto counts, marginal deltas, and recommendation metadata. | Diagnostic summary | To compare the marginal benefit of larger panels. | Yes |
| `outputs/tables/panel_summary.csv` | Stage 05 | Compatibility alias of `panel_size_comparison.csv`. | Compatibility summary | For legacy downstream scripts and report generation. | Yes |
| `outputs/tables/panel_search_summary.csv` | Stage 05 | Candidate counts before and after filtering/reduction, per-size search method summary, and recommendation metadata. | Diagnostic summary | To check whether exact or approximate search was used. | Yes |
| `outputs/figures/recommended_provisional_panel_panels.png` | Stage 05 | Recommended panel layout plus a compact best-by-size summary. | Diagnostic figure | To inspect panel layout and marginal improvement at a glance. | Yes |
| `outputs/figures/selected_panel_panels.png` | Stage 05 | Compatibility alias of `recommended_provisional_panel_panels.png`. | Compatibility figure | For legacy downstream scripts. | Yes |
| `outputs/tables/validation_panel_summary.csv` | Stage 06 | Complete-genome validation summaries by panel size. | Legacy diagnostic | To compare parameter sets. | Yes |
| `outputs/tables/validation_tip_summary.csv` | Stage 06 | Per-tip retained signal summaries. | Legacy diagnostic | To find genomes with poor retained signal. | Yes |
| `outputs/tables/validation_fragment_summary.csv` | Stage 06 | Random fixed-fragment signal summaries. | Legacy diagnostic baseline | To compare selected panels to generic fragments. | Yes |
| `outputs/tables/validation_baseline_summary.csv` | Stage 06 | Random matched-window baseline summaries. | Legacy diagnostic baseline | To check selected panel benefit over random windows. | Yes |
| `outputs/models/wssv_classifier.rds` | Stage 07 | Saved classifier object. | Legacy model for current config | Before partial classification; archive with experiment outputs if needed. | Yes, but rerun stage 07 before classifying. |
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

Stage 04 now reports candidate amplicon regions rather than a final assay panel. Stage 05 now compares candidate amplicon panels explicitly and retains both absolute and target-relative evidence vectors. Later stages keep compatibility outputs where needed, but the primary stage-05 products are the new panel comparison tables above.
