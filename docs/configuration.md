# Configuration

The main configuration file is `config/config.yml`. Run scripts read this file directly. The table below documents the current parameters and defaults in this repository.

Suggested ranges below are starting points for cautious experiments, not validated biological limits. Entries marked `TODO: Michael review` are especially dependent on domain judgement or the current dataset.

Tuning labels:

- **Safe for student tuning**: reasonable to adjust during parameter experiments.
- **Advanced**: change only with a clear reason and careful output checks.
- **Do not change unless you understand the code**: these affect file layout, coordinate mapping assumptions, or computational behavior.

## Paths

| Parameter | Current value | Tuning label | Meaning | Increase/decrease effect | Safe range | Risks and checks |
|---|---:|---|---|---|---|---|
| `paths.alignment` | `data/raw/alignment/41seqscollagen_Edited.fasta` | Do not change unless you understand the code | Whole-genome alignment FASTA. | Not numeric. | Existing FASTA path. | Wrong file gives mismatched tips or invalid coordinates. Check `precomputed.rds` messages and `dropped_alignment`. |
| `paths.tree` | `data/raw/tree/41SEQSUPDATEDEDITEDMCC` | Do not change unless you understand the code | MCC tree for the same genomes. | Not numeric. | Existing Newick/Nexus path. | Wrong tree changes all node labels and clades. Check retained tips and eligible nodes. |
| `paths.partials` | `data/raw/partials/` | Safe for student tuning | Directory of optional partial FASTA files. | Not numeric. | Existing directory. | Empty directory creates empty partial tables. Check `partial_classifications.csv`. |
| `paths.processed` | `data/processed/` | Do not change unless you understand the code | Generated processed-object directory. | Not numeric. | Directory under `data/processed/`. | Changing it can make stages miss upstream files. |
| `paths.outputs` | `outputs/` | Advanced | Generated output root. | Not numeric. | Directory under `outputs/` or a named experiment output directory. | If changed, all downstream scripts must use the same config. Useful for run comparisons. |

## Diagnostic Clades And SNP Scoring

| Parameter | Current value | Tuning label | Plain-language meaning | If increased | If decreased | Suggested range | Inspect after changing |
|---|---:|---|---|---|---|---|---|
| `analysis.min_clade_size` | `2` | Advanced | Smallest eligible diagnostic clade size. | Excludes smaller clades; fewer eligible nodes. | Includes smaller clades; may score fragile clades. | TODO: Michael review; try 2 to 5 only with node-count checks. | Stage 01 messages, `node_summary_all.csv`. |
| `analysis.max_clade_frac` | `0.95` | Advanced | Largest eligible diagnostic clade as a fraction of retained tips. | Allows very broad clades if closer to 1. | Excludes broad/root-like clades. | TODO: Michael review; try 0.8 to 0.98 only with node-count checks. | Eligible node count, `node_summary_all.csv`. |
| `analysis.min_target_posterior_support` | `0.90` | Advanced | Minimum tree posterior support required for an explicitly selected target clade. | Tightens the target set to higher-confidence clades. | Allows lower-support targets. | TODO: Michael review; 0.90 to 0.99 is a cautious range. | Stage 01 target-clade messages, `target_node_metadata$posterior_support`. |
| `analysis.target_clades` | `4 node_id entries` | Advanced | Explicitly selected clades used for stage 01 target-clade summarisation and the legacy assay/classifier path. Stage 04 now uses `analysis.amplicons.primary_target_nodes` for its own filtering. | Not numeric. | Not numeric. | Keep the configured target set to the intended major clades. | Stage 01 target-clade messages, `target_clade_mask`, `target_node_metadata`. |
| `analysis.min_total_obs` | `30` | Safe for student tuning | Minimum observed non-missing bases required for a site-node test. | Fewer SNP rules, more conservative. | More SNP rules, more missingness risk. | TODO: Michael review; start near 20 to retained tip count. | `site_node_scores_all.rds`, `site_node_scores_targets.rds`, `site_summary_all.csv`. |
| `analysis.min_side_obs` | `2` | Safe for student tuning | Minimum observed bases required on both inside and outside sides of a clade. | Excludes poorly represented splits. | Allows fragile comparisons. | 2 to 5. | `node_summary_all.csv`, eligible scored nodes. |
| `analysis.min_site_maf` | `2` | Safe for student tuning | Minimum count for an allele before it can define a SNP rule. | Removes rare allele rules. | Includes rare allele rules. | 2 to 5. | `site_summary_all.csv`, `site_node_scores_all.rds`, `site_node_scores_targets.rds`. |
| `analysis.min_gain_norm` | `0.5` | Safe for student tuning | Minimum normalized improvement over the clade-size baseline. | Fewer, stronger SNP rules. | More, weaker SNP rules. | TODO: Michael review; 0.3 to 0.8 is a cautious exploration range. | `site_summary_all.csv`, `window_summary.csv`. |
| `analysis.chunk_size` | `5000` | Do not change unless you understand the code | Progress-message interval during SNP scoring. | Fewer progress messages. | More progress messages. | 1000 to 10000. | Runtime only. |
| `analysis.candidate_window_sizes` | `300, 500, 1000, 2000` | Advanced | Legacy fallback window widths if `analysis.windows` is absent. | Larger windows capture more sites. | Smaller windows are more targeted. | Use `analysis.windows.widths` instead. | `window_summary.csv`. |

## Amplicon Candidate Scoring

| Parameter | Current value | Tuning label | Plain-language meaning | If increased | If decreased | Suggested range | Inspect after changing |
|---|---:|---|---|---|---|---|---|
| `analysis.amplicons.primary_target_nodes` | `45, 66, 72, 80` | Do not change unless you understand the code | Primary assay targets used by stage 04 to build candidate amplicons. | Not numeric. | Not numeric. | Keep the intended major clades only. | `candidate_amplicons.csv`, `candidate_amplicon_target_scores.csv`. |
| `analysis.amplicons.diagnostic_target_nodes` | `67, 70` | Advanced | Extra diagnostic nodes that may remain in the broader target catalogue but should not drive primary amplicon selection. | Not numeric. | Not numeric. | Keep any optional diagnostics separate from the primary targets. | `site_node_scores_targets.rds`, stage 04 filtering messages. |
| `analysis.amplicons.min_length` | `500` | Safe for student tuning | Minimum candidate amplicon length. | Larger, broader candidate intervals. | Smaller, tighter candidate intervals. | 500 to 3000. | `candidate_amplicons.csv`. |
| `analysis.amplicons.max_length` | `3000` | Safe for student tuning | Maximum candidate amplicon length. | Allows wider candidate regions. | Forces tighter clusters. | 500 to 3000. | `candidate_amplicons.csv`. |
| `analysis.amplicons.flank_search_width` | `100` | Safe for student tuning | How far left and right stage 04 inspects for simple flank quality checks. | Wider flank diagnostics. | Narrower flank diagnostics. | 50 to 250. | `candidate_amplicons.csv`. |
| `analysis.amplicons.minimum_gain_threshold` | `0.5` | Safe for student tuning | Normalized-gain threshold counted as any signal for a target. | Fewer targets counted as signal-bearing. | More targets counted as signal-bearing. | 0.3 to 0.8. | `targets_with_any_signal`. |
| `analysis.amplicons.strong_gain_threshold` | `0.8` | Safe for student tuning | Normalized-gain threshold counted as strong signal for a target. | Fewer strong-signal targets. | More strong-signal targets. | 0.6 to 0.95. | `targets_with_strong_signal`. |
| `analysis.amplicons.flank_min_non_missing_fraction` | `0.9` | Safe for student tuning | Minimum observed completeness for a flank stretch to be considered promising. | Stricter conserved flank calls. | More permissive conserved flank calls. | 0.8 to 0.99. | `left_flank_conserved_stretch`, `right_flank_conserved_stretch`. |
| `analysis.amplicons.flank_max_variable_fraction` | `0.25` | Safe for student tuning | Maximum polymorphic-site fraction for a flank stretch to be considered conserved. | Fewer conserved flanks. | More conserved flanks. | 0.05 to 0.3. | Flank diagnostics in `candidate_amplicons.csv`. |
| `analysis.amplicons.flank_min_gc_fraction` | `0.35` | Safe for student tuning | Minimum GC fraction for a flank stretch to look primer-friendly. | Biases toward GC-richer flanks. | Allows more AT-rich flanks. | 0.3 to 0.6. | Flank diagnostics in `candidate_amplicons.csv`. |
| `analysis.amplicons.flank_min_conserved_fraction` | `0.9` | Safe for student tuning | Minimum complete-sequence fraction for a flank stretch to count as conserved. | Requires more complete flank sequence. | Allows more partial flank sequence. | 0.8 to 0.99. | Flank diagnostics in `candidate_amplicons.csv`. |
| `analysis.amplicons.flank_complete_fraction` | `0.95` | Safe for student tuning | Completeness threshold used for reporting near-complete candidate regions. | More sequences count as near-complete. | Fewer sequences count as near-complete. | 0.9 to 0.99. | `complete_fraction`, `near_complete_fraction`. |

## Legacy Window Scoring

| Parameter | Current value | Tuning label | Plain-language meaning | If increased | If decreased | Suggested range | Inspect after changing |
|---|---:|---|---|---|---|---|---|
| `analysis.windows.widths` | `300, 500, 1000, 2000` | Legacy | Fixed window sizes to tile across the alignment. This is no longer the primary stage-04 path. | Larger windows may capture more SNPs but are less compact. | Smaller windows are more specific but may miss signal. | TODO: Michael review; 300 to 5000 is an exploratory range. | Legacy `candidate_windows.csv`, `window_summary.csv`. |
| `analysis.windows.overlap` | `0.5` | Legacy | Fractional overlap between fixed windows. | More overlapping windows, slower and more redundant. | Fewer windows, faster but less fine-grained. | 0 to 0.75. | Number of candidate windows. |
| `analysis.windows.min_informative_snps` | `1` | Legacy | Minimum informative SNPs for a window to be retained in scoring. | Fewer windows, stronger minimum signal. | More windows, including sparse windows. | 1 to 5. | `window_summary.csv`. |
| `analysis.windows.snp_centered.enabled` | `true` | Legacy | Adds windows centred on informative SNPs. | Not numeric; `true` increases candidate windows. | `false` uses fixed windows only. | `true` for discovery, `false` for simple tiling. | Legacy `candidate_windows.csv`, runtime. |
| `analysis.windows.snp_centered.widths` | `300, 500` | Legacy | Widths for SNP-centred windows. | Wider SNP-centred regions. | Narrower SNP-centred regions. | 300 to 1000. | Legacy `candidate_windows.csv`, `window_summary.csv`. |

## Panel Selection

| Parameter | Current value | Tuning label | Plain-language meaning | If increased | If decreased | Suggested range | Inspect after changing |
|---|---:|---|---|---|---|---|---|
| `analysis.panel_selection.min_panel_size` | `2` | Safe for student tuning | Smallest panel size to evaluate. | Includes larger minimum panels only. | Allows smaller panels. | 2 to 4. | `panel_candidates.csv`, `panel_search_summary.csv`. |
| `analysis.panel_selection.max_panel_size` | `6` | Safe for student tuning | Largest panel size to evaluate. | Allows larger comparison panels. | Narrows the search. | 3 to 8. | `panel_candidates.csv`, `panel_size_comparison.csv`. |
| `analysis.panel_selection.top_n_panels` | `20` | Safe for student tuning | Number of top-ranked panels retained for convenience. | Keeps more alternatives. | Keeps fewer alternatives. | 10 to 50. | `panel_candidates.csv`. |
| `analysis.panel_selection.exact_enumeration_limit` | `100000` | Advanced | Maximum combination count before falling back to a deterministic approximate search. | Makes exact search more likely. | Makes approximate search more likely. | Raise only if the reduced candidate set is small. | `panel_search_summary.csv`. |
| `analysis.panel_selection.approximate_search_budget` | `5000` | Advanced | Maximum number of approximate-search panels retained per size during deterministic beam search. | Explores more candidate panels for large sizes. | Speeds up large-size search but can miss alternatives. | 1000 to 20000. | `panel_search_summary.csv`, `panel_size_comparison.csv`. |
| `analysis.panel_selection.max_reduced_candidates` | `30` | Advanced | Cap on the reduced candidate pool before panel enumeration. | Retains more candidate diversity. | Forces more aggressive reduction. | 20 to 40. | `panel_search_summary.csv`. |
| `analysis.panel_selection.evidence_threshold` | `0.5` | Safe for student tuning | Threshold used when counting how many amplicons support a target. | Counts fewer amplicons as supporting evidence. | Counts more amplicons as supporting evidence. | 0.3 to 0.8. | `panel_target_scores.csv`. |
| `analysis.panel_selection.strong_evidence_threshold` | `0.8` | Safe for student tuning | Higher evidence threshold used for “strong” target support counts. | Fewer amplicons count as strong support. | More amplicons count as strong support. | 0.6 to 0.95. | `panel_candidates.csv`, `panel_target_scores.csv`. |
| `analysis.panel_selection.material_gain_improvement` | `0.01` | Safe for student tuning | Minimum improvement in minimum absolute gain required to keep extending the recommended panel size. | Makes the recommendation more willing to add amplicons. | Makes the recommendation stop earlier. | 0.005 to 0.05. | `panel_size_comparison.csv`, `panel_search_summary.csv`. |
| `analysis.panel_selection.material_relative_gain_improvement` | `0.01` | Safe for student tuning | Minimum improvement in minimum relative gain required to keep extending the recommended panel size. | Makes the recommendation more willing to add amplicons. | Makes the recommendation stop earlier. | 0.005 to 0.05. | `panel_size_comparison.csv`, `panel_search_summary.csv`. |
| `analysis.panel_selection.redundancy_cap_per_target` | `2` | Advanced | Maximum per-target redundancy rewarded in within-size ranking. | Rewards extra supporting amplicons more strongly. | Saturates redundancy sooner. | 1 to 3. | `panel_candidates.csv`, `panel_size_comparison.csv`. |
| `analysis.panel_selection.max_candidate_length` | `3000` | Safe for student tuning | Maximum amplicon length allowed into the panel search. | Allows longer amplicons. | Excludes longer amplicons. | 500 to 3000. | `panel_candidates.csv`. |
| `analysis.panel_selection.min_complete_fraction` | `0.9` | Safe for student tuning | Minimum candidate-region completeness allowed into the search. | Keeps only more complete amplicons. | Allows more incomplete regions. | 0.8 to 0.99. | `panel_candidates.csv`. |
| `analysis.panel_selection.max_missing_fraction` | `null` | Safe for student tuning | Optional maximum missing-cell fraction for a candidate amplicon. | Excludes more gappy regions. | Allows more gappy regions. | `null` or 0.05 to 0.3. | `panel_candidates.csv`. |
| `analysis.panel_selection.require_left_flank_conserved` | `false` | Advanced | Require the left flank to pass the stage-04 flank feasibility check. | Restricts candidates. | Keeps more candidates. | Usually `false` until flank review is complete. | `panel_candidates.csv`. |
| `analysis.panel_selection.require_right_flank_conserved` | `false` | Advanced | Require the right flank to pass the stage-04 flank feasibility check. | Restricts candidates. | Keeps more candidates. | Usually `false` until flank review is complete. | `panel_candidates.csv`. |
| `analysis.panel_selection.require_both_flanks_conserved` | `false` | Advanced | Require both flanks to pass the stage-04 flank feasibility check. | Most restrictive. | Least restrictive. | `false` unless the search is already well behaved. | `panel_candidates.csv`. |
| `analysis.panel_selection.min_target_signal` | `null` | Advanced | Optional minimum strongest target gain for a candidate amplicon. | Filters weaker candidates. | Keeps weaker candidates. | `null` or 0.4 to 0.8. | `panel_candidates.csv`. |
| `analysis.panel_selection.min_targets_with_signal` | `null` | Advanced | Optional minimum number of primary targets with any signal in a candidate amplicon. | Makes candidates more broadly informative. | Allows more target-specific candidates. | `null` or 1 to 4. | `panel_candidates.csv`. |
| `analysis.panel_selection.min_targets_with_strong_signal` | `null` | Advanced | Optional minimum number of primary targets with strong signal in a candidate amplicon. | Makes candidates stricter. | Makes candidates more permissive. | `null` or 1 to 3. | `panel_candidates.csv`. |
| `analysis.panel_selection.max_allowed_overlap` | `0.5` | Safe for student tuning | Maximum overlap allowed between panel amplicons. | Permits more overlap. | Forbids more overlap. | 0 to 0.5. | `recommended_provisional_panel.csv`, `selected_panel.csv`, `panel_candidates.csv`. |
| `analysis.panel_selection.min_distance` | `null` | Safe for student tuning | Optional minimum base spacing between panel amplicons. | Spreads panel members apart. | Allows nearby amplicons. | `null` or 20 to 1000. | `panel_candidates.csv`. |
| `analysis.panel_selection.dominance_overlap_threshold` | `0.85` | Advanced | Overlap threshold used when removing dominated local duplicates. | Removes only near-identical local candidates. | Removes more overlapping candidates. | 0.75 to 0.95. | `panel_search_summary.csv`. |
| `analysis.panel_selection.dominance_min_distance` | `100` | Advanced | Distance threshold used when removing dominated local duplicates. | Keeps more nearby candidates. | Collapses local duplicates more aggressively. | 50 to 200. | `panel_search_summary.csv`. |
| `analysis.panel_selection.local_redundancy_overlap_threshold` | `0.9` | Advanced | Overlap threshold used to identify locally redundant candidates. | Keeps more nearby variants. | Collapses more local duplicates. | 0.8 to 0.95. | `panel_search_summary.csv`. |
| `analysis.panel_selection.local_redundancy_min_distance` | `100` | Advanced | Minimum distance used when identifying locally redundant candidates. | Keeps more nearby variants. | Collapses more local duplicates. | 50 to 200. | `panel_search_summary.csv`. |
| `analysis.panel_selection.target_gain_round_digits` | `3` | Safe for student tuning | Decimal precision used when comparing and reporting candidate target-gain vectors. | Preserves more precision. | Makes reduction more coarse. | 2 to 4. | `panel_candidates.csv`, `panel_target_scores.csv`. |

## Validation

| Parameter | Current value | Tuning label | Plain-language meaning | If increased | If decreased | Suggested range | Inspect after changing |
|---|---:|---|---|---|---|---|---|
| `analysis.validation.panel_sizes` | `1, 2, 3, 5, 10` | Safe for student tuning | Panel sizes evaluated in complete genomes. | Evaluates larger panels. | Evaluates fewer panels. | Match panel selection sizes. | `validation_panel_summary.csv`. |
| `analysis.validation.fragment_lengths` | `300, 500, 1000, 2000, 5000` | Safe for student tuning | Random fragment lengths for signal comparison. | Tests longer fragments. | Tests shorter fragments. | 300 to 5000. | `validation_fragment_summary.csv`. |
| `analysis.validation.random_baseline_replicates` | `10` | Safe for student tuning | Number of random matched-window baselines. | More stable baseline, slower. | Faster, noisier baseline. | 10 to 100 for exploratory work; more if runtime allows. | `validation_baseline_summary.csv`. |
| `analysis.validation.min_informative_sites` | `1` | Safe for student tuning | Minimum observed informative sites for validation signal. | Stricter resolved-signal summaries. | More permissive summaries. | 1 to 5. | `validation_tip_summary.csv`. |
| `analysis.validation.min_informative_nodes` | `1` | Safe for student tuning | Minimum informative nodes for validation signal. | Stricter summaries. | More permissive summaries. | 1 to 3. | `validation_panel_summary.csv`. |
| `analysis.validation.random_seed` | `1` | Safe for student tuning | Seed for random baseline reproducibility. | Different random samples. | Different random samples. | Any integer. | Baseline reproducibility. |

## Classifier

| Parameter | Current value | Tuning label | Plain-language meaning | If increased | If decreased | Suggested range | Inspect after changing |
|---|---:|---|---|---|---|---|---|
| `analysis.classifier.min_sites_per_node` | `1` | Safe for student tuning | Minimum rule sites required for a node to be classifiable. | Fewer classifiable nodes, stronger evidence. | More classifiable nodes, sparse evidence risk. | 1 to 5. | `classifier_training_summary.csv`, `classifier_tip_predictions.csv`. |
| `analysis.classifier.min_total_informative_sites` | `1` | Safe for student tuning | Minimum observed informative sites in a query. | More weak/no-call outcomes. | More permissive calls. | 1 to 5 for panels; higher for larger panels. | Status counts in classifier and partial outputs. |
| `analysis.classifier.min_support` | `0.8` | Safe for student tuning | Minimum classifier evidence support fraction for a node to be strong. This is not tree posterior support. | Fewer resolved calls, more weak calls. | More resolved calls, possible over-permissive calls. | TODO: Michael review; 0.7 to 0.95 is an exploratory range. | `status`, `assigned_support`, `strongest_off_path_support`. |
| `analysis.classifier.max_conflict` | `0.2` | Safe for student tuning | Maximum classifier evidence conflict fraction allowed for a strong node. | More permissive conflict tolerance. | More conflicting/weak outcomes. | TODO: Michael review; 0.05 to 0.3 is an exploratory range. | `conflict`, `conflict_reason`, conflicting counts. |
| `analysis.classifier.support_margin` | `0.05` | Safe for student tuning | Required gap between best supported node and incompatible competitors in classifier evidence. | More conservative; more conflict/weak calls. | More permissive; close calls may resolve. | TODO: Michael review; 0 to 0.2 is an exploratory range. | `support_margin`, conflict counts. |
| `analysis.classifier.use_selected_panel` | `true` | Advanced | Train classifier only from selected panel sites. | Not numeric. `false` uses all scored informative SNPs. | Not numeric. | Usually `true` for designed panels. | `classifier_training_summary.csv`. |
| `analysis.classifier.panel_size` | `null` | Safe for student tuning | Number of selected panel windows used; `null` means all selected windows. | Larger panel if available. | Smaller panel. | `null` or one selected panel size. | `classifier_training_summary.csv`, `validation_panel_summary.csv`. |
## Partials

| Parameter | Current value | Tuning label | Plain-language meaning | If increased | If decreased | Suggested range | Inspect after changing |
|---|---:|---|---|---|---|---|---|
| `analysis.partials.input_dir` | `data/raw/partials/` | Safe for student tuning | Directory of partial FASTA records. | Not numeric. | Not numeric. | Existing directory. | Partial records read message. |
| `analysis.partials.reference_id` | `CN01_1994` | Do not change unless you understand the code | Alignment row used as the reference sequence for exact-substring and pairwise-local partial mapping. This was previously hard-coded in `scripts/08_classify_partials.R`. | Not numeric. | Not numeric. | Keep `CN01_1994` unless Michael confirms a different reference row should define coordinates. | `partial_region_diagnostics.csv`, mapping columns. |
| `analysis.partials.fasta_extensions` | `fa, fas, fasta, fna` | Safe for student tuning | File extensions read as FASTA. | More file types included. | Fewer files included. | Add only real FASTA extensions. | Partial records read message. |
| `analysis.partials.mapping_mode` | `auto` | Safe for student tuning | How partial sequences are mapped to alignment coordinates. | Not numeric. | Not numeric. | `auto`, `aligned_full_length`, `exact_substring`, `pairwise_local`. | `partial_region_diagnostics.csv`, mapping columns. |
| `analysis.partials.classification_mode` | `auto` | Safe for student tuning | Whether to use panel, opportunistic, or auto fallback classification. | Not numeric. | Not numeric. | `panel`, `opportunistic`, `auto`. | `classification_source`, separate panel/opportunistic tables. |
| `analysis.partials.pairwise_local.min_identity` | `0.85` | Safe for student tuning | Minimum local-alignment identity for mapping. | Fewer mapped partials, safer. | More mapped partials, more false mapping risk. | TODO: Michael review; 0.8 to 0.95 is an exploratory range. | `mapping_identity`, unmapped count. |
| `analysis.partials.pairwise_local.min_query_coverage` | `0.6` | Safe for student tuning | Minimum fraction of query aligned in local mapping. | Fewer mappings. | More partial mappings. | 0.5 to 0.9. | `mapping_query_coverage`. |
| `analysis.partials.pairwise_local.min_aligned_length` | `20` | Safe for student tuning | Minimum aligned bases for local mapping. | Rejects short matches. | Allows short matches. | 20 to 100. | `mapping_aligned_length`. |
| `analysis.partials.pairwise_local.gap_opening` | `-10` | Advanced | Pairwise local alignment gap-opening score. | More negative penalizes gaps more. | Less negative allows gaps. | Leave unless mapping is reviewed. | Mapping diagnostics. |
| `analysis.partials.pairwise_local.gap_extension` | `-0.5` | Advanced | Pairwise local alignment gap-extension score. | More negative penalizes long gaps more. | Less negative allows longer gaps. | Leave unless mapping is reviewed. | Mapping diagnostics. |
| `analysis.partials.min_observed_informative_sites` | `1` | Safe for student tuning | Minimum observed informative SNPs for partial classification. | More weak/no-call partials. | More permissive partial calls. | 1 to 5. | `observed_informative_sites`, status counts. |
| `analysis.partials.min_support` | `0.8` | Safe for student tuning | Opportunistic classifier evidence support threshold. This is not tree posterior support. | More conservative opportunistic calls. | More permissive opportunistic calls. | TODO: Michael review; 0.7 to 0.95 is an exploratory range. | `partial_classifications_opportunistic.csv`. |
| `analysis.partials.min_support_margin` | `0.05` | Safe for student tuning | Opportunistic classifier margin over competitors. | Fewer close resolved calls. | More close resolved calls. | TODO: Michael review; 0 to 0.2 is an exploratory range. | `support_margin`. |
| `analysis.partials.conflict_support_threshold` | `0.2` | Safe for student tuning | Opportunistic classifier conflict threshold. | More conflict tolerated. | More conflicting outcomes. | TODO: Michael review; 0.05 to 0.3 is an exploratory range. | `conflict_sites`, `status`. |
| `analysis.partials.use_weighted_support` | `true` | Advanced | Uses SNP rule weights for opportunistic support. | Not numeric. `true` weights stronger SNPs more. | `false` treats rules more equally. | Usually `true`. | Compare opportunistic status counts. |
| `analysis.partials.write_unresolved_evidence` | `false` | Safe for student tuning | Whether to write evidence rows for unresolved panel classifications. | `true` writes more diagnostic rows. | `false` keeps files smaller. | `false` normally; `true` for debugging. | `partial_node_evidence.csv`. |

## General Tuning Advice

Always compare a changed run against a baseline. Look at at least:

- `site_summary_all.csv`
- `candidate_amplicons.csv`
- `candidate_amplicon_target_scores.csv`
- `recommended_provisional_panel.csv`
- `panel_summary.csv`
- `validation_panel_summary.csv`
- `classifier_training_summary.csv`
- `classifier_tip_predictions.csv`
- `partial_classifications.csv` if partials are used

Do not call a parameter set better only because it produces more resolved classifications.
