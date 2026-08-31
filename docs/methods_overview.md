# Methods Overview

This page explains the method in plain language. It avoids mathematical detail and points to the output files that show each result.

## Core Idea

The repository asks: which SNPs and genomic windows help identify where a WSSV genome belongs on an existing tree?

The workflow starts with a complete-genome alignment and an MCC tree. Each internal tree node represents a clade: a group of related tips. The pipeline first defines a broad diagnostic set of eligible clades, then separately selects a smaller target set for initial assay and classifier design. SNPs are searched against the eligible diagnostic set, grouped into windows, selected into marker panels, and then used by a classifier to make conservative tree-path assignments.

## Marker Scoring

A marker is an informative SNP rule. For a given site and node, the code asks whether one allele is more consistent with the clade or with genomes outside the clade.

For example, a SNP may help if most genomes inside a clade have `A` while genomes outside the clade usually do not. The SNP is kept only if enough genomes are observed and the improvement over a simple clade-size baseline is large enough.

Relevant outputs:

- `outputs/tables/site_node_scores_all.rds`
- `outputs/tables/site_node_scores_targets.rds`
- `outputs/tables/site_summary_all.csv`
- `outputs/tables/node_summary_all.csv`
- `outputs/tables/site_summary_targets.csv`
- `outputs/tables/node_summary_targets.csv`

## Window Scoring

Single SNPs may be too narrow for practical marker design. The workflow groups informative SNPs into candidate genomic windows.

Two window types can be created:

- Fixed windows that tile across the alignment.
- SNP-centred windows around informative SNPs when enabled.

Each window receives summary scores based on the informative SNPs and nodes it covers.

Relevant outputs:

- `outputs/tables/candidate_windows.csv`
- `outputs/tables/window_summary.csv`
- `outputs/tables/window_node_summary.csv`

## Panel Selection

A good panel should contain complementary windows, not just several windows that all support the same clade. The greedy panel selector adds one window at a time. At each step it prefers a compatible window that adds the most new weighted node coverage under the configured overlap and distance rules.

Relevant outputs:

- `outputs/tables/selected_panel.csv`
- `outputs/tables/selected_panel_steps.csv`
- `outputs/tables/panel_summary.csv`
- `outputs/tables/panel_node_coverage.csv`

## Validation

Validation here checks how much complete-genome signal remains when only selected panel regions or simulated fragments are considered. It is a signal-retention check, not proof that future unknown samples will always classify correctly.

Relevant outputs:

- `outputs/tables/validation_panel_summary.csv`
- `outputs/tables/validation_tip_summary.csv`
- `outputs/tables/validation_fragment_summary.csv`
- `outputs/tables/validation_baseline_summary.csv`

## Classifier

The classifier stores SNP rules for eligible nodes. For a query genome, it counts observed informative sites and calculates support and conflict for each node. It then tries to choose a supported nested tree path.

It may return:

- `resolved` when support, conflict, and margin thresholds pass.
- `weak_support` when evidence exists but is not strong enough.
- `conflicting` when incompatible nodes have strong competing evidence.
- `no_informative_sites` when no classifier sites are observed.
- `unresolved` when no conservative assignment can be made.

Relevant outputs:

- `outputs/models/wssv_classifier.rds`
- `outputs/tables/classifier_training_summary.csv`
- `outputs/tables/classifier_tip_predictions.csv`

## Partial-Genome Classification

Partial genomes must first be mapped to alignment coordinates. A partial record may cover only a few informative SNPs, so the classifier must be interpreted more cautiously than for complete genomes.

The partial workflow can use:

- Panel mode: only the trained selected-panel classifier.
- Opportunistic mode: all scored informative SNPs that happen to fall inside the mapped partial interval.
- Auto mode: panel first, then opportunistic fallback when panel sites are not observed.

Relevant outputs:

- `outputs/tables/partial_classifications.csv`
- `outputs/tables/partial_region_diagnostics.csv`
- `outputs/tables/partial_opportunistic_site_summary.csv`

## Interpretation Principle

More resolved calls are not automatically better. A useful parameter set should preserve signal, avoid unsupported or conflicting calls, and explain non-resolved records clearly.
