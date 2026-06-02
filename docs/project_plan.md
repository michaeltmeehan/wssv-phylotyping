# WSSV Phylotyping Project Plan

This project will turn the current exploratory SNP-scoring workflow into a reproducible, tree-informed phylotyping pipeline for WSSV.

The linked products are:

- An informative SNP catalogue.
- An optimal marker/window panel.
- A validated classifier.
- A partial-genome phylotyping workflow.

## Stage 0: Repository Setup

Create the repository skeleton, documentation, configuration, data/output conventions, placeholder scripts, placeholder function files, and initial tests.

Primary output: version-controlled project structure.

## Stage 1: Clean Preprocessing and SNP Scoring

Preprocess the whole-genome alignment and MCC tree into reusable objects:

- MCC tree pruned to alignment tips.
- Integer-encoded alignment matrix.
- Site coordinate mapping.
- Node and clade metadata.
- Node-by-tip target masks.
- Optional tip metadata.

Then score polymorphic sites by how well they distinguish eligible MCC-tree clades from the rest of the tree.

Primary outputs:

- `data/processed/precomputed.rds`
- `data/processed/site_node_scores.rds`
- `outputs/tables/site_summary.csv`

## Stage 2: Window Scoring

Aggregate site-level signal into candidate genomic windows. Candidate windows may include fixed sliding windows, SNP-centred windows, known marker regions, GenBank-covered regions, and future primer-compatible amplicons.

Primary outputs:

- `outputs/tables/candidate_windows.csv`
- `outputs/tables/window_summary.csv`
- Genome-wide window informativeness figures.

## Stage 3: Greedy Marker Panel Selection

Select complementary non-redundant windows that jointly resolve MCC-tree nodes. The panel search should reward marginal coverage of unresolved nodes and avoid repeatedly selecting windows that explain the same split.

Primary outputs:

- `outputs/tables/selected_panel.csv`
- Node coverage tables.
- Panel comparison summaries.

## Stage 4: Validation Using Complete Genomes

Validate the selected marker panel with complete genomes using leave-one-out classification, artificial partial-genome experiments, random-window baselines, and window-only tree reconstruction.

Primary outputs:

- `outputs/tables/validation_results.csv`
- Fragment-length performance summaries.
- Window-only tree recovery summaries.

## Stage 5: Classifier Construction

Train interpretable classifiers that use observed informative sites only. Candidate approaches include a naive probabilistic clade classifier and a conservative tree-path classifier that stops when support becomes weak.

Primary outputs:

- `outputs/models/wssv_classifier.rds`
- Training-set sanity-check prediction tables.
- Per-node evidence diagnostics when compact enough.

Training-set predictions are internal checks only; independent validation is
handled by the Stage 4 validation scaffolding and future external partial-genome
work.

## Stage 6: Partial-Genome Classification

Map external partial WSSV sequences to reference coordinates, extract overlapping informative sites, classify them where possible, and report uncertainty or conflicts when signal is insufficient.

Primary outputs:

- `outputs/tables/partial_classifications.csv`
- Partial sequence coverage and informativeness summaries.

## Analysis Products

The final project should produce:

- A genome-wide informativeness map.
- A top SNP catalogue.
- An optimal marker/window panel.
- A validation report.
- A partial-genome classifier.
- A retrospective classification of informative partial WSSV sequences.
