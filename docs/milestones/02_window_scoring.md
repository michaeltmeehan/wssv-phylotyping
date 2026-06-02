# Milestone 2: Window Scoring

## Goal

Aggregate informative SNPs into candidate genomic windows and rank regions by tree-informed signal.

## Expected Inputs

- Site summaries.
- Per-site/per-node SNP scores.
- Candidate window-size parameters.
- Optional known marker or GenBank-covered regions.

## Expected Outputs

- Candidate window table.
- Window scoring table.
- Genome-wide window informativeness figures.

## Success Criteria

- Windows are generated reproducibly from configuration.
- Window scores summarize total signal, weighted signal, node coverage, missingness, and gap burden.
- Outputs distinguish high-scoring regions from redundant or low-quality regions.

