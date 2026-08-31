# Milestone 1: Preprocess Alignment, Tree, and Score SNPs

## Goal

Create reproducible preprocessing objects and score SNPs by their ability to distinguish eligible MCC-tree clades.

## Expected Inputs

- Whole-genome alignment FASTA.
- MCC tree file.
- Optional sample metadata.
- Analysis parameters from `config/config.yml`.

## Expected Outputs

- `data/processed/precomputed.rds`
- `outputs/tables/site_node_scores_all.rds`
- `outputs/tables/site_node_scores_targets.rds`
- `outputs/tables/site_node_scores.rds` (legacy alias for the full diagnostic set)
- `outputs/tables/site_summary.csv`
- `outputs/tables/node_summary.csv`
- `outputs/tables/target_clade_diagnostics.csv`
- `outputs/tables/target_clade_strongest_snps.csv`

The long scoring table includes node index/id, site, gain, normalized gain,
best allele, allele direction, observed inside/outside counts, total observed
count, and node metadata, including MCC posterior support. Stage 01 now keeps
the full eligible diagnostic clade set separate from the smaller configured
target clade set. The processed object stores the broad masks and metadata as
`all_clade_mask` and `all_node_metadata`, plus the restricted target versions
as `target_clade_mask` and `target_node_metadata`. Legacy aliases are retained
for downstream stages that have not yet been migrated.

Stage 03 adds a compact target-clade checkpoint on top of those summaries. It
is the review point for deciding whether the four configured target clades are
already discriminated well enough to justify proceeding to window and primer
design. If a target clade shows no informative SNPs, weak maximum normalized
gain, or heavy missingness at its best site, that is a sign to revisit the
clade simplification before moving on.

## Success Criteria

- Alignment and tree labels are matched and validated.
- Eligible clades are consistently defined.
- SNP scoring outputs include enough information for later interpretation and classifier training.
- The target-clade checkpoint makes it obvious whether the four-clade simplification is viable before window or primer design.
- Initial tests cover label matching, encoding, and basic scoring expectations.
