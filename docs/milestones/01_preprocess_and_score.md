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
- `outputs/tables/site_node_scores.rds`
- `outputs/tables/site_summary.csv`
- `outputs/tables/node_summary.csv`

The long scoring table includes node index/id, site, gain, normalized gain,
best allele, allele direction, observed inside/outside counts, total observed
count, and node metadata. Eligible clades are consistently represented as
`target_mask` rows, fixing the prototype's clade/split naming ambiguity.

## Success Criteria

- Alignment and tree labels are matched and validated.
- Eligible clades are consistently defined.
- SNP scoring outputs include enough information for later interpretation and classifier training.
- Initial tests cover label matching, encoding, and basic scoring expectations.
