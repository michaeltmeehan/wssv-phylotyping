# Milestone 1: Preprocess Alignment, Tree, and Score SNPs

## Goal

Create reproducible preprocessing objects and score SNPs by their ability to distinguish eligible MCC-tree clades.

## Expected Inputs

- Whole-genome alignment FASTA.
- MCC tree file.
- Optional sample metadata.
- Analysis parameters from `config/config.yml`.

## Expected Outputs

- Precomputed alignment/tree object.
- Node metadata and target masks.
- Per-site/per-node SNP scoring table.
- Site-level summary table.

## Success Criteria

- Alignment and tree labels are matched and validated.
- Eligible clades are consistently defined.
- SNP scoring outputs include enough information for later interpretation and classifier training.
- Initial tests cover label matching, encoding, and basic scoring expectations.

