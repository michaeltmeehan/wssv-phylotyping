# WSSV Phylotyping

This repository contains a staged tree-informed phylotyping pipeline for white
spot syndrome virus (WSSV). The aim is to identify compact genomic marker
regions that preserve useful phylogenetic signal from whole-genome alignments
and an MCC tree, then use those regions to classify complete and partial genomes
conservatively.

Implemented workflow:

1. Preprocess a whole-genome alignment and MCC tree.
2. Score SNPs by their ability to distinguish MCC-tree clades.
3. Aggregate SNP scores into candidate genomic windows.
4. Select complementary marker windows.
5. Validate reduced-region placement using complete genomes.
6. Train a conservative tree-path classifier.
7. Classify partial genomes where sufficient signal is available.

Milestones 0-6 are implemented. The runnable entry points are the numbered
scripts in `scripts/`, from `01_preprocess_alignment_tree.R` through
`08_classify_partials.R`; paths, thresholds, window settings, marker-panel
selection settings, validation settings, classifier thresholds, and partial
FASTA settings are read from `config/config.yml`.

Generated processed data, tables, figures, models, and local raw inputs are
ignored by git. Keep user-supplied raw inputs under `data/raw/`, and write
analysis outputs under `data/processed/` or `outputs/`.
