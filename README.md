# WSSV Phylotyping

This repository will develop a tree-informed phylotyping pipeline for white spot syndrome virus (WSSV). The aim is to identify compact genomic marker regions that preserve useful phylogenetic signal from whole-genome alignments and an MCC tree, then use those regions to classify complete and partial genomes conservatively.

Planned workflow:

1. Preprocess a whole-genome alignment and MCC tree.
2. Score SNPs by their ability to distinguish MCC-tree clades.
3. Aggregate SNP scores into candidate genomic windows.
4. Select complementary marker windows.
5. Validate reduced-region placement using complete genomes.
6. Train a probabilistic/tree-path classifier.
7. Classify partial genomes where sufficient signal is available.

Milestone 1 implements preprocessing and SNP-node scoring. The runnable entry
points are `scripts/01_preprocess_alignment_tree.R`,
`scripts/02_score_snps.R`, and `scripts/03_summarise_sites.R`; paths and
thresholds are read from `config/config.yml`.
