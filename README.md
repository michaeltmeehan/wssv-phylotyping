# WSSV Phylotyping

This repository contains a staged R workflow for WSSV phylotyping. It uses a whole-genome alignment and an MCC tree to find informative SNPs, group them into candidate marker windows, select compact marker panels, train a conservative tree-path classifier, validate retained signal, and classify complete or partial genomes when enough evidence is available.

The pipeline is designed to be conservative. A resolved classifier call means the observed SNP evidence passed the configured thresholds; it does not automatically prove that the biological classification is correct.

## Main Workflow

Plain-language summary:

1. Match the raw alignment and tree so they contain the same genomes.
2. Find polymorphic sites and score SNPs that help distinguish tree clades.
3. Review the stage-03 target-clade checkpoint to see whether the four-clade simplification is viable before window or primer-panel design.
4. Combine useful SNPs into genomic windows.
5. Select a small panel of complementary marker windows.
6. Check how much signal the selected panel retains in complete genomes.
7. Train a tree-path classifier from the selected panel or scored SNPs.
8. Classify partial genomes if they can be mapped and contain enough observed informative sites.
9. Generate a lightweight Markdown report from the outputs.

The numbered scripts in `scripts/` are the canonical workflow stages. The simplest baseline command is:

```sh
Rscript scripts/run_all.R
```

You can also run stages one by one:

```sh
Rscript scripts/01_preprocess_alignment_tree.R
Rscript scripts/02_score_snps.R
Rscript scripts/03_summarise_sites.R
Rscript scripts/04_score_windows.R
Rscript scripts/05_select_marker_panel.R
Rscript scripts/06_validate_panel.R
Rscript scripts/07_train_classifier.R
Rscript scripts/08_classify_partials.R
Rscript scripts/09_generate_report.R
```

## Required Inputs

The default paths are set in `config/config.yml`:

- Whole-genome alignment FASTA: `data/raw/alignment/41seqscollagen_Edited.fasta`
- MCC tree file: `data/raw/tree/41SEQSUPDATEDEDITEDMCC`
- Optional partial-genome FASTA files: `data/raw/partials/`
- Optional metadata files: `data/raw/metadata/`

Raw data are user-supplied and should stay under `data/raw/`. Do not edit raw input files as part of routine analysis.

## Minimal Setup

Use R from the repository root. Required or commonly used packages are:

- `yaml`
- `ape`
- `phangorn`
- `Biostrings`
- `pwalign` or a Biostrings version with pairwise-alignment functions, for pairwise local mapping of partial genomes

Install missing CRAN packages with `install.packages()`. Bioconductor packages such as `Biostrings` and `pwalign` are usually installed with `BiocManager::install()`.

## For Students Tuning Parameters

Start with the current `config/config.yml` as the baseline. Change one parameter group at a time, such as only `analysis.panel_selection`, only `analysis.classifier`, or only `analysis.partials`. After each run, record the config, command, output directory or copied output snapshot, and a short interpretation.

Good tuning is not just maximizing resolved classifications. Prefer settings that keep validation signal high, avoid many conflicts, avoid many weak calls, retain enough informative sites, and do not overstate partial-genome evidence.

## Detailed Documentation

- [Workflow](docs/workflow.md)
- [Configuration](docs/configuration.md)
- [Tuning guide](docs/tuning_guide.md)
- [Classifier outputs](docs/classifier_outputs.md)
- [Outputs catalogue](docs/outputs.md)
- [Data inputs](docs/data_inputs.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Methods overview](docs/methods_overview.md)
- [Glossary](docs/glossary.md)
- [Tuning experiment template](docs/tuning_experiment_template.md)

Generated processed data, tables, figures, models, and local raw inputs are ignored by git. Write generated analysis outputs under `data/processed/` or `outputs/`.
