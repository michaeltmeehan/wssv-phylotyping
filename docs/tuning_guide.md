# Tuning Guide

Use this guide when comparing parameter settings in `config/config.yml`. The goal is to find settings that are reproducible, conservative, and useful, not simply settings that produce the most resolved classifications.

## Basic Protocol

1. Save the current `config/config.yml` as your baseline record, for example `config/examples/baseline.yml` or a dated copy outside git.
2. Run the baseline:

   ```sh
   Rscript scripts/run_all.R
   ```

3. Copy or rename the generated `outputs/` directory if you want to preserve the baseline before the next run.
4. Change one parameter group at a time.
5. Rerun the workflow from the earliest affected stage.
6. Record the exact change, outputs inspected, and interpretation in an experiment log.

## Change One Group At A Time

Recommended groups:

| Group | Config section | Earliest stage to rerun |
|---|---|---|
| Clade eligibility | `analysis.min_clade_size`, `analysis.max_clade_frac` | `01_preprocess_alignment_tree.R` |
| SNP scoring | `analysis.min_total_obs`, `min_side_obs`, `min_site_maf`, `min_gain_norm` | `02_score_snps.R` after rerunning stage 01 if clade eligibility changed |
| Window design | `analysis.windows` | `04_score_windows.R` |
| Panel selection | `analysis.panel_selection` | `05_select_marker_panel.R` |
| Validation summaries | `analysis.validation` | `06_validate_panel.R` |
| Complete-genome classifier | `analysis.classifier` | `07_train_classifier.R` |
| Partial-genome classification | `analysis.partials` | `08_classify_partials.R` |

When in doubt, run `Rscript scripts/run_all.R` from a clean config so all derived files match the same parameter set.

## Naming Experiment Configs

Use names that encode the purpose, not just a number:

```text
config/examples/baseline.yml
config/examples/conservative.yml
config/examples/sensitive.yml
config/examples/partial_genomes.yml
config/experiments/2026-06-22_min-support-090.yml
```

The scripts currently read `config/config.yml`. To run an experiment, copy the experiment YAML into `config/config.yml`, run the workflow, and record the copy you used. Do not edit raw input files for an experiment.

## Comparing Runs

Compare each run against baseline using the same set of outputs:

| Question | Output to inspect | What to prefer |
|---|---|---|
| Did SNP scoring become too sparse? | `site_summary.csv`, `node_summary.csv` | Enough scored SNPs across multiple nodes, not just one strong site. |
| Did windows become too sparse or too broad? | `window_summary.csv` | Windows with useful signal and acceptable missingness. |
| Does the panel add complementary signal? | `selected_panel_steps.csv`, `panel_summary.csv` | Positive marginal gain across early panel steps. |
| Does validation retain signal? | `validation_panel_summary.csv` | Higher observed informative sites/nodes and resolved fraction, without relying only on random luck. |
| Are classifier calls conservative? | `classifier_tip_predictions.csv` | Few conflicts, reasonable weak calls, assigned nodes compatible with known training paths. |
| Are partial calls trustworthy? | `partial_classifications.csv`, `partial_region_diagnostics.csv` | Mapped records with enough observed informative sites, support, margin, and low conflict. |

## Multi-Criterion Scorecard

Do not reduce tuning to one metric. Score each parameter set as better, similar, or worse than baseline for each criterion:

| Criterion | Better means |
|---|---|
| SNP availability | More useful scored SNPs without relying on rare or poorly observed alleles. |
| Node coverage | More eligible nodes helped by selected panel windows. |
| Panel compactness | Similar or better signal with equal or fewer/lower-width windows. |
| Missingness | Lower missing fraction in selected windows. |
| Validation signal | Higher `mean_observed_informative_sites`, `mean_observed_informative_nodes`, and `resolved_fraction`. |
| Conflict control | Fewer `conflicting` outcomes and lower off-path support. |
| Weak/no-call control | Weak/unclassified calls are explainable by low evidence, not caused by overly strict settings everywhere. |
| Partial-genome caution | Partial calls are resolved only when mapping and observed informative sites are adequate. |
| Reproducibility | Same config and seed reproduce the same outputs. |

## Experiment Log Template

Copy [tuning_experiment_template.md](tuning_experiment_template.md) for each experiment. At minimum record:

- Config file used.
- Exact parameter changes.
- Command run.
- Output snapshot location.
- Scorecard comparison against baseline.
- Decision: keep, reject, or needs review.
