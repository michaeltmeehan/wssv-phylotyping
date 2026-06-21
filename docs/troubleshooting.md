# Troubleshooting

Start with the first failing stage and the most recent parameter change. If outputs look surprising, check whether they were regenerated after the config changed.

| Symptom | Likely cause | Suggested check |
|---|---|---|
| Missing input file error | Config path points to a file that does not exist. | Check `paths.alignment`, `paths.tree`, and `analysis.partials.input_dir` in `config/config.yml`. |
| Many `NA` values or unclassified partial outputs | Partial records are unmapped, missing informative sites, or have ambiguous/gapped states at informative sites. | Inspect `partial_region_diagnostics.csv` and `observed_informative_sites`. |
| Many `weak_support` classifications | `min_support`, `support_margin`, or minimum observed-site thresholds may be too strict, or evidence is genuinely sparse. | Compare `assigned_support`, `support_margin`, `observed_informative_sites`, and validation summaries. |
| Many `conflicting` classifications | Markers support incompatible branches, conflict threshold is too permissive, or input sequence/mapping may be problematic. | Inspect `conflict_reason`, `conflict_sites`, `off_path_nodes_supported`, and mapping diagnostics. |
| Unexpectedly few eligible nodes | `min_clade_size` too high, `max_clade_frac` too low, or many tips dropped during alignment/tree matching. | Check stage 01 messages and retained tip count. |
| Unexpectedly few scored SNPs | SNP filters are too strict or alignment has little usable variation. | Lower one of `min_total_obs`, `min_side_obs`, `min_site_maf`, or `min_gain_norm` cautiously; inspect `site_summary.csv`. |
| Unexpectedly few candidate windows | Window minimum informative SNPs too strict or SNP scoring produced few informative sites. | Inspect `site_summary.csv`, `window_summary.csv`, and `analysis.windows.min_informative_snps`. |
| No panel satisfying constraints | Panel filters, overlap, distance, or marginal-gain settings are too strict. | Check `window_summary.csv`, relax `min_total_weighted_gain`, `max_allowed_overlap`, `min_distance`, or `min_marginal_gain`. |
| Selected panel stops early | No remaining compatible window has positive/allowed marginal gain. | Inspect `selected_panel_steps.csv` and `panel_summary.csv`. |
| Partial FASTA files present but records read as zero | File extensions are not included or files are not valid FASTA. | Check `analysis.partials.fasta_extensions` and FASTA headers. |
| Many partial records unmapped | Mapping mode or pairwise thresholds are too strict, sequences are unrelated, or orientation/format is wrong. | Inspect `mapping_mode`, `mapping_identity`, `mapping_query_coverage`, and `mapping_note`. |
| Stale outputs from a previous run | Only later stages were rerun after changing an upstream parameter. | Rerun from the earliest affected stage or use `Rscript scripts/run_all.R`. |
| Report says outputs are missing | The report only summarises existing files and does not rerun missing stages. | Run the missing upstream script named in the report or script error. |
| R package error for `yaml`, `ape`, or `phangorn` | Required package is not installed in the active R library. | Install the package and restart R if needed. |
| Pairwise local mapping package error | Neither `pwalign` nor suitable `Biostrings` pairwise alignment functions are available. | Install Bioconductor packages `pwalign` and `Biostrings`. |
| Output CSV cannot be overwritten | File may be open in Excel or another program. | Close the file. Stage 08 has fallback timestamped writes for some locked partial outputs. |

## Quick Recovery Steps

1. Close any open output CSV files.
2. Confirm raw input paths in `config/config.yml`.
3. Rerun the earliest affected stage.
4. If unsure, rerun the full workflow:

   ```sh
   Rscript scripts/run_all.R
   ```

5. Compare the new `outputs/reports/wssv_phylotyping_report.md` with your experiment log.
