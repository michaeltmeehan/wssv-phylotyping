# Outputs

Generated analysis products should be written under `outputs/`.

- `outputs/tables/`: CSV/TSV summaries and reports.
- `outputs/figures/`: generated plots and visual diagnostics.
- `outputs/models/`: trained classifiers and model objects.

Generated outputs are ignored by Git by default. Keep only lightweight placeholders such as `.gitkeep` files under version control unless a future milestone explicitly identifies a result that should be committed.

The stage 04-09 outputs are intentionally kept for compatibility with the legacy assay-oriented workflow while the stage 01-03 cleanup settles the new names and summaries.
