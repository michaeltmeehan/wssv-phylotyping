# Codex Instructions

This repository is for a staged WSSV phylotyping pipeline. Keep changes milestone-scoped and avoid implementing future analysis logic before the relevant milestone.

Project conventions:

- Use R for analysis scripts and reusable functions.
- Keep workflow scripts numbered in `scripts/`.
- Keep reusable functions in `R/`.
- Treat `data/raw/` as user-supplied input storage. Do not modify or generate raw data files.
- Write generated analysis outputs under `data/processed/` or `outputs/`.
- Do not commit generated tables, figures, models, or R serialized objects unless a future milestone explicitly changes that policy.
- Prefer small, testable functions once implementation begins.
- Keep documentation updated when workflow files, parameters, or outputs change.

