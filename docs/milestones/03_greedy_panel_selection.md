# Milestone 3: Candidate Amplicon Panel Selection

## Goal

Compare small candidate amplicon panels from the stage-04 outputs using transparent multi-objective criteria.

## Expected Inputs

- Candidate amplicon table from stage 04.
- Candidate-amplicon target-score table from stage 04.
- Optional candidate-amplicon SNP table from stage 04.
- Selection constraints such as maximum panel size, overlap limits, and candidate-feasibility filters.

## Expected Outputs

- Candidate panel comparison table.
- Panel-by-target score table.
- Pareto-optimal panel table.
- Top-ranked selected panel table.
- Panel-size comparison summary.

Implemented outputs are written by `scripts/05_select_panel.R` under
`outputs/tables/`:

- `panel_candidates.csv`
- `panel_target_scores.csv`
- `pareto_panels.csv`
- `selected_panel.csv`
- `panel_size_comparison.csv`
- `panel_summary.csv`

The script reads `analysis.panel_selection` from `config/config.yml`, including
panel-size limits, candidate filtering thresholds, overlap and spacing
constraints, local-redundancy reduction settings, and the evidence thresholds
used for support summaries.

## Success Criteria

- Candidate panels are compared deterministically.
- Pareto-optimal panels are retained and reported.
- Redundant neighbouring amplicons are avoided when configured.
- Panel summaries explain which primary targets each selected panel supports.
