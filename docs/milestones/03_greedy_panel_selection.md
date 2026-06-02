# Milestone 3: Greedy Panel Selection

## Goal

Select complementary marker windows that jointly cover informative MCC-tree nodes while limiting redundancy.

## Expected Inputs

- Window scoring table.
- Window-by-node informativeness matrix.
- Selection constraints such as maximum panel size, window spacing, and minimum marginal gain.

## Expected Outputs

- Selected marker panel table.
- Marginal gain and cumulative score summaries.
- Node coverage table for each panel size.

Implemented outputs are written by `scripts/05_select_marker_panel.R` under
`outputs/tables/`:

- `selected_panel.csv`
- `selected_panel_steps.csv`
- `panel_node_coverage.csv`
- `panel_summary.csv`

The script reads `analysis.panel_selection` from `config/config.yml`, including
requested panel sizes, minimum marginal gain, retained-window SNP and score
thresholds, optional missingness and width filters, selected-window overlap or
spacing constraints, and whether repeated node coverage is capped or partially
credited.

## Success Criteria

- The selected panel improves cumulative node coverage at each step.
- Redundant neighbouring windows are avoided when configured.
- Panel summaries explain which MCC-tree nodes are newly covered by each selected window.
