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

## Success Criteria

- The selected panel improves cumulative node coverage at each step.
- Redundant neighbouring windows are avoided when configured.
- Panel summaries explain which MCC-tree nodes are newly covered by each selected window.

