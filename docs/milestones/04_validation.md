# Milestone 4: Validation

## Goal

Validate whether selected marker windows recover useful tree placement from complete genomes and synthetic partial genomes.

## Expected Inputs

- Selected marker panel.
- Complete-genome alignment.
- MCC tree.
- Trained or provisional classification rules.

## Expected Outputs

- Selected-panel signal summaries by panel size and tip.
- Artificial fragment signal summaries for configured random fixed-length fragments.
- Top-weighted and random matched-window baseline comparisons.
- Conservative true-path support summaries where existing site-node rules apply.

This milestone reports observed informative-site and informative-node signal.
It does not train the final classifier or claim final classification accuracy.

## Success Criteria

- Selected windows outperform random-window baselines.
- Validation reports major-clade accuracy, deepest supported node, overclassification, unresolved calls, and confidence calibration.
- Limitations are reported clearly for short or low-signal fragments.
