# Milestone 6: Partial-Genome Classification

## Goal

Classify external partial WSSV sequences where enough informative signal overlaps the selected marker or SNP set.

## Expected Inputs

- Partial WSSV sequences.
- Partial sequence metadata.
- Reference coordinate mapping.
- Trained classifier.

## Expected Outputs

- Partial classification table.
- Mapping and coverage summaries.
- Informativeness report for partial sequences.

## Success Criteria

- Partial sequences are classified only when sufficient evidence is present.
- Outputs distinguish strong, weak, unresolved, conflicting, and out-of-reference calls.
- Reports include informative sites used, mapped coordinates, assignment, support, and status.

