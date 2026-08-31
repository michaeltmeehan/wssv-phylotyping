# Milestone 6: Partial-Genome Classification

## Goal

Classify external partial WSSV sequences where enough informative signal overlaps
the selected marker panel or any scored informative SNPs in the mapped interval.
The workflow distinguishes designed-panel classification from opportunistic,
region-dependent assessment of arbitrary partial genomes.

## Expected Inputs

- Partial WSSV sequences.
- Partial sequence metadata.
- Reference coordinate mapping.
- Trained classifier.
- Full scored SNP tables (`outputs/tables/site_node_scores_all.rds` and the
  transitional alias `outputs/tables/site_node_scores.rds`) for opportunistic
  classification. Future primer-panel work should draw from
  `outputs/tables/site_node_scores_targets.rds`.

## Expected Outputs

- Partial classification table with explicit `classification_mode` and
  `classification_source`.
- Separate panel and opportunistic classification tables.
- Opportunistic node-evidence and mapped-region diagnostic tables.
- Mapping and coverage summaries.
- Informativeness report for partial sequences.
- Optional local-alignment mapping for partial sequences that do not match
  exactly to the training reference coordinates.

## Success Criteria

- Partial sequences are classified only when sufficient evidence is present.
- Panel mode preserves the trained selected-panel classifier behavior.
- Opportunistic mode uses all scored informative SNPs in the mapped interval,
  not only selected panel sites.
- Auto mode tries panel classification first and falls back to opportunistic
  classification only when no panel-informative sites are observed.
- Outputs distinguish strong, weak, unresolved, conflicting, and out-of-reference calls.
- Reports include informative sites used, mapped coordinates, assignment, classifier support, and status.
- Partial mapping reports strand, identity, query coverage, aligned length, and
  reference/alignment coordinate spans when local alignment is used.
- Reports warn that opportunistic calls are region-dependent and are not
  equivalent to designed-panel validation.
