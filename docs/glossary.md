# Glossary

| Term | Definition |
|---|---|
| Node | A branching point in the tree. In this repository, eligible internal nodes represent clades that may be scored and classified. |
| Tip | A terminal genome in the tree, usually one complete WSSV genome from the input alignment. |
| Clade | A group containing a node and all descendant tips below that node. |
| Marker | A SNP or genomic window used as evidence for classification. |
| SNP | Single nucleotide polymorphism; an alignment site where observed genomes vary. |
| Informative site | A SNP that passes filters and helps distinguish at least one eligible tree node. |
| Window | A genomic interval with start and end alignment coordinates. Windows collect one or more informative SNPs. |
| Panel | A selected set of marker windows intended to provide complementary classification evidence. |
| Compatible evidence | Evidence supporting nodes that can exist together on one nested tree path. |
| Conflicting evidence | Evidence supporting incompatible nodes on different branches. |
| Off-path support | Support for nodes outside the expected or assigned path. In training-tip diagnostics this means outside the known MCC path for that tip. |
| Resolved classification | A call that passes configured support, conflict, and margin thresholds. It still needs cautious interpretation. |
| Weak classification | Evidence was present, but not strong enough to pass thresholds. Output status is usually `weak_support`. |
| Conflicting classification | Strong evidence supports incompatible nodes or close competitors. Output status is `conflicting`. |
| Unclassified | Plain-language label for any record without a resolved assignment. Check the exact `status` and `reason`. |
| Missing | A gap, ambiguous base, absent base, or unobserved coordinate. Missing values do not contribute support. |
| Partial genome | An incomplete sequence record that covers only part of the whole-genome alignment. |
| MCC tree | Maximum clade credibility tree used as the reference tree for clade definitions. |
| Support | Fraction of observed rule weight or count that agrees with a node. |
| Conflict | Fraction of observed rule weight or count that disagrees with a node. |
| Support margin | Gap between the best supported assignment and competing incompatible evidence. |
| Opportunistic classification | Partial-genome classification using all scored SNPs observed in the mapped interval rather than only the selected marker panel. |
| Baseline run | The unchanged reference run used to compare tuning experiments. |
