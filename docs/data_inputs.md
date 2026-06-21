# Data Inputs

Raw input files belong under `data/raw/`. These files are user-supplied and should not be modified by pipeline scripts.

## Whole-Genome Alignment

Default path:

```text
data/raw/alignment/41seqscollagen_Edited.fasta
```

Expected format:

- FASTA alignment.
- All sequences should be aligned to the same coordinate system and have the same length.
- Sequence names should match tree tip labels exactly, or enough should match after direct name comparison.
- Bases should mainly be `A`, `C`, `G`, `T`, gaps, or standard ambiguous characters. Ambiguous/gap states are treated as missing for scoring/classification.

Common mistakes:

- Unaligned FASTA used instead of aligned FASTA.
- Sequence names differ from tree labels because of spaces, prefixes, suffixes, or version numbers.
- Mixed coordinate systems or sequences with different lengths.
- Accidental editing of raw FASTA after outputs were generated, leaving stale outputs.

## MCC Tree

Default path:

```text
data/raw/tree/41SEQSUPDATEDEDITEDMCC
```

Expected format:

- Newick or Nexus tree readable by `ape::read.tree()` or `ape::read.nexus()`.
- Tip labels should correspond to alignment sequence names.
- The tree should represent the complete-genome set used for scoring.

The preprocessing stage prunes alignment and tree entries to the shared names. Check stage 01 messages for dropped alignment sequences and dropped tree tips.

## Partial Genomes

Default directory:

```text
data/raw/partials/
```

Accepted extensions by default:

```text
fa, fas, fasta, fna
```

Partial records may be:

- Full-length aligned records already in the alignment coordinate system.
- Ungapped subsequences that can be found exactly in the reference sequence.
- Sequences that can be mapped by pairwise local alignment when mapping thresholds pass.

For exact-substring and pairwise-local mapping, the reference sequence is the alignment row named by `analysis.partials.reference_id` in `config/config.yml`. The current default is `CN01_1994`. Treat this as a deliberate coordinate-reference choice; change it only if the project lead confirms a different alignment row should define partial-genome mapping coordinates.

Common mistakes:

- Reverse-complement or unrelated records that fail mapping.
- Very short sequences that map ambiguously or contain no informative sites.
- FASTA files with unsupported extensions.
- Headers that are not unique; the script uses unique sequence IDs internally, but clear headers make review easier.

## Metadata

Directory:

```text
data/raw/metadata/
```

Metadata files are reserved for sample, accession, collection, or source information. The current main workflow does not require metadata for scoring or classification. Keep metadata filenames and column names clear for future domain review.

## Naming Assumptions

- Alignment sequence names and tree tip labels are the critical names for complete-genome stages.
- Partial FASTA headers are used to create `sequence_id` and accession-like fields where possible.
- Output tables refer to tree nodes by numeric/string `node_id` values derived from the tree object, not by biological clade names.
