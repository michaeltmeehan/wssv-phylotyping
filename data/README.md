# Data Directory

`data/raw/` is for user-supplied inputs such as the whole-genome alignment, MCC tree, partial sequences, and metadata. Raw data should be treated as immutable and should not be generated or modified by workflow scripts.

`data/processed/` is for derived intermediate files created by the pipeline, such as precomputed alignment/tree objects and scoring tables. These files are generated outputs and are ignored by Git by default.

