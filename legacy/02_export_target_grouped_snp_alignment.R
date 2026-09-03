#!/usr/bin/env Rscript

# Stage 02: export a target-grouped SNP-only alignment for manual inspection.
#
# Input:
# - data/processed/precomputed.rds from stage 01
#
# Outputs:
# - data/processed/target_grouped_snps.fasta
# - data/processed/target_grouped_snp_positions.csv
# - data/processed/target_grouped_taxa.csv
#
# The SNP-only alignment contains all polymorphic sites identified during
# preprocessing. Taxa are grouped by configured target-clade membership and
# retain their MCC-tree order within each target clade. Taxa not belonging to
# any target clade, if present, are placed last.
#
# No SNP gain, missingness, MAF, or other discriminatory filtering is applied
# here beyond the polymorphic-site definition used in stage 01. The purpose of
# this output is descriptive/manual exploration of SNP and multi-SNP patterns.
#
# Run directly from the repository root after stage 01.

source("R/encoding.R")

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

precomputed_path <- file.path(processed_dir, "precomputed.rds")
if (!file.exists(precomputed_path)) {
  stop(
    "Stage-01 output does not exist: ",
    precomputed_path,
    "\nRun scripts/01_preprocess_alignment_tree.R first.",
    call. = FALSE
  )
}

pre <- readRDS(precomputed_path)

required <- c(
  "tree",
  "alignment_names",
  "aln_int",
  "site_map",
  "target_clade_mask",
  "target_node_metadata",
  "polymorphic_sites"
)

missing_fields <- setdiff(required, names(pre))
if (length(missing_fields) > 0L) {
  stop(
    "precomputed.rds is missing required fields: ",
    paste(missing_fields, collapse = ", "),
    call. = FALSE
  )
}

aln_int <- pre$aln_int
polymorphic_sites <- as.integer(pre$polymorphic_sites)
target_mask <- pre$target_clade_mask
target_metadata <- pre$target_node_metadata

if (length(polymorphic_sites) == 0L) {
  stop("No polymorphic sites were found in the stage-01 alignment.", call. = FALSE)
}

if (!identical(colnames(target_mask), rownames(aln_int))) {
  stop(
    "Target-clade mask columns do not match alignment taxa/order.",
    call. = FALSE
  )
}

if (!all(target_metadata$node_id %in% as.integer(rownames(target_mask)))) {
  stop(
    "Target-node metadata and target-clade mask are inconsistent.",
    call. = FALSE
  )
}

# -------------------------------------------------------------------------
# Assign each taxon to its target clade
# -------------------------------------------------------------------------

membership_count <- colSums(target_mask)

if (any(membership_count > 1L)) {
  overlapping_taxa <- colnames(target_mask)[membership_count > 1L]
  stop(
    "Target clades overlap. The following taxa belong to more than one ",
    "configured target clade: ",
    paste(overlapping_taxa, collapse = ", "),
    call. = FALSE
  )
}

target_nodes <- as.integer(target_metadata$node_id)

target_membership <- rep(NA_integer_, ncol(target_mask))
names(target_membership) <- colnames(target_mask)

for (node in target_nodes) {
  mask_row <- match(as.character(node), rownames(target_mask))
  target_membership[target_mask[mask_row, ]] <- node
}

# The alignment was reordered to MCC-tree tip order during stage 01.
# Group taxa by target clade while retaining this tree order within each group.
grouped_taxa <- unlist(
  lapply(target_nodes, function(node) {
    pre$alignment_names[target_membership[pre$alignment_names] == node]
  }),
  use.names = FALSE
)

unassigned_taxa <- pre$alignment_names[
  is.na(target_membership[pre$alignment_names])
]

grouped_taxa <- c(grouped_taxa, unassigned_taxa)

if (!setequal(grouped_taxa, rownames(aln_int))) {
  stop(
    "Failed to construct a complete target-grouped taxon ordering.",
    call. = FALSE
  )
}

# -------------------------------------------------------------------------
# Extract the SNP-only alignment
# -------------------------------------------------------------------------

snp_aln_int <- aln_int[grouped_taxa, polymorphic_sites, drop = FALSE]

# Integer encoding:
#   A = 1, C = 2, G = 3, T = 4
#   anything else = 0
#
# Preserve non-ACGT states visibly as N in the diagnostic FASTA.
decode_sequence <- function(x) {
  alleles <- base_code_to_allele(x)
  alleles[is.na(alleles)] <- "N"
  paste0(alleles, collapse = "")
}

snp_sequences <- apply(snp_aln_int, 1L, decode_sequence)
names(snp_sequences) <- rownames(snp_aln_int)

# -------------------------------------------------------------------------
# Write SNP-only FASTA
# -------------------------------------------------------------------------

fasta_path <- file.path(processed_dir, "target_grouped_snps.fasta")

fasta_lines <- unlist(
  Map(
    function(name, sequence) {
      c(paste0(">", name), sequence)
    },
    names(snp_sequences),
    snp_sequences
  ),
  use.names = FALSE
)

writeLines(fasta_lines, fasta_path)


# -------------------------------------------------------------------------
# Write target-grouped SNP-only alignment as CSV
# -------------------------------------------------------------------------

decode_alleles <- function(x) {
  alleles <- base_code_to_allele(x)
  alleles[is.na(alleles)] <- "N"
  alleles
}

snp_aln_chr <- t(
  apply(
    snp_aln_int,
    1L,
    decode_alleles
  )
)

colnames(snp_aln_chr) <- as.character(polymorphic_sites)
rownames(snp_aln_chr) <- grouped_taxa

snp_alignment <- data.frame(
  taxon = grouped_taxa,
  target_node = unname(target_membership[grouped_taxa]),
  snp_aln_chr,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

snp_alignment$target_node[
  is.na(snp_alignment$target_node)
] <- "unassigned"

snp_alignment_path <- file.path(
  processed_dir,
  "target_grouped_snps.csv"
)

write.csv(
  snp_alignment,
  snp_alignment_path,
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------------------------------------------------
# Write SNP-column -> whole-genome alignment-coordinate map
# -------------------------------------------------------------------------

site_map_index <- match(polymorphic_sites, pre$site_map$site)

if (anyNA(site_map_index)) {
  stop(
    "One or more polymorphic sites could not be matched to site_map.",
    call. = FALSE
  )
}

snp_positions <- data.frame(
  snp_column = seq_along(polymorphic_sites),
  site = polymorphic_sites,
  alignment_position = pre$site_map$alignment_position[site_map_index],
  stringsAsFactors = FALSE
)

snp_positions_path <- file.path(
  processed_dir,
  "target_grouped_snp_positions.csv"
)

write.csv(
  snp_positions,
  snp_positions_path,
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Write taxon ordering and target-clade membership
# -------------------------------------------------------------------------

taxon_groups <- data.frame(
  alignment_order = seq_along(grouped_taxa),
  taxon = grouped_taxa,
  target_node = unname(target_membership[grouped_taxa]),
  target_clade = ifelse(
    is.na(target_membership[grouped_taxa]),
    "unassigned",
    paste0("node_", target_membership[grouped_taxa])
  ),
  stringsAsFactors = FALSE
)

taxon_groups_path <- file.path(
  processed_dir,
  "target_grouped_taxa.csv"
)

write.csv(
  taxon_groups,
  taxon_groups_path,
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

message("Target-grouped SNP alignment export complete")
message("  taxa: ", nrow(snp_aln_int))
message("  polymorphic sites: ", ncol(snp_aln_int))
message("  target clades: ", length(target_nodes))
message("  target node IDs: ", paste(target_nodes, collapse = ", "))
message(
  "  target clade sizes: ",
  paste(
    vapply(
      target_nodes,
      function(node) sum(target_membership == node, na.rm = TRUE),
      integer(1L)
    ),
    collapse = ", "
  )
)
message("  unassigned taxa: ", length(unassigned_taxa))
message("  wrote: ", fasta_path)
message("  wrote: ", snp_positions_path)
message("  wrote: ", taxon_groups_path)