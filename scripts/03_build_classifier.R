#!/usr/bin/env Rscript

# Stage 03: build the phylogeny-wide classifier structure.
#
# Input:
#   data/processed/precomputed.rds
#
# Output:
#   data/processed/classifier.rds
#
# The classifier uses every reference taxon as a mutually exclusive terminal
# state. Posterior probabilities over these terminal states therefore sum to 1.
#
# A complete tree-node x terminal-state membership matrix is constructed so
# that posterior probability for any tip or internal node can subsequently be
# obtained by summing posterior mass over its descendant terminal states.
#
# This stage deliberately does not:
#   - specify target clades;
#   - select SNPs or amplicons;
#   - assign likelihood/error parameters;
#   - classify query sequences.
#
# Those operations are handled downstream.

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required for tree operations.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

precomputed_path <- file.path(processed_dir, "precomputed.rds")

if (!file.exists(precomputed_path)) {
  stop(
    "Precomputed data not found: ", precomputed_path,
    "\nRun stage 01 first.",
    call. = FALSE
  )
}

precomputed <- readRDS(precomputed_path)

tree <- precomputed$tree
aln_int <- precomputed$aln_int
polymorphic_sites <- precomputed$polymorphic_sites
site_map <- precomputed$site_map


# -------------------------------------------------------------------------
# Validate tree/alignment correspondence
# -------------------------------------------------------------------------

if (!identical(tree$tip.label, rownames(aln_int))) {
  stop(
    "Tree-tip ordering does not match alignment ordering.",
    call. = FALSE
  )
}

n_tips <- ape::Ntip(tree)
n_internal <- ape::Nnode(tree)
n_nodes <- n_tips + n_internal

if (n_tips != nrow(aln_int)) {
  stop(
    "Number of tree tips does not match number of alignment sequences.",
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Define mutually exclusive terminal states
# -------------------------------------------------------------------------

terminal_states <- data.frame(
  state_id = seq_len(n_tips),
  node_id = seq_len(n_tips),
  taxon = tree$tip.label,
  prior = rep(1 / n_tips, n_tips),
  stringsAsFactors = FALSE
)

# Equal tip priors are used as a neutral initial choice.
#
# NOTE:
# Equal tip priors imply that an internal clade's prior probability is
# proportional to the number of sampled reference taxa it contains.
# This is convenient for the baseline implementation but should not
# automatically be interpreted as an epidemiological prior.


# -------------------------------------------------------------------------
# Store the observed alleles for each terminal state
# -------------------------------------------------------------------------

tip_alleles <- aln_int[
  ,
  polymorphic_sites,
  drop = FALSE
]

rownames(tip_alleles) <- tree$tip.label
colnames(tip_alleles) <- as.character(polymorphic_sites)


# -------------------------------------------------------------------------
# Construct complete node x terminal-state membership matrix
# -------------------------------------------------------------------------
#
# Rows correspond to every node in the tree:
#
#   1:Ntip(tree)                    = terminal tips
#   (Ntip(tree) + 1):...            = internal nodes
#
# Columns correspond to terminal states/reference taxa.
#
# Entry [v, i] is TRUE iff terminal state i descends from node v.
#
# This allows posterior probabilities for all nodes to be calculated as:
#
#   node_posterior = node_tip_mask %*% tip_posterior
#

node_tip_mask <- matrix(
  FALSE,
  nrow = n_nodes,
  ncol = n_tips,
  dimnames = list(
    as.character(seq_len(n_nodes)),
    tree$tip.label
  )
)

# Each tip contains only itself.
node_tip_mask[
  seq_len(n_tips),
  seq_len(n_tips)
] <- diag(n_tips) == 1L

# Internal nodes contain all descendant tips.
internal_nodes <- n_tips + seq_len(n_internal)

for (node in internal_nodes) {
  
  descendant_tips <- ape::extract.clade(
    tree,
    node = node
  )$tip.label
  
  tip_indices <- match(
    descendant_tips,
    tree$tip.label
  )
  
  if (anyNA(tip_indices)) {
    stop(
      "Could not match descendant tips for node ",
      node,
      ".",
      call. = FALSE
    )
  }
  
  node_tip_mask[
    as.character(node),
    tip_indices
  ] <- TRUE
}


# -------------------------------------------------------------------------
# Construct node metadata
# -------------------------------------------------------------------------

node_metadata <- data.frame(
  node_id = seq_len(n_nodes),
  node_type = c(
    rep("tip", n_tips),
    rep("internal", n_internal)
  ),
  clade_size = rowSums(node_tip_mask),
  stringsAsFactors = FALSE
)

node_metadata$tip_label <- NA_character_
node_metadata$tip_label[seq_len(n_tips)] <- tree$tip.label

node_metadata$posterior_support <- NA_real_

if (!is.null(tree$node.label)) {
  
  support <- suppressWarnings(
    as.numeric(tree$node.label)
  )
  
  if (length(support) == n_internal) {
    node_metadata$posterior_support[internal_nodes] <- support
  }
}


# -------------------------------------------------------------------------
# Validate probability-propagation structure
# -------------------------------------------------------------------------

# Every terminal state must belong to itself and only itself at tip level.
tip_mask <- node_tip_mask[seq_len(n_tips), , drop = FALSE]

if (!all(rowSums(tip_mask) == 1L)) {
  stop(
    "Terminal-state membership matrix is invalid.",
    call. = FALSE
  )
}

# Identify the root and verify that it contains every terminal state.
root_node <- setdiff(
  tree$edge[, 1],
  tree$edge[, 2]
)

if (length(root_node) != 1L) {
  stop(
    "Could not identify a unique tree root.",
    call. = FALSE
  )
}

if (!all(node_tip_mask[as.character(root_node), ])) {
  stop(
    "Root node does not contain every terminal state.",
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Save classifier structure
# -------------------------------------------------------------------------

classifier <- list(
  tree = tree,
  
  # Mutually exclusive probability states
  terminal_states = terminal_states,
  
  # Reference alleles available to the likelihood model
  polymorphic_sites = polymorphic_sites,
  site_map = site_map,
  tip_alleles = tip_alleles,
  
  # Mapping from terminal probabilities to arbitrary tree-node probabilities
  node_tip_mask = node_tip_mask,
  node_metadata = node_metadata,
  
  root_node = root_node
)

output_path <- file.path(
  processed_dir,
  "classifier.rds"
)

saveRDS(
  classifier,
  output_path
)


# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

message("Classifier structure complete")
message("  terminal states: ", n_tips)
message("  internal nodes: ", n_internal)
message("  total tree nodes: ", n_nodes)
message("  polymorphic sites: ", length(polymorphic_sites))
message("  root node: ", root_node)
message("  wrote: ", output_path)