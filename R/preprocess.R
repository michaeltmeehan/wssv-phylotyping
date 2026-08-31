read_alignment_dnastringset <- function(path) {
  if (!file.exists(path)) {
    stop("Alignment file does not exist: ", path, call. = FALSE)
  }

  if (requireNamespace("Biostrings", quietly = TRUE)) {
    alignment <- Biostrings::readDNAStringSet(path)
    validate_alignment_lengths(alignment)
    return(alignment)
  }

  warning("Package 'Biostrings' is not installed; using a minimal FASTA reader.", call. = FALSE)
  alignment <- read_fasta_character(path)
  validate_alignment_lengths(alignment)
  alignment
}

read_fasta_character <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0L) {
    stop("No FASTA headers found in alignment file: ", path, call. = FALSE)
  }

  seqs <- character(length(header_idx))
  names(seqs) <- sub("^>\\s*", "", lines[header_idx])
  for (i in seq_along(header_idx)) {
    start <- header_idx[[i]] + 1L
    end <- if (i == length(header_idx)) length(lines) else header_idx[[i + 1L]] - 1L
    if (start > end) {
      seqs[[i]] <- ""
    } else {
      seqs[[i]] <- paste0(gsub("\\s+", "", lines[start:end]), collapse = "")
    }
  }
  seqs
}

read_tree_any <- function(path) {
  if (!file.exists(path)) {
    stop("Tree file does not exist: ", path, call. = FALSE)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required to read tree files.", call. = FALSE)
  }

  has_beast_support <- file_contains_beast_posterior_support(path)

  tree <- tryCatch(ape::read.tree(path), error = function(e) NULL)
  if (!is.null(tree) && length(tree$tip.label) > 0L && (!has_beast_support || tree_has_internal_node_labels(tree))) {
    return(tree)
  }

  tree <- tryCatch(ape::read.nexus(path), error = function(e) NULL)
  if (!is.null(tree) && length(tree$tip.label) > 0L && (!has_beast_support || tree_has_internal_node_labels(tree))) {
    return(tree)
  }

  if (has_beast_support) {
    if (!requireNamespace("treeio", quietly = TRUE)) {
      stop(
        "Tree file contains BEAST posterior annotations but package 'treeio' is required to read them: ",
        path,
        call. = FALSE
      )
    }

    beast_tree <- tryCatch(treeio::read.beast(path), error = function(e) e)
    if (inherits(beast_tree, "error")) {
      stop("Failed to read BEAST-annotated tree: ", beast_tree$message, call. = FALSE)
    }

    tree <- treeio::as.phylo(beast_tree)
    tree$node.label <- beast_internal_posterior_labels(beast_tree, tree)
    return(tree)
  }

  stop("Failed to read tree as Newick or Nexus: ", path, call. = FALSE)
}

file_contains_beast_posterior_support <- function(path) {
  lines <- readLines(path, warn = FALSE)
  any(grepl("posterior=", lines, fixed = TRUE))
}

tree_has_internal_node_labels <- function(tree) {
  if (is.null(tree$node.label) || tree$Nnode <= 0L || length(tree$node.label) != tree$Nnode) {
    return(FALSE)
  }
  !anyNA(suppressWarnings(as.numeric(tree$node.label)))
}

beast_internal_posterior_labels <- function(beast_tree, phylo_tree) {
  if (is.null(beast_tree) || is.null(phylo_tree)) {
    stop("Internal tree conversion failed while attaching posterior support.", call. = FALSE)
  }

  data <- as.data.frame(beast_tree@data)
  if (!("posterior" %in% names(data))) {
    stop("BEAST tree data do not contain a posterior column.", call. = FALSE)
  }
  if (!("node" %in% names(data))) {
    stop("BEAST tree data do not contain node identifiers.", call. = FALSE)
  }
  if (is.null(phylo_tree$Nnode) || phylo_tree$Nnode < 1L) {
    stop("Tree must contain at least one internal node.", call. = FALSE)
  }

  node_ids <- suppressWarnings(as.integer(as.character(data$node)))
  internal <- data[!is.na(node_ids) & node_ids > length(phylo_tree$tip.label), c("node", "posterior")]
  if (nrow(internal) != phylo_tree$Nnode) {
    stop(
      "Expected posterior support for ",
      phylo_tree$Nnode,
      " internal nodes but found ",
      nrow(internal),
      " in the BEAST tree.",
      call. = FALSE
    )
  }

  internal <- internal[order(as.integer(as.character(internal$node))), , drop = FALSE]
  posterior <- suppressWarnings(as.numeric(internal$posterior))
  if (anyNA(posterior)) {
    bad_nodes <- internal$node[is.na(posterior)]
    stop(
      "Malformed posterior support values in BEAST tree for node(s): ",
      paste(bad_nodes, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(!is.finite(posterior))) {
    stop("Posterior support values must be finite.", call. = FALSE)
  }
  if (any(posterior < 0 | posterior > 1)) {
    warning("Posterior support values fall outside [0, 1] in the supplied tree.", call. = FALSE)
  }

  sprintf("%.15g", posterior)
}

match_prune_reorder_alignment_tree <- function(alignment, tree) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for tree pruning.", call. = FALSE)
  }

  seqs <- alignment_to_character(alignment)
  validate_alignment_lengths(seqs)
  if (is.null(tree$tip.label) || length(tree$tip.label) < 2L) {
    stop("Tree contains fewer than 2 tips.", call. = FALSE)
  }

  common <- intersect(tree$tip.label, names(seqs))
  if (length(common) < 2L) {
    stop("Too few shared labels between tree tips and alignment names.", call. = FALSE)
  }

  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common))
  reordered <- alignment[pruned_tree$tip.label]
  stopifnot(identical(names(reordered), pruned_tree$tip.label))

  list(
    alignment = reordered,
    tree = pruned_tree,
    matched_tips = pruned_tree$tip.label,
    dropped_alignment = setdiff(names(seqs), common),
    dropped_tree = setdiff(tree$tip.label, common)
  )
}

make_site_map <- function(aln_int) {
  data.frame(
    site = seq_len(ncol(aln_int)),
    alignment_position = seq_len(ncol(aln_int))
  )
}
