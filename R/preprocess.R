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

  tree <- tryCatch(ape::read.tree(path), error = function(e) NULL)
  if (!is.null(tree) && length(tree$tip.label) > 0L) {
    return(tree)
  }

  tree <- tryCatch(ape::read.nexus(path), error = function(e) NULL)
  if (!is.null(tree) && length(tree$tip.label) > 0L) {
    return(tree)
  }

  stop("Failed to read tree as Newick or Nexus: ", path, call. = FALSE)
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
