`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_tip_descendants <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node <= n_tip) {
    return(node)
  }
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    stop("Package 'phangorn' is required to derive clade tips.", call. = FALSE)
  }
  phangorn::Descendants(tree, node, type = "tips")[[1L]]
}

generate_clade_masks <- function(tree, min_clade_size = 2L, max_clade_frac = 0.95) {
  n_tip <- length(tree$tip.label)
  internal_nodes <- seq.int(n_tip + 1L, n_tip + tree$Nnode)
  root_node <- n_tip + 1L
  candidate_nodes <- setdiff(internal_nodes, root_node)

  all_clade_mask <- matrix(FALSE, nrow = length(candidate_nodes), ncol = n_tip)
  all_node_id <- integer(length(candidate_nodes))
  all_clade_tips <- vector("list", length(candidate_nodes))
  keep <- logical(length(candidate_nodes))

  max_clade_size <- floor(max_clade_frac * n_tip)
  for (i in seq_along(candidate_nodes)) {
    node <- candidate_nodes[[i]]
    tips <- get_tip_descendants(tree, node)
    clade_size <- length(tips)

    if (clade_size >= min_clade_size && clade_size <= max_clade_size) {
      all_clade_mask[i, tips] <- TRUE
      all_node_id[i] <- node
      all_clade_tips[[i]] <- tips
      keep[i] <- TRUE
    }
  }

  all_clade_mask <- all_clade_mask[keep, , drop = FALSE]
  colnames(all_clade_mask) <- tree$tip.label
  rownames(all_clade_mask) <- as.character(all_node_id[keep])

  list(
    all_node_id = all_node_id[keep],
    all_clade_tips = all_clade_tips[keep],
    all_clade_mask = all_clade_mask
  )
}

normalize_target_clade_specs <- function(target_clades) {
  if (is.null(target_clades)) {
    return(list())
  }
  if (!is.list(target_clades)) {
    stop("analysis.target_clades must be a list of target clade definitions.", call. = FALSE)
  }
  if (length(target_clades) == 0L) {
    return(list())
  }

  if (is.list(target_clades[[1L]])) {
    return(target_clades)
  }

  lapply(target_clades, function(node_id) {
    list(node_id = node_id)
  })
}

resolve_target_clade_node_id <- function(tree, target_spec) {
  if (!is.list(target_spec)) {
    stop("Each target clade definition must be a list.", call. = FALSE)
  }

  if (!is.null(target_spec$node_id)) {
    node_id <- suppressWarnings(as.integer(target_spec$node_id)[1L])
    if (is.na(node_id)) {
      stop("Target clade node_id must be an integer.", call. = FALSE)
    }
    return(list(node_id = node_id, resolution = "node_id"))
  }

  tip_labels <- target_spec$tip_labels %||% target_spec$mrca_tips %||% target_spec$anchor_tips
  if (!is.null(tip_labels)) {
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package 'ape' is required to resolve tip-based target clades.", call. = FALSE)
    }
    if (length(tip_labels) < 2L) {
      stop("Tip-based target clades require at least two tip labels.", call. = FALSE)
    }
    node_id <- ape::getMRCA(tree, tip_labels)
    if (is.na(node_id)) {
      stop(
        "Could not resolve a unique MRCA for target clade tips: ",
        paste(tip_labels, collapse = ", "),
        call. = FALSE
      )
    }
    return(list(node_id = as.integer(node_id), resolution = "mrca_tips"))
  }

  stop(
    "Target clade definition must provide node_id for this stage 01 implementation.",
    call. = FALSE
  )
}

resolve_target_clades <- function(tree, all_clades, target_clades,
                                  min_posterior_support = 0.90,
                                  min_clade_size = 2L,
                                  all_node_metadata = NULL) {
  specs <- normalize_target_clade_specs(target_clades)
  if (length(specs) == 0L) {
    stop("analysis.target_clades must define at least one target clade.", call. = FALSE)
  }
  if (is.null(all_node_metadata)) {
    all_node_metadata <- calculate_node_metadata(tree, all_clades)
  }
  if (!all(c("node_id", "clade_size", "posterior_support") %in% names(all_node_metadata))) {
    stop("all_node_metadata must include node_id, clade_size, and posterior_support.", call. = FALSE)
  }

  eligible_node_ids <- all_node_metadata$node_id
  resolved_node_id <- integer(length(specs))
  resolved_method <- character(length(specs))
  for (i in seq_along(specs)) {
    resolved <- resolve_target_clade_node_id(tree, specs[[i]])
    node_id <- resolved$node_id
    if (length(node_id) != 1L || is.na(node_id)) {
      stop("Each target clade must resolve to exactly one internal node.", call. = FALSE)
    }
    if (!(node_id %in% eligible_node_ids)) {
      stop("Requested target clade node ", node_id, " is not an eligible diagnostic clade.", call. = FALSE)
    }

    row <- match(node_id, all_node_metadata$node_id)
    clade_size <- all_node_metadata$clade_size[[row]]
    posterior_support <- all_node_metadata$posterior_support[[row]]
    if (clade_size < min_clade_size) {
      stop(
        "Requested target clade node ",
        node_id,
        " has only ",
        clade_size,
        " tips, which is below the configured minimum target clade size of ",
        min_clade_size,
        ".",
        call. = FALSE
      )
    }
    if (posterior_support < min_posterior_support) {
      stop(
        "Requested target clade node ",
        node_id,
        " has posterior support ",
        sprintf("%.6f", posterior_support),
        ", which is below the configured minimum of ",
        sprintf("%.6f", min_posterior_support),
        ".",
        call. = FALSE
      )
    }

    resolved_node_id[[i]] <- node_id
    resolved_method[[i]] <- resolved$resolution
  }

  if (anyDuplicated(resolved_node_id)) {
    stop("Target clades must resolve to unique node IDs.", call. = FALSE)
  }

  all_index <- match(resolved_node_id, all_clades$all_node_id)
  if (anyNA(all_index)) {
    stop("Target clade resolution failed for one or more node IDs.", call. = FALSE)
  }

  target_clade_mask <- all_clades$all_clade_mask[all_index, , drop = FALSE]
  target_clade_tips <- all_clades$all_clade_tips[all_index]

  rownames(target_clade_mask) <- as.character(resolved_node_id)

  list(
    target_clade_spec = specs,
    target_node_id = resolved_node_id,
    target_resolution = resolved_method,
    target_clade_tips = target_clade_tips,
    target_clade_mask = target_clade_mask
  )
}

extract_node_posterior_support <- function(tree) {
  n_node <- tree$Nnode
  if (is.null(n_node) || n_node < 1L) {
    stop("Tree must contain internal nodes to extract posterior support.", call. = FALSE)
  }

  labels <- tree$node.label
  if (is.null(labels) || length(labels) == 0L) {
    stop(
      "Tree is missing internal-node posterior support. Expected numeric values in tree$node.label.",
      call. = FALSE
    )
  }
  if (length(labels) != n_node) {
    stop(
      "Tree node labels must contain one posterior support value per internal node; expected ",
      n_node,
      " values but found ",
      length(labels),
      ".",
      call. = FALSE
    )
  }

  support <- suppressWarnings(as.numeric(labels))
  if (anyNA(support)) {
    stop(
      "Malformed posterior support values in tree$node.label. Expected numeric posterior support for every internal node.",
      call. = FALSE
    )
  }
  if (any(!is.finite(support))) {
    stop("Posterior support values must be finite.", call. = FALSE)
  }
  if (any(support < 0 | support > 1)) {
    warning("Posterior support values fall outside [0, 1] in the supplied tree.", call. = FALSE)
  }

  support
}

calculate_node_metadata <- function(tree, clades) {
  n_tip <- length(tree$tip.label)
  posterior_support <- extract_node_posterior_support(tree)
  node_id <- clades$all_node_id %||% clades$target_node_id %||% clades$node_id
  clade_mask <- clades$all_clade_mask %||% clades$target_clade_mask %||% clades$target_mask

  if (is.null(node_id) || is.null(clade_mask)) {
    stop("clades must contain node ids and a clade mask.", call. = FALSE)
  }
  if (nrow(clade_mask) == 0L) {
    return(data.frame(
      node_index = integer(),
      node_id = integer(),
      posterior_support = numeric(),
      clade_size = integer(),
      complement_size = integer(),
      depth = numeric(),
      balance = numeric(),
      weight = numeric()
    ))
  }

  depths <- ape::node.depth.edgelength(tree)
  if (is.null(tree$edge.length)) {
    depths <- ape::node.depth(tree) - 1
  }

  clade_size <- rowSums(clade_mask)
  complement_size <- n_tip - clade_size
  balance <- (2 * pmin(clade_size, complement_size)) / n_tip
  node_depth <- depths[node_id]
  posterior_by_node <- posterior_support[node_id - n_tip]
  if (anyNA(posterior_by_node)) {
    stop("Posterior support could not be aligned to one or more internal node ids.", call. = FALSE)
  }

  data.frame(
    node_index = seq_along(node_id),
    node_id = node_id,
    posterior_support = as.numeric(posterior_by_node),
    clade_size = as.integer(clade_size),
    complement_size = as.integer(complement_size),
    depth = as.numeric(node_depth),
    balance = as.numeric(balance),
    # Legacy/deprecated: retained for backward compatibility only.
    weight = as.numeric((1 / (node_depth + 1)) * balance)
  )
}
