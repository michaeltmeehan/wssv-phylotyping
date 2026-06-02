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

  target_mask <- matrix(FALSE, nrow = length(candidate_nodes), ncol = n_tip)
  node_id <- integer(length(candidate_nodes))
  clade_tips <- vector("list", length(candidate_nodes))
  keep <- logical(length(candidate_nodes))

  max_clade_size <- floor(max_clade_frac * n_tip)
  for (i in seq_along(candidate_nodes)) {
    node <- candidate_nodes[[i]]
    tips <- get_tip_descendants(tree, node)
    clade_size <- length(tips)

    if (clade_size >= min_clade_size && clade_size <= max_clade_size) {
      target_mask[i, tips] <- TRUE
      node_id[i] <- node
      clade_tips[[i]] <- tips
      keep[i] <- TRUE
    }
  }

  target_mask <- target_mask[keep, , drop = FALSE]
  colnames(target_mask) <- tree$tip.label
  rownames(target_mask) <- as.character(node_id[keep])

  list(
    node_id = node_id[keep],
    clade_tips = clade_tips[keep],
    target_mask = target_mask
  )
}

calculate_node_metadata <- function(tree, clades) {
  n_tip <- length(tree$tip.label)
  if (nrow(clades$target_mask) == 0L) {
    return(data.frame(
      node_index = integer(),
      node_id = integer(),
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

  clade_size <- rowSums(clades$target_mask)
  complement_size <- n_tip - clade_size
  balance <- (2 * pmin(clade_size, complement_size)) / n_tip
  node_depth <- depths[clades$node_id]
  depth_weight <- 1 / (node_depth + 1)

  data.frame(
    node_index = seq_along(clades$node_id),
    node_id = clades$node_id,
    clade_size = as.integer(clade_size),
    complement_size = as.integer(complement_size),
    depth = as.numeric(node_depth),
    balance = as.numeric(balance),
    weight = as.numeric(depth_weight * balance)
  )
}
