classifier_allele_codes <- c(A = 1L, C = 2L, G = 3L, T = 4L)

#' Train the conservative tree-path classifier.
#'
#' Builds node-specific SNP rules from the selected panel or full informative
#' catalogue, applies minimum rule-count filters, and stores classification
#' thresholds used by complete and partial genome prediction.
train_classifier <- function(aln_int, target_mask, node_metadata, site_node_scores,
                             selected_panel = NULL, use_selected_panel = TRUE,
                             panel_size = NULL, min_sites_per_node = 1L,
                             min_total_informative_sites = 1L,
                             min_support = 0.8, max_conflict = 0.2,
                             support_margin = 0.05) {
  if (!is.matrix(aln_int)) {
    stop("aln_int must be an integer matrix.", call. = FALSE)
  }
  if (is.null(target_mask) || !is.matrix(target_mask)) {
    stop("target_mask must be a logical node-by-tip matrix.", call. = FALSE)
  }
  if (is.null(rownames(target_mask))) {
    rownames(target_mask) <- as.character(node_metadata$node_id)
  }

  rules <- make_classifier_rules(site_node_scores, selected_panel, use_selected_panel, panel_size)
  if (nrow(rules) > 0L) {
    rules <- merge(
      rules,
      node_metadata[, intersect(c("node_id", "node_index", "clade_size", "complement_size", "depth", "weight"), names(node_metadata)), drop = FALSE],
      by = "node_id",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("", "_node")
    )
    site_counts <- aggregate(site ~ node_id, rules, function(x) length(unique(x)))
    keep_nodes <- site_counts$node_id[site_counts$site >= as.integer(min_sites_per_node)]
    rules <- rules[rules$node_id %in% keep_nodes, , drop = FALSE]
  }

  node_table <- make_classifier_node_table(target_mask, node_metadata, rules)
  structure(
    list(
      version = 1L,
      rules = rules,
      node_table = node_table,
      target_mask = target_mask[match(node_table$node_id, rownames(target_mask)), , drop = FALSE],
      informative_sites = sort(unique(as.integer(rules$site))),
      tip_names = rownames(aln_int),
      alignment_length = ncol(aln_int),
      settings = list(
        use_selected_panel = isTRUE(use_selected_panel),
        panel_size = panel_size,
        min_sites_per_node = as.integer(min_sites_per_node),
        min_total_informative_sites = as.integer(min_total_informative_sites),
        min_support = as.numeric(min_support),
        max_conflict = as.numeric(max_conflict),
        support_margin = as.numeric(support_margin)
      )
    ),
    class = "wssv_classifier"
  )
}

#' Convert scored SNPs into classifier rules.
#'
#' Restricts rules to selected panel windows when requested and records allele,
#' direction, coordinate, and rule weight for downstream evidence calculations.
make_classifier_rules <- function(site_node_scores, selected_panel = NULL,
                                  use_selected_panel = TRUE, panel_size = NULL) {
  required <- c("node_id", "site", "best_allele", "direction")
  missing <- setdiff(required, names(site_node_scores))
  if (length(missing) > 0L) {
    stop("site_node_scores is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  rules <- site_node_scores
  if (isTRUE(use_selected_panel) && !is.null(selected_panel) && nrow(selected_panel) > 0L) {
    panel <- selected_panel
    if (!is.null(panel_size)) {
      panel <- panel[seq_len(min(as.integer(panel_size), nrow(panel))), , drop = FALSE]
    }
    sites <- sort(unique(unlist(Map(seq.int, as.integer(panel$start), as.integer(panel$end)), use.names = FALSE)))
    rules <- rules[rules$site %in% sites, , drop = FALSE]
  }
  rules$node_id <- as.character(rules$node_id)
  rules$site <- as.integer(rules$site)
  rules$allele_code <- unname(classifier_allele_codes[as.character(rules$best_allele)])
  rules$rule_weight <- if ("normalized_gain" %in% names(rules)) rules$normalized_gain else 1
  rules$rule_weight[is.na(rules$rule_weight) | rules$rule_weight <= 0] <- 1
  rules[!is.na(rules$allele_code) & rules$direction %in% c("clade", "outside"), , drop = FALSE]
}

make_classifier_node_table <- function(target_mask, node_metadata, rules) {
  node_ids <- unique(as.character(rules$node_id))
  if (length(node_ids) == 0L) {
    node_ids <- character()
  }
  metadata <- node_metadata[as.character(node_metadata$node_id) %in% node_ids, , drop = FALSE]
  metadata$node_id <- as.character(metadata$node_id)
  metadata$parent_node_id <- find_classifier_parents(target_mask, metadata$node_id)
  metadata$rule_sites <- vapply(metadata$node_id, function(id) length(unique(rules$site[rules$node_id == id])), integer(1L))
  metadata[order(metadata$depth, metadata$clade_size, metadata$node_id), , drop = FALSE]
}

find_classifier_parents <- function(target_mask, node_ids) {
  mask <- target_mask[match(node_ids, rownames(target_mask)), , drop = FALSE]
  sizes <- rowSums(mask)
  vapply(seq_along(node_ids), function(i) {
    contains_node <- apply(mask, 1L, function(row) all(mask[i, ] <= row))
    ancestors <- which(sizes > sizes[[i]] & contains_node)
    if (length(ancestors) == 0L) {
      return(NA_character_)
    }
    node_ids[ancestors[which.min(sizes[ancestors])]]
  }, character(1L))
}

extract_observed_informative_sites <- function(query, classifier) {
  x <- encode_classifier_query(query)
  sites <- classifier$informative_sites
  sites <- sites[sites >= 1L & sites <= length(x)]
  observed <- sites[x[sites] != 0L]
  data.frame(site = observed, allele_code = as.integer(x[observed]))
}

encode_classifier_query <- function(query) {
  if (is.matrix(query)) {
    if (nrow(query) != 1L) {
      stop("query matrix must have one row.", call. = FALSE)
    }
    return(as.integer(query[1L, ]))
  }
  if (is.numeric(query) || is.integer(query)) {
    return(as.integer(query))
  }
  if (is.character(query) && length(query) == 1L && nchar(query) > 1L) {
    chars <- strsplit(query, "", fixed = TRUE)[[1L]]
  } else if (is.character(query)) {
    chars <- query
  } else {
    stop("query must be an encoded vector, one-row matrix, or sequence string.", call. = FALSE)
  }
  codes <- classifier_allele_codes[toupper(chars)]
  codes[is.na(codes)] <- 0L
  as.integer(codes)
}

#' Calculate support and conflict evidence for every classifier node.
#'
#' Evidence compares observed query alleles with each node's rules and returns
#' weighted support/conflict fractions used by classify_tree_path().
calculate_node_evidence <- function(query, classifier) {
  x <- encode_classifier_query(query)
  rules <- classifier$rules
  if (nrow(rules) == 0L) {
    return(empty_node_evidence(classifier$node_table$node_id))
  }
  rules <- rules[rules$site <= length(x), , drop = FALSE]
  observed <- x[rules$site]
  rules <- rules[observed != 0L, , drop = FALSE]
  observed <- observed[observed != 0L]
  if (nrow(rules) == 0L) {
    return(empty_node_evidence(classifier$node_table$node_id))
  }
  supports <- ifelse(rules$direction == "clade", observed == rules$allele_code, observed != rules$allele_code)
  conflicts <- !supports
  rows <- split(data.frame(rules, supports = supports, conflicts = conflicts), rules$node_id)
  evidence <- do.call(rbind, lapply(rows, function(r) {
    support_weight <- sum(r$rule_weight[r$supports])
    conflict_weight <- sum(r$rule_weight[r$conflicts])
    total_weight <- support_weight + conflict_weight
    data.frame(
      node_id = r$node_id[[1L]],
      observed_sites = length(unique(r$site)),
      support_sites = sum(r$supports),
      conflict_sites = sum(r$conflicts),
      support_weight = support_weight,
      conflict_weight = conflict_weight,
      support = if (total_weight > 0) support_weight / total_weight else NA_real_,
      conflict = if (total_weight > 0) conflict_weight / total_weight else NA_real_
    )
  }))
  merge(evidence, classifier$node_table, by = "node_id", all.y = TRUE, sort = FALSE)
}

#' Classify one complete or aligned query sequence along the tree path.
#'
#' Returns resolved, weak, conflicting, no-informative-sites, or unresolved
#' status without assuming that a resolved training-set call is independent
#' biological validation.
classify_tree_path <- function(query, classifier, min_total_informative_sites = NULL,
                               min_support = NULL, max_conflict = NULL,
                               support_margin = NULL) {
  settings <- classifier$settings
  min_total_informative_sites <- setting_or_default(min_total_informative_sites, settings$min_total_informative_sites)
  min_support <- setting_or_default(min_support, settings$min_support)
  max_conflict <- setting_or_default(max_conflict, settings$max_conflict)
  support_margin <- setting_or_default(support_margin, settings$support_margin)

  observed <- extract_observed_informative_sites(query, classifier)
  evidence <- calculate_node_evidence(query, classifier)
  evaluated <- evidence[!is.na(evidence$observed_sites) & evidence$observed_sites > 0L, , drop = FALSE]
  status <- "unresolved"
  assigned <- NA_character_
  conflict_reason <- NA_character_
  strong <- evaluated[0, , drop = FALSE]

  if (nrow(observed) == 0L || length(unique(observed$site)) < as.integer(min_total_informative_sites)) {
    status <- "no_informative_sites"
    conflict_reason <- "no_informative_evidence"
  } else {
    strong <- supported_evidence(evaluated, min_support, max_conflict)
    if (nrow(strong) == 0L) {
      status <- "weak_support"
      conflict_reason <- "weak_evidence"
    } else if (!supported_nodes_form_nested_path(strong, classifier$target_mask)) {
      status <- "conflicting"
      conflict_reason <- "incompatible_strong_off_path_support"
    } else {
      strong <- strong[order(-strong$depth, -strong$support, strong$conflict), , drop = FALSE]
      assigned <- strong$node_id[[1L]]
      assigned_support <- strong$support[[1L]]
      competitors <- evaluated[evaluated$node_id != assigned, , drop = FALSE]
      close_competitor <- nrow(competitors) > 0L &&
        any(!is.na(competitors$support) &
          competitors$support >= assigned_support - support_margin &
          competitors$observed_sites > 0L &
          !nodes_nested(assigned, competitors$node_id, classifier$target_mask))
      if (close_competitor) {
        status <- "conflicting"
        conflict_reason <- "close_incompatible_competitor"
      } else {
        status <- "resolved"
        conflict_reason <- if (nrow(strong) > 1L) "nested_true_path_compatible_support" else "single_supported_node"
      }
    }
  }

  if (!is.na(assigned) && length(assigned) == 1L) {
    assigned_row <- evidence[evidence$node_id == assigned, , drop = FALSE]
  } else {
    assigned_row <- evidence[0, , drop = FALSE]
  }
  if (nrow(evaluated) > 0L && "depth" %in% names(evaluated)) {
    best_overall <- evaluated[order(-evaluated$support, evaluated$conflict, -evaluated$depth), , drop = FALSE]
  } else {
    best_overall <- evaluated[0, , drop = FALSE]
  }
  competing <- best_overall[if (!is.na(assigned)) best_overall$node_id != assigned else rep(TRUE, nrow(best_overall)), , drop = FALSE]
  list(
    assigned_node = assigned,
    assigned_depth = if (nrow(assigned_row) > 0L) assigned_row$depth[[1L]] else NA_real_,
    observed_informative_sites = length(unique(observed$site)),
    informative_nodes_evaluated = nrow(evaluated),
    support_score = if (nrow(assigned_row) > 0L) assigned_row$support[[1L]] else NA_real_,
    competing_node = if (nrow(competing) > 0L) competing$node_id[[1L]] else NA_character_,
    competing_support = if (nrow(competing) > 0L) competing$support[[1L]] else NA_real_,
    strong_supported_nodes = nrow(strong),
    conflict_reason = conflict_reason,
    status = status,
    evidence = evidence
  )
}

classifier_prediction_row <- function(prediction, query_id = NA_character_) {
  data.frame(
    query_id = query_id,
    assigned_node = prediction$assigned_node,
    assigned_depth = prediction$assigned_depth,
    observed_informative_sites = prediction$observed_informative_sites,
    informative_nodes_evaluated = prediction$informative_nodes_evaluated,
    support_score = prediction$support_score,
    competing_node = prediction$competing_node,
    competing_support = prediction$competing_support,
    strong_supported_nodes = prediction$strong_supported_nodes,
    conflict_reason = prediction$conflict_reason,
    status = prediction$status
  )
}

supported_evidence <- function(evidence, min_support, max_conflict) {
  evidence[!is.na(evidence$support) &
    evidence$support >= min_support &
    !is.na(evidence$conflict) &
    evidence$conflict <= max_conflict, , drop = FALSE]
}

nodes_nested <- function(node_id, other_node_ids, target_mask) {
  node_id <- as.character(node_id)
  other_node_ids <- as.character(other_node_ids)
  mask <- target_mask
  rownames(mask) <- as.character(rownames(mask))
  if (!node_id %in% rownames(mask)) {
    return(rep(FALSE, length(other_node_ids)))
  }
  node_mask <- mask[node_id, , drop = TRUE]
  vapply(other_node_ids, function(other_id) {
    if (!other_id %in% rownames(mask)) {
      return(FALSE)
    }
    other_mask <- mask[other_id, , drop = TRUE]
    all(node_mask <= other_mask) || all(other_mask <= node_mask)
  }, logical(1L))
}

supported_nodes_form_nested_path <- function(evidence, target_mask) {
  ids <- as.character(evidence$node_id)
  if (length(ids) <= 1L) {
    return(TRUE)
  }
  pairs <- utils::combn(ids, 2L)
  all(vapply(seq_len(ncol(pairs)), function(i) {
    nodes_nested(pairs[1L, i], pairs[2L, i], target_mask)
  }, logical(1L)))
}

save_classifier <- function(classifier, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(classifier, path)
  invisible(path)
}

load_classifier <- function(path) {
  readRDS(path)
}

direct_child_of <- function(node_id, parent_node_id, current) {
  if (is.na(current)) {
    return(is.na(parent_node_id))
  }
  !is.na(parent_node_id) & parent_node_id == current
}

setting_or_default <- function(value, default) {
  if (is.null(value)) default else value
}

empty_node_evidence <- function(node_ids) {
  data.frame(
    node_id = as.character(node_ids),
    observed_sites = 0L,
    support_sites = 0L,
    conflict_sites = 0L,
    support_weight = 0,
    conflict_weight = 0,
    support = NA_real_,
    conflict = NA_real_
  )
}
