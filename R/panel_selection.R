# Panel selection helpers for stage 05.
#
# These helpers consume stage-04 candidate amplicon tables and compare small
# candidate panels with explicit target evidence vectors. The implementation is
# deterministic, keeps absolute and target-relative gain separate, and avoids
# the legacy whole-tree weighted-window objective.

resolve_panel_selection_config <- function(config) {
  panel_cfg <- value_or(value_or(config$analysis, list())$panel_selection, list())
  amp_cfg <- value_or(value_or(config$analysis, list())$amplicons, list())

  primary_target_nodes <- unique(as.integer(value_or(
    panel_cfg$primary_target_nodes,
    value_or(amp_cfg$primary_target_nodes, c(45L, 66L, 72L, 80L))
  )))
  primary_target_nodes <- primary_target_nodes[!is.na(primary_target_nodes)]
  if (length(primary_target_nodes) == 0L) {
    stop("analysis.amplicons.primary_target_nodes must list at least one target.", call. = FALSE)
  }

  min_panel_size <- max(1L, as.integer(value_or(panel_cfg$min_panel_size, 2L)))
  max_panel_size <- max(min_panel_size, as.integer(value_or(panel_cfg$max_panel_size, 6L)))

  list(
    primary_target_nodes = sort(primary_target_nodes),
    min_panel_size = min_panel_size,
    max_panel_size = max_panel_size,
    top_n_panels = max(1L, as.integer(value_or(panel_cfg$top_n_panels, 20L))),
    exact_enumeration_limit = max(1L, as.integer(value_or(panel_cfg$exact_enumeration_limit, 100000L))),
    approximate_search_budget = max(1L, as.integer(value_or(panel_cfg$approximate_search_budget, 5000L))),
    max_reduced_candidates = max(max_panel_size, as.integer(value_or(panel_cfg$max_reduced_candidates, 30L))),
    evidence_threshold = as.numeric(value_or(panel_cfg$evidence_threshold, value_or(amp_cfg$minimum_gain_threshold, 0.5))),
    strong_evidence_threshold = as.numeric(value_or(panel_cfg$strong_evidence_threshold, value_or(amp_cfg$strong_gain_threshold, 0.8))),
    material_gain_improvement = as.numeric(value_or(panel_cfg$material_gain_improvement, 0.01)),
    material_relative_gain_improvement = as.numeric(value_or(panel_cfg$material_relative_gain_improvement, 0.01)),
    redundancy_cap_per_target = max(1L, as.integer(value_or(panel_cfg$redundancy_cap_per_target, 2L))),
    max_candidate_length = value_or(panel_cfg$max_candidate_length, value_or(amp_cfg$max_length, NULL)),
    min_complete_fraction = panel_cfg$min_complete_fraction,
    max_missing_fraction = panel_cfg$max_missing_fraction,
    require_left_flank_conserved = isTRUE(panel_cfg$require_left_flank_conserved),
    require_right_flank_conserved = isTRUE(panel_cfg$require_right_flank_conserved),
    require_both_flanks_conserved = isTRUE(panel_cfg$require_both_flanks_conserved),
    min_target_signal = panel_cfg$min_target_signal,
    min_targets_with_signal = panel_cfg$min_targets_with_signal,
    min_targets_with_strong_signal = panel_cfg$min_targets_with_strong_signal,
    max_ambiguity_fraction = panel_cfg$max_ambiguity_fraction,
    max_gap_fraction = panel_cfg$max_gap_fraction,
    max_allowed_overlap = as.numeric(value_or(panel_cfg$max_allowed_overlap, 0.5)),
    min_distance = panel_cfg$min_distance,
    dominance_overlap_threshold = as.numeric(value_or(panel_cfg$dominance_overlap_threshold, 0.85)),
    dominance_min_distance = panel_cfg$dominance_min_distance,
    local_redundancy_overlap_threshold = as.numeric(value_or(panel_cfg$local_redundancy_overlap_threshold, 0.9)),
    local_redundancy_min_distance = as.integer(value_or(panel_cfg$local_redundancy_min_distance, 100L)),
    target_gain_round_digits = as.integer(value_or(panel_cfg$target_gain_round_digits, 3L))
  )
}

safe_max <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_min <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

objective_value <- function(x, direction) {
  if (is.na(x)) {
    return(if (identical(direction, "min")) Inf else -Inf)
  }
  as.numeric(x)
}

dominates_objectives <- function(a, b, objectives, directions) {
  lhs <- vapply(seq_along(objectives), function(i) {
    objective_value(a[[objectives[[i]]]][[1L]], directions[[i]])
  }, numeric(1L))
  rhs <- vapply(seq_along(objectives), function(i) {
    objective_value(b[[objectives[[i]]]][[1L]], directions[[i]])
  }, numeric(1L))

  transformed_lhs <- lhs
  transformed_rhs <- rhs
  min_indices <- which(directions == "min")
  if (length(min_indices) > 0L) {
    transformed_lhs[min_indices] <- -transformed_lhs[min_indices]
    transformed_rhs[min_indices] <- -transformed_rhs[min_indices]
  }

  no_worse <- all(transformed_lhs >= transformed_rhs)
  strictly_better <- any(transformed_lhs > transformed_rhs)
  no_worse && strictly_better
}

panel_overlap_fraction <- function(a_start, a_end, b_start, b_end) {
  overlap <- max(0L, min(a_end, b_end) - max(a_start, b_start) + 1L)
  overlap / max(1L, min(a_end - a_start + 1L, b_end - b_start + 1L))
}

panel_distance <- function(a_start, a_end, b_start, b_end) {
  if (a_end < b_start) {
    b_start - a_end - 1L
  } else if (b_end < a_start) {
    a_start - b_end - 1L
  } else {
    0L
  }
}

candidate_snp_positions <- function(candidate_amplicon_snps) {
  if (is.null(candidate_amplicon_snps) || nrow(candidate_amplicon_snps) == 0L ||
      !all(c("candidate_id", "site") %in% names(candidate_amplicon_snps))) {
    return(data.frame(candidate_id = character(), candidate_snp_positions = I(list())))
  }

  rows <- lapply(split(candidate_amplicon_snps, candidate_amplicon_snps$candidate_id), function(x) {
    data.frame(
      candidate_id = x$candidate_id[[1L]],
      candidate_snp_positions = I(list(sort(unique(as.integer(x$site))))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

candidate_target_wide <- function(candidate_amplicon_target_scores, primary_target_nodes) {
  primary_target_nodes <- sort(unique(as.integer(primary_target_nodes)))
  candidate_ids <- sort(unique(as.character(candidate_amplicon_target_scores$candidate_id)))
  if (length(candidate_ids) == 0L) {
    out <- data.frame(candidate_id = character())
    for (target in primary_target_nodes) {
      out[[paste0("gain_", target)]] <- numeric()
      out[[paste0("target_supporting_snps_", target)]] <- integer()
      out[[paste0("best_snp_position_", target)]] <- integer()
    }
    return(out)
  }

  rows <- lapply(candidate_ids, function(candidate_id) {
    candidate_rows <- candidate_amplicon_target_scores[
      candidate_amplicon_target_scores$candidate_id == candidate_id &
        candidate_amplicon_target_scores$target_node %in% primary_target_nodes,
      , drop = FALSE
    ]
    row <- data.frame(candidate_id = candidate_id, stringsAsFactors = FALSE)
    for (target in primary_target_nodes) {
      row[[paste0("gain_", target)]] <- 0
      row[[paste0("target_supporting_snps_", target)]] <- 0L
      row[[paste0("best_snp_position_", target)]] <- NA_integer_
    }
    if (nrow(candidate_rows) > 0L) {
      for (i in seq_len(nrow(candidate_rows))) {
        target <- candidate_rows$target_node[[i]]
        row[[paste0("gain_", target)]] <- as.numeric(candidate_rows$max_normalized_gain[[i]])
        if ("n_informative_snps" %in% names(candidate_rows)) {
          row[[paste0("target_supporting_snps_", target)]] <- as.integer(candidate_rows$n_informative_snps[[i]])
        }
        if ("best_snp_position" %in% names(candidate_rows)) {
          row[[paste0("best_snp_position_", target)]] <- as.integer(candidate_rows$best_snp_position[[i]])
        }
      }
    }
    row
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

build_candidate_universe <- function(candidate_amplicons, candidate_amplicon_target_scores,
                                     candidate_amplicon_snps = NULL, config) {
  cfg <- resolve_panel_selection_config(config)
  primary_targets <- cfg$primary_target_nodes
  missing_amplicon_cols <- setdiff(c("candidate_id", "start", "end", "length"), names(candidate_amplicons))
  if (length(missing_amplicon_cols) > 0L) {
    stop("candidate_amplicons is missing required columns: ", paste(missing_amplicon_cols, collapse = ", "), call. = FALSE)
  }
  missing_target_cols <- setdiff(c("candidate_id", "target_node", "max_normalized_gain"), names(candidate_amplicon_target_scores))
  if (length(missing_target_cols) > 0L) {
    stop("candidate_amplicon_target_scores is missing required columns: ", paste(missing_target_cols, collapse = ", "), call. = FALSE)
  }

  table <- candidate_amplicons
  if (!"window_id" %in% names(table)) {
    table$window_id <- table$candidate_id
  }

  target_wide <- candidate_target_wide(candidate_amplicon_target_scores, primary_targets)
  table <- merge(table, target_wide, by = "candidate_id", all.x = TRUE, sort = FALSE)
  table <- table[match(candidate_amplicons$candidate_id, table$candidate_id), , drop = FALSE]
  row.names(table) <- NULL

  for (target in primary_targets) {
    gain_col <- paste0("gain_", target)
    if (!gain_col %in% names(table)) {
      table[[gain_col]] <- 0
    }
    table[[gain_col]][is.na(table[[gain_col]])] <- 0
  }

  for (target in primary_targets) {
    gain_col <- paste0("gain_", target)
    rel_col <- paste0("relative_gain_", target)
    max_gain <- safe_max(table[[gain_col]])
    table[[rel_col]] <- if (is.na(max_gain) || max_gain <= 0) 0 else table[[gain_col]] / max_gain
  }

  gain_cols <- paste0("gain_", primary_targets)
  rel_cols <- paste0("relative_gain_", primary_targets)
  table$max_target_gain <- vapply(seq_len(nrow(table)), function(i) safe_max(table[i, gain_cols]), numeric(1L))
  table$min_target_gain <- vapply(seq_len(nrow(table)), function(i) safe_min(table[i, gain_cols]), numeric(1L))
  table$max_target_relative_gain <- vapply(seq_len(nrow(table)), function(i) safe_max(table[i, rel_cols]), numeric(1L))
  table$min_target_relative_gain <- vapply(seq_len(nrow(table)), function(i) safe_min(table[i, rel_cols]), numeric(1L))
  table$targets_with_any_signal <- rowSums(table[gain_cols] > 0, na.rm = TRUE)
  table$targets_with_strong_signal <- rowSums(table[gain_cols] >= cfg$strong_evidence_threshold, na.rm = TRUE)
  table$target_gain_vector <- vapply(seq_len(nrow(table)), function(i) {
    paste(formatC(as.numeric(table[i, gain_cols]), format = "f", digits = cfg$target_gain_round_digits), collapse = "|")
  }, character(1L))
  table$target_relative_vector <- vapply(seq_len(nrow(table)), function(i) {
    paste(formatC(as.numeric(table[i, rel_cols]), format = "f", digits = cfg$target_gain_round_digits), collapse = "|")
  }, character(1L))

  snp_positions <- candidate_snp_positions(candidate_amplicon_snps)
  if (nrow(snp_positions) > 0L) {
    table <- merge(table, snp_positions, by = "candidate_id", all.x = TRUE, sort = FALSE)
    table <- table[match(candidate_amplicons$candidate_id, table$candidate_id), , drop = FALSE]
    row.names(table) <- NULL
  } else {
    table$candidate_snp_positions <- I(vector("list", nrow(table)))
  }
  if (nrow(table) > 0L) {
    missing_rows <- which(vapply(table$candidate_snp_positions, function(x) length(x) == 0L || all(is.na(x)), logical(1L)))
    for (i in missing_rows) {
      table$candidate_snp_positions[[i]] <- integer()
    }
  }
  table$candidate_snp_positions_key <- vapply(table$candidate_snp_positions, function(x) paste(unique(as.integer(x)), collapse = ","), character(1L))
  table
}

filter_candidate_universe <- function(candidate_table, cfg) {
  out <- candidate_table
  if (nrow(out) == 0L) {
    return(out)
  }

  if (!is.null(cfg$max_candidate_length)) {
    out <- out[!is.na(out$length) & out$length <= as.integer(cfg$max_candidate_length), , drop = FALSE]
  }
  if (!is.null(cfg$min_complete_fraction) && "complete_fraction" %in% names(out)) {
    out <- out[is.na(out$complete_fraction) | out$complete_fraction >= as.numeric(cfg$min_complete_fraction), , drop = FALSE]
  }
  if (!is.null(cfg$max_missing_fraction) && "missing_fraction" %in% names(out)) {
    out <- out[is.na(out$missing_fraction) | out$missing_fraction <= as.numeric(cfg$max_missing_fraction), , drop = FALSE]
  }
  if (cfg$require_left_flank_conserved && "left_flank_conserved_stretch" %in% names(out)) {
    out <- out[!is.na(out$left_flank_conserved_stretch) & out$left_flank_conserved_stretch, , drop = FALSE]
  }
  if (cfg$require_right_flank_conserved && "right_flank_conserved_stretch" %in% names(out)) {
    out <- out[!is.na(out$right_flank_conserved_stretch) & out$right_flank_conserved_stretch, , drop = FALSE]
  }
  if (cfg$require_both_flanks_conserved && all(c("left_flank_conserved_stretch", "right_flank_conserved_stretch") %in% names(out))) {
    out <- out[
      !is.na(out$left_flank_conserved_stretch) &
        !is.na(out$right_flank_conserved_stretch) &
        out$left_flank_conserved_stretch &
        out$right_flank_conserved_stretch,
      , drop = FALSE
    ]
  }
  if (!is.null(cfg$min_target_signal)) {
    out <- out[out$max_target_gain >= as.numeric(cfg$min_target_signal), , drop = FALSE]
  }
  if (!is.null(cfg$min_targets_with_signal)) {
    out <- out[out$targets_with_any_signal >= as.integer(cfg$min_targets_with_signal), , drop = FALSE]
  }
  if (!is.null(cfg$min_targets_with_strong_signal)) {
    out <- out[out$targets_with_strong_signal >= as.integer(cfg$min_targets_with_strong_signal), , drop = FALSE]
  }
  if (!is.null(cfg$max_ambiguity_fraction) && "ambiguity_fraction" %in% names(out)) {
    out <- out[is.na(out$ambiguity_fraction) | out$ambiguity_fraction <= as.numeric(cfg$max_ambiguity_fraction), , drop = FALSE]
  }
  if (!is.null(cfg$max_gap_fraction) && "gap_fraction" %in% names(out)) {
    out <- out[is.na(out$gap_fraction) | out$gap_fraction <= as.numeric(cfg$max_gap_fraction), , drop = FALSE]
  }

  row.names(out) <- NULL
  out
}

candidate_redundant <- function(a, b, cfg) {
  overlap_fraction <- panel_overlap_fraction(a$start[[1L]], a$end[[1L]], b$start[[1L]], b$end[[1L]])
  distance <- panel_distance(a$start[[1L]], a$end[[1L]], b$start[[1L]], b$end[[1L]])
  same_evidence <- identical(a$target_gain_vector[[1L]], b$target_gain_vector[[1L]]) &&
    identical(a$target_relative_vector[[1L]], b$target_relative_vector[[1L]])
  same_snps <- identical(a$candidate_snp_positions_key[[1L]], b$candidate_snp_positions_key[[1L]])
  same_evidence && same_snps && (overlap_fraction >= cfg$local_redundancy_overlap_threshold || distance <= cfg$local_redundancy_min_distance)
}

candidate_dominates <- function(a, b, cfg) {
  if (!candidate_redundant(a, b, cfg)) {
    return(FALSE)
  }
  gain_cols <- grep("^gain_[0-9]+$", names(a), value = TRUE)
  better_or_equal_targets <- all(as.numeric(a[gain_cols]) >= as.numeric(b[gain_cols]))
  better_or_equal_complete <- !("complete_fraction" %in% names(a)) ||
    (safe_max(a$complete_fraction) >= safe_max(b$complete_fraction))
  better_or_equal_flanks <- TRUE
  if (all(c("left_flank_conserved_stretch", "right_flank_conserved_stretch") %in% names(a))) {
    better_or_equal_flanks <- (!isTRUE(b$left_flank_conserved_stretch[[1L]]) || isTRUE(a$left_flank_conserved_stretch[[1L]])) &&
      (!isTRUE(b$right_flank_conserved_stretch[[1L]]) || isTRUE(a$right_flank_conserved_stretch[[1L]]))
  }
  not_longer <- as.integer(a$length[[1L]]) <= as.integer(b$length[[1L]])
  better_or_equal_targets && better_or_equal_complete && better_or_equal_flanks && not_longer
}

reduce_candidate_universe <- function(candidate_table, cfg) {
  if (nrow(candidate_table) <= 1L) {
    return(candidate_table)
  }

  ordered <- candidate_table[order(
    -candidate_table$max_target_gain,
    -candidate_table$max_target_relative_gain,
    -candidate_table$targets_with_strong_signal,
    -candidate_table$targets_with_any_signal,
    -candidate_table$complete_fraction,
    candidate_table$length,
    candidate_table$start,
    candidate_table$end,
    candidate_table$candidate_id
  ), , drop = FALSE]

  keep <- rep(TRUE, nrow(ordered))
  for (i in seq_len(nrow(ordered) - 1L)) {
    if (!keep[[i]]) {
      next
    }
    for (j in seq.int(i + 1L, nrow(ordered))) {
      if (!keep[[j]]) {
        next
      }
      if (candidate_dominates(ordered[i, , drop = FALSE], ordered[j, , drop = FALSE], cfg)) {
        keep[[j]] <- FALSE
      } else if (candidate_dominates(ordered[j, , drop = FALSE], ordered[i, , drop = FALSE], cfg)) {
        keep[[i]] <- FALSE
        break
      }
    }
  }

  reduced <- ordered[keep, , drop = FALSE]
  if (nrow(reduced) > cfg$max_reduced_candidates) {
    reduced <- reduced[seq_len(cfg$max_reduced_candidates), , drop = FALSE]
  }
  row.names(reduced) <- NULL
  reduced
}

panel_compatible <- function(panel, cfg) {
  if (nrow(panel) <= 1L) {
    return(TRUE)
  }
  for (i in seq_len(nrow(panel) - 1L)) {
    for (j in seq.int(i + 1L, nrow(panel))) {
      overlap <- panel_overlap_fraction(panel$start[[i]], panel$end[[i]], panel$start[[j]], panel$end[[j]])
      distance <- panel_distance(panel$start[[i]], panel$end[[i]], panel$start[[j]], panel$end[[j]])
      if (overlap > cfg$max_allowed_overlap) {
        return(FALSE)
      }
      if (!is.null(cfg$min_distance) && distance < as.integer(cfg$min_distance)) {
        return(FALSE)
      }
    }
  }
  TRUE
}

panel_compatible_with_candidate <- function(panel, candidate, cfg) {
  if (nrow(panel) == 0L) {
    return(TRUE)
  }
  for (i in seq_len(nrow(panel))) {
    overlap <- panel_overlap_fraction(panel$start[[i]], panel$end[[i]], candidate$start[[1L]], candidate$end[[1L]])
    distance <- panel_distance(panel$start[[i]], panel$end[[i]], candidate$start[[1L]], candidate$end[[1L]])
    if (overlap > cfg$max_allowed_overlap) {
      return(FALSE)
    }
    if (!is.null(cfg$min_distance) && distance < as.integer(cfg$min_distance)) {
      return(FALSE)
    }
  }
  TRUE
}

panel_metrics <- function(panel, cfg) {
  targets <- cfg$primary_target_nodes
  gain_cols <- paste0("gain_", targets)
  rel_cols <- paste0("relative_gain_", targets)

  panel <- panel[order(panel$start, panel$end, panel$candidate_id), , drop = FALSE]
  abs_vec <- vapply(gain_cols, function(col) safe_max(panel[[col]]), numeric(1L))
  rel_vec <- vapply(rel_cols, function(col) safe_max(panel[[col]]), numeric(1L))
  supporting <- vapply(targets, function(target) sum(panel[[paste0("gain_", target)]] > 0, na.rm = TRUE), integer(1L))
  strong_supporting <- vapply(targets, function(target) sum(panel[[paste0("gain_", target)]] >= cfg$evidence_threshold, na.rm = TRUE), integer(1L))
  capped_supporting <- pmin(supporting, as.integer(cfg$redundancy_cap_per_target))
  snp_positions <- if ("candidate_snp_positions" %in% names(panel)) unique(unlist(panel$candidate_snp_positions, use.names = FALSE)) else integer()

  data.frame(
    panel_size = nrow(panel),
    candidate_ids = paste(panel$candidate_id, collapse = ";"),
    candidate_id_key = paste(sort(panel$candidate_id), collapse = "|"),
    total_length = sum(panel$length, na.rm = TRUE),
    mean_length = mean(panel$length, na.rm = TRUE),
    min_regional_complete_fraction = if ("complete_fraction" %in% names(panel)) safe_min(panel$complete_fraction) else NA_real_,
    mean_regional_complete_fraction = if ("complete_fraction" %in% names(panel)) mean(panel$complete_fraction, na.rm = TRUE) else NA_real_,
    left_flank_acceptable = if ("left_flank_conserved_stretch" %in% names(panel)) sum(!is.na(panel$left_flank_conserved_stretch) & panel$left_flank_conserved_stretch) else NA_integer_,
    right_flank_acceptable = if ("right_flank_conserved_stretch" %in% names(panel)) sum(!is.na(panel$right_flank_conserved_stretch) & panel$right_flank_conserved_stretch) else NA_integer_,
    both_flanks_acceptable = if (all(c("left_flank_conserved_stretch", "right_flank_conserved_stretch") %in% names(panel))) {
      sum(!is.na(panel$left_flank_conserved_stretch) & !is.na(panel$right_flank_conserved_stretch) & panel$left_flank_conserved_stretch & panel$right_flank_conserved_stretch)
    } else {
      NA_integer_
    },
    unique_informative_snp_positions = length(unique(as.integer(snp_positions))),
    min_absolute_target_gain = safe_min(abs_vec),
    min_target_relative_gain = safe_min(rel_vec),
    number_of_targets_with_two_supporting_amplicons = sum(supporting >= 2L),
    number_of_targets_with_signal_threshold_support = sum(strong_supporting >= 1L),
    capped_target_redundancy = sum(capped_supporting, na.rm = TRUE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

panel_row_for_indices <- function(candidate_table, indices, cfg,
                                  panel_id = NA_character_,
                                  search_method = NA_character_) {
  if (length(indices) == 0L) {
    panel <- candidate_table[0, , drop = FALSE]
  } else {
    panel <- candidate_table[indices, , drop = FALSE]
  }
  
  metrics <- panel_metrics(panel, cfg)
  
  for (target in cfg$primary_target_nodes) {
    col_gain <- paste0("gain_", target)
    col_rel <- paste0("relative_gain_", target)
    
    metrics[[col_gain]] <- if (nrow(panel) == 0L) 0 else safe_max(panel[[col_gain]])
    metrics[[col_rel]] <- if (nrow(panel) == 0L) 0 else safe_max(panel[[col_rel]])
  }
  
  metrics$panel_id <- panel_id
  metrics$search_method <- search_method
  metrics$panel_key <- if (nrow(panel) == 0L) {
    ""
  } else {
    paste(sort(panel$candidate_id), collapse = "|")
  }
  
  metrics$indices_key <- paste(sort(indices), collapse = ",")
  
  metrics
}

panel_target_rows <- function(panel_id, panel, cfg) {
  rows <- lapply(cfg$primary_target_nodes, function(target) {
    gain_col <- paste0("gain_", target)
    rel_col <- paste0("relative_gain_", target)
    gains <- as.numeric(panel[[gain_col]])
    rels <- as.numeric(panel[[rel_col]])
    ordered <- panel[order(-gains, panel$candidate_id), , drop = FALSE]
    best <- ordered[1L, , drop = FALSE]
    second_gain <- {
      uniq <- sort(unique(gains[!is.na(gains)]), decreasing = TRUE)
      if (length(uniq) >= 2L) uniq[[2L]] else NA_real_
    }
    data.frame(
      panel_id = panel_id,
      target_node = target,
      best_amplicon = best$candidate_id[[1L]],
      best_normalized_gain = best[[gain_col]][[1L]],
      target_relative_gain = best[[rel_col]][[1L]],
      second_best_normalized_gain = second_gain,
      number_of_supporting_amplicons = sum(gains > 0, na.rm = TRUE),
      number_above_threshold = sum(gains >= cfg$evidence_threshold, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

panel_objectives <- function(within_size = TRUE) {
  if (within_size) {
    list(
      columns = c(
        "min_absolute_target_gain",
        "min_target_relative_gain",
        "capped_target_redundancy",
        "min_regional_complete_fraction",
        "both_flanks_acceptable",
        "total_length"
      ),
      directions = c("max", "max", "max", "max", "max", "min")
    )
  } else {
    list(
      columns = c(
        "min_absolute_target_gain",
        "min_target_relative_gain",
        "capped_target_redundancy",
        "min_regional_complete_fraction",
        "both_flanks_acceptable",
        "panel_size",
        "total_length"
      ),
      directions = c("max", "max", "max", "max", "max", "min", "min")
    )
  }
}

identify_pareto_panels <- function(panel_candidates, within_size = TRUE) {
  if (nrow(panel_candidates) == 0L) {
    return(logical())
  }
  objectives <- panel_objectives(within_size = within_size)
  dominated <- rep(FALSE, nrow(panel_candidates))
  for (i in seq_len(nrow(panel_candidates) - 1L)) {
    if (dominated[[i]]) {
      next
    }
    for (j in seq.int(i + 1L, nrow(panel_candidates))) {
      if (dominated[[j]]) {
        next
      }
      if (dominates_objectives(panel_candidates[i, , drop = FALSE], panel_candidates[j, , drop = FALSE], objectives$columns, objectives$directions)) {
        dominated[[j]] <- TRUE
      } else if (dominates_objectives(panel_candidates[j, , drop = FALSE], panel_candidates[i, , drop = FALSE], objectives$columns, objectives$directions)) {
        dominated[[i]] <- TRUE
        break
      }
    }
  }
  !dominated
}

sort_panels_within_size <- function(panel_candidates) {
  if (nrow(panel_candidates) == 0L) {
    return(panel_candidates)
  }
  panel_candidates <- panel_candidates[order(
    -panel_candidates$min_absolute_target_gain,
    -panel_candidates$min_target_relative_gain,
    -panel_candidates$capped_target_redundancy,
    -panel_candidates$min_regional_complete_fraction,
    -panel_candidates$both_flanks_acceptable,
    panel_candidates$total_length,
    panel_candidates$candidate_id_key
  ), , drop = FALSE]
  row.names(panel_candidates) <- NULL
  panel_candidates$deterministic_rank <- seq_len(nrow(panel_candidates))
  panel_candidates$panel_rank <- panel_candidates$deterministic_rank
  panel_candidates$panel_id <- sprintf("panel_k%d_%03d", panel_candidates$panel_size, panel_candidates$deterministic_rank)
  panel_candidates
}

sort_panels_lexicographically <- function(panel_candidates) {
  sort_panels_within_size(panel_candidates)
}

enumerate_panels_exact_size <- function(candidate_table, cfg, panel_size) {
  if (panel_size > nrow(candidate_table) || panel_size <= 0L) {
    return(data.frame())
  }

  combos <- combn(seq_len(nrow(candidate_table)), panel_size, simplify = FALSE)
  results <- lapply(combos, function(indices) {
    panel <- candidate_table[indices, , drop = FALSE]
    if (!panel_compatible(panel, cfg)) {
      return(NULL)
    }
    panel_row_for_indices(candidate_table, indices, cfg, search_method = "exact")
  })
  results <- results[!vapply(results, is.null, logical(1L))]
  if (length(results) == 0L) {
    return(data.frame())
  }
  do.call(rbind, results)
}

score_state <- function(indices, candidate_table, cfg, search_method) {
  panel_row_for_indices(
    candidate_table,
    indices,
    cfg,
    search_method = search_method
  )
}

enumerate_panels_beam_size <- function(candidate_table, cfg, panel_size) {
  if (panel_size > nrow(candidate_table) || panel_size <= 0L) {
    return(data.frame())
  }

  beam_width <- max(50L, as.integer(cfg$approximate_search_budget))
  complete_order <- if ("complete_fraction" %in% names(candidate_table)) -candidate_table$complete_fraction else 0
  flank_order <- if (all(c("left_flank_conserved_stretch", "right_flank_conserved_stretch") %in% names(candidate_table))) {
    -as.integer(candidate_table$left_flank_conserved_stretch %in% TRUE & candidate_table$right_flank_conserved_stretch %in% TRUE)
  } else {
    0
  }
  ordered_candidates <- candidate_table[order(
    -candidate_table$max_target_gain,
    -candidate_table$max_target_relative_gain,
    -candidate_table$targets_with_strong_signal,
    -candidate_table$targets_with_any_signal,
    complete_order,
    flank_order,
    candidate_table$length,
    candidate_table$start,
    candidate_table$end,
    candidate_table$candidate_id
  ), , drop = FALSE]

  states_by_size <- vector("list", panel_size + 1L)
  states_by_size[[1L]] <- list(list(indices = integer(), score = score_state(integer(), candidate_table, cfg, "approximate")))

  for (i in seq_len(nrow(ordered_candidates))) {
    candidate_index <- match(ordered_candidates$candidate_id[[i]], candidate_table$candidate_id)
    candidate <- candidate_table[candidate_index, , drop = FALSE]
    for (current_size in seq.int(min(i - 1L, panel_size - 1L), 0L, by = -1L)) {
      current_states <- states_by_size[[current_size + 1L]]
      if (length(current_states) == 0L) {
        next
      }
      expanded <- list()
      for (state in current_states) {
        indices <- state$indices
        if (length(indices) >= panel_size) {
          next
        }
        panel <- if (length(indices) == 0L) candidate_table[0, , drop = FALSE] else candidate_table[indices, , drop = FALSE]
        if (!panel_compatible_with_candidate(panel, candidate, cfg)) {
          next
        }
        new_indices <- sort(unique(c(indices, candidate_index)))
        if (length(new_indices) > panel_size) {
          next
        }
        expanded[[length(expanded) + 1L]] <- list(
          indices = new_indices,
          score = score_state(new_indices, candidate_table, cfg, "approximate")
        )
      }
      if (length(expanded) > 0L) {
        expanded <- c(states_by_size[[current_size + 2L]], expanded)
        score_keys <- vapply(expanded, function(x) x$score$panel_key[[1L]], character(1L))
        expanded <- expanded[!duplicated(score_keys)]
        score_rows <- do.call(rbind, lapply(expanded, function(x) x$score))
        order_cols <- order(
          -score_rows$min_absolute_target_gain,
          -score_rows$min_target_relative_gain,
          -score_rows$capped_target_redundancy,
          -score_rows$min_regional_complete_fraction,
          -score_rows$both_flanks_acceptable,
          score_rows$total_length,
          score_rows$panel_key
        )
        if (length(order_cols) > beam_width) {
          order_cols <- order_cols[seq_len(beam_width)]
        }
        states_by_size[[current_size + 2L]] <- expanded[order_cols]
      }
    }
  }

  full_states <- states_by_size[[panel_size + 1L]]
  if (length(full_states) == 0L) {
    return(data.frame())
  }
  full_scores <- lapply(full_states, function(x) x$score)
  full_scores <- full_scores[!vapply(full_scores, is.null, logical(1L))]
  if (length(full_scores) == 0L) {
    return(data.frame())
  }
  do.call(rbind, full_scores)
}

evaluate_panel_size <- function(candidate_table, cfg, panel_size) {
  if (nrow(candidate_table) == 0L || panel_size <= 0L || panel_size > nrow(candidate_table)) {
    return(list(
      panel_candidates = data.frame(),
      search_method = "none",
      panels_evaluated = 0L,
      possible_combinations = 0,
      panel_size = panel_size
    ))
  }

  possible_combinations <- choose(nrow(candidate_table), panel_size)
  if (is.finite(possible_combinations) && possible_combinations <= cfg$exact_enumeration_limit) {
    panels <- enumerate_panels_exact_size(candidate_table, cfg, panel_size)
    method <- "exact"
    evaluated <- as.integer(possible_combinations)
  } else {
    panels <- enumerate_panels_beam_size(candidate_table, cfg, panel_size)
    method <- "approximate"
    evaluated <- nrow(panels)
  }

  if (nrow(panels) == 0L) {
    return(list(
      panel_candidates = data.frame(),
      search_method = method,
      panels_evaluated = evaluated,
      possible_combinations = possible_combinations,
      panel_size = panel_size
    ))
  }

  panels <- panels[!duplicated(panels$panel_key), , drop = FALSE]
  panels <- sort_panels_within_size(panels)
  panels$within_size_pareto <- identify_pareto_panels(panels, within_size = TRUE)
  panels$pareto_status <- ifelse(panels$within_size_pareto, "pareto", "dominated")
  panels$search_method <- method

  list(
    panel_candidates = panels,
    search_method = method,
    panels_evaluated = evaluated,
    possible_combinations = possible_combinations,
    panel_size = panel_size
  )
}

evaluate_panels <- function(candidate_table, cfg) {
  evaluate_panel_sizes(candidate_table, cfg)
}

best_panel_by_size <- function(panel_candidates) {
  if (nrow(panel_candidates) == 0L) {
    return(data.frame())
  }
  rows <- lapply(split(panel_candidates, panel_candidates$panel_size), function(x) sort_panels_within_size(x)[1L, , drop = FALSE])
  out <- do.call(rbind, rows)
  out[order(out$panel_size), , drop = FALSE]
}

evaluate_panel_sizes <- function(candidate_table, cfg) {
  sizes <- seq.int(cfg$min_panel_size, cfg$max_panel_size)
  size_results <- lapply(sizes, function(size) evaluate_panel_size(candidate_table, cfg, size))
  names(size_results) <- as.character(sizes)

  candidate_rows <- lapply(size_results, function(result) result$panel_candidates)
  candidate_rows <- candidate_rows[vapply(candidate_rows, function(x) nrow(x) > 0L, logical(1L))]
  panel_candidates <- if (length(candidate_rows) == 0L) data.frame() else do.call(rbind, candidate_rows)
  if (nrow(panel_candidates) > 0L) {
    panel_candidates <- panel_candidates[order(panel_candidates$panel_size, panel_candidates$deterministic_rank, panel_candidates$candidate_id_key), , drop = FALSE]
    row.names(panel_candidates) <- NULL
    panel_candidates$cross_size_pareto <- identify_pareto_panels(panel_candidates, within_size = FALSE)
  } else {
    panel_candidates$cross_size_pareto <- logical()
  }

  size_summary_rows <- lapply(size_results, function(result) {
    size_panels <- result$panel_candidates
    best <- if (nrow(size_panels) > 0L) size_panels[1L, , drop = FALSE] else data.frame()
    data.frame(
      panel_size = as.integer(result$panel_size),
      possible_combinations = as.numeric(result$possible_combinations),
      panels_evaluated = as.integer(result$panels_evaluated),
      search_method = as.character(result$search_method),
      within_size_pareto_count = if (nrow(size_panels) > 0L) as.integer(sum(size_panels$within_size_pareto, na.rm = TRUE)) else 0L,
      number_of_pareto_panels = if (nrow(size_panels) > 0L) as.integer(sum(size_panels$within_size_pareto, na.rm = TRUE)) else 0L,
      best_panel_id = if (nrow(best) > 0L) best$panel_id[[1L]] else NA_character_,
      candidate_ids = if (nrow(best) > 0L) best$candidate_ids[[1L]] else NA_character_,
      candidate_id_key = if (nrow(best) > 0L) best$candidate_id_key[[1L]] else NA_character_,
      gain_45 = if (nrow(best) > 0L && "gain_45" %in% names(best)) best$gain_45[[1L]] else NA_real_,
      gain_66 = if (nrow(best) > 0L && "gain_66" %in% names(best)) best$gain_66[[1L]] else NA_real_,
      gain_72 = if (nrow(best) > 0L && "gain_72" %in% names(best)) best$gain_72[[1L]] else NA_real_,
      gain_80 = if (nrow(best) > 0L && "gain_80" %in% names(best)) best$gain_80[[1L]] else NA_real_,
      min_absolute_target_gain = if (nrow(best) > 0L) best$min_absolute_target_gain[[1L]] else NA_real_,
      relative_gain_45 = if (nrow(best) > 0L && "relative_gain_45" %in% names(best)) best$relative_gain_45[[1L]] else NA_real_,
      relative_gain_66 = if (nrow(best) > 0L && "relative_gain_66" %in% names(best)) best$relative_gain_66[[1L]] else NA_real_,
      relative_gain_72 = if (nrow(best) > 0L && "relative_gain_72" %in% names(best)) best$relative_gain_72[[1L]] else NA_real_,
      relative_gain_80 = if (nrow(best) > 0L && "relative_gain_80" %in% names(best)) best$relative_gain_80[[1L]] else NA_real_,
      min_target_relative_gain = if (nrow(best) > 0L) best$min_target_relative_gain[[1L]] else NA_real_,
      capped_target_redundancy = if (nrow(best) > 0L) best$capped_target_redundancy[[1L]] else NA_real_,
      number_of_targets_with_two_supporting_amplicons = if (nrow(best) > 0L) best$number_of_targets_with_two_supporting_amplicons[[1L]] else NA_real_,
      number_of_targets_with_signal_threshold_support = if (nrow(best) > 0L) best$number_of_targets_with_signal_threshold_support[[1L]] else NA_real_,
      total_length = if (nrow(best) > 0L) best$total_length[[1L]] else NA_real_,
      mean_length = if (nrow(best) > 0L) best$mean_length[[1L]] else NA_real_,
      min_regional_complete_fraction = if (nrow(best) > 0L) best$min_regional_complete_fraction[[1L]] else NA_real_,
      mean_regional_complete_fraction = if (nrow(best) > 0L) best$mean_regional_complete_fraction[[1L]] else NA_real_,
      both_flanks_acceptable = if (nrow(best) > 0L) best$both_flanks_acceptable[[1L]] else NA_real_,
      unique_informative_snp_positions = if (nrow(best) > 0L) best$unique_informative_snp_positions[[1L]] else NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  size_summary <- do.call(rbind, size_summary_rows)
  size_summary <- size_summary[order(size_summary$panel_size), , drop = FALSE]
  row.names(size_summary) <- NULL
  if (nrow(size_summary) > 0L) {
    size_summary$delta_min_absolute_target_gain <- c(NA_real_, diff(size_summary$min_absolute_target_gain))
    size_summary$delta_min_target_relative_gain <- c(NA_real_, diff(size_summary$min_target_relative_gain))
    size_summary$delta_capped_target_redundancy <- c(NA_real_, diff(size_summary$capped_target_redundancy))
    size_summary$delta_number_of_targets_with_two_supporting_amplicons <- c(NA_real_, diff(size_summary$number_of_targets_with_two_supporting_amplicons))
    size_summary$delta_total_length <- c(NA_real_, diff(size_summary$total_length))
    size_summary$delta_min_regional_complete_fraction <- c(NA_real_, diff(size_summary$min_regional_complete_fraction))
  }

  recommendation <- recommend_provisional_panel_size(size_summary, cfg)
  if (is.na(recommendation$panel_size)) {
    size_summary$recommended_provisional_size <- FALSE
  } else {
    size_summary$recommended_provisional_size <- size_summary$panel_size == recommendation$panel_size
  }
  size_summary$recommendation_reason <- NA_character_
  if (!is.na(recommendation$panel_size)) {
    idx <- match(recommendation$panel_size, size_summary$panel_size)
    size_summary$recommendation_reason[[idx]] <- recommendation$reason
  }

  list(
    panel_candidates = panel_candidates,
    panel_size_comparison = size_summary,
    recommended_panel_size = recommendation$panel_size,
    recommended_reason = recommendation$reason,
    recommended_panel_id = recommendation$panel_id,
    size_results = size_results
  )
}

panel_target_table <- function(panel_candidates, candidate_table, cfg) {
  if (nrow(panel_candidates) == 0L) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(panel_candidates)), function(i) {
    panel_row <- panel_candidates[i, , drop = FALSE]
    ids <- unlist(strsplit(panel_row$candidate_ids[[1L]], ";", fixed = TRUE))
    ids <- ids[nzchar(ids)]
    panel <- candidate_table[match(ids, candidate_table$candidate_id), , drop = FALSE]
    panel_target_rows(panel_row$panel_id[[1L]], panel, cfg)
  })
  do.call(rbind, rows)
}

panel_members_table <- function(panel_row, candidate_table) {
  if (nrow(panel_row) == 0L) {
    return(data.frame())
  }
  ids <- unlist(strsplit(panel_row$candidate_ids[[1L]], ";", fixed = TRUE))
  ids <- ids[nzchar(ids)]
  selected <- candidate_table[match(ids, candidate_table$candidate_id), , drop = FALSE]
  selected <- selected[!is.na(selected$candidate_id), , drop = FALSE]
  selected <- selected[order(selected$start, selected$end, selected$candidate_id), , drop = FALSE]
  row.names(selected) <- NULL
  selected$panel_id <- panel_row$panel_id[[1L]]
  selected$panel_rank <- panel_row$panel_rank[[1L]]
  selected$panel_candidate_order <- seq_len(nrow(selected))
  selected$window_id <- selected$candidate_id
  selected
}

selected_panel_table <- function(panel_candidates, candidate_table) {
  if (nrow(panel_candidates) == 0L) {
    return(data.frame())
  }
  panel_members_table(panel_candidates[1L, , drop = FALSE], candidate_table)
}

recommend_provisional_panel_size <- function(panel_size_comparison, cfg) {
  if (nrow(panel_size_comparison) == 0L) {
    return(list(panel_size = NA_integer_, panel_id = NA_character_, reason = "No feasible panels were evaluated."))
  }

  valid <- panel_size_comparison[!is.na(panel_size_comparison$best_panel_id), , drop = FALSE]
  if (nrow(valid) == 0L) {
    return(list(panel_size = NA_integer_, panel_id = NA_character_, reason = "No feasible panel was found at any allowed size."))
  }

  chosen <- valid[1L, , drop = FALSE]
  reason <- "Smallest feasible size with a panel to evaluate."
  if (nrow(valid) > 1L) {
    for (i in seq_len(nrow(valid) - 1L)) {
      current <- valid[i, , drop = FALSE]
      next_row <- valid[i + 1L, , drop = FALSE]
      improves_gain <- !is.na(next_row$delta_min_absolute_target_gain[[1L]]) &&
        next_row$delta_min_absolute_target_gain[[1L]] >= cfg$material_gain_improvement
      improves_relative <- !is.na(next_row$delta_min_target_relative_gain[[1L]]) &&
        next_row$delta_min_target_relative_gain[[1L]] >= cfg$material_relative_gain_improvement
      improves_redundancy <- !is.na(next_row$delta_capped_target_redundancy[[1L]]) &&
        next_row$delta_capped_target_redundancy[[1L]] > 0
      if (improves_gain || improves_relative || improves_redundancy) {
        chosen <- next_row
        reason <- paste0(
          "Size ", current$panel_size[[1L]], " -> ", next_row$panel_size[[1L]],
          " materially improved discrimination/redundancy."
        )
      } else {
        reason <- paste0(
          "Size ", current$panel_size[[1L]], " -> ", next_row$panel_size[[1L]],
          " did not materially improve the principal discriminatory metrics."
        )
        break
      }
    }
  }

  list(
    panel_size = as.integer(chosen$panel_size[[1L]]),
    panel_id = as.character(chosen$best_panel_id[[1L]]),
    reason = reason
  )
}

panel_best_alternatives <- function(panel_candidates, recommended_size, top_n_panels) {
  if (nrow(panel_candidates) == 0L || is.na(recommended_size)) {
    return(data.frame())
  }
  size_panels <- panel_candidates[panel_candidates$panel_size == recommended_size, , drop = FALSE]
  if (nrow(size_panels) == 0L) {
    return(data.frame())
  }
  size_panels <- sort_panels_within_size(size_panels)
  size_panels[seq_len(min(top_n_panels, nrow(size_panels))), , drop = FALSE]
}

panel_cross_size_pareto <- function(panel_candidates) {
  if (nrow(panel_candidates) == 0L) {
    return(panel_candidates)
  }
  panel_candidates[panel_candidates$cross_size_pareto, , drop = FALSE]
}

panel_within_size_pareto <- function(panel_candidates) {
  if (nrow(panel_candidates) == 0L) {
    return(panel_candidates)
  }
  panel_candidates[panel_candidates$within_size_pareto, , drop = FALSE]
}

panel_search_summary <- function(before_count, after_filter_count, after_reduction_count, panel_size_comparison, cfg) {
  search_by_size <- if (nrow(panel_size_comparison) == 0L) {
    "none"
  } else {
    paste(panel_size_comparison$panel_size, panel_size_comparison$search_method, sep = ":", collapse = ";")
  }
  recommended_row <- if (nrow(panel_size_comparison) == 0L || !any(panel_size_comparison$recommended_provisional_size, na.rm = TRUE)) {
    data.frame()
  } else {
    panel_size_comparison[which(panel_size_comparison$recommended_provisional_size)[1L], , drop = FALSE]
  }
  data.frame(
    candidate_count_before = as.integer(before_count),
    candidate_count_after_filter = as.integer(after_filter_count),
    candidate_count_after_reduction = as.integer(after_reduction_count),
    reduced_candidate_count = as.integer(after_reduction_count),
    panels_evaluated = if (nrow(panel_size_comparison) == 0L) 0L else as.integer(sum(panel_size_comparison$panels_evaluated, na.rm = TRUE)),
    search_method = search_by_size,
    exact_enumeration_limit = as.integer(cfg$exact_enumeration_limit),
    approximate_search_budget = as.integer(cfg$approximate_search_budget),
    min_panel_size = as.integer(cfg$min_panel_size),
    max_panel_size = as.integer(cfg$max_panel_size),
    top_n_panels = as.integer(cfg$top_n_panels),
    recommended_provisional_size = if (nrow(recommended_row) == 0L) NA_integer_ else as.integer(recommended_row$panel_size[[1L]]),
    recommended_provisional_panel_id = if (nrow(recommended_row) == 0L) NA_character_ else as.character(recommended_row$best_panel_id[[1L]]),
    recommended_reason = if (nrow(recommended_row) == 0L) NA_character_ else as.character(recommended_row$recommendation_reason[[1L]]),
    stringsAsFactors = FALSE
  )
}

select_amplicon_panels <- function(candidate_amplicons, candidate_amplicon_target_scores,
                                   candidate_amplicon_snps = NULL, config) {
  cfg <- resolve_panel_selection_config(config)
  candidate_table <- build_candidate_universe(
    candidate_amplicons = candidate_amplicons,
    candidate_amplicon_target_scores = candidate_amplicon_target_scores,
    candidate_amplicon_snps = candidate_amplicon_snps,
    config = config
  )
  before_count <- nrow(candidate_table)
  filtered <- filter_candidate_universe(candidate_table, cfg)
  after_filter_count <- nrow(filtered)
  reduced <- reduce_candidate_universe(filtered, cfg)
  after_reduction_count <- nrow(reduced)

  evaluated <- evaluate_panel_sizes(reduced, cfg)
  panel_candidates <- evaluated$panel_candidates
  panel_size_comparison <- evaluated$panel_size_comparison
  recommended_size <- evaluated$recommended_panel_size
  recommended_reason <- evaluated$recommended_reason
  recommended_panel_id <- evaluated$recommended_panel_id

  if (nrow(panel_candidates) > 0L) {
    panel_target_scores <- panel_target_table(panel_candidates, reduced, cfg)
    panel_candidates$panel_rank <- panel_candidates$deterministic_rank
    panel_target_scores$panel_rank <- match(panel_target_scores$panel_id, panel_candidates$panel_id)
  } else {
    panel_target_scores <- data.frame()
  }

  pareto_panels_by_size <- panel_within_size_pareto(panel_candidates)
  cross_size_pareto_panels <- panel_cross_size_pareto(panel_candidates)
  recommended_provisional_panel <- if (nrow(panel_candidates) > 0L && !is.na(recommended_panel_id)) {
    panel_members_table(panel_candidates[panel_candidates$panel_id == recommended_panel_id, , drop = FALSE], reduced)
  } else {
    data.frame()
  }
  if (nrow(recommended_provisional_panel) > 0L) {
    recommended_provisional_panel$recommended_provisional_size <- recommended_size
    recommended_provisional_panel$recommendation_reason <- recommended_reason
  }
  recommended_panel_alternatives <- panel_best_alternatives(panel_candidates, recommended_size, cfg$top_n_panels)
  selected_panel <- recommended_provisional_panel
  top_panels <- recommended_panel_alternatives

  list(
    candidate_table = reduced,
    panel_candidates = panel_candidates,
    panel_target_scores = panel_target_scores,
    pareto_panels_by_size = pareto_panels_by_size,
    cross_size_pareto_panels = cross_size_pareto_panels,
    pareto_panels = cross_size_pareto_panels,
    top_panels = top_panels,
    recommended_panel_alternatives = recommended_panel_alternatives,
    panel_size_comparison = panel_size_comparison,
    recommended_provisional_panel = recommended_provisional_panel,
    selected_panel = selected_panel,
    search_method = if (nrow(panel_size_comparison) == 0L) "none" else paste(panel_size_comparison$search_method, collapse = ";"),
    panels_evaluated = if (nrow(panel_size_comparison) == 0L) 0L else sum(panel_size_comparison$panels_evaluated, na.rm = TRUE),
    candidate_count_before = before_count,
    candidate_count_after_filter = after_filter_count,
    candidate_count_after_reduction = after_reduction_count,
    recommended_panel_size = recommended_size,
    recommended_reason = recommended_reason,
    recommended_panel_id = recommended_panel_id,
    config = cfg
  )
}
