#' Filter candidate windows before greedy panel selection.
#'
#' Applies student-tunable constraints such as minimum informative SNPs, missing
#' fraction, permitted widths, and optional near-duplicate removal.
filter_candidate_windows <- function(window_summary,
                                     min_informative_snps = 1L,
                                     min_total_weighted_gain = 0,
                                     max_missing_fraction = NULL,
                                     permitted_widths = NULL,
                                     exclude_near_duplicates = FALSE,
                                     near_duplicate_max_overlap = 1,
                                     near_duplicate_min_distance = NULL) {
  required <- c("window_id", "start", "end", "width", "informative_snps", "total_weighted_gain")
  missing <- setdiff(required, names(window_summary))
  if (length(missing) > 0L) {
    stop("window_summary is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  out <- window_summary
  out <- out[!is.na(out$informative_snps) & out$informative_snps >= min_informative_snps, , drop = FALSE]
  out <- out[!is.na(out$total_weighted_gain) & out$total_weighted_gain >= min_total_weighted_gain, , drop = FALSE]

  if (!is.null(max_missing_fraction) && "missing_fraction" %in% names(out)) {
    out <- out[is.na(out$missing_fraction) | out$missing_fraction <= max_missing_fraction, , drop = FALSE]
  }

  if (!is.null(permitted_widths)) {
    permitted_widths <- as.integer(permitted_widths)
    out <- out[out$width %in% permitted_widths, , drop = FALSE]
  }

  out <- order_windows_for_selection(out)
  row.names(out) <- NULL

  if (isTRUE(exclude_near_duplicates)) {
    out <- retain_nonredundant_windows(
      out,
      max_overlap_fraction = near_duplicate_max_overlap,
      min_distance = near_duplicate_min_distance
    )
  }

  out
}

#' Calculate each candidate window's node-level coverage score.
#'
#' Converts window-node summaries into a compact table used to estimate marginal
#' gain while selecting complementary panel windows.
calculate_window_node_coverage <- function(window_node_summary,
                                           retained_window_ids = NULL,
                                           score_col = "total_weighted_gain") {
  required <- c("window_id", "node_id")
  missing <- setdiff(required, names(window_node_summary))
  if (length(missing) > 0L) {
    stop("window_node_summary is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!score_col %in% names(window_node_summary)) {
    stop("window_node_summary is missing score column: ", score_col, call. = FALSE)
  }

  out <- window_node_summary
  if (!is.null(retained_window_ids)) {
    out <- out[out$window_id %in% retained_window_ids, , drop = FALSE]
  }
  if (nrow(out) == 0L) {
    return(data.frame(window_id = character(), node_id = character(), node_score = numeric()))
  }

  grouped <- split(out, paste(out$window_id, out$node_id, sep = "\r"))
  coverage <- do.call(rbind, lapply(grouped, function(x) {
    data.frame(
      window_id = as.character(x$window_id[[1L]]),
      node_id = as.character(x$node_id[[1L]]),
      node_score = sum(x[[score_col]], na.rm = TRUE)
    )
  }))

  coverage[order(coverage$window_id, coverage$node_id), , drop = FALSE]
}

calculate_marginal_gain <- function(candidate_window_ids,
                                    node_coverage,
                                    selected_window_ids = character(),
                                    repeated_node_credit = 0,
                                    cap_repeated_coverage = TRUE) {
  candidate_window_ids <- as.character(candidate_window_ids)
  selected_window_ids <- as.character(selected_window_ids)
  coverage_by_window <- split(node_coverage, node_coverage$window_id)
  selected_coverage <- node_coverage[node_coverage$window_id %in% selected_window_ids, , drop = FALSE]
  covered_nodes <- unique(selected_coverage$node_id)

  rows <- lapply(candidate_window_ids, function(window_id) {
    node_rows <- coverage_by_window[[window_id]]
    if (is.null(node_rows) || nrow(node_rows) == 0L) {
      return(data.frame(
        window_id = window_id,
        marginal_gain = 0,
        raw_gain = 0,
        newly_covered_nodes = 0L,
        covered_nodes_after = length(covered_nodes),
        new_node_ids = ""
      ))
    }

    is_new <- !node_rows$node_id %in% covered_nodes
    repeated_credit <- if (isTRUE(cap_repeated_coverage)) 0 else repeated_node_credit
    credit <- ifelse(is_new, 1, repeated_credit)
    new_nodes <- sort(unique(node_rows$node_id[is_new]))

    data.frame(
      window_id = window_id,
      marginal_gain = sum(node_rows$node_score * credit, na.rm = TRUE),
      raw_gain = sum(node_rows$node_score, na.rm = TRUE),
      newly_covered_nodes = length(new_nodes),
      covered_nodes_after = length(unique(c(covered_nodes, new_nodes))),
      new_node_ids = paste(new_nodes, collapse = ";")
    )
  })

  do.call(rbind, rows)
}

#' Greedily select complementary marker windows.
#'
#' At each step, choose the compatible remaining window with the largest
#' marginal gain under overlap, distance, and repeated-node-credit settings.
greedy_select_panel <- function(window_summary,
                                node_coverage,
                                max_panel_size = 10L,
                                min_marginal_gain = 0,
                                max_allowed_overlap = 1,
                                min_distance = NULL,
                                repeated_node_credit = 0,
                                cap_repeated_coverage = TRUE) {
  if (nrow(window_summary) == 0L || max_panel_size < 1L) {
    empty <- empty_panel_steps()
    return(list(selected_panel = empty, steps = empty, node_coverage = empty_panel_node_coverage()))
  }

  candidates <- order_windows_for_selection(window_summary)
  selected_ids <- character()
  steps <- list()
  cumulative_gain <- 0

  for (step in seq_len(max_panel_size)) {
    available <- candidates[!candidates$window_id %in% selected_ids, , drop = FALSE]
    if (length(selected_ids) > 0L) {
      available <- available[is_compatible_with_panel(
        available,
        candidates[candidates$window_id %in% selected_ids, , drop = FALSE],
        max_allowed_overlap = max_allowed_overlap,
        min_distance = min_distance
      ), , drop = FALSE]
    }
    if (nrow(available) == 0L) {
      break
    }

    gains <- calculate_marginal_gain(
      candidate_window_ids = available$window_id,
      node_coverage = node_coverage,
      selected_window_ids = selected_ids,
      repeated_node_credit = repeated_node_credit,
      cap_repeated_coverage = cap_repeated_coverage
    )
    ranked <- merge(available, gains, by = "window_id", sort = FALSE)
    ranked <- ranked[order(
      -ranked$marginal_gain,
      -ranked$total_weighted_gain,
      -ranked$informative_snps,
      ranked$start,
      ranked$end,
      ranked$window_id
    ), , drop = FALSE]

    best <- ranked[1L, , drop = FALSE]
    if (is.na(best$marginal_gain) || best$marginal_gain < min_marginal_gain) {
      break
    }

    cumulative_gain <- cumulative_gain + best$marginal_gain
    selected_ids <- c(selected_ids, best$window_id)
    best$selection_step <- step
    best$cumulative_gain <- cumulative_gain
    best$total_covered_nodes <- best$covered_nodes_after
    steps[[step]] <- best
  }

  if (length(steps) == 0L) {
    empty <- empty_panel_steps()
    return(list(selected_panel = empty, steps = empty, node_coverage = empty_panel_node_coverage()))
  }

  panel_steps <- do.call(rbind, steps)
  panel_steps <- panel_steps[order(panel_steps$selection_step), , drop = FALSE]
  row.names(panel_steps) <- NULL

  list(
    selected_panel = panel_steps,
    steps = panel_steps,
    node_coverage = summarise_panel_node_coverage(panel_steps$window_id, node_coverage)
  )
}

summarise_selected_panels <- function(selected_panel,
                                      node_coverage,
                                      panel_sizes = c(1L, 2L, 3L, 5L, 10L),
                                      label = "greedy") {
  panel_sizes <- sort(unique(as.integer(panel_sizes)))
  panel_sizes <- panel_sizes[!is.na(panel_sizes) & panel_sizes > 0L]
  if (nrow(selected_panel) == 0L || length(panel_sizes) == 0L) {
    return(empty_panel_summary())
  }

  rows <- lapply(panel_sizes[panel_sizes <= nrow(selected_panel)], function(size) {
    ids <- selected_panel$window_id[seq_len(size)]
    covered <- summarise_panel_node_coverage(ids, node_coverage)
    data.frame(
      method = label,
      panel_size = size,
      windows_selected = length(ids),
      cumulative_gain = sum(selected_panel$marginal_gain[seq_len(size)], na.rm = TRUE),
      total_window_score = sum(selected_panel$total_weighted_gain[seq_len(size)], na.rm = TRUE),
      nodes_covered = length(unique(covered$node_id)),
      newly_covered_nodes_last_step = selected_panel$newly_covered_nodes[[size]]
    )
  })

  do.call(rbind, rows)
}

summarise_panel_node_coverage_by_size <- function(selected_panel,
                                                  node_coverage,
                                                  panel_sizes = c(1L, 2L, 3L, 5L, 10L),
                                                  label = "greedy") {
  panel_sizes <- sort(unique(as.integer(panel_sizes)))
  panel_sizes <- panel_sizes[!is.na(panel_sizes) & panel_sizes > 0L]
  if (nrow(selected_panel) == 0L || length(panel_sizes) == 0L) {
    out <- empty_panel_node_coverage()
    out$method <- character()
    out$panel_size <- integer()
    return(out[c("method", "panel_size", names(empty_panel_node_coverage()))])
  }

  rows <- lapply(panel_sizes[panel_sizes <= nrow(selected_panel)], function(size) {
    ids <- selected_panel$window_id[seq_len(size)]
    covered <- summarise_panel_node_coverage(ids, node_coverage)
    if (nrow(covered) == 0L) {
      return(NULL)
    }
    covered$method <- label
    covered$panel_size <- size
    covered[c("method", "panel_size", setdiff(names(covered), c("method", "panel_size")))]
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L) {
    out <- empty_panel_node_coverage()
    out$method <- character()
    out$panel_size <- integer()
    return(out[c("method", "panel_size", names(empty_panel_node_coverage()))])
  }
  do.call(rbind, rows)
}

compare_panel_alternatives <- function(window_summary,
                                       node_coverage,
                                       panel_sizes = c(1L, 2L, 3L, 5L, 10L),
                                       max_allowed_overlap = 1,
                                       min_distance = NULL) {
  ordered <- order_windows_for_selection(window_summary)
  top_rows <- summarise_ranked_alternative("top_weighted_gain", ordered, node_coverage, panel_sizes)
  nonoverlap <- retain_nonredundant_windows(
    ordered,
    max_overlap_fraction = max_allowed_overlap,
    min_distance = min_distance
  )
  nonoverlap_rows <- summarise_ranked_alternative(
    "nonoverlapping_top_weighted_gain",
    nonoverlap,
    node_coverage,
    panel_sizes
  )
  rbind(top_rows, nonoverlap_rows)
}

summarise_panel_node_coverage <- function(window_ids, node_coverage) {
  window_ids <- as.character(window_ids)
  selected <- node_coverage[node_coverage$window_id %in% window_ids, , drop = FALSE]
  if (nrow(selected) == 0L) {
    return(empty_panel_node_coverage())
  }

  selected$selection_rank <- match(selected$window_id, window_ids)
  grouped <- split(selected, selected$node_id)
  out <- do.call(rbind, lapply(grouped, function(x) {
    x <- x[order(x$selection_rank), , drop = FALSE]
    data.frame(
      node_id = x$node_id[[1L]],
      coverage_count = length(unique(x$window_id)),
      total_node_score = sum(x$node_score, na.rm = TRUE),
      first_covered_step = min(x$selection_rank, na.rm = TRUE),
      covering_windows = paste(unique(x$window_id), collapse = ";")
    )
  }))

  out[order(out$first_covered_step, out$node_id), , drop = FALSE]
}

order_windows_for_selection <- function(windows) {
  if (nrow(windows) == 0L) {
    return(windows)
  }
  windows[order(
    -windows$total_weighted_gain,
    -windows$informative_snps,
    windows$start,
    windows$end,
    windows$window_id
  ), , drop = FALSE]
}

retain_nonredundant_windows <- function(windows, max_overlap_fraction = 1, min_distance = NULL) {
  if (nrow(windows) <= 1L) {
    return(windows)
  }
  ordered <- order_windows_for_selection(windows)
  keep <- logical(nrow(ordered))
  kept_rows <- ordered[FALSE, , drop = FALSE]
  for (i in seq_len(nrow(ordered))) {
    current <- ordered[i, , drop = FALSE]
    if (nrow(kept_rows) == 0L || is_compatible_with_panel(
      current,
      kept_rows,
      max_allowed_overlap = max_overlap_fraction,
      min_distance = min_distance
    )) {
      keep[[i]] <- TRUE
      kept_rows <- rbind(kept_rows, current)
    }
  }
  out <- ordered[keep, , drop = FALSE]
  row.names(out) <- NULL
  out
}

is_compatible_with_panel <- function(candidates, selected, max_allowed_overlap = 1, min_distance = NULL) {
  if (nrow(candidates) == 0L || nrow(selected) == 0L) {
    return(rep(TRUE, nrow(candidates)))
  }
  vapply(seq_len(nrow(candidates)), function(i) {
    overlaps <- window_overlap_fraction(
      candidates$start[[i]],
      candidates$end[[i]],
      selected$start,
      selected$end
    )
    distances <- window_distances(candidates$start[[i]], candidates$end[[i]], selected$start, selected$end)
    overlap_ok <- all(overlaps <= max_allowed_overlap)
    distance_ok <- is.null(min_distance) || all(distances >= min_distance)
    overlap_ok && distance_ok
  }, logical(1L))
}

window_overlap_fraction <- function(start, end, other_start, other_end) {
  overlap <- pmax(0L, pmin(end, other_end) - pmax(start, other_start) + 1L)
  width <- end - start + 1L
  overlap / width
}

window_distances <- function(start, end, other_start, other_end) {
  ifelse(end < other_start, other_start - end - 1L, ifelse(other_end < start, start - other_end - 1L, 0L))
}

summarise_ranked_alternative <- function(method, ranked_windows, node_coverage, panel_sizes) {
  panel_sizes <- sort(unique(as.integer(panel_sizes)))
  panel_sizes <- panel_sizes[!is.na(panel_sizes) & panel_sizes > 0L & panel_sizes <= nrow(ranked_windows)]
  if (length(panel_sizes) == 0L) {
    return(empty_panel_summary())
  }

  rows <- lapply(panel_sizes, function(size) {
    ids <- ranked_windows$window_id[seq_len(size)]
    covered <- summarise_panel_node_coverage(ids, node_coverage)
    data.frame(
      method = method,
      panel_size = size,
      windows_selected = length(ids),
      cumulative_gain = NA_real_,
      total_window_score = sum(ranked_windows$total_weighted_gain[seq_len(size)], na.rm = TRUE),
      nodes_covered = length(unique(covered$node_id)),
      newly_covered_nodes_last_step = NA_integer_
    )
  })
  do.call(rbind, rows)
}

empty_panel_steps <- function() {
  data.frame(
    window_id = character(),
    selection_step = integer(),
    start = integer(),
    end = integer(),
    width = integer(),
    total_weighted_gain = numeric(),
    marginal_gain = numeric(),
    cumulative_gain = numeric(),
    newly_covered_nodes = integer(),
    total_covered_nodes = integer(),
    new_node_ids = character()
  )
}

empty_panel_node_coverage <- function() {
  data.frame(
    node_id = character(),
    coverage_count = integer(),
    total_node_score = numeric(),
    first_covered_step = integer(),
    covering_windows = character()
  )
}

empty_panel_summary <- function() {
  data.frame(
    method = character(),
    panel_size = integer(),
    windows_selected = integer(),
    cumulative_gain = numeric(),
    total_window_score = numeric(),
    nodes_covered = integer(),
    newly_covered_nodes_last_step = integer()
  )
}
