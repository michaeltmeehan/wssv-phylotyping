#' Validate and normalise marker-window coordinate tables.
validate_windows <- function(windows, alignment_length = NULL) {
  required <- c("window_id", "start", "end")
  missing <- setdiff(required, names(windows))
  if (length(missing) > 0L) {
    stop("windows is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- windows
  out$start <- as.integer(out$start)
  out$end <- as.integer(out$end)
  bad <- is.na(out$start) | is.na(out$end) | out$start < 1L | out$end < out$start
  if (any(bad)) {
    stop("windows contains invalid start/end coordinates.", call. = FALSE)
  }
  if (!is.null(alignment_length) && any(out$end > alignment_length)) {
    stop("windows extend beyond the alignment length.", call. = FALSE)
  }
  out$width <- out$end - out$start + 1L
  out
}

extract_window_columns <- function(aln_int, windows, keep_order = TRUE) {
  windows <- validate_windows(windows, ncol(aln_int))
  cols <- unlist(Map(seq.int, windows$start, windows$end), use.names = FALSE)
  if (!isTRUE(keep_order)) {
    cols <- sort(unique(cols))
  }
  aln_int[, cols, drop = FALSE]
}

window_column_indices <- function(windows, alignment_length = NULL) {
  windows <- validate_windows(windows, alignment_length)
  sort(unique(unlist(Map(seq.int, windows$start, windows$end), use.names = FALSE)))
}

mask_alignment_to_windows <- function(aln_int, windows, mask_value = 0L) {
  cols <- window_column_indices(windows, ncol(aln_int))
  out <- matrix(mask_value, nrow = nrow(aln_int), ncol = ncol(aln_int), dimnames = dimnames(aln_int))
  out[, cols] <- aln_int[, cols, drop = FALSE]
  out
}

make_random_fragments <- function(alignment_length, fragment_lengths, n_per_length = 1L,
                                  seed = NULL, label = "random_fixed") {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  fragment_lengths <- as.integer(fragment_lengths)
  fragment_lengths <- fragment_lengths[!is.na(fragment_lengths) & fragment_lengths > 0L]
  if (length(fragment_lengths) == 0L) {
    return(empty_fragment_windows())
  }
  rows <- list()
  k <- 0L
  for (len in fragment_lengths) {
    width <- min(len, alignment_length)
    max_start <- max(1L, alignment_length - width + 1L)
    starts <- sample(seq_len(max_start), size = as.integer(n_per_length), replace = TRUE)
    for (start in starts) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        fragment_id = sprintf("%s_%s_%03d", label, width, k),
        fragment_type = label,
        start = as.integer(start),
        end = as.integer(start + width - 1L),
        width = as.integer(width)
      )
    }
  }
  do.call(rbind, rows)
}

make_panel_fragments <- function(panel_windows, panel_size = NULL, label = "selected_panel") {
  windows <- validate_windows(panel_windows)
  if (!is.null(panel_size)) {
    windows <- windows[seq_len(min(as.integer(panel_size), nrow(windows))), , drop = FALSE]
  }
  data.frame(
    fragment_id = paste(label, seq_len(nrow(windows)), windows$window_id, sep = "_"),
    fragment_type = label,
    window_id = windows$window_id,
    start = windows$start,
    end = windows$end,
    width = windows$width
  )
}

make_matched_random_windows <- function(panel_windows, alignment_length, n_replicates = 1L,
                                        seed = NULL, label = "random_matched") {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  panel_windows <- validate_windows(panel_windows, alignment_length)
  rows <- list()
  k <- 0L
  for (replicate in seq_len(as.integer(n_replicates))) {
    for (i in seq_len(nrow(panel_windows))) {
      width <- min(panel_windows$width[[i]], alignment_length)
      start <- sample(seq_len(max(1L, alignment_length - width + 1L)), 1L)
      k <- k + 1L
      rows[[k]] <- data.frame(
        replicate = replicate,
        window_id = sprintf("%s_%03d_%03d", label, replicate, i),
        window_type = label,
        start = as.integer(start),
        end = as.integer(start + width - 1L),
        width = as.integer(width),
        requested_width = as.integer(width),
        centre_site = NA_integer_
      )
    }
  }
  do.call(rbind, rows)
}

#' Count informative signal observed in masked or partial alignments.
#'
#' Used for complete-genome panel validation and synthetic fragment summaries.
count_observed_signal <- function(partial_aln_int, site_node_scores, target_mask = NULL,
                                  min_informative_sites = 1L, min_informative_nodes = 1L) {
  tips <- rownames(partial_aln_int)
  if (is.null(tips)) {
    tips <- as.character(seq_len(nrow(partial_aln_int)))
  }
  if (nrow(site_node_scores) == 0L) {
    return(empty_tip_signal_summary(tips))
  }
  sites <- sort(unique(as.integer(site_node_scores$site)))
  sites <- sites[sites >= 1L & sites <= ncol(partial_aln_int)]
  score_rows <- site_node_scores[site_node_scores$site %in% sites, , drop = FALSE]
  by_tip <- lapply(seq_len(nrow(partial_aln_int)), function(tip_i) {
    observed_sites <- sites[partial_aln_int[tip_i, sites] != 0L]
    observed_rows <- score_rows[score_rows$site %in% observed_sites, , drop = FALSE]
    observed_nodes <- unique(observed_rows$node_id)
    true_rows <- observed_rows
    if (!is.null(target_mask) && nrow(observed_rows) > 0L) {
      tip_col <- if (!is.null(colnames(target_mask)) && tips[[tip_i]] %in% colnames(target_mask)) {
        match(tips[[tip_i]], colnames(target_mask))
      } else {
        tip_i
      }
      true_node <- target_mask[match(observed_rows$node_id, rownames(target_mask)), tip_col]
      true_rows <- observed_rows[!is.na(true_node) & true_node, , drop = FALSE]
    }
    supported <- count_supported_true_nodes(partial_aln_int[tip_i, , drop = TRUE], true_rows)
    data.frame(
      tip = tips[[tip_i]],
      observed_informative_sites = length(observed_sites),
      observed_informative_nodes = length(observed_nodes),
      true_path_informative_nodes = length(unique(true_rows$node_id)),
      supported_true_nodes = supported$supported_nodes,
      deepest_supported_depth = supported$deepest_supported_depth,
      resolved_signal = length(observed_sites) >= min_informative_sites &&
        length(observed_nodes) >= min_informative_nodes,
      insufficient_information = length(observed_sites) < min_informative_sites ||
        length(observed_nodes) < min_informative_nodes
    )
  })
  do.call(rbind, by_tip)
}

count_supported_true_nodes <- function(x, true_site_node_rows) {
  if (nrow(true_site_node_rows) == 0L ||
      !"best_allele" %in% names(true_site_node_rows) ||
      !"direction" %in% names(true_site_node_rows)) {
    return(list(supported_nodes = 0L, deepest_supported_depth = NA_real_))
  }
  allele_code <- c(A = 1L, C = 2L, G = 3L, T = 4L)
  rows <- true_site_node_rows
  expected <- allele_code[as.character(rows$best_allele)]
  observed <- x[as.integer(rows$site)]
  supports <- ifelse(
    rows$direction == "clade",
    observed == expected,
    observed != 0L & observed != expected
  )
  supported_rows <- rows[!is.na(supports) & supports, , drop = FALSE]
  if (nrow(supported_rows) == 0L) {
    return(list(supported_nodes = 0L, deepest_supported_depth = NA_real_))
  }
  depth <- if ("depth" %in% names(supported_rows)) max(supported_rows$depth, na.rm = TRUE) else NA_real_
  if (!is.finite(depth)) {
    depth <- NA_real_
  }
  list(supported_nodes = length(unique(supported_rows$node_id)), deepest_supported_depth = depth)
}

#' Evaluate retained signal for a selected panel size across complete genomes.
evaluate_panel_signal <- function(aln_int, windows, site_node_scores, target_mask = NULL,
                                  panel_size = nrow(windows), method = "selected_panel",
                                  replicate = NA_integer_, min_informative_sites = 1L,
                                  min_informative_nodes = 1L) {
  selected <- validate_windows(windows)[seq_len(min(panel_size, nrow(windows))), , drop = FALSE]
  masked <- mask_alignment_to_windows(aln_int, selected)
  out <- count_observed_signal(
    masked, site_node_scores, target_mask,
    min_informative_sites = min_informative_sites,
    min_informative_nodes = min_informative_nodes
  )
  out$method <- method
  out$replicate <- replicate
  out$panel_size <- nrow(selected)
  out$panel_width <- sum(selected$width)
  out[c("method", "replicate", "panel_size", "panel_width", setdiff(names(out), c("method", "replicate", "panel_size", "panel_width")))]
}

summarise_validation_by_panel <- function(tip_summary) {
  if (nrow(tip_summary) == 0L) {
    return(data.frame())
  }
  grouped <- split(tip_summary, paste(tip_summary$method, tip_summary$replicate, tip_summary$panel_size, sep = "\r"))
  out <- do.call(rbind, lapply(grouped, function(x) {
    data.frame(
      method = x$method[[1L]],
      replicate = x$replicate[[1L]],
      panel_size = x$panel_size[[1L]],
      panel_width = x$panel_width[[1L]],
      tips_assessed = nrow(x),
      mean_observed_informative_sites = mean(x$observed_informative_sites),
      median_observed_informative_sites = stats::median(x$observed_informative_sites),
      mean_observed_informative_nodes = mean(x$observed_informative_nodes),
      median_observed_informative_nodes = stats::median(x$observed_informative_nodes),
      mean_supported_true_nodes = mean(x$supported_true_nodes),
      resolved_tips = sum(x$resolved_signal),
      unresolved_tips = sum(x$insufficient_information),
      resolved_fraction = mean(x$resolved_signal)
    )
  }))
  out[order(out$method, out$replicate, out$panel_size), , drop = FALSE]
}

summarise_fragment_signal <- function(aln_int, fragments, site_node_scores,
                                      min_informative_sites = 1L, min_informative_nodes = 1L) {
  if (nrow(fragments) == 0L) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(fragments)), function(i) {
    fragment <- fragments[i, , drop = FALSE]
    fragment$window_id <- fragment$fragment_id
    masked <- mask_alignment_to_windows(aln_int, fragment)
    tip_rows <- count_observed_signal(masked, site_node_scores, NULL, min_informative_sites, min_informative_nodes)
    data.frame(
      fragment_id = fragment$fragment_id,
      fragment_type = fragment$fragment_type,
      width = fragment$width,
      start = fragment$start,
      end = fragment$end,
      tips_assessed = nrow(tip_rows),
      mean_observed_informative_sites = mean(tip_rows$observed_informative_sites),
      mean_observed_informative_nodes = mean(tip_rows$observed_informative_nodes),
      resolved_fraction = mean(tip_rows$resolved_signal)
    )
  })
  do.call(rbind, rows)
}

empty_fragment_windows <- function() {
  data.frame(fragment_id = character(), fragment_type = character(), start = integer(), end = integer(), width = integer())
}

empty_tip_signal_summary <- function(tips) {
  data.frame(
    tip = tips,
    observed_informative_sites = integer(length(tips)),
    observed_informative_nodes = integer(length(tips)),
    true_path_informative_nodes = integer(length(tips)),
    supported_true_nodes = integer(length(tips)),
    deepest_supported_depth = rep(NA_real_, length(tips)),
    resolved_signal = rep(FALSE, length(tips)),
    insufficient_information = rep(TRUE, length(tips))
  )
}
