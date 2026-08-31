#' Generate fixed-width windows across the alignment.
#'
#' Windows are created for each requested width using either an explicit step,
#' an overlap fraction, or non-overlapping defaults.
generate_fixed_windows <- function(alignment_length, widths, step = NULL, overlap = NULL) {
  if (alignment_length < 1L) {
    stop("alignment_length must be positive.", call. = FALSE)
  }
  widths <- unique(as.integer(widths))
  widths <- widths[!is.na(widths) & widths > 0L]
  if (length(widths) == 0L) {
    stop("At least one positive window width is required.", call. = FALSE)
  }

  rows <- vector("list", length(widths))
  for (i in seq_along(widths)) {
    width <- widths[[i]]
    this_step <- resolve_window_step(width, step, overlap)
    starts <- seq.int(1L, alignment_length, by = this_step)
    ends <- pmin(starts + width - 1L, alignment_length)
    keep <- ends >= starts
    rows[[i]] <- data.frame(
      window_type = "fixed",
      width = as.integer(ends[keep] - starts[keep] + 1L),
      requested_width = width,
      start = as.integer(starts[keep]),
      end = as.integer(ends[keep]),
      centre_site = NA_integer_
    )
  }

  add_window_ids(do.call(rbind, rows))
}

resolve_window_step <- function(width, step = NULL, overlap = NULL) {
  if (!is.null(step)) {
    if (length(step) == 1L) {
      out <- as.integer(step)
    } else {
      names_step <- names(step)
      out <- if (!is.null(names_step) && as.character(width) %in% names_step) {
        as.integer(step[[as.character(width)]])
      } else {
        NA_integer_
      }
    }
    if (!is.na(out) && out > 0L) {
      return(out)
    }
  }

  if (!is.null(overlap)) {
    overlap <- as.numeric(overlap)
    if (length(overlap) != 1L || is.na(overlap) || overlap < 0 || overlap >= 1) {
      stop("overlap must be a single value in [0, 1).", call. = FALSE)
    }
    return(max(1L, as.integer(round(width * (1 - overlap)))))
  }

  width
}

#' Generate windows centred on informative SNP coordinates.
#'
#' These windows supplement fixed tiling when the config enables
#' analysis.windows.snp_centered.
generate_snp_centered_windows <- function(sites, alignment_length, widths) {
  sites <- sort(unique(as.integer(sites)))
  sites <- sites[!is.na(sites) & sites >= 1L & sites <= alignment_length]
  widths <- unique(as.integer(widths))
  widths <- widths[!is.na(widths) & widths > 0L]
  if (length(sites) == 0L || length(widths) == 0L) {
    return(empty_windows())
  }

  rows <- vector("list", length(widths))
  for (i in seq_along(widths)) {
    width <- widths[[i]]
    left <- floor((width - 1L) / 2L)
    starts <- pmax(1L, sites - left)
    ends <- pmin(alignment_length, starts + width - 1L)
    starts <- pmax(1L, ends - width + 1L)
    rows[[i]] <- data.frame(
      window_type = "snp_centered",
      width = as.integer(ends - starts + 1L),
      requested_width = width,
      start = as.integer(starts),
      end = as.integer(ends),
      centre_site = as.integer(sites)
    )
  }

  out <- unique(do.call(rbind, rows))
  add_window_ids(out[order(out$requested_width, out$start, out$end, out$centre_site), ])
}

assign_snps_to_windows <- function(sites, windows) {
  sites <- sort(unique(as.integer(sites)))
  if (length(sites) == 0L || nrow(windows) == 0L) {
    return(data.frame(window_id = character(), site = integer()))
  }

  rows <- vector("list", nrow(windows))
  for (i in seq_len(nrow(windows))) {
    in_window <- sites[sites >= windows$start[[i]] & sites <= windows$end[[i]]]
    rows[[i]] <- data.frame(
      window_id = rep(windows$window_id[[i]], length(in_window)),
      site = in_window
    )
  }
  do.call(rbind, rows)
}

#' Aggregate informative SNP scores by window and node.
#'
#' The output feeds panel selection by showing which nodes each window helps and
#' how much weighted signal the window contributes. The `weight` field is the
#' legacy node-weight compatibility value from stage 01 metadata, while the SNP
#' score comes from stage-02 `normalized_gain`; neither should be confused with
#' classifier support thresholds.
aggregate_window_node_scores <- function(windows, site_node_scores, informative_assignments = NULL) {
  if (is.null(informative_assignments)) {
    informative_assignments <- assign_snps_to_windows(site_node_scores$site, windows)
  }
  if (nrow(informative_assignments) == 0L || nrow(site_node_scores) == 0L) {
    return(empty_window_node_summary(site_node_scores))
  }

  scored <- merge(informative_assignments, site_node_scores, by = "site", sort = FALSE)
  if (nrow(scored) == 0L) {
    return(empty_window_node_summary(site_node_scores))
  }

  rows <- split(scored, paste(scored$window_id, scored$node_id, sep = "\r"))
  out <- do.call(rbind, lapply(rows, function(x) {
    data.frame(
      window_id = x$window_id[[1L]],
      node_index = x$node_index[[1L]],
      node_id = x$node_id[[1L]],
      informative_snps = length(unique(x$site)),
      node_site_score_entries = nrow(x),
      total_gain = sum(x$gain),
      total_weighted_gain = sum(x$weight * x$normalized_gain),
      mean_gain = mean(x$gain),
      max_gain = max(x$gain)
    )
  }))

  metadata_cols <- setdiff(names(site_node_scores), c(
    "site", "gain", "normalized_gain", "best_allele", "direction",
    "allele_count_inside", "allele_count_outside", "observed_inside",
    "observed_outside", "total_observed"
  ))
  metadata_cols <- setdiff(metadata_cols, names(out))
  if (length(metadata_cols) > 0L) {
    metadata <- unique(site_node_scores[c("node_index", "node_id", metadata_cols)])
    out <- merge(out, metadata, by = c("node_index", "node_id"), sort = FALSE)
  }

  out[order(out$window_id, out$node_index), ]
}

summarise_windows <- function(windows, polymorphic_sites, site_node_scores,
                              window_node_summary, aln_int = NULL,
                              min_informative_snps = 1L, deep_depth_quantile = 0.5) {
  polymorphic_assignments <- assign_snps_to_windows(polymorphic_sites, windows)
  informative_sites <- unique(site_node_scores$site)
  informative_assignments <- assign_snps_to_windows(informative_sites, windows)
  grouped_nodes <- split(window_node_summary, window_node_summary$window_id)
  grouped_poly <- split(polymorphic_assignments$site, polymorphic_assignments$window_id)
  grouped_info <- split(informative_assignments$site, informative_assignments$window_id)

  deep_cutoff <- NA_real_
  if ("depth" %in% names(window_node_summary) && any(!is.na(window_node_summary$depth))) {
    deep_cutoff <- as.numeric(stats::quantile(window_node_summary$depth, deep_depth_quantile, na.rm = TRUE))
  }

  rows <- lapply(seq_len(nrow(windows)), function(i) {
    window_id <- windows$window_id[[i]]
    node_rows <- grouped_nodes[[window_id]]
    n_poly <- length(unique(grouped_poly[[window_id]]))
    n_info <- length(unique(grouped_info[[window_id]]))

    out <- data.frame(
      window_id = window_id,
      window_type = windows$window_type[[i]],
      start = windows$start[[i]],
      end = windows$end[[i]],
      width = windows$width[[i]],
      requested_width = windows$requested_width[[i]],
      polymorphic_snps = n_poly,
      informative_snps = n_info,
      node_site_score_entries = if (is.null(node_rows)) 0L else sum(node_rows$node_site_score_entries),
      total_gain = if (is.null(node_rows)) 0 else sum(node_rows$total_gain),
      total_weighted_gain = if (is.null(node_rows)) 0 else sum(node_rows$total_weighted_gain),
      mean_gain = if (is.null(node_rows)) NA_real_ else sum(node_rows$total_gain) / sum(node_rows$node_site_score_entries),
      max_gain = if (is.null(node_rows)) NA_real_ else max(node_rows$max_gain),
      nodes_covered = if (is.null(node_rows)) 0L else length(unique(node_rows$node_id)),
      deep_nodes_covered = count_deep_nodes(node_rows, deep_cutoff)
    )
    cbind(out, window_missingness_summary(aln_int, windows$start[[i]], windows$end[[i]]))
  })

  out <- do.call(rbind, rows)
  out <- out[out$informative_snps >= min_informative_snps, ]
  out[order(-out$total_weighted_gain, -out$total_gain, out$start), ]
}

score_candidate_windows <- function(windows, polymorphic_sites, site_node_scores, aln_int = NULL,
                                    min_informative_snps = 1L) {
  informative_assignments <- assign_snps_to_windows(unique(site_node_scores$site), windows)
  window_node_summary <- aggregate_window_node_scores(windows, site_node_scores, informative_assignments)
  window_summary <- summarise_windows(
    windows = windows,
    polymorphic_sites = polymorphic_sites,
    site_node_scores = site_node_scores,
    window_node_summary = window_node_summary,
    aln_int = aln_int,
    min_informative_snps = min_informative_snps
  )

  list(
    candidate_windows = windows[windows$window_id %in% window_summary$window_id, ],
    window_summary = window_summary,
    window_node_summary = window_node_summary[window_node_summary$window_id %in% window_summary$window_id, ]
  )
}

window_missingness_summary <- function(aln_int, start, end) {
  if (is.null(aln_int)) {
    return(data.frame(missing_sites_mean = NA_real_, missing_cells = NA_integer_, missing_fraction = NA_real_))
  }
  slice <- aln_int[, start:end, drop = FALSE]
  missing_by_site <- colMeans(slice == 0L)
  data.frame(
    missing_sites_mean = mean(missing_by_site),
    missing_cells = sum(slice == 0L),
    missing_fraction = mean(slice == 0L)
  )
}

count_deep_nodes <- function(node_rows, deep_cutoff) {
  if (is.null(node_rows) || is.na(deep_cutoff) || !"depth" %in% names(node_rows)) {
    return(NA_integer_)
  }
  length(unique(node_rows$node_id[!is.na(node_rows$depth) & node_rows$depth >= deep_cutoff]))
}

add_window_ids <- function(windows) {
  if (nrow(windows) == 0L) {
    return(empty_windows())
  }
  windows <- windows[order(windows$window_type, windows$requested_width, windows$start, windows$end, windows$centre_site), ]
  windows$window_id <- sprintf("win_%05d", seq_len(nrow(windows)))
  windows[c("window_id", "window_type", "start", "end", "width", "requested_width", "centre_site")]
}

empty_windows <- function() {
  data.frame(
    window_id = character(), window_type = character(), start = integer(),
    end = integer(), width = integer(), requested_width = integer(),
    centre_site = integer()
  )
}

empty_window_node_summary <- function(site_node_scores) {
  cols <- c("window_id", "node_index", "node_id", "informative_snps",
            "node_site_score_entries", "total_gain", "total_weighted_gain",
            "mean_gain", "max_gain")
  out <- data.frame(
    window_id = character(), node_index = integer(), node_id = integer(),
    informative_snps = integer(), node_site_score_entries = integer(),
    total_gain = numeric(), total_weighted_gain = numeric(),
    mean_gain = numeric(), max_gain = numeric()
  )
  metadata_cols <- setdiff(names(site_node_scores), c(
    "site", "gain", "normalized_gain", "best_allele", "direction",
    "allele_count_inside", "allele_count_outside", "observed_inside",
    "observed_outside", "total_observed", cols
  ))
  for (col in metadata_cols) {
    out[[col]] <- site_node_scores[[col]][0]
  }
  out
}
