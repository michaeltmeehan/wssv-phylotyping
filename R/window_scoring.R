# New stage-04 candidate amplicon helpers.
#
# These helpers build candidate amplicons from informative SNP geography using
# the target-specific stage-02/03 SNP scores. The legacy fixed-window helpers
# remain below for backward compatibility, but the stage-04 script now calls the
# amplicon-specific functions in this section.

resolve_amplicon_config <- function(config) {
  amp <- value_or(config$analysis$amplicons, list())
  min_length <- as.integer(value_or(amp$min_length, 500L))
  max_length <- as.integer(value_or(amp$max_length, 3000L))
  flank_search_width <- as.integer(value_or(amp$flank_search_width, 100L))
  strong_gain_threshold <- as.numeric(value_or(amp$strong_gain_threshold, 0.8))
  minimum_gain_threshold <- as.numeric(value_or(amp$minimum_gain_threshold, 0.5))
  flank_min_non_missing_fraction <- as.numeric(value_or(amp$flank_min_non_missing_fraction, 0.9))
  flank_max_variable_fraction <- as.numeric(value_or(amp$flank_max_variable_fraction, 0.25))
  flank_min_gc_fraction <- as.numeric(value_or(amp$flank_min_gc_fraction, 0.35))
  flank_min_conserved_fraction <- as.numeric(value_or(amp$flank_min_conserved_fraction, 0.9))
  flank_complete_fraction <- as.numeric(value_or(amp$flank_complete_fraction, 0.95))

  if (is.na(min_length) || min_length < 1L) {
    stop("analysis.amplicons.min_length must be a positive integer.", call. = FALSE)
  }
  if (is.na(max_length) || max_length < min_length) {
    stop("analysis.amplicons.max_length must be >= min_length.", call. = FALSE)
  }
  if (is.na(flank_search_width) || flank_search_width < 0L) {
    stop("analysis.amplicons.flank_search_width must be a non-negative integer.", call. = FALSE)
  }

  list(
    min_length = min_length,
    max_length = max_length,
    flank_search_width = flank_search_width,
    strong_gain_threshold = strong_gain_threshold,
    minimum_gain_threshold = minimum_gain_threshold,
    flank_min_non_missing_fraction = flank_min_non_missing_fraction,
    flank_max_variable_fraction = flank_max_variable_fraction,
    flank_min_gc_fraction = flank_min_gc_fraction,
    flank_min_conserved_fraction = flank_min_conserved_fraction,
    flank_complete_fraction = flank_complete_fraction
  )
}

resolve_primary_target_nodes <- function(config) {
  amp <- value_or(config$analysis$amplicons, list())
  primary <- unique(as.integer(value_or(amp$primary_target_nodes, integer())))
  primary <- primary[!is.na(primary)]
  if (length(primary) == 0L) {
    stop("analysis.amplicons.primary_target_nodes must list the primary assay targets.", call. = FALSE)
  }
  sort(primary)
}

filter_primary_target_scores <- function(site_node_scores, primary_target_nodes) {
  if (nrow(site_node_scores) == 0L) {
    return(site_node_scores[0, , drop = FALSE])
  }
  if (!"node_id" %in% names(site_node_scores)) {
    stop("site_node_scores must contain node_id.", call. = FALSE)
  }
  out <- site_node_scores[site_node_scores$node_id %in% primary_target_nodes, , drop = FALSE]
  out[order(out$site, out$node_id, out$gain, out$normalized_gain), , drop = FALSE]
}

site_position_lookup <- function(site_node_scores, site_map = NULL) {
  sites <- sort(unique(as.integer(site_node_scores$site)))
  if (length(sites) == 0L) {
    return(data.frame(site = integer(), alignment_position = integer(), stringsAsFactors = FALSE))
  }

  if (is.null(site_map) || !all(c("site", "alignment_position") %in% names(site_map))) {
    return(data.frame(site = sites, alignment_position = sites, stringsAsFactors = FALSE))
  }

  lookup <- unique(site_map[, c("site", "alignment_position"), drop = FALSE])
  lookup <- lookup[match(sites, lookup$site), , drop = FALSE]
  if (anyNA(lookup$alignment_position)) {
    lookup$alignment_position[is.na(lookup$alignment_position)] <- lookup$site[is.na(lookup$alignment_position)]
  }
  lookup[order(lookup$alignment_position, lookup$site), , drop = FALSE]
}

expand_interval_to_length <- function(core_start, core_end, target_length, alignment_length) {
  core_length <- as.integer(core_end - core_start + 1L)
  if (target_length < core_length) {
    stop("target_length must be at least the core interval length.", call. = FALSE)
  }
  if (target_length > alignment_length) {
    return(NULL)
  }

  start <- as.integer(core_start - floor((target_length - core_length) / 2))
  end <- start + target_length - 1L
  if (start < 1L) {
    end <- end + (1L - start)
    start <- 1L
  }
  if (end > alignment_length) {
    start <- start - (end - alignment_length)
    end <- alignment_length
  }
  if (start < 1L || end > alignment_length) {
    return(NULL)
  }
  if (end - start + 1L != target_length) {
    return(NULL)
  }

  c(start = as.integer(start), end = as.integer(end))
}

enumerate_candidate_cores <- function(site_positions, max_length) {
  site_positions <- sort(unique(as.integer(site_positions)))
  site_positions <- site_positions[!is.na(site_positions)]
  if (length(site_positions) == 0L) {
    return(data.frame(core_start = integer(), core_end = integer(), stringsAsFactors = FALSE))
  }

  rows <- list()
  index <- 1L
  for (i in seq_along(site_positions)) {
    if (index < i) {
      index <- i
    }
    while (index <= length(site_positions) && site_positions[[index]] - site_positions[[i]] + 1L <= max_length) {
      index <- index + 1L
    }
    if (index - i < 1L) {
      next
    }
    for (j in seq.int(i, index - 1L)) {
      rows[[length(rows) + 1L]] <- data.frame(
        core_start = as.integer(site_positions[[i]]),
        core_end = as.integer(site_positions[[j]]),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

deduplicate_candidate_intervals <- function(candidates) {
  if (nrow(candidates) == 0L) {
    return(candidates)
  }
  candidates$key <- as.character(candidates$covered_sites)
  candidates <- candidates[order(candidates$key, candidates$length, candidates$start, candidates$end), , drop = FALSE]
  candidates <- candidates[!duplicated(candidates$key), , drop = FALSE]
  candidates$key <- NULL
  candidates
}

generate_candidate_amplicons <- function(site_node_scores, site_map = NULL, alignment_length, config) {
  amp_cfg <- resolve_amplicon_config(config)
  primary_nodes <- resolve_primary_target_nodes(config)
  primary_scores <- filter_primary_target_scores(site_node_scores, primary_nodes)
  if (nrow(primary_scores) == 0L) {
    return(data.frame(
      candidate_id = character(), start = integer(), end = integer(), length = integer(),
      core_start = integer(), core_end = integer(), core_length = integer(),
      covered_sites = integer(), covered_site_count = integer(),
      stringsAsFactors = FALSE
    ))
  }

  positions <- site_position_lookup(primary_scores, site_map)
  candidate_cores <- enumerate_candidate_cores(positions$alignment_position, amp_cfg$max_length)
  if (nrow(candidate_cores) == 0L) {
    return(data.frame(
      candidate_id = character(), start = integer(), end = integer(), length = integer(),
      core_start = integer(), core_end = integer(), core_length = integer(),
      covered_sites = integer(), covered_site_count = integer(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- vector("list", nrow(candidate_cores))
  kept <- 0L
  for (i in seq_len(nrow(candidate_cores))) {
    core_start <- candidate_cores$core_start[[i]]
    core_end <- candidate_cores$core_end[[i]]
    core_length <- as.integer(core_end - core_start + 1L)
    target_length <- max(core_length, amp_cfg$min_length)
    interval <- expand_interval_to_length(core_start, core_end, target_length, alignment_length)
    if (is.null(interval)) {
      next
    }
    covered_sites <- positions$site[positions$alignment_position >= interval[[1L]] & positions$alignment_position <= interval[[2L]]]
    if (length(covered_sites) == 0L) {
      next
    }
    kept <- kept + 1L
    rows[[kept]] <- data.frame(
      start = as.integer(interval[[1L]]),
      end = as.integer(interval[[2L]]),
      length = as.integer(interval[[2L]] - interval[[1L]] + 1L),
      core_start = as.integer(core_start),
      core_end = as.integer(core_end),
      core_length = as.integer(core_length),
      covered_sites = I(list(as.integer(sort(unique(covered_sites))))),
      covered_site_count = as.integer(length(unique(covered_sites))),
      stringsAsFactors = FALSE
    )
  }

  candidates <- if (kept == 0L) {
    data.frame(
      start = integer(), end = integer(), length = integer(),
      core_start = integer(), core_end = integer(), core_length = integer(),
      covered_sites = I(list()), covered_site_count = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, rows[seq_len(kept)])
  }

  if (nrow(candidates) == 0L) {
    return(data.frame(
      candidate_id = character(), start = integer(), end = integer(), length = integer(),
      core_start = integer(), core_end = integer(), core_length = integer(),
      covered_sites = integer(), covered_site_count = integer(),
      stringsAsFactors = FALSE
    ))
  }

  candidates$covered_sites <- vapply(candidates$covered_sites, function(x) {
    paste(as.integer(x), collapse = ",")
  }, character(1L))
  candidates <- candidates[!duplicated(candidates$covered_sites), , drop = FALSE]
  candidates <- candidates[order(candidates$length, candidates$start, candidates$core_start, candidates$core_end), , drop = FALSE]
  candidates$candidate_id <- sprintf("amplicon_%03d", seq_len(nrow(candidates)))
  candidates[ , c("candidate_id", "start", "end", "length", "core_start", "core_end", "core_length", "covered_sites", "covered_site_count")]
}

summarise_alignment_region <- function(aln_int, start, end, polymorphic_sites = NULL, complete_fraction = 0.95) {
  if (is.null(aln_int)) {
    return(data.frame(
      region_length = as.integer(end - start + 1L),
      non_missing_fraction = NA_real_,
      missing_fraction = NA_real_,
      missing_cells = NA_integer_,
      complete_sequences = NA_integer_,
      complete_fraction = NA_real_,
      near_complete_sequences = NA_integer_,
      near_complete_fraction = NA_real_,
      variable_sites = NA_integer_,
      variable_site_fraction = NA_real_,
      gc_fraction = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  slice <- aln_int[, start:end, drop = FALSE]
  non_missing <- slice != 0L
  observed_cells <- sum(non_missing)
  total_cells <- length(slice)
  row_completeness <- rowSums(non_missing) / ncol(slice)
  variable_sites <- if (is.null(polymorphic_sites)) {
    NA_integer_
  } else {
    sum(polymorphic_sites >= start & polymorphic_sites <= end)
  }
  observed_gc <- sum(slice %in% c(2L, 3L))
  out <- data.frame(
    region_length = as.integer(end - start + 1L),
    non_missing_fraction = if (total_cells == 0L) NA_real_ else observed_cells / total_cells,
    missing_fraction = if (total_cells == 0L) NA_real_ else 1 - (observed_cells / total_cells),
    missing_cells = as.integer(total_cells - observed_cells),
    complete_sequences = as.integer(sum(row_completeness == 1)),
    complete_fraction = if (nrow(slice) == 0L) NA_real_ else mean(row_completeness == 1),
    near_complete_sequences = as.integer(sum(row_completeness >= complete_fraction)),
    near_complete_fraction = if (nrow(slice) == 0L) NA_real_ else mean(row_completeness >= complete_fraction),
    variable_sites = variable_sites,
    variable_site_fraction = if (is.na(variable_sites)) NA_real_ else variable_sites / as.integer(end - start + 1L),
    gc_fraction = if (observed_cells == 0L) NA_real_ else observed_gc / observed_cells,
    stringsAsFactors = FALSE
  )
  out
}

summarise_flank_region <- function(aln_int, start, end, polymorphic_sites = NULL, config = NULL) {
  if (start > end) {
    return(data.frame(
      length = 0L,
      non_missing_fraction = NA_real_,
      missing_fraction = NA_real_,
      variable_sites = NA_integer_,
      variable_site_fraction = NA_real_,
      gc_fraction = NA_real_,
      conserved_stretch = NA,
      stringsAsFactors = FALSE
    ))
  }
  amp_cfg <- if (is.null(config)) resolve_amplicon_config(list(analysis = list(amplicons = list()))) else resolve_amplicon_config(config)
  region <- summarise_alignment_region(aln_int, start, end, polymorphic_sites = polymorphic_sites, complete_fraction = amp_cfg$flank_complete_fraction)
  data.frame(
    length = region$region_length,
    non_missing_fraction = region$non_missing_fraction,
    missing_fraction = region$missing_fraction,
    variable_sites = region$variable_sites,
    variable_site_fraction = region$variable_site_fraction,
    gc_fraction = region$gc_fraction,
    conserved_stretch = isTRUE(!is.na(region$non_missing_fraction) &&
      region$non_missing_fraction >= amp_cfg$flank_min_non_missing_fraction &&
      !is.na(region$complete_fraction) &&
      region$complete_fraction >= amp_cfg$flank_min_conserved_fraction &&
      !is.na(region$variable_site_fraction) &&
      region$variable_site_fraction <= amp_cfg$flank_max_variable_fraction &&
      !is.na(region$gc_fraction) &&
      region$gc_fraction >= amp_cfg$flank_min_gc_fraction),
    stringsAsFactors = FALSE
  )
}

score_candidate_amplicons <- function(site_node_scores, polymorphic_sites, aln_int, site_map = NULL, config) {
  amp_cfg <- resolve_amplicon_config(config)
  primary_nodes <- resolve_primary_target_nodes(config)
  site_node_scores <- filter_primary_target_scores(site_node_scores, primary_nodes)
  alignment_length <- ncol(aln_int)
  candidates <- generate_candidate_amplicons(
    site_node_scores = site_node_scores,
    site_map = site_map,
    alignment_length = alignment_length,
    config = config
  )
  if (nrow(candidates) == 0L) {
    empty_candidate_summary <- data.frame(
      candidate_id = character(), start = integer(), end = integer(), length = integer(),
      core_start = integer(), core_end = integer(), core_length = integer(),
      informative_snp_associations = integer(), unique_informative_snp_positions = integer(),
      primary_targets_represented = integer(), target_nodes_represented = character(),
      max_normalized_gain = numeric(), targets_with_any_signal = integer(),
      targets_with_strong_signal = integer(),
      region_length = integer(), non_missing_fraction = numeric(), missing_fraction = numeric(),
      missing_cells = integer(), complete_sequences = integer(), complete_fraction = numeric(),
      near_complete_sequences = integer(), near_complete_fraction = numeric(),
      variable_sites = integer(), variable_site_fraction = numeric(), gc_fraction = numeric(),
      left_flank_length = integer(), left_flank_non_missing_fraction = numeric(),
      left_flank_missing_fraction = numeric(), left_flank_variable_sites = integer(),
      left_flank_variable_site_fraction = numeric(), left_flank_gc_fraction = numeric(),
      left_flank_conserved_stretch = logical(),
      right_flank_length = integer(), right_flank_non_missing_fraction = numeric(),
      right_flank_missing_fraction = numeric(), right_flank_variable_sites = integer(),
      right_flank_variable_site_fraction = numeric(), right_flank_gc_fraction = numeric(),
      right_flank_conserved_stretch = logical(),
      stringsAsFactors = FALSE
    )
    return(list(
      candidate_amplicons = empty_candidate_summary,
      candidate_amplicon_target_scores = data.frame(),
      candidate_amplicon_snps = data.frame()
    ))
  }

  candidate_site_pool <- sort(unique(as.integer(unlist(strsplit(candidates$covered_sites, ",", fixed = TRUE)))))
  site_scores <- site_node_scores[site_node_scores$site %in% candidate_site_pool, , drop = FALSE]
  if (nrow(site_scores) == 0L) {
    return(list(
      candidate_amplicons = candidates[0, , drop = FALSE],
      candidate_amplicon_target_scores = data.frame(),
      candidate_amplicon_snps = data.frame()
    ))
  }

  candidate_snp_rows <- vector("list", nrow(candidates))
  candidate_target_rows <- vector("list", nrow(candidates))
  summary_rows <- vector("list", nrow(candidates))

  for (i in seq_len(nrow(candidates))) {
    candidate_id <- candidates$candidate_id[[i]]
    start <- candidates$start[[i]]
    end <- candidates$end[[i]]
    candidate_sites <- unique(site_scores$site[site_scores$site >= start & site_scores$site <= end])
    candidate_rows <- site_scores[site_scores$site %in% candidate_sites, , drop = FALSE]
    if (nrow(candidate_rows) == 0L) {
      next
    }

    candidate_snp_rows[[i]] <- within(candidate_rows, {
      candidate_id <- candidate_id
      start <- start
      end <- end
      length <- as.integer(end - start + 1L)
      raw_gain <- gain
      best_snp_rule <- best_allele
      target_node <- node_id
    })[c(
      "candidate_id", "start", "end", "length", "target_node", "node_index", "node_id", "site",
      "raw_gain", "normalized_gain", "best_snp_rule", "direction", "allele_count_inside",
      "allele_count_outside", "observed_inside", "observed_outside", "total_observed"
    )]

    grouped <- split(candidate_rows, candidate_rows$node_id)
    target_rows <- lapply(grouped, function(x) {
      best_index <- order(-x$normalized_gain, -x$gain, x$site)[1L]
      best_row <- x[best_index, , drop = FALSE]
      data.frame(
        candidate_id = candidate_id,
        start = start,
        end = end,
        length = as.integer(end - start + 1L),
        target_node = best_row$node_id[[1L]],
        n_informative_snps = length(unique(x$site)),
        max_raw_gain = max(x$gain),
        max_normalized_gain = max(x$normalized_gain),
        mean_normalized_gain = mean(x$normalized_gain),
        median_normalized_gain = stats::median(x$normalized_gain),
        best_snp_position = best_row$site[[1L]],
        best_snp_rule = best_row$best_allele[[1L]],
        stringsAsFactors = FALSE
      )
    })
    candidate_target_rows[[i]] <- do.call(rbind, target_rows)

    target_signal <- candidate_target_rows[[i]]
    summary_rows[[i]] <- data.frame(
      candidate_id = candidate_id,
      start = start,
      end = end,
      length = as.integer(end - start + 1L),
      core_start = candidates$core_start[[i]],
      core_end = candidates$core_end[[i]],
      core_length = candidates$core_length[[i]],
      informative_snp_associations = nrow(candidate_rows),
      unique_informative_snp_positions = length(unique(candidate_rows$site)),
      primary_targets_represented = length(unique(candidate_rows$node_id)),
      target_nodes_represented = paste(sort(unique(candidate_rows$node_id)), collapse = ","),
      max_normalized_gain = max(candidate_rows$normalized_gain),
      targets_with_any_signal = sum(tapply(target_signal$max_normalized_gain >= amp_cfg$minimum_gain_threshold, target_signal$target_node, any)),
      targets_with_strong_signal = sum(tapply(target_signal$max_normalized_gain >= amp_cfg$strong_gain_threshold, target_signal$target_node, any)),
      stringsAsFactors = FALSE
    )
    summary_rows[[i]] <- cbind(
      summary_rows[[i]],
      summarise_alignment_region(aln_int, start, end, polymorphic_sites = polymorphic_sites, complete_fraction = amp_cfg$flank_complete_fraction)
    )
    left_start <- max(1L, start - amp_cfg$flank_search_width)
    left_end <- start - 1L
    right_start <- end + 1L
    right_end <- min(ncol(aln_int), end + amp_cfg$flank_search_width)
    summary_rows[[i]] <- cbind(
      summary_rows[[i]],
      {
        left_flank <- summarise_flank_region(aln_int, left_start, left_end, polymorphic_sites = polymorphic_sites, config = config)
        setNames(left_flank, paste0("left_flank_", names(left_flank)))
      },
      {
        right_flank <- summarise_flank_region(aln_int, right_start, right_end, polymorphic_sites = polymorphic_sites, config = config)
        setNames(right_flank, paste0("right_flank_", names(right_flank)))
      }
    )
  }

  candidate_amplicon_snps <- do.call(rbind, candidate_snp_rows)
  candidate_amplicon_target_scores <- do.call(rbind, candidate_target_rows)
  candidate_amplicons <- do.call(rbind, summary_rows)

  candidate_amplicons <- candidate_amplicons[order(
    -candidate_amplicons$primary_targets_represented,
    -candidate_amplicons$targets_with_strong_signal,
    -candidate_amplicons$complete_fraction,
    -candidate_amplicons$unique_informative_snp_positions,
    candidate_amplicons$start,
    candidate_amplicons$end
  ), , drop = FALSE]
  candidate_amplicons$candidate_id <- sprintf("amplicon_%03d", seq_len(nrow(candidate_amplicons)))
  candidate_amplicons$start <- as.integer(candidate_amplicons$start)
  candidate_amplicons$end <- as.integer(candidate_amplicons$end)
  candidate_amplicons$length <- as.integer(candidate_amplicons$length)

  candidate_amplicon_target_scores$candidate_id <- NULL
  candidate_amplicon_target_scores <- merge(
    candidate_amplicons[c("start", "end", "length", "candidate_id")],
    candidate_amplicon_target_scores,
    by = c("start", "end", "length"),
    sort = FALSE
  )
  candidate_amplicon_target_scores <- candidate_amplicon_target_scores[order(
    candidate_amplicon_target_scores$candidate_id,
    candidate_amplicon_target_scores$target_node
  ), , drop = FALSE]

  candidate_amplicon_snps$candidate_id <- NULL
  candidate_amplicon_snps <- merge(
    candidate_amplicons[c("start", "end", "length", "candidate_id")],
    candidate_amplicon_snps,
    by = c("start", "end", "length"),
    sort = FALSE
  )
  candidate_amplicon_snps <- candidate_amplicon_snps[order(
    candidate_amplicon_snps$candidate_id,
    candidate_amplicon_snps$target_node,
    candidate_amplicon_snps$site
  ), , drop = FALSE]

  list(
    candidate_amplicons = candidate_amplicons,
    candidate_amplicon_target_scores = candidate_amplicon_target_scores,
    candidate_amplicon_snps = candidate_amplicon_snps
  )
}

plot_candidate_amplicons_spatial <- function(candidate_amplicons, candidate_amplicon_snps, primary_target_nodes, output_path, alignment_length) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(output_path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  y_levels <- seq_along(primary_target_nodes)
  names(y_levels) <- as.character(primary_target_nodes)
  colors <- stats::setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"), as.character(primary_target_nodes))

  graphics::plot(NA,
    xlim = c(1, alignment_length),
    ylim = c(0.5, length(primary_target_nodes) + 0.5),
    xlab = "Alignment position",
    ylab = "Primary target clade",
    yaxt = "n",
    main = "Informative SNP geography and candidate amplicons"
  )
  graphics::axis(2, at = y_levels, labels = primary_target_nodes, las = 1)

  if (nrow(candidate_amplicons) > 0L) {
    for (i in seq_len(nrow(candidate_amplicons))) {
      graphics::rect(
        xleft = candidate_amplicons$start[[i]],
        ybottom = 0.5,
        xright = candidate_amplicons$end[[i]],
        ytop = length(primary_target_nodes) + 0.5,
        col = grDevices::adjustcolor("#999999", alpha.f = 0.12),
        border = NA
      )
    }
  }

  if (nrow(candidate_amplicon_snps) > 0L) {
    for (node in primary_target_nodes) {
      rows <- candidate_amplicon_snps[candidate_amplicon_snps$target_node == node, , drop = FALSE]
      if (nrow(rows) == 0L) {
        next
      }
      graphics::points(rows$site, rep(y_levels[[as.character(node)]], nrow(rows)),
        pch = 16, cex = 0.7, col = colors[[as.character(node)]]
      )
    }
  }
  graphics::legend("topright", legend = as.character(primary_target_nodes), col = colors, pch = 16, bty = "n", title = "Target nodes")
}

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
