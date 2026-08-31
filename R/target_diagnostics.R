#' Build compact diagnostics for the configured target clades.
#'
#' The checkpoint tables below are intentionally descriptive rather than scored
#' summaries. They help decide whether the configured target-clade simplification
#' is strong enough to justify proceeding to window and primer-panel design.

configured_target_name <- function(target_spec) {
  if (!is.list(target_spec)) {
    return(NA_character_)
  }

  for (field in c("target_name", "name", "label", "clade_name")) {
    value <- target_spec[[field]]
    if (!is.null(value)) {
      value <- as.character(value)[1L]
      if (!is.na(value) && nzchar(value)) {
        return(value)
      }
    }
  }

  NA_character_
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

target_display_label <- function(target_name, node_id) {
  if (!is.na(target_name) && nzchar(target_name)) {
    target_name
  } else {
    paste0("node ", node_id)
  }
}

site_to_alignment_position <- function(site, site_map = NULL) {
  if (is.null(site_map) || !all(c("site", "alignment_position") %in% names(site_map))) {
    return(site)
  }

  mapped <- site_map$alignment_position[match(site, site_map$site)]
  ifelse(is.na(mapped), site, mapped)
}

make_target_clade_checkpoint <- function(pre, site_node_scores_targets, top_n = 3L) {
  if (is.null(pre$target_node_metadata)) {
    stop("pre must include target_node_metadata.", call. = FALSE)
  }
  if (is.null(pre$aln_int)) {
    stop("pre must include aln_int.", call. = FALSE)
  }
  if (is.null(site_node_scores_targets)) {
    site_node_scores_targets <- data.frame()
  }

  target_meta <- pre$target_node_metadata
  n_tip <- nrow(pre$aln_int)
  n_polymorphic_sites <- length(pre$polymorphic_sites)
  target_specs <- pre$target_clades$target_clade_spec %||% vector("list", nrow(target_meta))
  if (length(target_specs) < nrow(target_meta)) {
    target_specs <- c(target_specs, vector("list", nrow(target_meta) - length(target_specs)))
  }

  target_meta$target_name <- vapply(
    seq_len(nrow(target_meta)),
    function(i) configured_target_name(target_specs[[i]]),
    character(1L)
  )
  target_meta$target <- vapply(
    seq_len(nrow(target_meta)),
    function(i) target_display_label(target_meta$target_name[[i]], target_meta$node_id[[i]]),
    character(1L)
  )

  summary <- target_meta[, c(
    "target",
    "target_name",
    "node_id",
    "posterior_support",
    "clade_size",
    "complement_size"
  ), drop = FALSE]
  names(summary)[names(summary) == "clade_size"] <- "descendant_sequences"
  names(summary)[names(summary) == "complement_size"] <- "non_descendant_sequences"
  summary$polymorphic_sites_assessed <- n_polymorphic_sites
  summary$informative_sites_passing_threshold <- 0L
  summary$max_normalized_gain <- NA_real_
  summary$mean_normalized_gain <- NA_real_
  summary$best_site <- NA_integer_
  summary$best_alignment_position <- NA_integer_
  summary$best_site_observed_sequences <- NA_integer_
  summary$best_site_missing_sequences <- NA_integer_
  summary$best_site_missing_fraction <- NA_real_

  if (nrow(site_node_scores_targets) > 0L) {
    node_summary <- make_node_summary(site_node_scores_targets, target_meta)
    node_idx <- match(node_summary$node_id, summary$node_id)
    summary$informative_sites_passing_threshold[node_idx] <- node_summary$informative_sites
    summary$max_normalized_gain[node_idx] <- node_summary$normalized_gain_max
    summary$mean_normalized_gain[node_idx] <- node_summary$normalized_gain_mean
    summary$best_site[node_idx] <- node_summary$best_site
    summary$best_alignment_position[node_idx] <- site_to_alignment_position(node_summary$best_site, pre$site_map)

    best_row_index <- match(
      paste(node_summary$node_id, node_summary$best_site, sep = "::"),
      paste(site_node_scores_targets$node_id, site_node_scores_targets$site, sep = "::")
    )
    best_rows <- site_node_scores_targets[best_row_index, , drop = FALSE]
    best_node_idx <- match(best_rows$node_id, summary$node_id)
    summary$best_site_observed_sequences[best_node_idx] <- best_rows$total_observed
    summary$best_site_missing_sequences[best_node_idx] <- n_tip - best_rows$total_observed
    summary$best_site_missing_fraction[best_node_idx] <- summary$best_site_missing_sequences[best_node_idx] / n_tip
  }

  summary <- summary[match(target_meta$node_id, summary$node_id), , drop = FALSE]
  summary <- summary[, c(
    "target",
    "target_name",
    "node_id",
    "posterior_support",
    "descendant_sequences",
    "non_descendant_sequences",
    "polymorphic_sites_assessed",
    "informative_sites_passing_threshold",
    "max_normalized_gain",
    "mean_normalized_gain",
    "best_site",
    "best_alignment_position",
    "best_site_observed_sequences",
    "best_site_missing_sequences",
    "best_site_missing_fraction"
  ), drop = FALSE]

  strongest <- data.frame()
  if (nrow(site_node_scores_targets) > 0L) {
    strongest <- site_node_scores_targets
    strongest$target_name <- target_meta$target_name[match(strongest$node_id, target_meta$node_id)]
    strongest$target <- target_meta$target[match(strongest$node_id, target_meta$node_id)]
    strongest$target_order <- match(strongest$node_id, target_meta$node_id)
    strongest$alignment_position <- site_to_alignment_position(strongest$site, pre$site_map)
    strongest$baseline_accuracy <- pmax(strongest$observed_inside, strongest$observed_outside) / strongest$total_observed
    strongest$accuracy <- strongest$baseline_accuracy + strongest$gain
    strongest$best_rule <- paste0(strongest$best_allele, " / ", strongest$direction)
    strongest$missing_sequences <- n_tip - strongest$total_observed
    strongest$missing_fraction <- strongest$missing_sequences / n_tip
    strongest <- strongest[order(
      strongest$target_order,
      -strongest$normalized_gain,
      -strongest$gain,
      strongest$site
    ), , drop = FALSE]
    strongest$rank_within_target <- ave(seq_len(nrow(strongest)), strongest$node_id, FUN = seq_along)
    strongest <- strongest[strongest$rank_within_target <= top_n, , drop = FALSE]
    strongest <- strongest[, c(
      "target",
      "target_name",
      "node_id",
      "rank_within_target",
      "site",
      "alignment_position",
      "best_rule",
      "best_allele",
      "direction",
      "observed_inside",
      "observed_outside",
      "accuracy",
      "baseline_accuracy",
      "gain",
      "normalized_gain",
      "missing_sequences",
      "missing_fraction"
    ), drop = FALSE]
    names(strongest)[names(strongest) == "direction"] <- "best_direction"
    names(strongest)[names(strongest) == "gain"] <- "raw_gain"
  } else {
    strongest <- data.frame(
      target = character(),
      target_name = character(),
      node_id = integer(),
      rank_within_target = integer(),
      site = integer(),
      alignment_position = integer(),
      best_rule = character(),
      best_allele = character(),
      best_direction = character(),
      observed_inside = integer(),
      observed_outside = integer(),
      accuracy = numeric(),
      baseline_accuracy = numeric(),
      raw_gain = numeric(),
      normalized_gain = numeric(),
      missing_sequences = integer(),
      missing_fraction = numeric(),
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = summary,
    strongest_snps = strongest
  )
}
