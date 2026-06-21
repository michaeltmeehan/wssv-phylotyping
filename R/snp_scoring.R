#' Find the best allele rule for one alignment column across eligible nodes.
#'
#' Used by SNP scoring to ask whether a nucleotide state is enriched inside or
#' outside each target clade after applying observation and minor-allele filters.
best_rule_for_site <- function(x, target_mask, min_total_obs, min_side_obs, min_site_maf) {
  observed <- x != 0L
  allele_totals <- tabulate(x[observed], nbins = 4L)

  obs_in <- as.vector(target_mask %*% as.integer(observed))
  obs_out <- as.vector((!target_mask) %*% as.integer(observed))
  total_obs <- obs_in + obs_out
  ok_node <- total_obs >= min_total_obs & obs_in >= min_side_obs & obs_out >= min_side_obs

  best_acc <- rep(NA_real_, nrow(target_mask))
  best_allele <- rep(NA_integer_, nrow(target_mask))
  best_dir <- rep(NA_integer_, nrow(target_mask))
  count_in <- rep(NA_integer_, nrow(target_mask))
  count_out <- rep(NA_integer_, nrow(target_mask))
  best_acc[ok_node] <- -Inf

  for (allele in 1:4) {
    if (allele_totals[[allele]] < min_site_maf) {
      next
    }

    has_allele <- x == allele
    allele_in <- as.vector(target_mask %*% as.integer(has_allele))
    allele_out <- as.vector((!target_mask) %*% as.integer(has_allele))

    acc_marks_clade <- (allele_in + (obs_out - allele_out)) / total_obs
    acc_marks_outside <- (allele_out + (obs_in - allele_in)) / total_obs
    acc <- pmax(acc_marks_clade, acc_marks_outside)
    direction <- ifelse(acc_marks_clade >= acc_marks_outside, 1L, 2L)

    improve <- ok_node & acc > best_acc
    if (any(improve)) {
      best_acc[improve] <- acc[improve]
      best_allele[improve] <- allele
      best_dir[improve] <- direction[improve]
      count_in[improve] <- allele_in[improve]
      count_out[improve] <- allele_out[improve]
    }
  }

  bad <- ok_node & !is.finite(best_acc)
  best_acc[bad] <- NA_real_

  list(
    best_acc = best_acc,
    best_allele = best_allele,
    best_dir = best_dir,
    allele_count_in = count_in,
    allele_count_out = count_out,
    obs_in = obs_in,
    obs_out = obs_out,
    total_obs = total_obs
  )
}

#' Score one polymorphic site against all eligible clades.
#'
#' Returns one row for each node where the best allele rule improves over the
#' clade-size baseline by at least the configured normalized-gain threshold.
score_site <- function(site, x, target_mask, node_metadata, min_total_obs, min_side_obs,
                       min_site_maf, min_gain_norm) {
  rule <- best_rule_for_site(x, target_mask, min_total_obs, min_side_obs, min_site_maf)
  ok <- !is.na(rule$best_acc)
  if (!any(ok)) {
    return(data.frame())
  }

  baseline <- pmax(rule$obs_in, rule$obs_out) / rule$total_obs
  gain <- rule$best_acc - baseline
  max_gain <- 1 - baseline
  gain_norm <- gain / max_gain
  gain_norm[!is.finite(gain_norm)] <- NA_real_

  helped <- ok & gain >= (1 / rule$total_obs) & gain_norm >= min_gain_norm & is.finite(gain_norm)
  if (!any(helped)) {
    return(data.frame())
  }

  node_index <- which(helped)
  out <- data.frame(
    node_index = node_index,
    node_id = node_metadata$node_id[node_index],
    site = site,
    gain = gain[node_index],
    normalized_gain = gain_norm[node_index],
    best_allele = base_code_to_allele(rule$best_allele[node_index]),
    direction = ifelse(rule$best_dir[node_index] == 1L, "clade", "outside"),
    allele_count_inside = rule$allele_count_in[node_index],
    allele_count_outside = rule$allele_count_out[node_index],
    observed_inside = rule$obs_in[node_index],
    observed_outside = rule$obs_out[node_index],
    total_observed = rule$total_obs[node_index]
  )

  merge(out, node_metadata, by = c("node_index", "node_id"), sort = FALSE)
}

#' Score many polymorphic SNP sites.
#'
#' This is the main stage-02 worker. It loops over candidate sites and combines
#' all informative site-node rules used by downstream window scoring and the
#' classifier.
score_snp_sites <- function(aln_int, target_mask, node_metadata, sites,
                            min_total_obs, min_side_obs, min_site_maf,
                            min_gain_norm = 0.5, chunk_size = 5000L,
                            progress = interactive()) {
  scores <- vector("list", length(sites))
  names(scores) <- as.character(sites)

  for (i in seq_along(sites)) {
    site <- sites[[i]]
    scores[[i]] <- score_site(
      site = site,
      x = aln_int[, site],
      target_mask = target_mask,
      node_metadata = node_metadata,
      min_total_obs = min_total_obs,
      min_side_obs = min_side_obs,
      min_site_maf = min_site_maf,
      min_gain_norm = min_gain_norm
    )

    if (progress && (i %% chunk_size == 0L || i == length(sites))) {
      message("  scored SNPs ", i, " of ", length(sites))
    }
  }

  do.call(rbind, scores)
}

summarise_site_scores <- function(site_node_scores) {
  if (nrow(site_node_scores) == 0L) {
    return(data.frame())
  }

  aggregate(
    cbind(
      nodes_helped = rep(1, nrow(site_node_scores)),
      gain_max = site_node_scores$gain,
      normalized_gain_max = site_node_scores$normalized_gain,
      weighted_gain_sum = site_node_scores$weight * site_node_scores$normalized_gain
    ) ~ site,
    data = site_node_scores,
    FUN = function(x) c(sum = sum(x), max = max(x), mean = mean(x))[1L]
  )
}

make_site_summary <- function(site_node_scores) {
  if (nrow(site_node_scores) == 0L) {
    return(data.frame(
      site = integer(), nodes_helped = integer(), gain_max = numeric(),
      gain_mean = numeric(), normalized_gain_max = numeric(),
      normalized_gain_mean = numeric(), weighted_gain_sum = numeric()
    ))
  }

  rows <- split(site_node_scores, site_node_scores$site)
  out <- do.call(rbind, lapply(rows, function(x) {
    data.frame(
      site = x$site[[1L]],
      nodes_helped = nrow(x),
      gain_max = max(x$gain),
      gain_mean = mean(x$gain),
      normalized_gain_max = max(x$normalized_gain),
      normalized_gain_mean = mean(x$normalized_gain),
      weighted_gain_sum = sum(x$weight * x$normalized_gain)
    )
  }))
  out[order(-out$weighted_gain_sum, -out$normalized_gain_max), ]
}

make_node_summary <- function(site_node_scores, node_metadata) {
  if (nrow(site_node_scores) == 0L) {
    return(node_metadata[0, ])
  }

  rows <- split(site_node_scores, site_node_scores$node_index)
  summary <- do.call(rbind, lapply(rows, function(x) {
    data.frame(
      node_index = x$node_index[[1L]],
      node_id = x$node_id[[1L]],
      sites_helpful = nrow(x),
      gain_max = max(x$gain),
      normalized_gain_max = max(x$normalized_gain),
      best_site = x$site[which.max(x$normalized_gain)]
    )
  }))
  merge(summary, node_metadata, by = c("node_index", "node_id"), sort = FALSE)
}
