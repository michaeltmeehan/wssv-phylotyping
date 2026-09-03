#!/usr/bin/env Rscript

# Stage 07: scan fixed-length candidate amplicons using leave-one-out
# target-clade classification.
#
# For every candidate window and every held-out reference taxon:
#   1. Rebuild the hierarchical allele model with that taxon removed.
#   2. Construct the query from polymorphic sites inside the window.
#   3. Classify using the available query sites.
#   4. If no query sites are available, return the leave-one-out prior.
#   5. Aggregate terminal posterior mass into target clades + Other.
#   6. Score the posterior relative to the corresponding leave-one-out prior.
#
# This means that every taxon contributes to every candidate window.
# A window that provides no evidence receives exactly zero prior-relative
# improvement rather than being omitted from the score.
#
# Outputs:
#   data/processed/amplicon_scan_cases.csv
#   data/processed/amplicon_scan_summary.csv
#   data/processed/amplicon_scan_by_class.csv
#   data/processed/amplicon_scan.rds

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

precomputed_path <- file.path(processed_dir, "precomputed.rds")
classifier_path <- file.path(processed_dir, "classifier.rds")

if (!file.exists(precomputed_path)) {
  stop("Precomputed data not found: ", precomputed_path, call. = FALSE)
}

if (!file.exists(classifier_path)) {
  stop("Classifier not found: ", classifier_path, call. = FALSE)
}

precomputed <- readRDS(precomputed_path)
classifier <- readRDS(classifier_path)

source("scripts/04_classify_query.R")


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

# Modify these as required.
window_lengths <- c(
  500L
)

window_step <- 500L

# Windows containing fewer than this many polymorphic sites are not evaluated.
# They remain in the summary table with NA scores.
min_snps <- 1L

# Classifier parameters.
lambda <- 20
epsilon <- 0.01
tau <- 0.1

# Recommended first-pass ranking criterion.
#
# This is the mean improvement in log probability assigned to the correct
# class, after first averaging within each true class and then giving those
# classes equal weight.
default_ranking <- "balanced_log_score_improvement"

probability_floor <- 1e-15


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0L) {
    return(0)
  }
  -sum(p * log(p))
}


calculate_tempered_posterior <- function(
    log_likelihood,
    excluded_index,
    tau
) {
  
  n_states <- length(log_likelihood)
  
  prior <- rep(
    1 / (n_states - 1L),
    n_states
  )
  
  prior[excluded_index] <- 0
  
  retained <- prior > 0
  
  log_posterior <- rep(
    -Inf,
    n_states
  )
  
  log_posterior[retained] <-
    log(prior[retained]) +
    tau * log_likelihood[retained]
  
  max_log_posterior <- max(
    log_posterior[retained]
  )
  
  posterior <- numeric(n_states)
  
  posterior[retained] <- exp(
    log_posterior[retained] -
      max_log_posterior
  )
  
  posterior / sum(posterior)
}


aggregate_target_posterior <- function(
    tip_posterior,
    target_tip_mask
) {
  
  target_posterior <- drop(
    target_tip_mask %*% tip_posterior
  )
  
  other_posterior <- 1 - sum(target_posterior)
  
  if (
    other_posterior < 0 &&
    abs(other_posterior) < 1e-12
  ) {
    other_posterior <- 0
  }
  
  if (other_posterior < -1e-12) {
    stop(
      "Target posterior probabilities exceed 1. ",
      "Check whether target clades overlap.",
      call. = FALSE
    )
  }
  
  posterior <- c(
    target_posterior,
    Other = other_posterior
  )
  
  posterior / sum(posterior)
}


multiclass_brier <- function(
    posterior,
    true_class
) {
  
  observed <- setNames(
    rep(0, length(posterior)),
    names(posterior)
  )
  
  observed[true_class] <- 1
  
  sum(
    (posterior - observed)^2
  )
}


summarise_vector <- function(x) {
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}


# -------------------------------------------------------------------------
# Target classes
# -------------------------------------------------------------------------

target_node_ids <- as.integer(
  precomputed$target_node_metadata$node_id
)

if (length(target_node_ids) == 0L) {
  stop("No target clades found in precomputed data.", call. = FALSE)
}

target_names <- paste0(
  "node_",
  target_node_ids
)

names(target_node_ids) <- target_names

target_tip_mask <- classifier$node_tip_mask[
  as.character(target_node_ids),
  ,
  drop = FALSE
]

rownames(target_tip_mask) <- target_names

n_target_memberships <- colSums(target_tip_mask)

if (any(n_target_memberships > 1L)) {
  stop(
    "Configured target clades overlap. ",
    "Target-clade scoring requires mutually exclusive targets.",
    call. = FALSE
  )
}

tip_names <- classifier$terminal_states$taxon

true_class <- setNames(
  rep("Other", length(tip_names)),
  tip_names
)

for (target in target_names) {
  members <- colnames(target_tip_mask)[
    target_tip_mask[target, ]
  ]
  
  true_class[members] <- target
}

class_names <- c(
  target_names,
  "Other"
)


# -------------------------------------------------------------------------
# Candidate windows
# -------------------------------------------------------------------------

alignment_length <- ncol(
  precomputed$aln_int
)

candidate_windows <- do.call(
  rbind,
  lapply(
    window_lengths,
    function(window_length) {
      
      if (window_length > alignment_length) {
        return(NULL)
      }
      
      starts <- seq.int(
        from = 1L,
        to = alignment_length - window_length + 1L,
        by = window_step
      )
      
      final_start <- alignment_length - window_length + 1L
      
      if (tail(starts, 1L) != final_start) {
        starts <- c(
          starts,
          final_start
        )
      }
      
      data.frame(
        window_length = window_length,
        window_start = starts,
        window_end = starts + window_length - 1L,
        stringsAsFactors = FALSE
      )
    }
  )
)

if (is.null(candidate_windows) || nrow(candidate_windows) == 0L) {
  stop("No candidate windows could be constructed.", call. = FALSE)
}

rownames(candidate_windows) <- NULL

candidate_windows$window_id <- sprintf(
  "L%d_%d_%d",
  candidate_windows$window_length,
  candidate_windows$window_start,
  candidate_windows$window_end
)

candidate_windows <- candidate_windows[
  ,
  c(
    "window_id",
    "window_length",
    "window_start",
    "window_end"
  )
]


# -------------------------------------------------------------------------
# Polymorphic sites within each window
# -------------------------------------------------------------------------

polymorphic_sites <- as.integer(
  classifier$polymorphic_sites
)

window_sites <- lapply(
  seq_len(nrow(candidate_windows)),
  function(i) {
    polymorphic_sites[
      polymorphic_sites >= candidate_windows$window_start[[i]] &
        polymorphic_sites <= candidate_windows$window_end[[i]]
    ]
  }
)

candidate_windows$n_polymorphic_sites <- vapply(
  window_sites,
  length,
  integer(1)
)

candidate_windows$sites <- vapply(
  window_sites,
  function(x) paste(x, collapse = ";"),
  character(1)
)

eligible_windows <- which(
  candidate_windows$n_polymorphic_sites >= min_snps
)

if (length(eligible_windows) == 0L) {
  stop(
    "No candidate window contains at least ",
    min_snps,
    " polymorphic site(s).",
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Leave-one-out scan
# -------------------------------------------------------------------------
#
# Iterate over held-out taxa first so the hierarchical allele model is rebuilt
# only once per taxon rather than once per taxon x window.

case_results <- vector(
  "list",
  length(eligible_windows) * length(tip_names)
)

case_counter <- 1L

for (i in seq_along(tip_names)) {
  
  held_out_taxon <- tip_names[[i]]
  
  message(
    "Preparing held-out taxon ",
    i,
    "/",
    length(tip_names),
    ": ",
    held_out_taxon
  )
  
  query_states <- classifier$tip_alleles[
    held_out_taxon,
    ,
    drop = TRUE
  ]
  
  loo_classifier <- classifier
  
  loo_classifier$tip_alleles[
    held_out_taxon,
  ] <- 0L
  
  allele_model <- build_hierarchical_allele_model(
    classifier = loo_classifier,
    lambda = lambda
  )
  
  excluded_index <- match(
    held_out_taxon,
    tip_names
  )
  
  true_target <- unname(
    true_class[held_out_taxon]
  )
  
  
  # -----------------------------------------------------------------------
  # Leave-one-out prior for this held-out taxon
  # -----------------------------------------------------------------------
  #
  # The terminal prior is uniform over all retained reference taxa.
  # Aggregating that terminal prior into the target classes produces the
  # appropriate baseline for this particular LOO case.
  
  prior_tip_posterior <- calculate_tempered_posterior(
    log_likelihood = rep(0, length(tip_names)),
    excluded_index = excluded_index,
    tau = tau
  )
  
  names(prior_tip_posterior) <- tip_names
  
  prior_target_posterior <- aggregate_target_posterior(
    tip_posterior = prior_tip_posterior,
    target_tip_mask = target_tip_mask
  )
  
  prior_true_class_posterior <- unname(
    prior_target_posterior[true_target]
  )
  
  prior_entropy <- entropy(
    prior_target_posterior
  )
  
  prior_log_score <- log(
    max(
      prior_true_class_posterior,
      probability_floor
    )
  )
  
  prior_brier_score <- multiclass_brier(
    posterior = prior_target_posterior,
    true_class = true_target
  )
  
  
  # -----------------------------------------------------------------------
  # Evaluate candidate windows
  # -----------------------------------------------------------------------
  
  for (window_index in eligible_windows) {
    
    sites <- window_sites[[window_index]]
    
    site_columns <- match(
      as.character(sites),
      names(query_states)
    )
    
    if (anyNA(site_columns)) {
      stop(
        "Could not match candidate-window SNPs to classifier site columns.",
        call. = FALSE
      )
    }
    
    query <- as.integer(
      query_states[site_columns]
    )
    
    names(query) <- as.character(sites)
    
    query <- query[
      query != 0L
    ]
    
    has_evidence <- length(query) > 0L
    
    
    # ---------------------------------------------------------------------
    # Posterior
    # ---------------------------------------------------------------------
    #
    # No evidence -> posterior equals the corresponding LOO prior exactly.
    
    if (has_evidence) {
      
      likelihood_result <- calculate_tip_log_likelihoods(
        query = query,
        classifier = loo_classifier,
        allele_model = allele_model,
        epsilon = epsilon
      )
      
      tip_posterior <- calculate_tempered_posterior(
        log_likelihood = likelihood_result$log_likelihood,
        excluded_index = excluded_index,
        tau = tau
      )
      
      names(tip_posterior) <- tip_names
      
      target_posterior <- aggregate_target_posterior(
        tip_posterior = tip_posterior,
        target_tip_mask = target_tip_mask
      )
      
    } else {
      
      target_posterior <- prior_target_posterior
    }
    
    
    # ---------------------------------------------------------------------
    # Absolute posterior scores
    # ---------------------------------------------------------------------
    
    predicted_target <- names(
      which.max(target_posterior)
    )
    
    true_target_posterior <- unname(
      target_posterior[true_target]
    )
    
    predicted_target_posterior <- max(
      target_posterior
    )
    
    posterior_entropy <- entropy(
      target_posterior
    )
    
    log_score <- log(
      max(
        true_target_posterior,
        probability_floor
      )
    )
    
    brier_score <- multiclass_brier(
      posterior = target_posterior,
      true_class = true_target
    )
    
    
    # ---------------------------------------------------------------------
    # Prior-relative scores
    # ---------------------------------------------------------------------
    #
    # Positive log-score improvement:
    #   more posterior probability assigned to the true class than under
    #   the LOO prior.
    #
    # Positive true-posterior improvement:
    #   posterior probability of the true class increased.
    #
    # Positive Brier improvement:
    #   Brier score decreased relative to the prior.
    #
    # Positive entropy reduction:
    #   posterior is more concentrated than the prior.
    #
    # With an empty query all four are exactly zero.
    
    true_posterior_improvement <-
      true_target_posterior -
      prior_true_class_posterior
    
    log_score_improvement <-
      log_score -
      prior_log_score
    
    brier_score_improvement <-
      prior_brier_score -
      brier_score
    
    entropy_reduction <-
      prior_entropy -
      posterior_entropy
    
    
    # ---------------------------------------------------------------------
    # Store case result
    # ---------------------------------------------------------------------
    
    row <- data.frame(
      window_id =
        candidate_windows$window_id[[window_index]],
      
      held_out_taxon =
        held_out_taxon,
      
      true_class =
        true_target,
      
      n_query_sites =
        length(query),
      
      has_evidence =
        has_evidence,
      
      predicted_class =
        predicted_target,
      
      correct =
        predicted_target == true_target,
      
      prior_true_class_posterior =
        prior_true_class_posterior,
      
      true_class_posterior =
        true_target_posterior,
      
      true_posterior_improvement =
        true_posterior_improvement,
      
      prior_log_score =
        prior_log_score,
      
      log_score =
        log_score,
      
      log_score_improvement =
        log_score_improvement,
      
      prior_brier_score =
        prior_brier_score,
      
      brier_score =
        brier_score,
      
      brier_score_improvement =
        brier_score_improvement,
      
      prior_entropy =
        prior_entropy,
      
      posterior_entropy =
        posterior_entropy,
      
      entropy_reduction =
        entropy_reduction,
      
      predicted_class_posterior =
        predicted_target_posterior,
      
      stringsAsFactors = FALSE
    )
    
    for (class_name in class_names) {
      row[[paste0("prior_p_", class_name)]] <-
        unname(prior_target_posterior[class_name])
      
      row[[paste0("p_", class_name)]] <-
        unname(target_posterior[class_name])
    }
    
    case_results[[case_counter]] <- row
    case_counter <- case_counter + 1L
  }
}

case_results <- case_results[
  !vapply(case_results, is.null, logical(1))
]

case_table <- do.call(
  rbind,
  case_results
)

rownames(case_table) <- NULL


# -------------------------------------------------------------------------
# Per-class summaries
# -------------------------------------------------------------------------
#
# These are particularly useful for later panel selection because they retain
# the utility profile of each amplicon across the target classes.

window_class_groups <- split(
  case_table,
  interaction(
    case_table$window_id,
    case_table$true_class,
    drop = TRUE,
    lex.order = TRUE
  )
)

by_class_results <- lapply(
  window_class_groups,
  function(x) {
    
    data.frame(
      window_id =
        x$window_id[[1]],
      
      true_class =
        x$true_class[[1]],
      
      n_taxa =
        nrow(x),
      
      n_with_evidence =
        sum(x$has_evidence),
      
      evidence_rate =
        mean(x$has_evidence),
      
      mean_query_sites =
        mean(x$n_query_sites),
      
      accuracy =
        mean(x$correct),
      
      mean_prior_true_class_posterior =
        mean(x$prior_true_class_posterior),
      
      mean_true_class_posterior =
        mean(x$true_class_posterior),
      
      mean_true_posterior_improvement =
        mean(x$true_posterior_improvement),
      
      mean_log_score =
        mean(x$log_score),
      
      mean_log_score_improvement =
        mean(x$log_score_improvement),
      
      mean_brier_score =
        mean(x$brier_score),
      
      mean_brier_score_improvement =
        mean(x$brier_score_improvement),
      
      mean_entropy_reduction =
        mean(x$entropy_reduction),
      
      stringsAsFactors = FALSE
    )
  }
)

by_class_table <- do.call(
  rbind,
  by_class_results
)

rownames(by_class_table) <- NULL


# -------------------------------------------------------------------------
# Window-level summaries
# -------------------------------------------------------------------------

window_groups <- split(
  case_table,
  case_table$window_id
)

summary_results <- lapply(
  window_groups,
  function(x) {
    
    window_id <- x$window_id[[1]]
    
    window_row <- candidate_windows[
      candidate_windows$window_id == window_id,
      ,
      drop = FALSE
    ]
    
    class_rows <- by_class_table[
      by_class_table$window_id == window_id,
      ,
      drop = FALSE
    ]
    
    
    # ---------------------------------------------------------------------
    # Empirical / tip-weighted scores
    # ---------------------------------------------------------------------
    
    empirical_accuracy <-
      mean(x$correct)
    
    empirical_mean_true_posterior <-
      mean(x$true_class_posterior)
    
    empirical_true_posterior_improvement <-
      mean(x$true_posterior_improvement)
    
    empirical_log_score <-
      mean(x$log_score)
    
    empirical_log_score_improvement <-
      mean(x$log_score_improvement)
    
    empirical_brier_score <-
      mean(x$brier_score)
    
    empirical_brier_score_improvement <-
      mean(x$brier_score_improvement)
    
    empirical_entropy_reduction <-
      mean(x$entropy_reduction)
    
    
    # ---------------------------------------------------------------------
    # Balanced class scores
    # ---------------------------------------------------------------------
    #
    # Every true class represented in the reference set contributes equally.
    # Importantly, a class with no query evidence is NOT omitted: its
    # prior-relative improvement is zero.
    
    balanced_accuracy <-
      mean(class_rows$accuracy)
    
    balanced_mean_true_posterior <-
      mean(class_rows$mean_true_class_posterior)
    
    balanced_true_posterior_improvement <-
      mean(class_rows$mean_true_posterior_improvement)
    
    balanced_log_score <-
      mean(class_rows$mean_log_score)
    
    balanced_log_score_improvement <-
      mean(class_rows$mean_log_score_improvement)
    
    balanced_brier_score <-
      mean(class_rows$mean_brier_score)
    
    balanced_brier_score_improvement <-
      mean(class_rows$mean_brier_score_improvement)
    
    balanced_entropy_reduction <-
      mean(class_rows$mean_entropy_reduction)
    
    
    # ---------------------------------------------------------------------
    # Utility breadth
    # ---------------------------------------------------------------------
    #
    # These are descriptive only. They do not require an amplicon to be useful
    # for every target. They indicate how broadly its useful information is
    # distributed across classes.
    
    n_classes_with_positive_log_gain <- sum(
      class_rows$mean_log_score_improvement > 0
    )
    
    n_classes_with_evidence <- sum(
      class_rows$n_with_evidence > 0
    )
    
    
    data.frame(
      window_id =
        window_id,
      
      window_length =
        window_row$window_length,
      
      window_start =
        window_row$window_start,
      
      window_end =
        window_row$window_end,
      
      n_polymorphic_sites =
        window_row$n_polymorphic_sites,
      
      sites =
        window_row$sites,
      
      n_taxa =
        nrow(x),
      
      n_with_evidence =
        sum(x$has_evidence),
      
      evidence_rate =
        mean(x$has_evidence),
      
      mean_query_sites =
        mean(x$n_query_sites),
      
      n_classes_with_evidence =
        n_classes_with_evidence,
      
      n_classes_with_positive_log_gain =
        n_classes_with_positive_log_gain,
      
      empirical_accuracy =
        empirical_accuracy,
      
      empirical_mean_true_posterior =
        empirical_mean_true_posterior,
      
      empirical_true_posterior_improvement =
        empirical_true_posterior_improvement,
      
      empirical_log_score =
        empirical_log_score,
      
      empirical_log_score_improvement =
        empirical_log_score_improvement,
      
      empirical_brier_score =
        empirical_brier_score,
      
      empirical_brier_score_improvement =
        empirical_brier_score_improvement,
      
      empirical_entropy_reduction =
        empirical_entropy_reduction,
      
      balanced_accuracy =
        balanced_accuracy,
      
      balanced_mean_true_posterior =
        balanced_mean_true_posterior,
      
      balanced_true_posterior_improvement =
        balanced_true_posterior_improvement,
      
      balanced_log_score =
        balanced_log_score,
      
      balanced_log_score_improvement =
        balanced_log_score_improvement,
      
      balanced_brier_score =
        balanced_brier_score,
      
      balanced_brier_score_improvement =
        balanced_brier_score_improvement,
      
      balanced_entropy_reduction =
        balanced_entropy_reduction,
      
      stringsAsFactors = FALSE
    )
  }
)

summary_table <- do.call(
  rbind,
  summary_results
)

rownames(summary_table) <- NULL


# -------------------------------------------------------------------------
# Add windows below min_snps
# -------------------------------------------------------------------------

ineligible_windows <- setdiff(
  seq_len(nrow(candidate_windows)),
  eligible_windows
)

if (length(ineligible_windows) > 0L) {
  
  skipped <- candidate_windows[
    ineligible_windows,
    ,
    drop = FALSE
  ]
  
  skipped$n_taxa <- length(tip_names)
  skipped$n_with_evidence <- 0L
  skipped$evidence_rate <- 0
  skipped$mean_query_sites <- 0
  
  skipped$n_classes_with_evidence <- 0L
  skipped$n_classes_with_positive_log_gain <- 0L
  
  score_columns <- setdiff(
    names(summary_table),
    names(skipped)
  )
  
  for (column in score_columns) {
    skipped[[column]] <- NA_real_
  }
  
  summary_table <- rbind(
    summary_table,
    skipped[
      ,
      names(summary_table)
    ]
  )
}


# -------------------------------------------------------------------------
# Rank candidate windows
# -------------------------------------------------------------------------

if (!default_ranking %in% names(summary_table)) {
  stop(
    "Unknown default ranking metric: ",
    default_ranking,
    call. = FALSE
  )
}

smaller_is_better <- grepl(
  "(^|_)brier_score$",
  default_ranking
)

ranking_value <- summary_table[[default_ranking]]

summary_table$rank <- if (smaller_is_better) {
  
  rank(
    ranking_value,
    ties.method = "min",
    na.last = "keep"
  )
  
} else {
  
  rank(
    -ranking_value,
    ties.method = "min",
    na.last = "keep"
  )
}

summary_table <- summary_table[
  order(
    summary_table$rank,
    summary_table$window_length,
    summary_table$window_start,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

rownames(summary_table) <- NULL


# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------

case_output_path <- file.path(
  processed_dir,
  "amplicon_scan_cases.csv"
)

summary_output_path <- file.path(
  processed_dir,
  "amplicon_scan_summary.csv"
)

by_class_output_path <- file.path(
  processed_dir,
  "amplicon_scan_by_class.csv"
)

rds_output_path <- file.path(
  processed_dir,
  "amplicon_scan.rds"
)

write.csv(
  case_table,
  case_output_path,
  row.names = FALSE
)

write.csv(
  summary_table,
  summary_output_path,
  row.names = FALSE
)

write.csv(
  by_class_table,
  by_class_output_path,
  row.names = FALSE
)

saveRDS(
  list(
    lambda = lambda,
    epsilon = epsilon,
    tau = tau,
    window_lengths = window_lengths,
    window_step = window_step,
    min_snps = min_snps,
    target_node_ids = target_node_ids,
    classes = class_names,
    default_ranking = default_ranking,
    windows = candidate_windows,
    cases = case_table,
    by_class = by_class_table,
    summary = summary_table
  ),
  rds_output_path
)


# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------

message("")
message("Amplicon scan complete")

message(
  "  lambda: ",
  lambda
)

message(
  "  epsilon: ",
  epsilon
)

message(
  "  tau: ",
  tau
)

message(
  "  window lengths: ",
  paste(window_lengths, collapse = ", ")
)

message(
  "  window step: ",
  window_step
)

message(
  "  candidate windows: ",
  nrow(candidate_windows)
)

message(
  "  eligible windows: ",
  length(eligible_windows)
)

message(
  "  targets: ",
  paste(target_names, collapse = ", ")
)

message("")
message(
  "Top windows by ",
  default_ranking,
  ":"
)

print(
  head(summary_table, 20L),
  row.names = FALSE
)

message("")
message("  wrote: ", case_output_path)
message("  wrote: ", summary_output_path)
message("  wrote: ", by_class_output_path)
message("  wrote: ", rds_output_path)
