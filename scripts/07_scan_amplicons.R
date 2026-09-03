#!/usr/bin/env Rscript

# Stage 07: scan fixed-length candidate amplicons using leave-one-out
# target-clade classification.
#
# Inputs:
#   data/processed/precomputed.rds
#   data/processed/classifier.rds
#   scripts/04_classify_query.R
#
# Outputs:
#   data/processed/amplicon_scan_cases.csv
#   data/processed/amplicon_scan_summary.csv
#   data/processed/amplicon_scan_by_class.csv
#   data/processed/amplicon_scan.rds
#
# For each fixed-length candidate window:
#   1. Identify polymorphic sites falling inside the window.
#   2. Leave each reference taxon out in turn.
#   3. Remove that taxon's sequence from the hierarchical allele model.
#   4. Classify the held-out taxon using only SNPs inside the window.
#   5. Aggregate terminal posterior probabilities into the configured target
#      clades plus an "Other" class.
#   6. Score the window using several posterior-based criteria.
#
# The target clades are used only for evaluation. The underlying classifier
# remains phylogeny-wide.
#
# IMPORTANT:
# Candidate windows are currently defined in alignment coordinates because
# classifier SNP indices are alignment-column indices. Reference-coordinate
# conversion can be added later for primer design.

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "Package 'yaml' is required to read config/config.yml.",
    call. = FALSE
  )
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed


# -------------------------------------------------------------------------
# Load objects and classifier functions
# -------------------------------------------------------------------------

precomputed_path <- file.path(
  processed_dir,
  "precomputed.rds"
)

classifier_path <- file.path(
  processed_dir,
  "classifier.rds"
)

if (!file.exists(precomputed_path)) {
  stop(
    "Precomputed data not found: ",
    precomputed_path,
    "\nRun stage 01 first.",
    call. = FALSE
  )
}

if (!file.exists(classifier_path)) {
  stop(
    "Classifier not found: ",
    classifier_path,
    "\nRun stage 03 first.",
    call. = FALSE
  )
}

precomputed <- readRDS(
  precomputed_path
)

classifier <- readRDS(
  classifier_path
)

source(
  "scripts/04_classify_query.R"
)


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------
#
# Start with one amplicon length if desired, e.g.:
#
#   window_lengths <- 500L
#
# Multiple lengths can be scanned in one run.
#
# A step of 1 evaluates every possible window and is exhaustive but may be
# slower. A larger step is useful for an initial coarse scan.

# window_lengths <- c(
#   250L,
#   500L,
#   750L,
#   1000L
# )

window_lengths <- c(500L)

# window_step <- 25L

window_step <- 500L

# Minimum number of polymorphic sites required for a candidate window to be
# evaluated. Windows below this threshold are retained in the summary table
# but receive NA classification scores.

min_snps <- 1L

# Classifier parameters.

lambda <- 20
epsilon <- 0.01
tau <- 0.1

# Prior used when averaging performance across held-out taxa.
#
# "empirical":
#   every held-out taxon contributes equally. Larger target clades therefore
#   contribute more observations to the mean.
#
# "balanced":
#   every target class, including Other, contributes equally to the final
#   window-level mean irrespective of its number of reference taxa.
#
# Both empirical and balanced metrics are actually calculated below, so this
# option only controls the default ranking column at the end.

default_ranking <- "balanced_log_score"

# Small floor used before taking log probabilities.

probability_floor <- 1e-15


# -------------------------------------------------------------------------
# Basic validation
# -------------------------------------------------------------------------

if (
  any(window_lengths <= 0L) ||
  any(window_lengths != as.integer(window_lengths))
) {
  stop(
    "'window_lengths' must contain positive integers.",
    call. = FALSE
  )
}

if (
  length(window_step) != 1L ||
  is.na(window_step) ||
  window_step <= 0L ||
  window_step != as.integer(window_step)
) {
  stop(
    "'window_step' must be a single positive integer.",
    call. = FALSE
  )
}

if (
  length(min_snps) != 1L ||
  is.na(min_snps) ||
  min_snps < 0L ||
  min_snps != as.integer(min_snps)
) {
  stop(
    "'min_snps' must be a single non-negative integer.",
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Identify configured target nodes
# -------------------------------------------------------------------------

target_node_ids <- as.integer(
  precomputed$target_node_metadata$node_id
)

if (length(target_node_ids) == 0L) {
  stop(
    "No target clades found in precomputed data.",
    call. = FALSE
  )
}

target_names <- paste0(
  "node_",
  target_node_ids
)

names(target_node_ids) <- target_names


# -------------------------------------------------------------------------
# Construct target x terminal-state membership matrix
# -------------------------------------------------------------------------

target_tip_mask <- classifier$node_tip_mask[
  as.character(target_node_ids),
  ,
  drop = FALSE
]

rownames(target_tip_mask) <- target_names

n_target_memberships <- colSums(
  target_tip_mask
)

if (any(n_target_memberships > 1L)) {
  overlapping_taxa <- names(
    n_target_memberships[
      n_target_memberships > 1L
    ]
  )

  stop(
    "Configured target clades overlap. Amplicon scoring at target-clade ",
    "resolution requires mutually exclusive target nodes.\nOverlapping taxa: ",
    paste(overlapping_taxa, collapse = ", "),
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Determine true target class for every reference taxon
# -------------------------------------------------------------------------

tip_names <- classifier$terminal_states$taxon

true_class <- rep(
  "Other",
  length(tip_names)
)

names(true_class) <- tip_names

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
# Helper: temper and normalise terminal log-likelihoods
# -------------------------------------------------------------------------

calculate_tempered_posterior <- function(
    log_likelihood,
    excluded_index,
    tau
) {

  n_states <- length(
    log_likelihood
  )

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

  posterior <- numeric(
    n_states
  )

  posterior[retained] <- exp(
    log_posterior[retained] -
      max_log_posterior
  )

  posterior <- posterior / sum(
    posterior
  )

  posterior
}


# -------------------------------------------------------------------------
# Helper: aggregate terminal posterior into targets + Other
# -------------------------------------------------------------------------

aggregate_target_posterior <- function(
    tip_posterior,
    target_tip_mask
) {

  target_posterior <- drop(
    target_tip_mask %*%
      tip_posterior
  )

  other_posterior <- 1 - sum(
    target_posterior
  )

  if (
    other_posterior < 0 &&
    abs(other_posterior) < 1e-12
  ) {
    other_posterior <- 0
  }

  if (other_posterior < -1e-12) {
    stop(
      "Target posterior probabilities exceed 1.",
      call. = FALSE
    )
  }

  posterior <- c(
    target_posterior,
    Other = other_posterior
  )

  posterior <- posterior / sum(
    posterior
  )

  posterior
}


# -------------------------------------------------------------------------
# Helper: entropy
# -------------------------------------------------------------------------

entropy <- function(p) {
  p <- p[
    is.finite(p) &
      p > 0
  ]

  if (length(p) == 0L) {
    return(0)
  }

  -sum(
    p * log(p)
  )
}


# -------------------------------------------------------------------------
# Construct candidate windows
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

      # Always include the final possible window even if the chosen step
      # does not land exactly on it.

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

if (
  is.null(candidate_windows) ||
  nrow(candidate_windows) == 0L
) {
  stop(
    "No candidate windows could be constructed.",
    call. = FALSE
  )
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
# Identify polymorphic sites in each window
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
  function(x) {
    paste(
      x,
      collapse = ";"
    )
  },
  character(1)
)


# -------------------------------------------------------------------------
# Prior target-class distribution
# -------------------------------------------------------------------------
#
# This is the target-class prior induced by equal prior weight over retained
# terminal states. It is useful for measuring posterior entropy reduction.

class_counts <- table(
  factor(
    true_class,
    levels = class_names
  )
)

empirical_class_prior <- as.numeric(
  class_counts
)

empirical_class_prior <- empirical_class_prior /
  sum(empirical_class_prior)

names(empirical_class_prior) <- class_names

prior_entropy <- entropy(
  empirical_class_prior
)


# -------------------------------------------------------------------------
# Preallocate case results
# -------------------------------------------------------------------------

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

case_results <- vector(
  "list",
  length(eligible_windows) * length(tip_names)
)

case_counter <- 1L


# -------------------------------------------------------------------------
# Leave-one-out scan
# -------------------------------------------------------------------------
#
# We iterate over held-out taxa first. The hierarchical allele model therefore
# needs to be rebuilt only once per taxon, not once per taxon x window.

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


  # -----------------------------------------------------------------------
  # Construct LOO classifier
  # -----------------------------------------------------------------------

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
    true_class[
      held_out_taxon
    ]
  )


  # -----------------------------------------------------------------------
  # Evaluate every eligible candidate window
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
      query_states[
        site_columns
      ]
    )

    names(query) <- as.character(
      sites
    )

    query <- query[
      query != 0L
    ]

    # If this held-out taxon has no called SNP within the window, there is no
    # query evidence. Record the case as unavailable rather than forcing the
    # classifier to operate on an empty query.

    if (length(query) == 0L) {

      row <- data.frame(
        window_id =
          candidate_windows$window_id[[window_index]],

        held_out_taxon =
          held_out_taxon,

        true_class =
          true_target,

        n_query_sites =
          0L,

        predicted_class =
          NA_character_,

        correct =
          NA,

        true_class_posterior =
          NA_real_,

        predicted_class_posterior =
          NA_real_,

        posterior_entropy =
          NA_real_,

        entropy_reduction =
          NA_real_,

        log_score =
          NA_real_,

        brier_score =
          NA_real_,

        stringsAsFactors = FALSE
      )

      for (class_name in class_names) {
        row[[paste0("p_", class_name)]] <- NA_real_
      }

      case_results[[case_counter]] <- row
      case_counter <- case_counter + 1L

      next
    }


    # ---------------------------------------------------------------------
    # Terminal likelihood and posterior
    # ---------------------------------------------------------------------

    likelihood_result <- calculate_tip_log_likelihoods(
      query = query,
      classifier = loo_classifier,
      allele_model = allele_model,
      epsilon = epsilon
    )

    tip_posterior <- calculate_tempered_posterior(
      log_likelihood =
        likelihood_result$log_likelihood,

      excluded_index =
        excluded_index,

      tau =
        tau
    )

    names(tip_posterior) <- tip_names


    # ---------------------------------------------------------------------
    # Aggregate to target-clade resolution
    # ---------------------------------------------------------------------

    target_posterior <- aggregate_target_posterior(
      tip_posterior =
        tip_posterior,

      target_tip_mask =
        target_tip_mask
    )

    predicted_target <- names(
      which.max(
        target_posterior
      )
    )

    true_target_posterior <- unname(
      target_posterior[
        true_target
      ]
    )

    predicted_target_posterior <- max(
      target_posterior
    )


    # ---------------------------------------------------------------------
    # Scoring metrics
    # ---------------------------------------------------------------------

    posterior_entropy <- entropy(
      target_posterior
    )

    # Entropy reduction relative to the empirical class prior.
    #
    # Positive values indicate that the posterior is more concentrated than
    # the prior. This is a case-level quantity; the window-level mean is an
    # empirical approximation to expected information gain.

    entropy_reduction <-
      prior_entropy -
      posterior_entropy

    # Log score: larger is better (best possible value is zero).

    log_score <- log(
      max(
        true_target_posterior,
        probability_floor
      )
    )

    # Multiclass Brier score: smaller is better.

    observed <- setNames(
      rep(
        0,
        length(class_names)
      ),
      class_names
    )

    observed[
      true_target
    ] <- 1

    brier_score <- sum(
      (
        target_posterior -
          observed
      )^2
    )


    # ---------------------------------------------------------------------
    # Store case-level result
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

      predicted_class =
        predicted_target,

      correct =
        predicted_target == true_target,

      true_class_posterior =
        true_target_posterior,

      predicted_class_posterior =
        predicted_target_posterior,

      posterior_entropy =
        posterior_entropy,

      entropy_reduction =
        entropy_reduction,

      log_score =
        log_score,

      brier_score =
        brier_score,

      stringsAsFactors = FALSE
    )

    for (class_name in class_names) {
      row[[paste0("p_", class_name)]] <- unname(
        target_posterior[
          class_name
        ]
      )
    }

    case_results[[case_counter]] <- row
    case_counter <- case_counter + 1L
  }
}


# -------------------------------------------------------------------------
# Assemble case table
# -------------------------------------------------------------------------

case_results <- case_results[
  !vapply(
    case_results,
    is.null,
    logical(1)
  )
]

case_table <- do.call(
  rbind,
  case_results
)

rownames(case_table) <- NULL


# -------------------------------------------------------------------------
# Class-specific window summaries
# -------------------------------------------------------------------------

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

    valid <- !is.na(
      x$true_class_posterior
    )

    x_valid <- x[
      valid,
      ,
      drop = FALSE
    ]

    data.frame(
      window_id =
        x$window_id[[1]],

      true_class =
        x$true_class[[1]],

      n_taxa =
        nrow(x),

      n_classifiable =
        nrow(x_valid),

      call_rate =
        mean(valid),

      mean_query_sites =
        if (nrow(x_valid) > 0L) {
          mean(x_valid$n_query_sites)
        } else {
          NA_real_
        },

      accuracy =
        if (nrow(x_valid) > 0L) {
          mean(x_valid$correct)
        } else {
          NA_real_
        },

      mean_true_class_posterior =
        if (nrow(x_valid) > 0L) {
          mean(x_valid$true_class_posterior)
        } else {
          NA_real_
        },

      mean_log_score =
        if (nrow(x_valid) > 0L) {
          mean(x_valid$log_score)
        } else {
          NA_real_
        },

      mean_brier_score =
        if (nrow(x_valid) > 0L) {
          mean(x_valid$brier_score)
        } else {
          NA_real_
        },

      mean_entropy_reduction =
        if (nrow(x_valid) > 0L) {
          mean(x_valid$entropy_reduction)
        } else {
          NA_real_
        },

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
# Window-level summary metrics
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

    valid <- !is.na(
      x$true_class_posterior
    )

    x_valid <- x[
      valid,
      ,
      drop = FALSE
    ]

    if (nrow(x_valid) == 0L) {

      return(
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

          n_classifiable =
            0L,

          call_rate =
            0,

          mean_query_sites =
            NA_real_,

          empirical_accuracy =
            NA_real_,

          empirical_mean_true_posterior =
            NA_real_,

          empirical_log_score =
            NA_real_,

          empirical_brier_score =
            NA_real_,

          empirical_entropy_reduction =
            NA_real_,

          balanced_accuracy =
            NA_real_,

          balanced_mean_true_posterior =
            NA_real_,

          balanced_log_score =
            NA_real_,

          balanced_brier_score =
            NA_real_,

          balanced_entropy_reduction =
            NA_real_,

          worst_class_accuracy =
            NA_real_,

          worst_class_mean_true_posterior =
            NA_real_,

          stringsAsFactors = FALSE
        )
      )
    }


    # ---------------------------------------------------------------------
    # Empirical/tip-weighted scores
    # ---------------------------------------------------------------------

    empirical_accuracy <- mean(
      x_valid$correct
    )

    empirical_mean_true_posterior <- mean(
      x_valid$true_class_posterior
    )

    empirical_log_score <- mean(
      x_valid$log_score
    )

    empirical_brier_score <- mean(
      x_valid$brier_score
    )

    empirical_entropy_reduction <- mean(
      x_valid$entropy_reduction
    )


    # ---------------------------------------------------------------------
    # Balanced target-class scores
    # ---------------------------------------------------------------------
    #
    # Calculate each score within each true class, then average the class
    # means. This gives each represented target class equal weight.

    class_rows <- by_class_table[
      by_class_table$window_id == window_id,
      ,
      drop = FALSE
    ]

    class_rows <- class_rows[
      class_rows$n_classifiable > 0L,
      ,
      drop = FALSE
    ]

    balanced_accuracy <- mean(
      class_rows$accuracy
    )

    balanced_mean_true_posterior <- mean(
      class_rows$mean_true_class_posterior
    )

    balanced_log_score <- mean(
      class_rows$mean_log_score
    )

    balanced_brier_score <- mean(
      class_rows$mean_brier_score
    )

    balanced_entropy_reduction <- mean(
      class_rows$mean_entropy_reduction
    )

    worst_class_accuracy <- min(
      class_rows$accuracy
    )

    worst_class_mean_true_posterior <- min(
      class_rows$mean_true_class_posterior
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

      n_classifiable =
        nrow(x_valid),

      call_rate =
        mean(valid),

      mean_query_sites =
        mean(
          x_valid$n_query_sites
        ),

      empirical_accuracy =
        empirical_accuracy,

      empirical_mean_true_posterior =
        empirical_mean_true_posterior,

      empirical_log_score =
        empirical_log_score,

      empirical_brier_score =
        empirical_brier_score,

      empirical_entropy_reduction =
        empirical_entropy_reduction,

      balanced_accuracy =
        balanced_accuracy,

      balanced_mean_true_posterior =
        balanced_mean_true_posterior,

      balanced_log_score =
        balanced_log_score,

      balanced_brier_score =
        balanced_brier_score,

      balanced_entropy_reduction =
        balanced_entropy_reduction,

      worst_class_accuracy =
        worst_class_accuracy,

      worst_class_mean_true_posterior =
        worst_class_mean_true_posterior,

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
# Add windows below min_snps to summary
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
  skipped$n_classifiable <- 0L
  skipped$call_rate <- 0
  skipped$mean_query_sites <- NA_real_
  skipped$empirical_accuracy <- NA_real_
  skipped$empirical_mean_true_posterior <- NA_real_
  skipped$empirical_log_score <- NA_real_
  skipped$empirical_brier_score <- NA_real_
  skipped$empirical_entropy_reduction <- NA_real_
  skipped$balanced_accuracy <- NA_real_
  skipped$balanced_mean_true_posterior <- NA_real_
  skipped$balanced_log_score <- NA_real_
  skipped$balanced_brier_score <- NA_real_
  skipped$balanced_entropy_reduction <- NA_real_
  skipped$worst_class_accuracy <- NA_real_
  skipped$worst_class_mean_true_posterior <- NA_real_

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
#
# Metrics where larger is better:
#   *_accuracy
#   *_mean_true_posterior
#   *_log_score
#   *_entropy_reduction
#   worst_class_*
#
# Metrics where smaller is better:
#   *_brier_score

larger_is_better <- !grepl(
  "brier",
  default_ranking,
  fixed = TRUE
)

if (!default_ranking %in% names(summary_table)) {
  stop(
    "Unknown default ranking metric: ",
    default_ranking,
    call. = FALSE
  )
}

ranking_value <- summary_table[[default_ranking]]

summary_table$rank <- if (larger_is_better) {
  rank(
    -ranking_value,
    ties.method = "min",
    na.last = "keep"
  )
} else {
  rank(
    ranking_value,
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
# Save results
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
    lambda =
      lambda,

    epsilon =
      epsilon,

    tau =
      tau,

    window_lengths =
      window_lengths,

    window_step =
      window_step,

    min_snps =
      min_snps,

    target_node_ids =
      target_node_ids,

    classes =
      class_names,

    empirical_class_prior =
      empirical_class_prior,

    prior_entropy =
      prior_entropy,

    default_ranking =
      default_ranking,

    windows =
      candidate_windows,

    cases =
      case_table,

    by_class =
      by_class_table,

    summary =
      summary_table
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
  paste(
    window_lengths,
    collapse = ", "
  )
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
  paste(
    target_names,
    collapse = ", "
  )
)

message("")
message(
  "Top windows by ",
  default_ranking,
  ":"
)

print(
  head(
    summary_table,
    20L
  ),
  row.names = FALSE
)

message("")
message(
  "  wrote: ",
  case_output_path
)

message(
  "  wrote: ",
  summary_output_path
)

message(
  "  wrote: ",
  by_class_output_path
)

message(
  "  wrote: ",
  rds_output_path
)
