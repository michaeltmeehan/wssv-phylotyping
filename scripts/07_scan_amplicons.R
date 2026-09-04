#!/usr/bin/env Rscript

# Stage 07: scan fixed-length candidate amplicons across a lambda x tau grid
# using leave-one-out target-clade classification.
#
# For every held-out taxon, candidate window, lambda and tau:
#   1. Remove the held-out taxon from the hierarchical allele model.
#   2. Build the allele model for the current lambda.
#   3. Construct the query from polymorphic sites within the window.
#   4. Calculate terminal log-likelihoods once for that lambda/window.
#   5. Reuse those log-likelihoods across all tau values.
#   6. Aggregate terminal posterior mass into target clades + Other.
#   7. Score posterior performance relative to the leave-one-out prior.
#
# Empty queries return the leave-one-out prior exactly, giving zero
# prior-relative improvement rather than being omitted.
#
# Outputs:
#   data/processed/amplicon_scan_cases.csv
#       One row per taxon x window x lambda x tau.
#
#   data/processed/amplicon_scan_by_class.csv
#       Per-target utility profiles for every window x lambda x tau.
#
#   data/processed/amplicon_scan_parameter_summary.csv
#       One row per window x lambda x tau.
#
#   data/processed/amplicon_scan_summary.csv
#       Robust window-level summary across the full lambda x tau grid.
#
#   data/processed/amplicon_scan.rds
#       All of the above plus configuration.
#
# The primary robust ranking is the mean balanced log-score improvement
# across the full lambda x tau grid. Worst-case and variability measures are
# also retained so that parameter-sensitive amplicons can be identified.

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "Package 'yaml' is required to read config/config.yml.",
    call. = FALSE
  )
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

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
    call. = FALSE
  )
}

if (!file.exists(classifier_path)) {
  stop(
    "Classifier not found: ",
    classifier_path,
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

# Candidate amplicon geometry.
#
# A coarse 500-bp tiling is retained here because the exhaustive scans were
# relatively slow. Reduce window_step later for finer local searches around
# promising regions.

window_lengths <- c(
  1000L
)

window_step <- 500L

min_snps <- 1L


# -------------------------------------------------------------------------
# Classifier sensitivity grid
# -------------------------------------------------------------------------
#
# lambda controls hierarchical/phylogenetic smoothing.
#
# The grid spans:
#   0     no hierarchical smoothing
#   1     weak smoothing
#   5     moderate smoothing
#   10    substantial smoothing
#   20    the previously used strong-smoothing setting
#
# tau controls likelihood tempering.
#
# The grid spans:
#   0.05  very strong tempering
#   0.10  previously used strong tempering
#   0.25  moderate tempering
#   0.50  mild tempering
#   1.00  untempered likelihood
#
# tau = 0 is deliberately omitted because it discards all query evidence and
# therefore produces the prior for every window.

lambda_values <- c(
  0,
  1,
  5,
  10,
  20
)

tau_values <- c(
  0.05,
  0.10,
  0.25,
  0.50,
  1.00
)

# Keep epsilon fixed for this first sensitivity scan.
epsilon <- 0.01

# Robust window-level ranking criterion.
#
# This averages balanced log-score improvement over the complete lambda x tau
# grid. It therefore rewards windows that perform well across plausible
# classifier settings rather than windows that have one exceptional setting.
default_ranking <- "mean_balanced_log_score_improvement"

probability_floor <- 1e-15


# -------------------------------------------------------------------------
# Validation
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

if (
  any(!is.finite(lambda_values)) ||
  any(lambda_values < 0)
) {
  stop(
    "'lambda_values' must contain finite non-negative values.",
    call. = FALSE
  )
}

if (
  any(!is.finite(tau_values)) ||
  any(tau_values <= 0)
) {
  stop(
    "'tau_values' must contain finite positive values.",
    call. = FALSE
  )
}

lambda_values <- sort(
  unique(lambda_values)
)

tau_values <- sort(
  unique(tau_values)
)


# -------------------------------------------------------------------------
# Helpers
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
  
  posterior / sum(
    posterior
  )
}


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
      "Target posterior probabilities exceed 1. ",
      "Check whether target clades overlap.",
      call. = FALSE
    )
  }
  
  posterior <- c(
    target_posterior,
    Other = other_posterior
  )
  
  posterior / sum(
    posterior
  )
}


multiclass_brier <- function(
    posterior,
    true_class
) {
  
  observed <- setNames(
    rep(
      0,
      length(posterior)
    ),
    names(posterior)
  )
  
  observed[
    true_class
  ] <- 1
  
  sum(
    (
      posterior -
        observed
    )^2
  )
}


safe_sd <- function(x) {
  
  if (length(x) <= 1L) {
    return(0)
  }
  
  sd(x)
}


# -------------------------------------------------------------------------
# Target classes
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

target_tip_mask <- classifier$node_tip_mask[
  as.character(target_node_ids),
  ,
  drop = FALSE
]

rownames(
  target_tip_mask
) <- target_names

n_target_memberships <- colSums(
  target_tip_mask
)

if (any(n_target_memberships > 1L)) {
  stop(
    "Configured target clades overlap. ",
    "Target-clade scoring requires mutually exclusive targets.",
    call. = FALSE
  )
}

tip_names <- classifier$terminal_states$taxon

true_class <- setNames(
  rep(
    "Other",
    length(tip_names)
  ),
  tip_names
)

for (target in target_names) {
  
  members <- colnames(
    target_tip_mask
  )[
    target_tip_mask[
      target,
    ]
  ]
  
  true_class[
    members
  ] <- target
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
      
      final_start <-
        alignment_length -
        window_length +
        1L
      
      if (
        tail(
          starts,
          1L
        ) != final_start
      ) {
        starts <- c(
          starts,
          final_start
        )
      }
      
      data.frame(
        window_length =
          window_length,
        
        window_start =
          starts,
        
        window_end =
          starts +
          window_length -
          1L,
        
        stringsAsFactors =
          FALSE
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

rownames(
  candidate_windows
) <- NULL

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
# Polymorphic sites in each window
# -------------------------------------------------------------------------

polymorphic_sites <- as.integer(
  classifier$polymorphic_sites
)

window_sites <- lapply(
  seq_len(
    nrow(candidate_windows)
  ),
  function(i) {
    
    polymorphic_sites[
      polymorphic_sites >=
        candidate_windows$window_start[[i]] &
        polymorphic_sites <=
        candidate_windows$window_end[[i]]
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

eligible_windows <- which(
  candidate_windows$n_polymorphic_sites >=
    min_snps
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
# Leave-one-out scan across lambda x tau
# -------------------------------------------------------------------------
#
# Computational structure:
#
#   held-out taxon
#       -> lambda
#           -> build hierarchical allele model once
#           -> window
#               -> calculate log-likelihood once
#               -> tau
#                   -> temper same log-likelihood
#
# Thus increasing the number of tau values is relatively cheap. Increasing
# the number of lambda values is the main computational cost.

n_case_rows <-
  length(tip_names) *
  length(eligible_windows) *
  length(lambda_values) *
  length(tau_values)

case_results <- vector(
  "list",
  n_case_rows
)

case_counter <- 1L

for (i in seq_along(tip_names)) {
  
  held_out_taxon <- tip_names[[i]]
  
  message(
    "Held-out taxon ",
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
  # Leave-one-out prior
  # -----------------------------------------------------------------------
  #
  # This prior does not depend on lambda or tau. With no likelihood evidence,
  # every retained terminal state has equal prior probability.
  
  prior_tip_posterior <- calculate_tempered_posterior(
    log_likelihood =
      rep(
        0,
        length(tip_names)
      ),
    
    excluded_index =
      excluded_index,
    
    tau =
      1
  )
  
  names(
    prior_tip_posterior
  ) <- tip_names
  
  prior_target_posterior <- aggregate_target_posterior(
    tip_posterior =
      prior_tip_posterior,
    
    target_tip_mask =
      target_tip_mask
  )
  
  prior_true_class_posterior <- unname(
    prior_target_posterior[
      true_target
    ]
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
    posterior =
      prior_target_posterior,
    
    true_class =
      true_target
  )
  
  
  # -----------------------------------------------------------------------
  # Preconstruct the query for every eligible window
  # -----------------------------------------------------------------------
  #
  # The observed query does not depend on lambda or tau.
  
  taxon_window_queries <- vector(
    "list",
    length(eligible_windows)
  )
  
  names(
    taxon_window_queries
  ) <- candidate_windows$window_id[
    eligible_windows
  ]
  
  for (w in seq_along(eligible_windows)) {
    
    window_index <- eligible_windows[[w]]
    
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
    
    names(
      query
    ) <- as.character(
      sites
    )
    
    query <- query[
      query != 0L
    ]
    
    taxon_window_queries[[w]] <- query
  }
  
  
  # -----------------------------------------------------------------------
  # lambda loop
  # -----------------------------------------------------------------------
  
  for (lambda in lambda_values) {
    
    message(
      "  lambda = ",
      lambda
    )
    
    allele_model <- build_hierarchical_allele_model(
      classifier =
        loo_classifier,
      
      lambda =
        lambda
    )
    
    
    # ---------------------------------------------------------------------
    # window loop
    # ---------------------------------------------------------------------
    
    for (w in seq_along(eligible_windows)) {
      
      window_index <- eligible_windows[[w]]
      
      query <- taxon_window_queries[[w]]
      
      has_evidence <- length(
        query
      ) > 0L
      
      
      # -------------------------------------------------------------------
      # Calculate likelihood once for this taxon/window/lambda
      # -------------------------------------------------------------------
      
      if (has_evidence) {
        
        likelihood_result <- calculate_tip_log_likelihoods(
          query =
            query,
          
          classifier =
            loo_classifier,
          
          allele_model =
            allele_model,
          
          epsilon =
            epsilon
        )
        
        log_likelihood <-
          likelihood_result$log_likelihood
        
      } else {
        
        log_likelihood <- NULL
      }
      
      
      # -------------------------------------------------------------------
      # tau loop
      # -------------------------------------------------------------------
      
      for (tau in tau_values) {
        
        if (has_evidence) {
          
          tip_posterior <- calculate_tempered_posterior(
            log_likelihood =
              log_likelihood,
            
            excluded_index =
              excluded_index,
            
            tau =
              tau
          )
          
          names(
            tip_posterior
          ) <- tip_names
          
          target_posterior <- aggregate_target_posterior(
            tip_posterior =
              tip_posterior,
            
            target_tip_mask =
              target_tip_mask
          )
          
        } else {
          
          target_posterior <-
            prior_target_posterior
        }
        
        
        # ---------------------------------------------------------------
        # Absolute posterior scores
        # ---------------------------------------------------------------
        
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
          posterior =
            target_posterior,
          
          true_class =
            true_target
        )
        
        
        # ---------------------------------------------------------------
        # Prior-relative utility scores
        # ---------------------------------------------------------------
        
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
        
        
        # ---------------------------------------------------------------
        # Store case result
        # ---------------------------------------------------------------
        
        row <- data.frame(
          window_id =
            candidate_windows$window_id[[window_index]],
          
          lambda =
            lambda,
          
          tau =
            tau,
          
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
            predicted_target ==
            true_target,
          
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
          
          stringsAsFactors =
            FALSE
        )
        
        for (class_name in class_names) {
          
          row[[paste0(
              "prior_p_",
              class_name
            )]] <- unname(
            prior_target_posterior[
              class_name
            ]
          )
          
          row[[paste0(
              "p_",
              class_name
            )]] <- unname(
            target_posterior[
              class_name
            ]
          )
        }
        
        case_results[[case_counter]] <- row
        
        case_counter <-
          case_counter +
          1L
      }
    }
  }
}

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

rownames(
  case_table
) <- NULL


# -------------------------------------------------------------------------
# Per-class utility profiles for each window x lambda x tau
# -------------------------------------------------------------------------

window_class_groups <- split(
  case_table,
  interaction(
    case_table$window_id,
    case_table$lambda,
    case_table$tau,
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
      
      lambda =
        x$lambda[[1]],
      
      tau =
        x$tau[[1]],
      
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
        mean(
          x$prior_true_class_posterior
        ),
      
      mean_true_class_posterior =
        mean(
          x$true_class_posterior
        ),
      
      mean_true_posterior_improvement =
        mean(
          x$true_posterior_improvement
        ),
      
      mean_log_score =
        mean(
          x$log_score
        ),
      
      mean_log_score_improvement =
        mean(
          x$log_score_improvement
        ),
      
      mean_brier_score =
        mean(
          x$brier_score
        ),
      
      mean_brier_score_improvement =
        mean(
          x$brier_score_improvement
        ),
      
      mean_entropy_reduction =
        mean(
          x$entropy_reduction
        ),
      
      stringsAsFactors =
        FALSE
    )
  }
)

by_class_table <- do.call(
  rbind,
  by_class_results
)

rownames(
  by_class_table
) <- NULL


# -------------------------------------------------------------------------
# Parameter-specific window summaries
# -------------------------------------------------------------------------
#
# One row per window x lambda x tau.

parameter_groups <- split(
  case_table,
  interaction(
    case_table$window_id,
    case_table$lambda,
    case_table$tau,
    drop = TRUE,
    lex.order = TRUE
  )
)

parameter_summary_results <- lapply(
  parameter_groups,
  function(x) {
    
    window_id <- x$window_id[[1]]
    lambda <- x$lambda[[1]]
    tau <- x$tau[[1]]
    
    window_row <- candidate_windows[
      candidate_windows$window_id ==
        window_id,
      ,
      drop = FALSE
    ]
    
    class_rows <- by_class_table[
      by_class_table$window_id ==
        window_id &
        by_class_table$lambda ==
        lambda &
        by_class_table$tau ==
        tau,
      ,
      drop = FALSE
    ]
    
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
      
      lambda =
        lambda,
      
      tau =
        tau,
      
      n_taxa =
        nrow(x),
      
      n_with_evidence =
        sum(x$has_evidence),
      
      evidence_rate =
        mean(x$has_evidence),
      
      mean_query_sites =
        mean(x$n_query_sites),
      
      n_classes_with_evidence =
        sum(
          class_rows$n_with_evidence > 0
        ),
      
      n_classes_with_positive_log_gain =
        sum(
          class_rows$mean_log_score_improvement >
            0
        ),
      
      empirical_accuracy =
        mean(x$correct),
      
      empirical_mean_true_posterior =
        mean(
          x$true_class_posterior
        ),
      
      empirical_true_posterior_improvement =
        mean(
          x$true_posterior_improvement
        ),
      
      empirical_log_score =
        mean(
          x$log_score
        ),
      
      empirical_log_score_improvement =
        mean(
          x$log_score_improvement
        ),
      
      empirical_brier_score =
        mean(
          x$brier_score
        ),
      
      empirical_brier_score_improvement =
        mean(
          x$brier_score_improvement
        ),
      
      empirical_entropy_reduction =
        mean(
          x$entropy_reduction
        ),
      
      balanced_accuracy =
        mean(
          class_rows$accuracy
        ),
      
      balanced_mean_true_posterior =
        mean(
          class_rows$mean_true_class_posterior
        ),
      
      balanced_true_posterior_improvement =
        mean(
          class_rows$mean_true_posterior_improvement
        ),
      
      balanced_log_score =
        mean(
          class_rows$mean_log_score
        ),
      
      balanced_log_score_improvement =
        mean(
          class_rows$mean_log_score_improvement
        ),
      
      balanced_brier_score =
        mean(
          class_rows$mean_brier_score
        ),
      
      balanced_brier_score_improvement =
        mean(
          class_rows$mean_brier_score_improvement
        ),
      
      balanced_entropy_reduction =
        mean(
          class_rows$mean_entropy_reduction
        ),
      
      stringsAsFactors =
        FALSE
    )
  }
)

parameter_summary_table <- do.call(
  rbind,
  parameter_summary_results
)

rownames(
  parameter_summary_table
) <- NULL


# -------------------------------------------------------------------------
# Robust window summary across lambda x tau
# -------------------------------------------------------------------------
#
# The grid is treated as a sensitivity analysis, not as a tuning exercise.
#
# For each window we therefore report:
#   mean   typical performance over the grid
#   min    worst-case performance
#   max    best-case performance
#   sd     parameter sensitivity
#   positive_fraction
#          proportion of grid settings for which log-score improvement > 0
#
# No window is ranked by its best-performing lambda/tau combination.

robust_groups <- split(
  parameter_summary_table,
  parameter_summary_table$window_id
)

summary_results <- lapply(
  robust_groups,
  function(x) {
    
    window_id <- x$window_id[[1]]
    
    # The following quantities do not depend on lambda/tau.
    invariant_columns <- c(
      "window_length",
      "window_start",
      "window_end",
      "n_polymorphic_sites",
      "sites",
      "n_taxa",
      "n_with_evidence",
      "evidence_rate",
      "mean_query_sites",
      "n_classes_with_evidence"
    )
    
    invariant <- x[
      1L,
      invariant_columns,
      drop = FALSE
    ]
    
    # Primary balanced log-score utility.
    balanced_log_gain <-
      x$balanced_log_score_improvement
    
    empirical_log_gain <-
      x$empirical_log_score_improvement
    
    balanced_brier_gain <-
      x$balanced_brier_score_improvement
    
    balanced_posterior_gain <-
      x$balanced_true_posterior_improvement
    
    balanced_entropy_gain <-
      x$balanced_entropy_reduction
    
    balanced_accuracy <-
      x$balanced_accuracy
    
    # Identify parameter settings corresponding to best and worst balanced
    # log-score improvement. These are descriptive only.
    best_index <- which.max(
      balanced_log_gain
    )
    
    worst_index <- which.min(
      balanced_log_gain
    )
    
    data.frame(
      window_id =
        window_id,
      
      invariant,
      
      n_parameter_combinations =
        nrow(x),
      
      mean_balanced_log_score_improvement =
        mean(
          balanced_log_gain
        ),
      
      min_balanced_log_score_improvement =
        min(
          balanced_log_gain
        ),
      
      max_balanced_log_score_improvement =
        max(
          balanced_log_gain
        ),
      
      sd_balanced_log_score_improvement =
        safe_sd(
          balanced_log_gain
        ),
      
      positive_balanced_log_gain_fraction =
        mean(
          balanced_log_gain > 0
        ),
      
      mean_empirical_log_score_improvement =
        mean(
          empirical_log_gain
        ),
      
      min_empirical_log_score_improvement =
        min(
          empirical_log_gain
        ),
      
      sd_empirical_log_score_improvement =
        safe_sd(
          empirical_log_gain
        ),
      
      mean_balanced_brier_score_improvement =
        mean(
          balanced_brier_gain
        ),
      
      min_balanced_brier_score_improvement =
        min(
          balanced_brier_gain
        ),
      
      sd_balanced_brier_score_improvement =
        safe_sd(
          balanced_brier_gain
        ),
      
      mean_balanced_true_posterior_improvement =
        mean(
          balanced_posterior_gain
        ),
      
      min_balanced_true_posterior_improvement =
        min(
          balanced_posterior_gain
        ),
      
      sd_balanced_true_posterior_improvement =
        safe_sd(
          balanced_posterior_gain
        ),
      
      mean_balanced_entropy_reduction =
        mean(
          balanced_entropy_gain
        ),
      
      min_balanced_entropy_reduction =
        min(
          balanced_entropy_gain
        ),
      
      sd_balanced_entropy_reduction =
        safe_sd(
          balanced_entropy_gain
        ),
      
      mean_balanced_accuracy =
        mean(
          balanced_accuracy
        ),
      
      min_balanced_accuracy =
        min(
          balanced_accuracy
        ),
      
      max_balanced_accuracy =
        max(
          balanced_accuracy
        ),
      
      sd_balanced_accuracy =
        safe_sd(
          balanced_accuracy
        ),
      
      best_lambda =
        x$lambda[
          best_index
        ],
      
      best_tau =
        x$tau[
          best_index
        ],
      
      best_balanced_log_score_improvement =
        balanced_log_gain[
          best_index
        ],
      
      worst_lambda =
        x$lambda[
          worst_index
        ],
      
      worst_tau =
        x$tau[
          worst_index
        ],
      
      worst_balanced_log_score_improvement =
        balanced_log_gain[
          worst_index
        ],
      
      stringsAsFactors =
        FALSE
    )
  }
)

summary_table <- do.call(
  rbind,
  summary_results
)

rownames(
  summary_table
) <- NULL


# -------------------------------------------------------------------------
# Add windows below min_snps
# -------------------------------------------------------------------------

ineligible_windows <- setdiff(
  seq_len(
    nrow(candidate_windows)
  ),
  eligible_windows
)

if (length(ineligible_windows) > 0L) {
  
  skipped <- candidate_windows[
    ineligible_windows,
    ,
    drop = FALSE
  ]
  
  skipped$n_taxa <- length(
    tip_names
  )
  
  skipped$n_with_evidence <- 0L
  skipped$evidence_rate <- 0
  skipped$mean_query_sites <- 0
  skipped$n_classes_with_evidence <- 0L
  skipped$n_parameter_combinations <-
    length(lambda_values) *
    length(tau_values)
  
  missing_columns <- setdiff(
    names(summary_table),
    names(skipped)
  )
  
  for (column in missing_columns) {
    skipped[[column]] <- NA_real_
  }
  
  skipped <- skipped[
    ,
    names(summary_table),
    drop = FALSE
  ]
  
  summary_table <- rbind(
    summary_table,
    skipped
  )
}


# -------------------------------------------------------------------------
# Rank robust window summary
# -------------------------------------------------------------------------

if (!default_ranking %in% names(summary_table)) {
  stop(
    "Unknown default ranking metric: ",
    default_ranking,
    call. = FALSE
  )
}

ranking_value <- summary_table[[default_ranking]]

summary_table$rank <- rank(
  -ranking_value,
  ties.method = "min",
  na.last = "keep"
)

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

rownames(
  summary_table
) <- NULL


# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------

case_output_path <- file.path(
  processed_dir,
  "amplicon_scan_cases.csv"
)

by_class_output_path <- file.path(
  processed_dir,
  "amplicon_scan_by_class.csv"
)

parameter_summary_output_path <- file.path(
  processed_dir,
  "amplicon_scan_parameter_summary.csv"
)

summary_output_path <- file.path(
  processed_dir,
  "amplicon_scan_summary.csv"
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
  by_class_table,
  by_class_output_path,
  row.names = FALSE
)

write.csv(
  parameter_summary_table,
  parameter_summary_output_path,
  row.names = FALSE
)

write.csv(
  summary_table,
  summary_output_path,
  row.names = FALSE
)

saveRDS(
  list(
    epsilon =
      epsilon,
    
    lambda_values =
      lambda_values,
    
    tau_values =
      tau_values,
    
    parameter_grid =
      expand.grid(
        lambda =
          lambda_values,
        
        tau =
          tau_values,
        
        KEEP.OUT.ATTRS =
          FALSE,
        
        stringsAsFactors =
          FALSE
      ),
    
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
    
    default_ranking =
      default_ranking,
    
    windows =
      candidate_windows,
    
    cases =
      case_table,
    
    by_class =
      by_class_table,
    
    parameter_summary =
      parameter_summary_table,
    
    summary =
      summary_table
  ),
  rds_output_path
)


# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------

message("")
message(
  "Amplicon sensitivity scan complete"
)

message(
  "  epsilon: ",
  epsilon
)

message(
  "  lambda values: ",
  paste(
    lambda_values,
    collapse = ", "
  )
)

message(
  "  tau values: ",
  paste(
    tau_values,
    collapse = ", "
  )
)

message(
  "  parameter combinations: ",
  length(lambda_values) *
    length(tau_values)
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
  by_class_output_path
)

message(
  "  wrote: ",
  parameter_summary_output_path
)

message(
  "  wrote: ",
  summary_output_path
)

message(
  "  wrote: ",
  rds_output_path
)
