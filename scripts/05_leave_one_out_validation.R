#!/usr/bin/env Rscript

# Stage 05: leave-one-out validation of the phylogeny-wide classifier.
#
# Inputs:
#   data/processed/classifier.rds
#   scripts/04_classify_query.R
#
# Outputs:
#   data/processed/leave_one_out_validation.rds
#   data/processed/leave_one_out_validation.csv
#
# For each reference taxon:
#   1. Treat that taxon's observed SNP states as the query.
#   2. Remove the held-out taxon from the candidate terminal states.
#   3. Calculate posterior probability across the remaining reference taxa.
#   4. Propagate terminal posterior mass to all internal nodes.
#   5. Record the posterior probability assigned to every ancestor of the
#      held-out taxon.
#
# Because the held-out taxon is removed from the candidate set, successful
# validation does not mean identifying the exact taxon. Instead, it measures
# whether posterior mass is concentrated within the correct ancestral clades.
#
# This stage uses all polymorphic SNPs available for each held-out taxon.
# Missing states in the held-out query contribute no evidence.

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

classifier_path <- file.path(
  processed_dir,
  "classifier.rds"
)

if (!file.exists(classifier_path)) {
  stop(
    "Classifier structure not found: ", classifier_path,
    "\nRun stage 03 first.",
    call. = FALSE
  )
}

classifier <- readRDS(classifier_path)


# -------------------------------------------------------------------------
# Load classification functions
# -------------------------------------------------------------------------

source("scripts/04_classify_query.R")


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

epsilon <- 0.01

n_tips <- nrow(classifier$terminal_states)
tip_names <- classifier$terminal_states$taxon


# -------------------------------------------------------------------------
# Helper: classify while excluding one terminal state
# -------------------------------------------------------------------------

classify_query_excluding_tip <- function(
    query,
    classifier,
    excluded_taxon,
    epsilon = 0.01
) {
  
  likelihood_result <- calculate_tip_log_likelihoods(
    query = query,
    classifier = classifier,
    epsilon = epsilon
  )
  
  log_likelihood <- likelihood_result$log_likelihood
  
  excluded_index <- match(
    excluded_taxon,
    classifier$terminal_states$taxon
  )
  
  if (is.na(excluded_index)) {
    stop(
      "Excluded taxon not found in classifier: ",
      excluded_taxon,
      call. = FALSE
    )
  }
  
  # Equal prior probability across all remaining reference taxa.
  prior <- rep(
    1 / (length(log_likelihood) - 1L),
    length(log_likelihood)
  )
  
  prior[excluded_index] <- 0
  
  log_posterior <- rep(
    -Inf,
    length(log_likelihood)
  )
  
  retained <- prior > 0
  
  log_posterior[retained] <-
    log(prior[retained]) +
    log_likelihood[retained]
  
  max_log_posterior <- max(
    log_posterior[retained]
  )
  
  posterior <- numeric(
    length(log_likelihood)
  )
  
  posterior[retained] <- exp(
    log_posterior[retained] -
      max_log_posterior
  )
  
  posterior <- posterior / sum(posterior)
  
  tip_posteriors <- data.frame(
    state_id = classifier$terminal_states$state_id,
    node_id = classifier$terminal_states$node_id,
    taxon = classifier$terminal_states$taxon,
    prior = prior,
    log_likelihood = unname(log_likelihood),
    posterior = posterior,
    n_informative = unname(
      likelihood_result$n_informative
    ),
    n_matches = unname(
      likelihood_result$n_matches
    ),
    n_mismatches = unname(
      likelihood_result$n_mismatches
    ),
    stringsAsFactors = FALSE
  )
  
  node_posteriors <- calculate_node_posteriors(
    tip_posteriors = tip_posteriors,
    classifier = classifier
  )
  
  list(
    tip_posteriors = tip_posteriors,
    node_posteriors = node_posteriors,
    query_sites_used = likelihood_result$query_sites_used
  )
}


# -------------------------------------------------------------------------
# Leave-one-out validation
# -------------------------------------------------------------------------

results <- vector(
  "list",
  n_tips
)

for (i in seq_len(n_tips)) {
  
  held_out_taxon <- tip_names[[i]]
  
  message(
    "Validating ",
    i,
    "/",
    n_tips,
    ": ",
    held_out_taxon
  )
  
  # Use every observed polymorphic site in the held-out sequence.
  query_states <- classifier$tip_alleles[
    held_out_taxon,
    ,
    drop = TRUE
  ]
  
  query <- as.integer(
    query_states
  )
  
  names(query) <- colnames(
    classifier$tip_alleles
  )
  
  # Remove missing states.
  query <- query[
    query != 0L
  ]
  
  if (length(query) == 0L) {
    stop(
      "Held-out taxon has no observed polymorphic SNPs: ",
      held_out_taxon,
      call. = FALSE
    )
  }
  
  classification <- classify_query_excluding_tip(
    query = query,
    classifier = classifier,
    excluded_taxon = held_out_taxon,
    epsilon = epsilon
  )
  
  
  # -----------------------------------------------------------------------
  # Best remaining reference taxon
  # -----------------------------------------------------------------------
  
  tip_result <- classification$tip_posteriors
  
  retained_tip_result <- tip_result[
    tip_result$taxon != held_out_taxon,
    ,
    drop = FALSE
  ]
  
  best_tip_index <- which.max(
    retained_tip_result$posterior
  )
  
  best_tip <- retained_tip_result[
    best_tip_index,
    ,
    drop = FALSE
  ]
  
  
  # -----------------------------------------------------------------------
  # Identify true ancestors of the held-out tip
  # -----------------------------------------------------------------------
  
  true_ancestor_mask <-
    classifier$node_tip_mask[
      ,
      held_out_taxon
    ]
  
  true_ancestor_nodes <- as.integer(
    rownames(
      classifier$node_tip_mask
    )[true_ancestor_mask]
  )
  
  # Remove the held-out terminal node itself because its posterior has
  # deliberately been set to zero.
  true_internal_ancestors <- setdiff(
    true_ancestor_nodes,
    i
  )
  
  node_result <- classification$node_posteriors
  
  ancestor_posteriors <- node_result[
    node_result$node_id %in% true_internal_ancestors,
    c(
      "node_id",
      "clade_size",
      "posterior"
    ),
    drop = FALSE
  ]
  
  ancestor_posteriors <- ancestor_posteriors[
    order(
      ancestor_posteriors$clade_size
    ),
    ,
    drop = FALSE
  ]
  
  
  # -----------------------------------------------------------------------
  # Find smallest ancestral node with strong posterior support
  # -----------------------------------------------------------------------
  
  supported_ancestors <- ancestor_posteriors[
    ancestor_posteriors$posterior >= 0.5,
    ,
    drop = FALSE
  ]
  
  if (nrow(supported_ancestors) > 0L) {
    
    smallest_supported <- supported_ancestors[
      which.min(
        supported_ancestors$clade_size
      ),
      ,
      drop = FALSE
    ]
    
    smallest_supported_node <-
      smallest_supported$node_id
    
    smallest_supported_size <-
      smallest_supported$clade_size
    
    smallest_supported_posterior <-
      smallest_supported$posterior
    
  } else {
    
    smallest_supported_node <- NA_integer_
    smallest_supported_size <- NA_integer_
    smallest_supported_posterior <- NA_real_
  }
  
  
  # -----------------------------------------------------------------------
  # Posterior at immediate parent
  # -----------------------------------------------------------------------
  
  parent_row <- classifier$tree$edge[
    classifier$tree$edge[, 2] == i,
    ,
    drop = FALSE
  ]
  
  if (nrow(parent_row) != 1L) {
    stop(
      "Could not identify unique parent for tip: ",
      held_out_taxon,
      call. = FALSE
    )
  }
  
  parent_node <- parent_row[1, 1]
  
  parent_posterior <- node_result$posterior[
    node_result$node_id == parent_node
  ]
  
  parent_size <- node_result$clade_size[
    node_result$node_id == parent_node
  ]
  
  
  # -----------------------------------------------------------------------
  # Store summary
  # -----------------------------------------------------------------------
  
  results[[i]] <- data.frame(
    held_out_taxon = held_out_taxon,
    
    n_query_sites = length(query),
    
    best_reference_taxon = best_tip$taxon,
    best_reference_posterior = best_tip$posterior,
    best_reference_matches = best_tip$n_matches,
    best_reference_mismatches = best_tip$n_mismatches,
    
    parent_node = parent_node,
    parent_clade_size = parent_size,
    parent_posterior = parent_posterior,
    
    smallest_supported_node = smallest_supported_node,
    smallest_supported_clade_size = smallest_supported_size,
    smallest_supported_posterior = smallest_supported_posterior,
    
    stringsAsFactors = FALSE
  )
  
  # Retain full posterior distributions as attributes for the RDS output.
  attr(
    results[[i]],
    "tip_posteriors"
  ) <- classification$tip_posteriors
  
  attr(
    results[[i]],
    "node_posteriors"
  ) <- classification$node_posteriors
  
  attr(
    results[[i]],
    "ancestor_posteriors"
  ) <- ancestor_posteriors
}


# -------------------------------------------------------------------------
# Assemble summary table
# -------------------------------------------------------------------------

summary_table <- do.call(
  rbind,
  lapply(
    results,
    function(x) {
      attributes(x)[
        c(
          "tip_posteriors",
          "node_posteriors",
          "ancestor_posteriors"
        )
      ] <- NULL
      
      x
    }
  )
)

rownames(summary_table) <- NULL


# -------------------------------------------------------------------------
# Construct richer validation object
# -------------------------------------------------------------------------

validation <- list(
  epsilon = epsilon,
  summary = summary_table,
  
  cases = lapply(
    seq_along(results),
    function(i) {
      
      list(
        held_out_taxon = summary_table$held_out_taxon[[i]],
        
        tip_posteriors = attr(
          results[[i]],
          "tip_posteriors"
        ),
        
        node_posteriors = attr(
          results[[i]],
          "node_posteriors"
        ),
        
        ancestor_posteriors = attr(
          results[[i]],
          "ancestor_posteriors"
        )
      )
    }
  )
)

names(validation$cases) <-
  summary_table$held_out_taxon


# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------

rds_output_path <- file.path(
  processed_dir,
  "leave_one_out_validation.rds"
)

csv_output_path <- file.path(
  processed_dir,
  "leave_one_out_validation.csv"
)

saveRDS(
  validation,
  rds_output_path
)

write.csv(
  summary_table,
  csv_output_path,
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

message("")
message("Leave-one-out validation complete")
message("  taxa tested: ", nrow(summary_table))
message("  epsilon: ", epsilon)

message(
  "  median parent posterior: ",
  sprintf(
    "%.4f",
    median(
      summary_table$parent_posterior,
      na.rm = TRUE
    )
  )
)

message(
  "  parent posterior >= 0.50: ",
  sum(
    summary_table$parent_posterior >= 0.50,
    na.rm = TRUE
  ),
  "/",
  nrow(summary_table)
)

message(
  "  parent posterior >= 0.90: ",
  sum(
    summary_table$parent_posterior >= 0.90,
    na.rm = TRUE
  ),
  "/",
  nrow(summary_table)
)

message(
  "  wrote: ",
  csv_output_path
)

message(
  "  wrote: ",
  rds_output_path
)