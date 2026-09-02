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
#   1. Retain its observed SNP states as the query.
#   2. Remove its sequence information from the training data.
#   3. Rebuild the hierarchical allele model.
#   4. Exclude the held-out taxon from the terminal posterior distribution.
#   5. Classify the query using all observed polymorphic sites.
#   6. Propagate terminal posterior mass to every internal node.
#   7. Record posterior support for the held-out taxon's true ancestors.
#
# The tree topology itself is retained during leave-one-out validation.
# Only the held-out taxon's nucleotide observations are removed.
#
# For lambda = 0, this reduces to the deterministic-tip classifier.
# For lambda > 0, terminal allele distributions borrow information
# hierarchically from ancestral nodes.

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "Package 'yaml' is required to read config/config.yml.",
    call. = FALSE
  )
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

classifier_path <- file.path(
  processed_dir,
  "classifier.rds"
)

if (!file.exists(classifier_path)) {
  stop(
    "Classifier structure not found: ",
    classifier_path,
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

lambda <- 20
epsilon <- 0.01

n_tips <- nrow(
  classifier$terminal_states
)

tip_names <- classifier$terminal_states$taxon


# -------------------------------------------------------------------------
# Helper: classify while excluding one terminal state
# -------------------------------------------------------------------------

classify_query_excluding_tip <- function(
    query,
    classifier,
    allele_model,
    excluded_taxon,
    epsilon = 0.01
) {
  
  likelihood_result <- calculate_tip_log_likelihoods(
    query = query,
    classifier = classifier,
    allele_model = allele_model,
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
  
  # Equal prior probability across all retained terminal states.
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
  
  
  # -----------------------------------------------------------------------
  # Construct query from the original held-out sequence
  # -----------------------------------------------------------------------
  
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
  
  # Missing query states provide no evidence.
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
  
  
  # -----------------------------------------------------------------------
  # Remove held-out nucleotide data from the training classifier
  # -----------------------------------------------------------------------
  #
  # The topology and terminal-state mapping remain unchanged.
  #
  # Encoding the held-out training row entirely as missing means it
  # contributes zero allele counts at every ancestral node.
  #
  
  loo_classifier <- classifier
  
  loo_classifier$tip_alleles[
    held_out_taxon,
  ] <- 0L
  
  
  # -----------------------------------------------------------------------
  # Rebuild hierarchical allele model without held-out sequence
  # -----------------------------------------------------------------------
  
  allele_model <- build_hierarchical_allele_model(
    classifier = loo_classifier,
    lambda = lambda
  )
  
  
  # -----------------------------------------------------------------------
  # Classify held-out sequence
  # -----------------------------------------------------------------------
  
  classification <- classify_query_excluding_tip(
    query = query,
    classifier = loo_classifier,
    allele_model = allele_model,
    excluded_taxon = held_out_taxon,
    epsilon = epsilon
  )
  
  
  # -----------------------------------------------------------------------
  # Best remaining reference terminal state
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
  # Identify true ancestors of held-out tip
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
  
  # Remove the terminal node itself because that state was deliberately
  # assigned prior probability zero.
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
  # Find smallest true ancestral node with posterior >= 0.5
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
  # Posterior assigned to the best true ancestral node
  # -----------------------------------------------------------------------
  
  best_true_ancestor_index <- which.max(
    ancestor_posteriors$posterior
  )
  
  best_true_ancestor <-
    ancestor_posteriors[
      best_true_ancestor_index,
      ,
      drop = FALSE
    ]
  
  
  # -----------------------------------------------------------------------
  # Store summary
  # -----------------------------------------------------------------------
  
  results[[i]] <- data.frame(
    held_out_taxon = held_out_taxon,
    
    n_query_sites = length(query),
    
    best_reference_taxon = best_tip$taxon,
    best_reference_posterior = best_tip$posterior,
    best_reference_log_likelihood = best_tip$log_likelihood,
    
    parent_node = parent_node,
    parent_clade_size = parent_size,
    parent_posterior = parent_posterior,
    
    smallest_supported_node = smallest_supported_node,
    smallest_supported_clade_size = smallest_supported_size,
    smallest_supported_posterior = smallest_supported_posterior,
    
    best_true_ancestor_node =
      best_true_ancestor$node_id,
    
    best_true_ancestor_clade_size =
      best_true_ancestor$clade_size,
    
    best_true_ancestor_posterior =
      best_true_ancestor$posterior,
    
    stringsAsFactors = FALSE
  )
  
  
  # -----------------------------------------------------------------------
  # Retain full results for RDS
  # -----------------------------------------------------------------------
  
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
      
      attr(x, "tip_posteriors") <- NULL
      attr(x, "node_posteriors") <- NULL
      attr(x, "ancestor_posteriors") <- NULL
      
      x
    }
  )
)

rownames(summary_table) <- NULL


# -------------------------------------------------------------------------
# Construct richer validation object
# -------------------------------------------------------------------------

validation <- list(
  lambda = lambda,
  epsilon = epsilon,
  summary = summary_table,
  
  cases = lapply(
    seq_along(results),
    function(i) {
      
      list(
        held_out_taxon =
          summary_table$held_out_taxon[[i]],
        
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
message("  lambda: ", lambda)
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
  "  median smallest supported clade size: ",
  median(
    summary_table$smallest_supported_clade_size,
    na.rm = TRUE
  )
)

message(
  "  wrote: ",
  csv_output_path
)

message(
  "  wrote: ",
  rds_output_path
)