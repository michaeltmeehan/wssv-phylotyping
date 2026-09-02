#!/usr/bin/env Rscript

# Stage 04: classify a query sequence from an arbitrary subset of SNP calls.
#
# Input:
#   data/processed/classifier.rds
#
# Output:
#   Functions for calculating:
#     - posterior probability for each reference tip/terminal state
#     - posterior probability for every node in the tree
#
# The classifier is target-agnostic. Any subset of polymorphic sites may be
# supplied. Missing/unobserved sites contribute no evidence.
#
# Baseline likelihood model:
#   For each observed query allele:
#
#     P(query allele | reference tip)
#
#   is determined by a symmetric nucleotide error model with error probability
#   epsilon:
#
#     match:       1 - epsilon
#     mismatch:    epsilon / 3
#
# Reference states encoded as 0 are treated as unavailable/missing and do not
# contribute evidence for that tip at that site.
#
# This is deliberately a simple baseline classifier. It does not yet model
# within-clade allele frequencies, SNP dependence, phylogenetic uncertainty,
# or calibrated likelihood tempering.

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
# Validate classifier
# -------------------------------------------------------------------------

required_components <- c(
  "terminal_states",
  "tip_alleles",
  "node_tip_mask",
  "node_metadata",
  "polymorphic_sites"
)

missing_components <- setdiff(
  required_components,
  names(classifier)
)

if (length(missing_components) > 0L) {
  stop(
    "Classifier object is missing required components: ",
    paste(missing_components, collapse = ", "),
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Core likelihood function
# -------------------------------------------------------------------------

calculate_tip_log_likelihoods <- function(
    query,
    classifier,
    epsilon = 0.01
) {
  
  if (
    !is.numeric(epsilon) ||
    length(epsilon) != 1L ||
    is.na(epsilon) ||
    epsilon <= 0 ||
    epsilon >= 0.75
  ) {
    stop(
      "'epsilon' must be a single number between 0 and 0.75.",
      call. = FALSE
    )
  }
  
  tip_alleles <- classifier$tip_alleles
  valid_sites <- classifier$polymorphic_sites
  
  # Query must be a named integer/numeric vector:
  #
  #   names(query) = genomic/alignment site indices
  #   values        = 1:A, 2:C, 3:G, 4:T, 0:missing
  #
  # Example:
  #
  #   query <- c(
  #     "6518" = 1L,
  #     "54397" = 3L,
  #     "108500" = 4L
  #   )
  
  if (is.null(names(query))) {
    stop(
      "'query' must be a named vector with site indices as names.",
      call. = FALSE
    )
  }
  
  query_sites <- suppressWarnings(
    as.integer(names(query))
  )
  
  if (anyNA(query_sites)) {
    stop(
      "All query names must be valid integer site indices.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(query_sites)) {
    stop(
      "Query contains duplicated site indices.",
      call. = FALSE
    )
  }
  
  if (any(!query %in% 0:4)) {
    stop(
      "Query states must use the alignment encoding 0, 1, 2, 3, 4.",
      call. = FALSE
    )
  }
  
  # Ignore query sites that are missing.
  keep <- query != 0L
  
  query <- query[keep]
  query_sites <- query_sites[keep]
  
  if (length(query) == 0L) {
    stop(
      "Query contains no observed A/C/G/T SNP calls.",
      call. = FALSE
    )
  }
  
  # Restrict to SNPs represented in the trained classifier.
  in_classifier <- query_sites %in% valid_sites
  
  if (!any(in_classifier)) {
    stop(
      "None of the observed query sites are present in the classifier.",
      call. = FALSE
    )
  }
  
  query <- query[in_classifier]
  query_sites <- query_sites[in_classifier]
  
  site_columns <- match(
    as.character(query_sites),
    colnames(tip_alleles)
  )
  
  if (anyNA(site_columns)) {
    stop(
      "Could not match one or more query sites to classifier columns.",
      call. = FALSE
    )
  }
  
  reference <- tip_alleles[
    ,
    site_columns,
    drop = FALSE
  ]
  
  n_tips <- nrow(reference)
  
  log_likelihood <- numeric(n_tips)
  n_informative <- integer(n_tips)
  n_matches <- integer(n_tips)
  n_mismatches <- integer(n_tips)
  
  log_match <- log(1 - epsilon)
  log_mismatch <- log(epsilon / 3)
  
  for (j in seq_along(query)) {
    
    query_state <- query[[j]]
    reference_state <- reference[, j]
    
    # A reference state of 0 is missing/unavailable and contributes no
    # evidence for that particular terminal state.
    called <- reference_state != 0L
    
    matched <- called & reference_state == query_state
    mismatched <- called & reference_state != query_state
    
    log_likelihood[matched] <-
      log_likelihood[matched] + log_match
    
    log_likelihood[mismatched] <-
      log_likelihood[mismatched] + log_mismatch
    
    n_informative <- n_informative + called
    n_matches <- n_matches + matched
    n_mismatches <- n_mismatches + mismatched
  }
  
  names(log_likelihood) <- rownames(reference)
  names(n_informative) <- rownames(reference)
  names(n_matches) <- rownames(reference)
  names(n_mismatches) <- rownames(reference)
  
  list(
    log_likelihood = log_likelihood,
    n_informative = n_informative,
    n_matches = n_matches,
    n_mismatches = n_mismatches,
    query_sites_used = query_sites
  )
}


# -------------------------------------------------------------------------
# Convert terminal likelihoods to posterior probabilities
# -------------------------------------------------------------------------

calculate_tip_posteriors <- function(
    query,
    classifier,
    epsilon = 0.01,
    prior = NULL
) {
  
  likelihood_result <- calculate_tip_log_likelihoods(
    query = query,
    classifier = classifier,
    epsilon = epsilon
  )
  
  log_likelihood <- likelihood_result$log_likelihood
  
  n_tips <- length(log_likelihood)
  
  if (is.null(prior)) {
    prior <- classifier$terminal_states$prior
  }
  
  if (
    !is.numeric(prior) ||
    length(prior) != n_tips ||
    anyNA(prior) ||
    any(prior < 0)
  ) {
    stop(
      "'prior' must contain one non-negative value per terminal state.",
      call. = FALSE
    )
  }
  
  if (sum(prior) <= 0) {
    stop(
      "'prior' must have positive total probability.",
      call. = FALSE
    )
  }
  
  prior <- prior / sum(prior)
  
  # log posterior up to a normalising constant
  log_posterior <- log(prior) + log_likelihood
  
  # Stable log-sum-exp normalisation
  max_log_posterior <- max(log_posterior)
  
  posterior <- exp(
    log_posterior - max_log_posterior
  )
  
  posterior <- posterior / sum(posterior)
  
  names(posterior) <- classifier$terminal_states$taxon
  
  data.frame(
    state_id = classifier$terminal_states$state_id,
    node_id = classifier$terminal_states$node_id,
    taxon = classifier$terminal_states$taxon,
    prior = prior,
    log_likelihood = unname(log_likelihood),
    posterior = unname(posterior),
    n_informative = unname(likelihood_result$n_informative),
    n_matches = unname(likelihood_result$n_matches),
    n_mismatches = unname(likelihood_result$n_mismatches),
    stringsAsFactors = FALSE
  )
}


# -------------------------------------------------------------------------
# Propagate terminal posterior probabilities to all nodes
# -------------------------------------------------------------------------

calculate_node_posteriors <- function(
    tip_posteriors,
    classifier
) {
  
  if (!all(c("taxon", "posterior") %in% names(tip_posteriors))) {
    stop(
      "'tip_posteriors' must contain columns 'taxon' and 'posterior'.",
      call. = FALSE
    )
  }
  
  tip_order <- colnames(
    classifier$node_tip_mask
  )
  
  posterior <- tip_posteriors$posterior[
    match(
      tip_order,
      tip_posteriors$taxon
    )
  ]
  
  if (anyNA(posterior)) {
    stop(
      "Could not match all terminal posterior probabilities to tree tips.",
      call. = FALSE
    )
  }
  
  node_posterior <- drop(
    classifier$node_tip_mask %*% posterior
  )
  
  metadata <- classifier$node_metadata
  
  metadata$posterior <- node_posterior[
    match(
      metadata$node_id,
      as.integer(names(node_posterior))
    )
  ]
  
  metadata
}


# -------------------------------------------------------------------------
# Full classification wrapper
# -------------------------------------------------------------------------

classify_query <- function(
    query,
    classifier,
    epsilon = 0.01,
    prior = NULL
) {
  
  tip_posteriors <- calculate_tip_posteriors(
    query = query,
    classifier = classifier,
    epsilon = epsilon,
    prior = prior
  )
  
  node_posteriors <- calculate_node_posteriors(
    tip_posteriors = tip_posteriors,
    classifier = classifier
  )
  
  # Sanity checks
  
  if (abs(sum(tip_posteriors$posterior) - 1) > 1e-12) {
    stop(
      "Terminal posterior probabilities do not sum to 1.",
      call. = FALSE
    )
  }
  
  root_posterior <- node_posteriors$posterior[
    node_posteriors$node_id == classifier$root_node
  ]
  
  if (
    length(root_posterior) != 1L ||
    abs(root_posterior - 1) > 1e-12
  ) {
    stop(
      "Root posterior probability is not equal to 1.",
      call. = FALSE
    )
  }
  
  list(
    tip_posteriors = tip_posteriors,
    node_posteriors = node_posteriors
  )
}


# -------------------------------------------------------------------------
# Example usage
# -------------------------------------------------------------------------
#
# Query states use:
#
#   A = 1
#   C = 2
#   G = 3
#   T = 4
#   missing = 0
#
# For example:
#
# query <- c(
#   "6518" = 1L,
#   "54397" = 3L,
#   "56235" = 4L,
#   "108500" = 4L
# )
#
# result <- classify_query(
#   query = query,
#   classifier = classifier,
#   epsilon = 0.01
# )
#
# head(
#   result$tip_posteriors[
#     order(-result$tip_posteriors$posterior),
#   ],
#   10
# )
#
# head(
#   result$node_posteriors[
#     order(-result$node_posteriors$posterior),
#   ],
#   20
# )


message("Stage 04 classifier functions loaded")
message("  terminal states: ", nrow(classifier$terminal_states))
message("  polymorphic sites available: ", length(classifier$polymorphic_sites))