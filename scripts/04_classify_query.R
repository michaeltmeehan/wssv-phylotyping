#!/usr/bin/env Rscript

# Stage 04: classify a query sequence using a hierarchical allele model.
#
# Input:
#   data/processed/classifier.rds
#
# The classifier uses all reference tips as mutually exclusive terminal states.
# Allele probabilities for each terminal state are estimated hierarchically
# through the reference tree:
#
#   theta[v,s,a] =
#       (n[v,s,a] + lambda * theta[parent(v),s,a]) /
#       (n[v,s]   + lambda)
#
# where:
#   n[v,s,a] = number of descendant taxa at node v carrying allele a at site s
#   n[v,s]   = number of descendant taxa at node v called at site s
#
# lambda controls shrinkage toward the parent:
#
#   lambda = 0:
#       reproduces the original deterministic-tip classifier for called
#       reference states.
#
#   lambda > 0:
#       progressively borrows allele-frequency information from ancestors.
#
# Query observations are linked to the latent reference allele distribution
# through a symmetric nucleotide-error model:
#
#   P(Y=b | theta) =
#       (1 - epsilon) * theta[b] +
#       (epsilon / 3) * (1 - theta[b])
#
# Missing query states contribute no evidence.
#
# For lambda > 0, a missing reference state at a tip inherits its allele
# distribution from its ancestors.
#
# For lambda = 0, a missing reference state contributes no evidence, matching
# the behaviour of the original Stage 04 classifier.

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
# Validate classifier
# -------------------------------------------------------------------------

required_components <- c(
  "tree",
  "terminal_states",
  "tip_alleles",
  "node_tip_mask",
  "node_metadata",
  "polymorphic_sites",
  "root_node"
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
# Build hierarchical allele probabilities
# -------------------------------------------------------------------------

build_hierarchical_allele_model <- function(
    classifier,
    lambda = 0
) {
  
  if (
    !is.numeric(lambda) ||
    length(lambda) != 1L ||
    is.na(lambda) ||
    lambda < 0
  ) {
    stop(
      "'lambda' must be a single non-negative number.",
      call. = FALSE
    )
  }
  
  tree <- classifier$tree
  tip_alleles <- classifier$tip_alleles
  node_tip_mask <- classifier$node_tip_mask
  
  n_tips <- nrow(tip_alleles)
  n_sites <- ncol(tip_alleles)
  n_nodes <- nrow(node_tip_mask)
  
  site_names <- colnames(tip_alleles)
  
  allele_names <- c("A", "C", "G", "T")
  
  
  # -----------------------------------------------------------------------
  # Count allele observations at every node
  # -----------------------------------------------------------------------
  
  allele_counts <- array(
    0,
    dim = c(
      n_nodes,
      n_sites,
      4L
    ),
    dimnames = list(
      node = rownames(node_tip_mask),
      site = site_names,
      allele = allele_names
    )
  )
  
  called_counts <- matrix(
    0L,
    nrow = n_nodes,
    ncol = n_sites,
    dimnames = list(
      rownames(node_tip_mask),
      site_names
    )
  )
  
  for (node_index in seq_len(n_nodes)) {
    
    descendant_tips <- node_tip_mask[
      node_index,
      ,
      drop = TRUE
    ]
    
    x <- tip_alleles[
      descendant_tips,
      ,
      drop = FALSE
    ]
    
    called_counts[node_index, ] <-
      colSums(x != 0L)
    
    allele_counts[node_index, , "A"] <-
      colSums(x == 1L)
    
    allele_counts[node_index, , "C"] <-
      colSums(x == 2L)
    
    allele_counts[node_index, , "G"] <-
      colSums(x == 3L)
    
    allele_counts[node_index, , "T"] <-
      colSums(x == 4L)
  }
  
  
  # -----------------------------------------------------------------------
  # Construct parent lookup
  # -----------------------------------------------------------------------
  
  parent <- rep(
    NA_integer_,
    n_nodes
  )
  
  names(parent) <- as.character(
    seq_len(n_nodes)
  )
  
  for (i in seq_len(nrow(tree$edge))) {
    
    parent[
      as.character(tree$edge[i, 2])
    ] <- tree$edge[i, 1]
  }
  
  
  # -----------------------------------------------------------------------
  # Order nodes from root toward tips
  # -----------------------------------------------------------------------
  #
  # Hierarchical probabilities must be calculated only after the parent's
  # probabilities are available.
  #
  
  children <- split(
    tree$edge[, 2],
    tree$edge[, 1]
  )
  
  traversal_order <- integer(0)
  
  queue <- classifier$root_node
  
  while (length(queue) > 0L) {
    
    node <- queue[[1]]
    
    queue <- queue[-1]
    
    traversal_order <- c(
      traversal_order,
      node
    )
    
    node_children <- children[[as.character(node)]]
    
    if (!is.null(node_children)) {
      queue <- c(
        queue,
        node_children
      )
    }
  }
  
  if (length(traversal_order) != n_nodes) {
    stop(
      "Could not construct a complete root-to-tip traversal.",
      call. = FALSE
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Allocate hierarchical allele probabilities
  # -----------------------------------------------------------------------
  
  theta <- array(
    NA_real_,
    dim = c(
      n_nodes,
      n_sites,
      4L
    ),
    dimnames = list(
      node = rownames(node_tip_mask),
      site = site_names,
      allele = allele_names
    )
  )
  
  
  # -----------------------------------------------------------------------
  # Root probabilities
  # -----------------------------------------------------------------------
  #
  # At the root there is no parent from which to borrow information.
  #
  # Use the empirical allele distribution among all called reference taxa.
  #
  # If an entire polymorphic site somehow has no called states, use a uniform
  # distribution. This should not ordinarily occur.
  #
  
  root <- classifier$root_node
  
  root_called <- called_counts[
    as.character(root),
  ]
  
  for (a in seq_len(4L)) {
    
    root_count <- allele_counts[
      as.character(root),
      ,
      a
    ]
    
    theta[
      as.character(root),
      ,
      a
    ] <- ifelse(
      root_called > 0L,
      root_count / root_called,
      0.25
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Recursively estimate child probabilities
  # -----------------------------------------------------------------------
  
  non_root_nodes <- traversal_order[
    traversal_order != root
  ]
  
  for (node in non_root_nodes) {
    
    node_name <- as.character(node)
    parent_name <- as.character(
      parent[node_name]
    )
    
    n_called <- called_counts[
      node_name,
    ]
    
    if (lambda == 0) {
      
      # Exact empirical distribution.
      #
      # Sites with no observations remain NA. At terminal tips this
      # reproduces the original deterministic reference model for called
      # states and the original "no evidence" treatment for missing states.
      
      observed <- n_called > 0L
      
      for (a in seq_len(4L)) {
        
        theta[
          node_name,
          observed,
          a
        ] <-
          allele_counts[
            node_name,
            observed,
            a
          ] /
          n_called[observed]
      }
      
    } else {
      
      # Hierarchically smoothed distribution.
      #
      # If n_called = 0:
      #
      #   theta_child = theta_parent
      #
      # exactly, so missing reference information is inherited naturally.
      
      denominator <- n_called + lambda
      
      for (a in seq_len(4L)) {
        
        theta[
          node_name,
          ,
          a
        ] <-
          (
            allele_counts[
              node_name,
              ,
              a
            ] +
              lambda *
              theta[
                parent_name,
                ,
                a
              ]
          ) /
          denominator
      }
    }
  }
  
  
  # -----------------------------------------------------------------------
  # Extract terminal-state distributions
  # -----------------------------------------------------------------------
  
  tip_theta <- theta[
    as.character(seq_len(n_tips)),
    ,
    ,
    drop = FALSE
  ]
  
  dimnames(tip_theta)[[1]] <-
    classifier$terminal_states$taxon
  
  
  # -----------------------------------------------------------------------
  # Validate probability distributions
  # -----------------------------------------------------------------------
  
  theta_sum <- apply(
    tip_theta,
    c(1, 2),
    sum
  )
  
  non_missing <- !is.na(
    theta_sum
  )
  
  if (
    any(
      abs(
        theta_sum[non_missing] - 1
      ) > 1e-10
    )
  ) {
    stop(
      "Hierarchical terminal allele probabilities do not sum to 1.",
      call. = FALSE
    )
  }
  
  
  list(
    lambda = lambda,
    theta = theta,
    tip_theta = tip_theta,
    allele_counts = allele_counts,
    called_counts = called_counts,
    parent = parent,
    traversal_order = traversal_order
  )
}


# -------------------------------------------------------------------------
# Core likelihood function
# -------------------------------------------------------------------------

calculate_tip_log_likelihoods <- function(
    query,
    classifier,
    allele_model,
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
  
  
  # -----------------------------------------------------------------------
  # Remove missing query states
  # -----------------------------------------------------------------------
  
  keep <- query != 0L
  
  query <- query[keep]
  query_sites <- query_sites[keep]
  
  if (length(query) == 0L) {
    stop(
      "Query contains no observed A/C/G/T SNP calls.",
      call. = FALSE
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Restrict to sites represented by the classifier
  # -----------------------------------------------------------------------
  
  valid_sites <- classifier$polymorphic_sites
  
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
    dimnames(allele_model$tip_theta)[[2]]
  )
  
  if (anyNA(site_columns)) {
    stop(
      "Could not match one or more query sites to classifier columns.",
      call. = FALSE
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Calculate likelihood
  # -----------------------------------------------------------------------
  
  tip_theta <- allele_model$tip_theta
  
  n_tips <- dim(tip_theta)[1]
  
  log_likelihood <- numeric(n_tips)
  n_informative <- integer(n_tips)
  
  for (j in seq_along(query)) {
    
    query_state <- query[[j]]
    
    theta_query <- tip_theta[
      ,
      site_columns[[j]],
      query_state
    ]
    
    # NA occurs only when lambda = 0 and that terminal reference has no
    # observed allele at this site.
    informative <- !is.na(
      theta_query
    )
    
    emission_probability <-
      (1 - epsilon) * theta_query +
      (epsilon / 3) * (1 - theta_query)
    
    log_likelihood[informative] <-
      log_likelihood[informative] +
      log(
        emission_probability[informative]
      )
    
    n_informative <-
      n_informative +
      informative
  }
  
  names(log_likelihood) <-
    classifier$terminal_states$taxon
  
  names(n_informative) <-
    classifier$terminal_states$taxon
  
  
  list(
    log_likelihood = log_likelihood,
    n_informative = n_informative,
    query_sites_used = query_sites
  )
}


# -------------------------------------------------------------------------
# Convert terminal likelihoods to posterior probabilities
# -------------------------------------------------------------------------

calculate_tip_posteriors <- function(
    query,
    classifier,
    allele_model,
    epsilon = 0.01,
    prior = NULL
) {
  
  likelihood_result <- calculate_tip_log_likelihoods(
    query = query,
    classifier = classifier,
    allele_model = allele_model,
    epsilon = epsilon
  )
  
  log_likelihood <-
    likelihood_result$log_likelihood
  
  n_tips <- length(
    log_likelihood
  )
  
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
  
  log_posterior <- rep(
    -Inf,
    n_tips
  )
  
  positive_prior <- prior > 0
  
  log_posterior[positive_prior] <-
    log(
      prior[positive_prior]
    ) +
    log_likelihood[positive_prior]
  
  max_log_posterior <- max(
    log_posterior
  )
  
  posterior <- exp(
    log_posterior -
      max_log_posterior
  )
  
  posterior <-
    posterior /
    sum(posterior)
  
  data.frame(
    state_id =
      classifier$terminal_states$state_id,
    
    node_id =
      classifier$terminal_states$node_id,
    
    taxon =
      classifier$terminal_states$taxon,
    
    prior =
      prior,
    
    log_likelihood =
      unname(log_likelihood),
    
    posterior =
      unname(posterior),
    
    n_informative =
      unname(
        likelihood_result$n_informative
      ),
    
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
  
  if (
    !all(
      c(
        "taxon",
        "posterior"
      ) %in%
      names(tip_posteriors)
    )
  ) {
    stop(
      "'tip_posteriors' must contain columns 'taxon' and 'posterior'.",
      call. = FALSE
    )
  }
  
  tip_order <- colnames(
    classifier$node_tip_mask
  )
  
  posterior <-
    tip_posteriors$posterior[
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
    classifier$node_tip_mask %*%
      posterior
  )
  
  metadata <-
    classifier$node_metadata
  
  metadata$posterior <-
    node_posterior[
      match(
        metadata$node_id,
        as.integer(
          names(node_posterior)
        )
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
    allele_model,
    epsilon = 0.01,
    prior = NULL
) {
  
  tip_posteriors <-
    calculate_tip_posteriors(
      query = query,
      classifier = classifier,
      allele_model = allele_model,
      epsilon = epsilon,
      prior = prior
    )
  
  node_posteriors <-
    calculate_node_posteriors(
      tip_posteriors = tip_posteriors,
      classifier = classifier
    )
  
  
  # -----------------------------------------------------------------------
  # Validate posterior conservation
  # -----------------------------------------------------------------------
  
  if (
    abs(
      sum(
        tip_posteriors$posterior
      ) - 1
    ) > 1e-12
  ) {
    stop(
      "Terminal posterior probabilities do not sum to 1.",
      call. = FALSE
    )
  }
  
  root_posterior <-
    node_posteriors$posterior[
      node_posteriors$node_id ==
        classifier$root_node
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
# Example model construction
# -------------------------------------------------------------------------
#
# Original classifier:
#
# allele_model <- build_hierarchical_allele_model(
#   classifier,
#   lambda = 0
# )
#
# Hierarchically smoothed classifier:
#
# allele_model <- build_hierarchical_allele_model(
#   classifier,
#   lambda = 1
# )
#
#
# Example query:
#
# query <- c(
#   "6518"   = 1L,
#   "54397"  = 3L,
#   "56235"  = 4L,
#   "108500" = 4L
# )
#
# result <- classify_query(
#   query = query,
#   classifier = classifier,
#   allele_model = allele_model,
#   epsilon = 0.01
# )


message("Stage 04 hierarchical classifier functions loaded")
message("  terminal states: ", nrow(classifier$terminal_states))
message(
  "  polymorphic sites available: ",
  length(classifier$polymorphic_sites)
)