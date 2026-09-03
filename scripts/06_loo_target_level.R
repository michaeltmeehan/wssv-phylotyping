#!/usr/bin/env Rscript

# Stage 06: leave-one-out validation at the configured target-clade resolution.
#
# Inputs:
#   data/processed/precomputed.rds
#   data/processed/classifier.rds
#   scripts/04_classify_query.R
#
# Outputs:
#   data/processed/target_validation_cases.csv
#   data/processed/target_validation_summary.csv
#   data/processed/target_validation_calibration.csv
#   data/processed/target_validation.rds
#
# For each reference taxon:
#   1. Hold its sequence out of the hierarchical allele model.
#   2. Calculate terminal-state log-likelihoods from all observed SNPs.
#   3. Apply several likelihood-tempering values tau.
#   4. Normalise posterior probability across retained terminal states.
#   5. Aggregate terminal posterior probabilities into the configured
#      target clades plus an "Other" class.
#   6. Compare the resulting probabilities with the held-out taxon's
#      known target-clade membership.
#
# The target clades play NO role in training the classifier. They are used
# only here to aggregate the already-computed phylogeny-wide posterior.

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "Package 'yaml' is required to read config/config.yml.",
    call. = FALSE
  )
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed


# -------------------------------------------------------------------------
# Load objects and functions
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

lambda <- 20

tau_values <- c(
  0.01,
  0.03,
  0.05,
  0.10,
  0.20,
  0.50,
  1.00
)

epsilon <- 0.01


# -------------------------------------------------------------------------
# Identify configured target nodes
# -------------------------------------------------------------------------

target_node_ids <- precomputed$target_node_metadata$node_id

if (length(target_node_ids) == 0L) {
  stop(
    "No target clades found in precomputed data.",
    call. = FALSE
  )
}

target_node_ids <- as.integer(
  target_node_ids
)

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


# Target clades must be mutually exclusive if they are to form a
# multiclass probability partition with "Other".

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
    "Configured target clades overlap. Multiclass target validation ",
    "requires mutually exclusive target nodes.\nOverlapping taxa: ",
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
  
  if (
    !is.numeric(tau) ||
    length(tau) != 1L ||
    is.na(tau) ||
    tau <= 0
  ) {
    stop(
      "'tau' must be a single positive number.",
      call. = FALSE
    )
  }
  
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
  
  # Protect against tiny floating-point errors.
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
# Leave-one-out validation
# -------------------------------------------------------------------------

case_results <- vector(
  "list",
  length(tip_names) * length(tau_values)
)

calibration_results <- vector(
  "list",
  length(tip_names) *
    length(tau_values) *
    length(class_names)
)

case_counter <- 1L
calibration_counter <- 1L


for (i in seq_along(tip_names)) {
  
  held_out_taxon <- tip_names[[i]]
  
  message(
    "Preparing ",
    i,
    "/",
    length(tip_names),
    ": ",
    held_out_taxon
  )
  
  
  # -----------------------------------------------------------------------
  # Query from original held-out sequence
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
  # Remove held-out sequence from training data
  # -----------------------------------------------------------------------
  
  loo_classifier <- classifier
  
  loo_classifier$tip_alleles[
    held_out_taxon,
  ] <- 0L
  
  
  # -----------------------------------------------------------------------
  # Build hierarchical allele model
  # -----------------------------------------------------------------------
  
  allele_model <- build_hierarchical_allele_model(
    classifier = loo_classifier,
    lambda = lambda
  )
  
  
  # -----------------------------------------------------------------------
  # Calculate terminal log-likelihoods once
  # -----------------------------------------------------------------------
  
  likelihood_result <- calculate_tip_log_likelihoods(
    query = query,
    classifier = loo_classifier,
    allele_model = allele_model,
    epsilon = epsilon
  )
  
  log_likelihood <- likelihood_result$log_likelihood
  
  excluded_index <- match(
    held_out_taxon,
    tip_names
  )
  
  
  # -----------------------------------------------------------------------
  # Evaluate each tempering value
  # -----------------------------------------------------------------------
  
  for (tau in tau_values) {
    
    tip_posterior <- calculate_tempered_posterior(
      log_likelihood = log_likelihood,
      excluded_index = excluded_index,
      tau = tau
    )
    
    names(tip_posterior) <- tip_names
    
    
    # ---------------------------------------------------------------------
    # Aggregate to target-clade resolution
    # ---------------------------------------------------------------------
    
    target_posterior <- aggregate_target_posterior(
      tip_posterior = tip_posterior,
      target_tip_mask = target_tip_mask
    )
    
    true_target <- true_class[
      held_out_taxon
    ]
    
    predicted_target <- names(
      which.max(
        target_posterior
      )
    )
    
    true_target_posterior <-
      target_posterior[
        true_target
      ]
    
    predicted_target_posterior <-
      max(
        target_posterior
      )
    
    
    # ---------------------------------------------------------------------
    # Store case-level result
    # ---------------------------------------------------------------------
    
    row <- data.frame(
      held_out_taxon = held_out_taxon,
      true_class = true_target,
      predicted_class = predicted_target,
      correct = predicted_target == true_target,
      
      lambda = lambda,
      tau = tau,
      epsilon = epsilon,
      
      n_query_sites = length(query),
      
      true_class_posterior =
        unname(
          true_target_posterior
        ),
      
      predicted_class_posterior =
        unname(
          predicted_target_posterior
        ),
      
      stringsAsFactors = FALSE
    )
    
    for (class_name in class_names) {
      row[[paste0("p_", class_name)]] <- target_posterior[class_name]
    }
    
    case_results[[case_counter]] <- row
    
    case_counter <- case_counter + 1L
    
    
    # ---------------------------------------------------------------------
    # Store long-form probabilities for calibration assessment
    # ---------------------------------------------------------------------
    
    for (class_name in class_names) {
      
      calibration_results[[calibration_counter]] <- data.frame(
        held_out_taxon =
          held_out_taxon,
        
        lambda =
          lambda,
        
        tau =
          tau,
        
        class =
          class_name,
        
        probability =
          unname(
            target_posterior[
              class_name
            ]
          ),
        
        observed =
          as.integer(
            true_target ==
              class_name
          ),
        
        stringsAsFactors = FALSE
      )
      
      calibration_counter <-
        calibration_counter + 1L
    }
  }
}


# -------------------------------------------------------------------------
# Assemble results
# -------------------------------------------------------------------------

case_table <- do.call(
  rbind,
  case_results
)

rownames(case_table) <- NULL

calibration_table <- do.call(
  rbind,
  calibration_results
)

rownames(calibration_table) <- NULL


# -------------------------------------------------------------------------
# Summary metrics by tau
# -------------------------------------------------------------------------

summary_results <- lapply(
  tau_values,
  function(tau) {
    
    x <- case_table[
      case_table$tau == tau,
      ,
      drop = FALSE
    ]
    
    # Multiclass Brier score:
    #
    # mean_i sum_k (p_ik - y_ik)^2
    
    probability_columns <- paste0(
      "p_",
      class_names
    )
    
    probability_matrix <- as.matrix(
      x[
        ,
        probability_columns,
        drop = FALSE
      ]
    )
    
    observed_matrix <- matrix(
      0,
      nrow = nrow(x),
      ncol = length(class_names)
    )
    
    colnames(observed_matrix) <-
      class_names
    
    for (i in seq_len(nrow(x))) {
      
      observed_matrix[
        i,
        x$true_class[[i]]
      ] <- 1
    }
    
    brier <- mean(
      rowSums(
        (
          probability_matrix -
            observed_matrix
        )^2
      )
    )
    
    # Multiclass log loss.
    true_probability <- pmax(
      x$true_class_posterior,
      .Machine$double.eps
    )
    
    log_loss <- -mean(
      log(
        true_probability
      )
    )
    
    data.frame(
      lambda = lambda,
      tau = tau,
      
      n = nrow(x),
      
      accuracy = mean(
        x$correct
      ),
      
      n_correct = sum(
        x$correct
      ),
      
      mean_true_class_posterior =
        mean(
          x$true_class_posterior
        ),
      
      median_true_class_posterior =
        median(
          x$true_class_posterior
        ),
      
      mean_predicted_class_posterior =
        mean(
          x$predicted_class_posterior
        ),
      
      log_loss = log_loss,
      
      brier_score = brier,
      
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
# Calibration bins
# -------------------------------------------------------------------------

breaks <- seq(
  0,
  1,
  by = 0.1
)

calibration_table$bin <- cut(
  calibration_table$probability,
  breaks = breaks,
  include.lowest = TRUE,
  right = TRUE
)

calibration_bins <- aggregate(
  cbind(
    probability,
    observed
  ) ~ tau + bin,
  data = calibration_table,
  FUN = mean
)

calibration_n <- aggregate(
  observed ~ tau + bin,
  data = calibration_table,
  FUN = length
)

names(calibration_n)[
  names(calibration_n) == "observed"
] <- "n"

calibration_bins <- merge(
  calibration_bins,
  calibration_n,
  by = c(
    "tau",
    "bin"
  ),
  all.x = TRUE,
  sort = FALSE
)

names(calibration_bins)[
  names(calibration_bins) == "probability"
] <- "mean_predicted_probability"

names(calibration_bins)[
  names(calibration_bins) == "observed"
] <- "observed_frequency"


# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------

case_output_path <- file.path(
  processed_dir,
  "target_validation_cases.csv"
)

summary_output_path <- file.path(
  processed_dir,
  "target_validation_summary.csv"
)

calibration_output_path <- file.path(
  processed_dir,
  "target_validation_calibration.csv"
)

rds_output_path <- file.path(
  processed_dir,
  "target_validation.rds"
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
  calibration_bins,
  calibration_output_path,
  row.names = FALSE
)

saveRDS(
  list(
    lambda = lambda,
    tau_values = tau_values,
    epsilon = epsilon,
    target_node_ids = target_node_ids,
    classes = class_names,
    cases = case_table,
    summary = summary_table,
    calibration_long = calibration_table,
    calibration_bins = calibration_bins
  ),
  rds_output_path
)


# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------

message("")
message("Target-clade validation complete")
message(
  "  lambda: ",
  lambda
)

message(
  "  targets: ",
  paste(
    target_names,
    collapse = ", "
  )
)

message(
  "  classes: ",
  paste(
    class_names,
    collapse = ", "
  )
)

message("")
print(
  summary_table,
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
  calibration_output_path
)

message(
  "  wrote: ",
  rds_output_path
)