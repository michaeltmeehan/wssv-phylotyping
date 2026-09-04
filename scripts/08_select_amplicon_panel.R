#!/usr/bin/env Rscript

# Stage 08: greedy forward selection of a panel of fixed-length amplicons.
#
# Fixed settings:
#   window length = 1000 bp
#   lambda = 1
#   tau = 1
#   epsilon = 0.01
#
# At each step, select the remaining non-overlapping window that maximizes the
# balanced mean log-score improvement across the configured target clades.
#
# The joint panel likelihood is formed by summing per-window Naive Bayes
# log-likelihood contributions. Because selected windows are not allowed to
# overlap, this is exactly equivalent to classifying on the union of SNPs
# across the selected amplicons.
#
# Outputs:
#   data/processed/panel_selection_history.csv
#   data/processed/panel_candidate_evaluations.csv
#   data/processed/panel_by_class.csv
#   data/processed/panel_cases.csv
#   data/processed/panel_selection.rds

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/config.yml.", call. = FALSE)
}

config <- yaml::read_yaml("config/config.yml")
processed_dir <- config$paths$processed

precomputed <- readRDS(file.path(processed_dir, "precomputed.rds"))
classifier <- readRDS(file.path(processed_dir, "classifier.rds"))

source("scripts/04_classify_query.R")


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

window_length <- 1000L
window_step <- 500L
min_snps <- 1L

lambda <- 1
tau <- 1
epsilon <- 0.01

max_panel_size <- 10L
min_marginal_gain <- 0

allow_overlap <- FALSE

# "Other" remains a classification outcome, but does not contribute to the
# primary panel-selection objective unless this is set to TRUE.
include_other_in_objective <- FALSE

probability_floor <- 1e-15


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0L) return(0)
  -sum(p * log(p))
}

calculate_tempered_posterior <- function(log_likelihood, excluded_index, tau) {
  n_states <- length(log_likelihood)

  prior <- rep(1 / (n_states - 1L), n_states)
  prior[excluded_index] <- 0

  retained <- prior > 0

  log_posterior <- rep(-Inf, n_states)
  log_posterior[retained] <-
    log(prior[retained]) +
    tau * log_likelihood[retained]

  max_log_posterior <- max(log_posterior[retained])

  posterior <- numeric(n_states)
  posterior[retained] <- exp(
    log_posterior[retained] -
      max_log_posterior
  )

  posterior / sum(posterior)
}

aggregate_target_posterior <- function(tip_posterior, target_tip_mask) {
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
      "Target posterior probabilities exceed 1. Check target overlap.",
      call. = FALSE
    )
  }

  posterior <- c(
    target_posterior,
    Other = other_posterior
  )

  posterior / sum(posterior)
}

multiclass_brier <- function(posterior, true_class) {
  observed <- setNames(
    rep(0, length(posterior)),
    names(posterior)
  )

  observed[true_class] <- 1

  sum(
    (posterior - observed)^2
  )
}

windows_overlap <- function(start_1, end_1, start_2, end_2) {
  max(start_1, start_2) <= min(end_1, end_2)
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

target_names <- paste0("node_", target_node_ids)
names(target_node_ids) <- target_names

target_tip_mask <- classifier$node_tip_mask[
  as.character(target_node_ids),
  ,
  drop = FALSE
]

rownames(target_tip_mask) <- target_names

if (any(colSums(target_tip_mask) > 1L)) {
  stop(
    "Configured target clades overlap. Panel scoring requires exclusive targets.",
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

objective_classes <- if (include_other_in_objective) {
  class_names
} else {
  target_names
}

message(
  "Objective classes: ",
  paste(objective_classes, collapse = ", ")
)


# -------------------------------------------------------------------------
# Candidate windows
# -------------------------------------------------------------------------

alignment_length <- ncol(
  precomputed$aln_int
)

starts <- seq.int(
  from = 1L,
  to = alignment_length - window_length + 1L,
  by = window_step
)

final_start <- alignment_length - window_length + 1L

if (tail(starts, 1L) != final_start) {
  starts <- c(starts, final_start)
}

candidate_windows <- data.frame(
  window_start = starts,
  window_end = starts + window_length - 1L,
  window_length = window_length,
  stringsAsFactors = FALSE
)

candidate_windows$window_id <- sprintf(
  "L%d_%d_%d",
  window_length,
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

polymorphic_sites <- as.integer(
  classifier$polymorphic_sites
)

all_window_sites <- lapply(
  seq_len(nrow(candidate_windows)),
  function(i) {
    polymorphic_sites[
      polymorphic_sites >= candidate_windows$window_start[[i]] &
        polymorphic_sites <= candidate_windows$window_end[[i]]
    ]
  }
)

candidate_windows$n_polymorphic_sites <- vapply(
  all_window_sites,
  length,
  integer(1)
)

candidate_windows$sites <- vapply(
  all_window_sites,
  function(x) paste(x, collapse = ";"),
  character(1)
)

keep <- candidate_windows$n_polymorphic_sites >= min_snps

candidate_windows <- candidate_windows[
  keep,
  ,
  drop = FALSE
]

window_sites <- all_window_sites[keep]

rownames(candidate_windows) <- NULL
names(window_sites) <- candidate_windows$window_id

if (nrow(candidate_windows) == 0L) {
  stop("No eligible candidate windows.", call. = FALSE)
}


# -------------------------------------------------------------------------
# Leave-one-out priors
# -------------------------------------------------------------------------

n_taxa <- length(tip_names)
n_states <- n_taxa

prior_target_posterior <- vector(
  "list",
  n_taxa
)

prior_true_class_posterior <- numeric(
  n_taxa
)

prior_log_score <- numeric(
  n_taxa
)

prior_brier_score <- numeric(
  n_taxa
)

prior_entropy <- numeric(
  n_taxa
)

for (i in seq_along(tip_names)) {
  prior_tip <- calculate_tempered_posterior(
    log_likelihood = rep(0, n_states),
    excluded_index = i,
    tau = 1
  )

  names(prior_tip) <- tip_names

  prior_target <- aggregate_target_posterior(
    tip_posterior = prior_tip,
    target_tip_mask = target_tip_mask
  )

  true_target <- unname(
    true_class[tip_names[[i]]]
  )

  prior_target_posterior[[i]] <- prior_target

  prior_true_class_posterior[i] <- unname(
    prior_target[true_target]
  )

  prior_log_score[i] <- log(
    max(
      prior_true_class_posterior[i],
      probability_floor
    )
  )

  prior_brier_score[i] <- multiclass_brier(
    posterior = prior_target,
    true_class = true_target
  )

  prior_entropy[i] <- entropy(
    prior_target
  )
}


# -------------------------------------------------------------------------
# Precompute per-window LOO log-likelihood contributions
# -------------------------------------------------------------------------

message("")
message("Precomputing window likelihood contributions...")

window_log_likelihoods <- vector(
  "list",
  n_taxa
)

window_query_site_counts <- matrix(
  0L,
  nrow = n_taxa,
  ncol = nrow(candidate_windows),
  dimnames = list(
    tip_names,
    candidate_windows$window_id
  )
)

for (i in seq_along(tip_names)) {
  held_out_taxon <- tip_names[[i]]

  message(
    "  taxon ",
    i,
    "/",
    n_taxa,
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

  taxon_ll <- vector(
    "list",
    nrow(candidate_windows)
  )

  names(taxon_ll) <- candidate_windows$window_id

  for (w in seq_len(nrow(candidate_windows))) {
    sites <- window_sites[[candidate_windows$window_id[[w]]]]

    site_columns <- match(
      as.character(sites),
      names(query_states)
    )

    if (anyNA(site_columns)) {
      stop(
        "Could not match window SNPs to classifier site columns.",
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

    window_query_site_counts[
      i,
      w
    ] <- length(query)

    if (length(query) == 0L) {
      taxon_ll[[w]] <- rep(
        0,
        n_states
      )
    } else {
      likelihood_result <- calculate_tip_log_likelihoods(
        query = query,
        classifier = loo_classifier,
        allele_model = allele_model,
        epsilon = epsilon
      )

      taxon_ll[[w]] <- likelihood_result$log_likelihood
    }
  }

  window_log_likelihoods[[i]] <- taxon_ll
}


# -------------------------------------------------------------------------
# Panel evaluation
# -------------------------------------------------------------------------

evaluate_panel <- function(
    current_log_likelihoods,
    candidate_index = NULL,
    panel_step,
    selected_window_ids
) {
  candidate_window_id <- if (is.null(candidate_index)) {
    NA_character_
  } else {
    candidate_windows$window_id[[candidate_index]]
  }

  proposed_window_ids <- selected_window_ids

  if (!is.null(candidate_index)) {
    proposed_window_ids <- c(
      proposed_window_ids,
      candidate_window_id
    )
  }

  case_rows <- vector(
    "list",
    n_taxa
  )

  for (i in seq_along(tip_names)) {
    panel_log_likelihood <- current_log_likelihoods[[i]]

    if (!is.null(candidate_index)) {
      panel_log_likelihood <-
        panel_log_likelihood +
        window_log_likelihoods[[i]][[candidate_index]]
    }

    tip_posterior <- calculate_tempered_posterior(
      log_likelihood = panel_log_likelihood,
      excluded_index = i,
      tau = tau
    )

    names(tip_posterior) <- tip_names

    target_posterior <- aggregate_target_posterior(
      tip_posterior = tip_posterior,
      target_tip_mask = target_tip_mask
    )

    true_target <- unname(
      true_class[tip_names[[i]]]
    )

    predicted_target <- names(
      which.max(target_posterior)
    )

    true_target_posterior <- unname(
      target_posterior[true_target]
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

    posterior_entropy <- entropy(
      target_posterior
    )

    n_query_sites <- if (length(proposed_window_ids) == 0L) {
      0L
    } else {
      sum(
        window_query_site_counts[
          i,
          proposed_window_ids,
          drop = TRUE
        ]
      )
    }

    row <- data.frame(
      panel_step = panel_step,
      held_out_taxon = tip_names[[i]],
      true_class = true_target,
      n_query_sites = n_query_sites,
      predicted_class = predicted_target,
      correct = predicted_target == true_target,
      prior_true_class_posterior =
        prior_true_class_posterior[i],
      true_class_posterior =
        true_target_posterior,
      true_posterior_improvement =
        true_target_posterior -
        prior_true_class_posterior[i],
      log_score = log_score,
      log_score_improvement =
        log_score -
        prior_log_score[i],
      brier_score = brier_score,
      brier_score_improvement =
        prior_brier_score[i] -
        brier_score,
      posterior_entropy = posterior_entropy,
      entropy_reduction =
        prior_entropy[i] -
        posterior_entropy,
      stringsAsFactors = FALSE
    )

    for (class_name in class_names) {
      row[[paste0("p_", class_name)]] <- unname(
        target_posterior[class_name]
      )
    }

    case_rows[[i]] <- row
  }

  cases <- do.call(
    rbind,
    case_rows
  )

  rownames(cases) <- NULL

  class_groups <- split(
    cases,
    cases$true_class
  )

  by_class_rows <- lapply(
    class_groups,
    function(x) {
      data.frame(
        panel_step = panel_step,
        true_class = x$true_class[[1]],
        n_taxa = nrow(x),
        mean_query_sites =
          mean(x$n_query_sites),
        accuracy =
          mean(x$correct),
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

  by_class <- do.call(
    rbind,
    by_class_rows
  )

  rownames(by_class) <- NULL

  objective_rows <- by_class[
    by_class$true_class %in% objective_classes,
    ,
    drop = FALSE
  ]

  if (nrow(objective_rows) != length(objective_classes)) {
    stop(
      "Not all objective classes are represented.",
      call. = FALSE
    )
  }

  summary <- data.frame(
    panel_step = panel_step,
    candidate_window_id = candidate_window_id,
    panel_size = length(proposed_window_ids),
    balanced_log_score_improvement =
      mean(objective_rows$mean_log_score_improvement),
    balanced_true_posterior_improvement =
      mean(objective_rows$mean_true_posterior_improvement),
    balanced_brier_score_improvement =
      mean(objective_rows$mean_brier_score_improvement),
    balanced_accuracy =
      mean(objective_rows$accuracy),
    balanced_entropy_reduction =
      mean(objective_rows$mean_entropy_reduction),
    empirical_log_score_improvement =
      mean(
        cases$log_score_improvement[
          cases$true_class %in% objective_classes
        ]
      ),
    empirical_accuracy =
      mean(
        cases$correct[
          cases$true_class %in% objective_classes
        ]
      ),
    mean_query_sites =
      mean(
        cases$n_query_sites[
          cases$true_class %in% objective_classes
        ]
      ),
    min_target_log_score_improvement =
      min(objective_rows$mean_log_score_improvement),
    max_target_log_score_improvement =
      max(objective_rows$mean_log_score_improvement),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    by_class = by_class,
    cases = cases
  )
}


# -------------------------------------------------------------------------
# Empty-panel baseline
# -------------------------------------------------------------------------

current_log_likelihoods <- lapply(
  seq_along(tip_names),
  function(i) rep(0, n_states)
)

selected_indices <- integer(0)
selected_window_ids <- character(0)

history_rows <- list()
candidate_evaluation_rows <- list()
selected_case_rows <- list()
selected_by_class_rows <- list()

candidate_eval_counter <- 1L

baseline <- evaluate_panel(
  current_log_likelihoods = current_log_likelihoods,
  candidate_index = NULL,
  panel_step = 0L,
  selected_window_ids = selected_window_ids
)

current_utility <-
  baseline$summary$balanced_log_score_improvement

message("")
message(
  "Empty-panel balanced log-score improvement: ",
  signif(current_utility, 5)
)


# -------------------------------------------------------------------------
# Greedy forward selection
# -------------------------------------------------------------------------

for (step in seq_len(max_panel_size)) {
  message("")
  message(
    "Greedy panel step ",
    step
  )

  remaining_indices <- setdiff(
    seq_len(nrow(candidate_windows)),
    selected_indices
  )

  if (
    !allow_overlap &&
    length(selected_indices) > 0L
  ) {
    non_overlapping <- vapply(
      remaining_indices,
      function(candidate_index) {
        !any(
          vapply(
            selected_indices,
            function(selected_index) {
              windows_overlap(
                candidate_windows$window_start[[candidate_index]],
                candidate_windows$window_end[[candidate_index]],
                candidate_windows$window_start[[selected_index]],
                candidate_windows$window_end[[selected_index]]
              )
            },
            logical(1)
          )
        )
      },
      logical(1)
    )

    remaining_indices <- remaining_indices[
      non_overlapping
    ]
  }

  if (length(remaining_indices) == 0L) {
    message("No eligible non-overlapping candidates remain.")
    break
  }

  step_results <- vector(
    "list",
    length(remaining_indices)
  )

  for (j in seq_along(remaining_indices)) {
    candidate_index <- remaining_indices[[j]]

    evaluation <- evaluate_panel(
      current_log_likelihoods = current_log_likelihoods,
      candidate_index = candidate_index,
      panel_step = step,
      selected_window_ids = selected_window_ids
    )

    candidate_summary <- evaluation$summary

    candidate_summary$window_start <-
      candidate_windows$window_start[[candidate_index]]

    candidate_summary$window_end <-
      candidate_windows$window_end[[candidate_index]]

    candidate_summary$n_polymorphic_sites <-
      candidate_windows$n_polymorphic_sites[[candidate_index]]

    candidate_summary$marginal_log_score_gain <-
      candidate_summary$balanced_log_score_improvement -
      current_utility

    target_rows <- evaluation$by_class[
      evaluation$by_class$true_class %in%
        objective_classes,
      ,
      drop = FALSE
    ]

    for (class_name in objective_classes) {
      class_row <- target_rows[
        target_rows$true_class == class_name,
        ,
        drop = FALSE
      ]

      candidate_summary[[paste0("utility_", class_name)]] <- class_row$mean_log_score_improvement
    }

    step_results[[j]] <- list(
      summary = candidate_summary,
      by_class = evaluation$by_class,
      cases = evaluation$cases,
      candidate_index = candidate_index
    )

    candidate_evaluation_rows[[candidate_eval_counter]] <- candidate_summary

    candidate_eval_counter <-
      candidate_eval_counter +
      1L
  }

  utilities <- vapply(
    step_results,
    function(x) {
      x$summary$balanced_log_score_improvement
    },
    numeric(1)
  )

  best_position <- which.max(
    utilities
  )

  best_result <- step_results[[best_position]]

  best_index <- best_result$candidate_index
  best_summary <- best_result$summary

  best_marginal_gain <-
    best_summary$marginal_log_score_gain

  if (
    !is.finite(best_marginal_gain) ||
    best_marginal_gain <= min_marginal_gain
  ) {
    message(
      "Stopping: best marginal gain = ",
      signif(best_marginal_gain, 5),
      " <= min_marginal_gain = ",
      min_marginal_gain
    )
    break
  }

  selected_indices <- c(
    selected_indices,
    best_index
  )

  selected_window_ids <- c(
    selected_window_ids,
    candidate_windows$window_id[[best_index]]
  )

  for (i in seq_along(tip_names)) {
    current_log_likelihoods[[i]] <-
      current_log_likelihoods[[i]] +
      window_log_likelihoods[[i]][[best_index]]
  }

  current_utility <-
    best_summary$balanced_log_score_improvement

  history_row <- best_summary

  history_row$selected_window_id <-
    candidate_windows$window_id[[best_index]]

  history_row$selected_window_start <-
    candidate_windows$window_start[[best_index]]

  history_row$selected_window_end <-
    candidate_windows$window_end[[best_index]]

  history_row$selected_n_polymorphic_sites <-
    candidate_windows$n_polymorphic_sites[[best_index]]

  history_row$panel_window_ids <- paste(
    selected_window_ids,
    collapse = ";"
  )

  history_rows[[length(history_rows) + 1L]] <- history_row

  selected_by_class <- best_result$by_class
  selected_by_class$panel_window_ids <- paste(
    selected_window_ids,
    collapse = ";"
  )

  selected_by_class_rows[[length(selected_by_class_rows) + 1L]] <- selected_by_class

  selected_cases <- best_result$cases
  selected_cases$panel_window_ids <- paste(
    selected_window_ids,
    collapse = ";"
  )

  selected_case_rows[[length(selected_case_rows) + 1L]] <- selected_cases

  message(
    "  selected: ",
    candidate_windows$window_id[[best_index]]
  )

  message(
    "  marginal balanced log-score gain: ",
    signif(best_marginal_gain, 5)
  )

  message(
    "  total balanced log-score improvement: ",
    signif(current_utility, 5)
  )

  message(
    "  balanced accuracy: ",
    signif(best_summary$balanced_accuracy, 4)
  )

  target_utility_text <- paste(
    vapply(
      objective_classes,
      function(class_name) {
        paste0(
          class_name,
          "=",
          signif(
            best_summary[[paste0("utility_", class_name)]],
            4
          )
        )
      },
      character(1)
    ),
    collapse = ", "
  )

  message(
    "  target utilities: ",
    target_utility_text
  )
}


# -------------------------------------------------------------------------
# Assemble outputs
# -------------------------------------------------------------------------

history_table <- if (length(history_rows) > 0L) {
  do.call(rbind, history_rows)
} else {
  data.frame()
}

candidate_evaluation_table <- if (
  length(candidate_evaluation_rows) > 0L
) {
  do.call(rbind, candidate_evaluation_rows)
} else {
  data.frame()
}

panel_by_class_table <- if (
  length(selected_by_class_rows) > 0L
) {
  do.call(rbind, selected_by_class_rows)
} else {
  data.frame()
}

panel_cases_table <- if (
  length(selected_case_rows) > 0L
) {
  do.call(rbind, selected_case_rows)
} else {
  data.frame()
}


# -------------------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------------------

history_output_path <- file.path(
  processed_dir,
  "panel_selection_history.csv"
)

candidate_output_path <- file.path(
  processed_dir,
  "panel_candidate_evaluations.csv"
)

by_class_output_path <- file.path(
  processed_dir,
  "panel_by_class.csv"
)

cases_output_path <- file.path(
  processed_dir,
  "panel_cases.csv"
)

rds_output_path <- file.path(
  processed_dir,
  "panel_selection.rds"
)

write.csv(
  history_table,
  history_output_path,
  row.names = FALSE
)

write.csv(
  candidate_evaluation_table,
  candidate_output_path,
  row.names = FALSE
)

write.csv(
  panel_by_class_table,
  by_class_output_path,
  row.names = FALSE
)

write.csv(
  panel_cases_table,
  cases_output_path,
  row.names = FALSE
)

saveRDS(
  list(
    window_length = window_length,
    window_step = window_step,
    min_snps = min_snps,
    lambda = lambda,
    tau = tau,
    epsilon = epsilon,
    max_panel_size = max_panel_size,
    min_marginal_gain = min_marginal_gain,
    allow_overlap = allow_overlap,
    include_other_in_objective = include_other_in_objective,
    target_node_ids = target_node_ids,
    objective_classes = objective_classes,
    candidate_windows = candidate_windows,
    selected_indices = selected_indices,
    selected_window_ids = selected_window_ids,
    history = history_table,
    candidate_evaluations = candidate_evaluation_table,
    by_class = panel_by_class_table,
    cases = panel_cases_table
  ),
  rds_output_path
)


# -------------------------------------------------------------------------
# Print final panel
# -------------------------------------------------------------------------

message("")
message("Greedy panel selection complete")
message(
  "  selected amplicons: ",
  length(selected_window_ids)
)

if (length(selected_window_ids) > 0L) {
  message("")
  message("Selected panel:")

  print(
    history_table[
      ,
      c(
        "panel_step",
        "selected_window_id",
        "selected_window_start",
        "selected_window_end",
        "selected_n_polymorphic_sites",
        "marginal_log_score_gain",
        "balanced_log_score_improvement",
        "balanced_accuracy"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}

message("")
message("  wrote: ", history_output_path)
message("  wrote: ", candidate_output_path)
message("  wrote: ", by_class_output_path)
message("  wrote: ", cases_output_path)
message("  wrote: ", rds_output_path)
