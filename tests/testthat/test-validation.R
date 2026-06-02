source(test_path("../../R/validation.R"))

toy_alignment <- function() {
  x <- matrix(c(
    1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L,
    1L, 2L, 0L, 4L, 1L, 2L, 0L, 4L,
    4L, 3L, 2L, 1L, 4L, 3L, 2L, 1L
  ), nrow = 3L, byrow = TRUE)
  rownames(x) <- c("tip_a", "tip_b", "tip_c")
  x
}

toy_scores <- function() {
  data.frame(
    node_id = c("n1", "n1", "n2"),
    site = c(2L, 5L, 7L),
    best_allele = c("C", "A", "G"),
    direction = c("clade", "clade", "outside"),
    depth = c(1, 1, 2)
  )
}

toy_target_mask <- function() {
  x <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE), nrow = 2L, byrow = TRUE)
  rownames(x) <- c("n1", "n2")
  colnames(x) <- c("tip_a", "tip_b", "tip_c")
  x
}

test_that("selected-window columns are extracted in window order", {
  windows <- data.frame(window_id = c("w2", "w1"), start = c(5L, 2L), end = c(6L, 3L))
  out <- extract_window_columns(toy_alignment(), windows)

  expect_equal(ncol(out), 4L)
  expect_equal(out[1L, ], c(1L, 2L, 2L, 3L))
})

test_that("alignment masking keeps only selected-window regions", {
  windows <- data.frame(window_id = "w1", start = 2L, end = 3L)
  masked <- mask_alignment_to_windows(toy_alignment(), windows)

  expect_equal(masked[1L, ], c(0L, 2L, 3L, 0L, 0L, 0L, 0L, 0L))
  expect_equal(masked[, 2:3], toy_alignment()[, 2:3])
})

test_that("artificial fragments are created within alignment bounds", {
  fragments <- make_random_fragments(100L, c(10L, 30L), n_per_length = 2L, seed = 2L)

  expect_equal(nrow(fragments), 4L)
  expect_true(all(fragments$start >= 1L))
  expect_true(all(fragments$end <= 100L))
  expect_equal(sort(unique(fragments$width)), c(10L, 30L))
})

test_that("informative sites and nodes are counted for partial sequences", {
  masked <- mask_alignment_to_windows(toy_alignment(), data.frame(window_id = "w", start = 2L, end = 5L))
  counts <- count_observed_signal(masked, toy_scores(), toy_target_mask(), 1L, 1L)

  expect_equal(counts$observed_informative_sites[counts$tip == "tip_a"], 2L)
  expect_equal(counts$observed_informative_nodes[counts$tip == "tip_a"], 1L)
  expect_equal(counts$supported_true_nodes[counts$tip == "tip_a"], 1L)
  expect_true(counts$resolved_signal[counts$tip == "tip_a"])
})

test_that("panel coverage is summarised per tip", {
  windows <- data.frame(window_id = c("w1", "w2"), start = c(2L, 5L), end = c(3L, 7L), width = c(2L, 3L))
  summary <- evaluate_panel_signal(toy_alignment(), windows, toy_scores(), toy_target_mask(), panel_size = 2L)

  expect_equal(nrow(summary), 3L)
  expect_equal(unique(summary$panel_size), 2L)
  expect_true("insufficient_information" %in% names(summary))
})

test_that("selected and baseline panels compare on toy data", {
  selected <- data.frame(window_id = "selected", start = 2L, end = 5L, width = 4L)
  baseline <- data.frame(window_id = "baseline", start = 1L, end = 1L, width = 1L)

  selected_summary <- summarise_validation_by_panel(evaluate_panel_signal(toy_alignment(), selected, toy_scores(), toy_target_mask()))
  baseline_summary <- summarise_validation_by_panel(evaluate_panel_signal(toy_alignment(), baseline, toy_scores(), toy_target_mask(), method = "baseline"))

  expect_gt(selected_summary$mean_observed_informative_sites, baseline_summary$mean_observed_informative_sites)
})
