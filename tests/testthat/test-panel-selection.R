source(test_path("../../R/panel_selection.R"))

synthetic_windows_for_panel <- function() {
  data.frame(
    window_id = c("w1", "w2", "w3", "w4"),
    window_type = "fixed",
    start = c(1L, 101L, 201L, 301L),
    end = c(100L, 200L, 300L, 400L),
    width = rep(100L, 4L),
    requested_width = rep(100L, 4L),
    informative_snps = c(3L, 3L, 2L, 1L),
    total_gain = c(3, 3, 2, 0.5),
    total_weighted_gain = c(10, 9, 6, 0.5),
    missing_fraction = c(0.01, 0.02, 0.03, 0.50),
    nodes_covered = c(2L, 2L, 2L, 1L)
  )
}

synthetic_window_node_summary_for_panel <- function() {
  data.frame(
    window_id = c("w1", "w1", "w2", "w2", "w3", "w3", "w4"),
    node_index = c(1L, 2L, 1L, 2L, 3L, 4L, 5L),
    node_id = c("n1", "n2", "n1", "n2", "n3", "n4", "n5"),
    informative_snps = 1L,
    node_site_score_entries = 1L,
    total_gain = c(6, 4, 5, 4, 3, 3, 0.5),
    total_weighted_gain = c(6, 4, 5, 4, 3, 3, 0.5),
    mean_gain = c(6, 4, 5, 4, 3, 3, 0.5),
    max_gain = c(6, 4, 5, 4, 3, 3, 0.5)
  )
}

test_that("candidate-window filtering applies configured criteria", {
  filtered <- filter_candidate_windows(
    synthetic_windows_for_panel(),
    min_informative_snps = 2L,
    min_total_weighted_gain = 5,
    max_missing_fraction = 0.1,
    permitted_widths = 100L
  )

  expect_equal(filtered$window_id, c("w1", "w2", "w3"))
})

test_that("window-node coverage aggregates node scores for retained windows", {
  node_summary <- rbind(
    synthetic_window_node_summary_for_panel(),
    synthetic_window_node_summary_for_panel()[1L, ]
  )

  coverage <- calculate_window_node_coverage(node_summary, retained_window_ids = c("w1", "w3"))

  w1_n1 <- coverage[coverage$window_id == "w1" & coverage$node_id == "n1", ]
  expect_equal(w1_n1$node_score, 12)
  expect_equal(coverage$window_id, c("w1", "w1", "w3", "w3"))
})

test_that("marginal gain credits unresolved nodes only when repeated coverage is capped", {
  coverage <- calculate_window_node_coverage(synthetic_window_node_summary_for_panel())

  gains <- calculate_marginal_gain(
    candidate_window_ids = c("w2", "w3"),
    node_coverage = coverage,
    selected_window_ids = "w1",
    cap_repeated_coverage = TRUE
  )

  expect_equal(gains$marginal_gain[gains$window_id == "w2"], 0)
  expect_equal(gains$marginal_gain[gains$window_id == "w3"], 6)
  expect_equal(gains$newly_covered_nodes[gains$window_id == "w3"], 2L)
})

test_that("greedy selection chooses complementary windows over redundant high-scoring windows", {
  windows <- synthetic_windows_for_panel()[1:3, ]
  coverage <- calculate_window_node_coverage(synthetic_window_node_summary_for_panel())

  selected <- greedy_select_panel(
    window_summary = windows,
    node_coverage = coverage,
    max_panel_size = 2L,
    cap_repeated_coverage = TRUE
  )

  expect_equal(selected$selected_panel$window_id, c("w1", "w3"))
  expect_equal(selected$selected_panel$marginal_gain, c(10, 6))
  expect_equal(selected$selected_panel$total_covered_nodes, c(2L, 4L))
})

test_that("greedy selection stops when marginal gain is below threshold", {
  windows <- synthetic_windows_for_panel()[1:3, ]
  coverage <- calculate_window_node_coverage(synthetic_window_node_summary_for_panel())

  selected <- greedy_select_panel(
    window_summary = windows,
    node_coverage = coverage,
    max_panel_size = 3L,
    min_marginal_gain = 7,
    cap_repeated_coverage = TRUE
  )

  expect_equal(selected$selected_panel$window_id, "w1")
})

test_that("tied windows are handled deterministically by score and coordinates", {
  windows <- data.frame(
    window_id = c("w_late", "w_early"),
    start = c(101L, 1L),
    end = c(200L, 100L),
    width = c(100L, 100L),
    informative_snps = c(2L, 2L),
    total_weighted_gain = c(5, 5)
  )
  coverage <- data.frame(
    window_id = c("w_late", "w_early"),
    node_id = c("n1", "n2"),
    node_score = c(5, 5)
  )

  selected <- greedy_select_panel(windows, coverage, max_panel_size = 1L)

  expect_equal(selected$selected_panel$window_id, "w_early")
})
