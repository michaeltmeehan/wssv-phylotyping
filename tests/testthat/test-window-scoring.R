source(test_path("../../R/window_scoring.R"))

synthetic_site_node_scores <- function() {
  data.frame(
    node_index = c(1L, 2L, 1L),
    node_id = c(11L, 12L, 11L),
    site = c(5L, 5L, 18L),
    gain = c(0.2, 0.4, 0.6),
    normalized_gain = c(0.5, 0.8, 1.0),
    best_allele = c("A", "C", "G"),
    direction = c("clade", "outside", "clade"),
    allele_count_inside = c(2L, 3L, 4L),
    allele_count_outside = c(0L, 1L, 0L),
    observed_inside = c(4L, 4L, 4L),
    observed_outside = c(4L, 4L, 4L),
    total_observed = c(8L, 8L, 8L),
    clade_size = c(2L, 3L, 2L),
    complement_size = c(6L, 5L, 6L),
    depth = c(0.2, 0.8, 0.2),
    balance = c(0.5, 0.75, 0.5),
    weight = c(1.0, 2.0, 1.0)
  )
}

test_that("fixed-width windows are generated with configured overlap", {
  windows <- generate_fixed_windows(alignment_length = 25L, widths = c(10L), overlap = 0.5)

  expect_equal(windows$start, c(1L, 6L, 11L, 16L, 21L))
  expect_equal(windows$end, c(10L, 15L, 20L, 25L, 25L))
  expect_equal(windows$width, c(10L, 10L, 10L, 10L, 5L))
  expect_equal(windows$window_type, rep("fixed", 5L))
})

test_that("SNPs are assigned to every overlapping window", {
  windows <- data.frame(
    window_id = c("w1", "w2"),
    window_type = "fixed",
    start = c(1L, 5L),
    end = c(10L, 15L),
    width = c(10L, 11L),
    requested_width = c(10L, 10L),
    centre_site = NA_integer_
  )

  assigned <- assign_snps_to_windows(c(3L, 5L, 12L), windows)

  expect_equal(assigned$window_id, c("w1", "w1", "w2", "w2"))
  expect_equal(assigned$site, c(3L, 5L, 5L, 12L))
})

test_that("window-node scores aggregate site-node rows and preserve metadata", {
  windows <- data.frame(
    window_id = "w1",
    window_type = "fixed",
    start = 1L,
    end = 20L,
    width = 20L,
    requested_width = 20L,
    centre_site = NA_integer_
  )
  scores <- synthetic_site_node_scores()

  node_summary <- aggregate_window_node_scores(windows, scores)

  node_11 <- node_summary[node_summary$node_id == 11L, ]
  expect_equal(node_11$informative_snps, 2L)
  expect_equal(node_11$node_site_score_entries, 2L)
  expect_equal(node_11$total_gain, 0.8)
  expect_equal(node_11$total_weighted_gain, 1.5)
  expect_equal(node_11$clade_size, 2L)
  expect_equal(node_11$depth, 0.2)
})

test_that("window summaries keep empty windows with zero signal when allowed", {
  windows <- data.frame(
    window_id = c("w1", "w2"),
    window_type = "fixed",
    start = c(1L, 21L),
    end = c(20L, 30L),
    width = c(20L, 10L),
    requested_width = c(20L, 10L),
    centre_site = NA_integer_
  )
  scores <- synthetic_site_node_scores()
  node_summary <- aggregate_window_node_scores(windows, scores)
  aln_int <- matrix(c(1L, 0L, 2L, 2L), nrow = 2L, ncol = 2L)

  summary <- summarise_windows(
    windows = windows,
    polymorphic_sites = c(5L, 18L, 25L),
    site_node_scores = scores,
    window_node_summary = node_summary,
    aln_int = NULL,
    min_informative_snps = 0L
  )

  empty <- summary[summary$window_id == "w2", ]
  expect_equal(empty$polymorphic_snps, 1L)
  expect_equal(empty$informative_snps, 0L)
  expect_equal(empty$node_site_score_entries, 0L)
  expect_equal(empty$total_gain, 0)
  expect_equal(empty$nodes_covered, 0L)
  expect_true(is.na(empty$mean_gain))
})

test_that("candidate scoring filters by minimum informative SNPs", {
  windows <- data.frame(
    window_id = c("w1", "w2"),
    window_type = "fixed",
    start = c(1L, 21L),
    end = c(20L, 30L),
    width = c(20L, 10L),
    requested_width = c(20L, 10L),
    centre_site = NA_integer_
  )

  scored <- score_candidate_windows(
    windows = windows,
    polymorphic_sites = c(5L, 18L, 25L),
    site_node_scores = synthetic_site_node_scores(),
    min_informative_snps = 1L
  )

  expect_equal(scored$candidate_windows$window_id, "w1")
  expect_true(all(scored$window_node_summary$window_id == "w1"))
})
