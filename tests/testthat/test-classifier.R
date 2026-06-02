source(test_path("../../R/classifier.R"))

toy_classifier_inputs <- function() {
  aln <- matrix(c(
    1L, 1L, 3L, 1L,
    1L, 1L, 3L, 1L,
    1L, 2L, 4L, 1L,
    1L, 2L, 4L, 1L
  ), nrow = 4L, byrow = TRUE)
  rownames(aln) <- c("a", "b", "c", "d")
  target_mask <- matrix(c(
    TRUE, TRUE, FALSE, FALSE,
    TRUE, FALSE, FALSE, FALSE,
    FALSE, FALSE, TRUE, TRUE
  ), nrow = 3L, byrow = TRUE)
  rownames(target_mask) <- c("10", "11", "12")
  colnames(target_mask) <- rownames(aln)
  node_metadata <- data.frame(
    node_index = 1:3,
    node_id = c("10", "11", "12"),
    clade_size = c(2L, 1L, 2L),
    complement_size = c(2L, 3L, 2L),
    depth = c(1, 2, 1),
    weight = c(1, 1, 1)
  )
  scores <- data.frame(
    node_id = c("10", "10", "11", "12", "12"),
    site = c(2L, 3L, 4L, 2L, 3L),
    best_allele = c("A", "G", "A", "C", "T"),
    direction = c("clade", "clade", "clade", "clade", "clade"),
    normalized_gain = c(1, 1, 1, 1, 1)
  )
  list(aln = aln, target_mask = target_mask, node_metadata = node_metadata, scores = scores)
}

make_toy_classifier <- function(min_support = 0.8, max_conflict = 0.2, ...) {
  x <- toy_classifier_inputs()
  train_classifier(
    x$aln, x$target_mask, x$node_metadata, x$scores,
    use_selected_panel = FALSE,
    min_support = min_support,
    max_conflict = max_conflict,
    ...
  )
}

test_that("node-level allele-rule models are trained", {
  classifier <- make_toy_classifier()

  expect_s3_class(classifier, "wssv_classifier")
  expect_equal(sort(unique(classifier$rules$node_id)), c("10", "11", "12"))
  expect_equal(classifier$node_table$parent_node_id[classifier$node_table$node_id == "11"], "10")
})

test_that("observed informative sites are extracted from partial vectors", {
  classifier <- make_toy_classifier()
  query <- c(0L, 1L, 0L, 1L)
  observed <- extract_observed_informative_sites(query, classifier)

  expect_equal(observed$site, c(2L, 4L))
  expect_equal(observed$allele_code, c(1L, 1L))
})

test_that("tree-path classifier resolves a query with enough evidence", {
  classifier <- make_toy_classifier()
  pred <- classify_tree_path(c(1L, 1L, 3L, 1L), classifier)

  expect_equal(pred$status, "resolved")
  expect_equal(pred$assigned_node, "11")
  expect_equal(pred$observed_informative_sites, 3L)
})

test_that("classifier returns no informative sites for all-missing queries", {
  classifier <- make_toy_classifier()
  pred <- classify_tree_path(c(0L, 0L, 0L, 0L), classifier)

  expect_equal(pred$status, "no_informative_sites")
  expect_true(is.na(pred$assigned_node))
})

test_that("classifier stops when support is below threshold", {
  classifier <- make_toy_classifier(min_support = 1, max_conflict = 1, support_margin = -1)
  pred <- classify_tree_path(c(1L, 1L, 4L, 0L), classifier)

  expect_equal(pred$status, "weak_support")
  expect_true(is.na(pred$assigned_node))
})

test_that("classifier detects conflicting evidence", {
  classifier <- make_toy_classifier(max_conflict = 0)
  pred <- classify_tree_path(c(1L, 2L, 4L, 1L), classifier)

  expect_equal(pred$status, "conflicting")
  expect_true(is.na(pred$assigned_node))
})

test_that("nested node support is not a conflict", {
  classifier <- make_toy_classifier(max_conflict = 0)
  pred <- classify_tree_path(c(1L, 1L, 3L, 1L), classifier)

  expect_equal(pred$status, "resolved")
  expect_equal(pred$assigned_node, "11")
  expect_equal(pred$conflict_reason, "nested_true_path_compatible_support")
})

test_that("incompatible off-path support is a conflict", {
  classifier <- make_toy_classifier(max_conflict = 0)
  pred <- classify_tree_path(c(1L, 2L, 4L, 1L), classifier)

  expect_equal(pred$status, "conflicting")
  expect_match(pred$conflict_reason, "incompatible")
})

test_that("weak off-path evidence does not automatically trigger conflict", {
  classifier <- make_toy_classifier(min_support = 0.8, max_conflict = 0.2, support_margin = -1)
  pred <- classify_tree_path(c(1L, 1L, 4L, 1L), classifier)

  expect_equal(pred$status, "resolved")
  expect_equal(pred$assigned_node, "11")
})

test_that("deepest supported node is selected correctly", {
  classifier <- make_toy_classifier()
  pred <- classify_tree_path(c(1L, 1L, 3L, 1L), classifier)

  expect_equal(pred$assigned_node, "11")
  expect_equal(pred$assigned_depth, 2)
})

test_that("unresolved weak support is returned when evidence is insufficient", {
  classifier <- make_toy_classifier(min_support = 1, max_conflict = 1, support_margin = -1)
  pred <- classify_tree_path(c(1L, 1L, 4L, 0L), classifier)

  expect_equal(pred$status, "weak_support")
  expect_true(is.na(pred$assigned_node))
  expect_equal(pred$conflict_reason, "weak_evidence")
})

test_that("classifier objects can be saved and loaded", {
  classifier <- make_toy_classifier()
  path <- tempfile(fileext = ".rds")
  save_classifier(classifier, path)
  loaded <- load_classifier(path)

  expect_s3_class(loaded, "wssv_classifier")
  expect_equal(loaded$informative_sites, classifier$informative_sites)
})
