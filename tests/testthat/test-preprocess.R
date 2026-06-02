source(test_path("../../R/encoding.R"))
source(test_path("../../R/preprocess.R"))
source(test_path("../../R/tree_utils.R"))

test_that("integer encoding maps canonical bases and missing states", {
  alignment <- c(s1 = "ACGTN-", s2 = "acgtRY")

  encoded <- encode_alignment_int(alignment)

  expect_equal(unname(encoded[1, ]), c(1L, 2L, 3L, 4L, 0L, 0L))
  expect_equal(unname(encoded[2, ]), c(1L, 2L, 3L, 4L, 0L, 0L))
})

test_that("tree and alignment are matched, pruned, and reordered", {
  tree <- ape::read.tree(text = "((b:1,a:1):1,c:1,d:1);")
  alignment <- c(a = "AC", b = "AC", c = "GT", extra = "GT")

  matched <- match_prune_reorder_alignment_tree(alignment, tree)

  expect_equal(matched$tree$tip.label, c("b", "a", "c"))
  expect_equal(names(matched$alignment), matched$tree$tip.label)
  expect_equal(matched$dropped_alignment, "extra")
  expect_equal(matched$dropped_tree, "d")
})

test_that("eligible clade masks and node metadata are generated", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")

  clades <- generate_clade_masks(tree, min_clade_size = 2L, max_clade_frac = 0.95)
  metadata <- calculate_node_metadata(tree, clades)

  expect_true(any(rowSums(clades$target_mask) == 2L))
  expect_named(metadata, c(
    "node_index", "node_id", "clade_size", "complement_size",
    "depth", "balance", "weight"
  ))
  expect_equal(metadata$node_index, seq_len(nrow(metadata)))
  expect_true(all(metadata$clade_size >= 2L))
})
