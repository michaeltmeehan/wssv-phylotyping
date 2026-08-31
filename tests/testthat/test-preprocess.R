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

test_that("the MCC tree posterior support is imported into node labels", {
  path <- test_path("../../data/raw/tree/41SEQSUPDATEDEDITEDMCC")
  tree <- read_tree_any(path)

  expect_length(tree$node.label, tree$Nnode)
  expect_true(all(!is.na(suppressWarnings(as.numeric(tree$node.label)))))

  if (requireNamespace("treeio", quietly = TRUE)) {
    beast <- treeio::read.beast(path)
    data <- as.data.frame(beast@data)
    expected <- data[!is.na(as.integer(data$node)) & as.integer(data$node) > length(tree$tip.label), c("node", "posterior")]
    expected <- expected[order(as.integer(expected$node)), , drop = FALSE]

    expect_equal(as.numeric(tree$node.label), expected$posterior, tolerance = 1e-7)
  }
})

test_that("eligible clade masks and node metadata are generated", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")
  tree$node.label <- c("0.9", "0.85", "0.75", "1.0")

  clades <- generate_clade_masks(tree, min_clade_size = 2L, max_clade_frac = 0.95)
  metadata <- calculate_node_metadata(tree, clades)

  expect_true(any(rowSums(clades$all_clade_mask) == 2L))
  expect_named(metadata, c(
    "node_index", "node_id", "posterior_support", "clade_size", "complement_size",
    "depth", "balance", "weight"
  ))
  expect_equal(metadata$node_index, seq_len(nrow(metadata)))
  expect_true(all(metadata$clade_size >= 2L))
  expect_equal(metadata$posterior_support, as.numeric(tree$node.label[clades$all_node_id - length(tree$tip.label)]))
})

test_that("configured target clades resolve against the shipped tree", {
  config <- yaml::read_yaml(test_path("../../config/config.yml"))
  tree <- read_tree_any(test_path("../../data/raw/tree/41SEQSUPDATEDEDITEDMCC"))
  all_clades <- generate_clade_masks(
    tree,
    min_clade_size = config$analysis$min_clade_size,
    max_clade_frac = config$analysis$max_clade_frac
  )
  all_metadata <- calculate_node_metadata(tree, all_clades)
  target_clades <- resolve_target_clades(
    tree,
    all_clades,
    config$analysis$target_clades,
    min_posterior_support = config$analysis$min_target_posterior_support,
    min_clade_size = config$analysis$min_clade_size,
    all_node_metadata = all_metadata
  )
  target_metadata <- calculate_node_metadata(tree, target_clades)

  expect_length(target_metadata$node_id, 4L)
  expect_equal(target_metadata$node_id, c(43L, 44L, 46L, 47L))
  expect_true(all(target_metadata$clade_size >= config$analysis$min_clade_size))
  expect_true(all(target_metadata$posterior_support >= config$analysis$min_target_posterior_support))
  expect_true(all(target_metadata$node_id %in% all_metadata$node_id))
})

test_that("node metadata validation fails when posterior support is missing or malformed", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")
  clades <- generate_clade_masks(tree, min_clade_size = 2L, max_clade_frac = 0.95)

  expect_error(
    calculate_node_metadata(tree, clades),
    "posterior support"
  )

  tree$node.label <- c("0.9", "bad", "0.75", "1.0")
  expect_error(
    calculate_node_metadata(tree, clades),
    "Malformed posterior support"
  )
})
