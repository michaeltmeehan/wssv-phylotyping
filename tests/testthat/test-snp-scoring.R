source(test_path("../../R/encoding.R"))
source(test_path("../../R/tree_utils.R"))
source(test_path("../../R/snp_scoring.R"))

test_that("polymorphic sites ignore ambiguous and missing bases", {
  alignment <- c(a = "AN", b = "A-", c = "CN")
  encoded <- encode_alignment_int(alignment)

  expect_equal(unname(find_polymorphic_sites(encoded)), 1L)
})

test_that("SNP scoring finds a clade-marking allele", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")
  tree$node.label <- c("0.9", "0.85", "0.75", "1.0")
  alignment <- c(a = "A", b = "A", c = "C", d = "C", e = "C")
  encoded <- encode_alignment_int(alignment)
  clades <- generate_clade_masks(tree, min_clade_size = 2L, max_clade_frac = 0.95)
  metadata <- calculate_node_metadata(tree, clades)

  scores <- score_snp_sites(
    aln_int = encoded,
    target_mask = clades$all_clade_mask,
    node_metadata = metadata,
    sites = find_polymorphic_sites(encoded),
    min_total_obs = 4L,
    min_side_obs = 1L,
    min_site_maf = 1L,
    min_gain_norm = 0.5,
    progress = FALSE
  )

  best <- scores[which.max(scores$normalized_gain), ]
  expect_equal(best$site, 1L)
  expect_equal(best$best_allele, "A")
  expect_equal(best$direction, "clade")
  expect_equal(best$allele_count_inside, 2L)
  expect_equal(best$allele_count_outside, 0L)
  expect_equal(best$total_observed, 5L)
  expect_equal(best$normalized_gain, 1)
})
