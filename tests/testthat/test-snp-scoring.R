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

test_that("stage 02 broad scores remain unchanged and target scores stay restricted", {
  pre <- readRDS(test_path("../../data/processed/precomputed.rds"))
  expected_all <- readRDS(test_path("../../outputs/tables/site_node_scores.rds"))

  actual_all <- score_snp_sites(
    aln_int = pre$aln_int,
    target_mask = pre$all_clade_mask,
    node_metadata = pre$all_node_metadata,
    sites = pre$polymorphic_sites,
    min_total_obs = pre$params$min_total_obs,
    min_side_obs = pre$params$min_side_obs,
    min_site_maf = pre$params$min_site_maf,
    min_gain_norm = pre$params$min_gain_norm,
    progress = FALSE
  )

  actual_targets <- score_snp_sites(
    aln_int = pre$aln_int,
    target_mask = pre$target_clade_mask,
    node_metadata = pre$target_node_metadata,
    sites = pre$polymorphic_sites,
    min_total_obs = pre$params$min_total_obs,
    min_side_obs = pre$params$min_side_obs,
    min_site_maf = pre$params$min_site_maf,
    min_gain_norm = pre$params$min_gain_norm,
    progress = FALSE
  )

  score_order <- order(actual_all$site, actual_all$node_index)
  actual_all <- actual_all[score_order, , drop = FALSE]
  expected_all <- expected_all[order(expected_all$site, expected_all$node_index), , drop = FALSE]

  expect_equal(actual_all[names(expected_all)], expected_all)
  expect_true("posterior_support" %in% names(actual_all))
  expect_true(all(actual_targets$node_id %in% pre$target_node_metadata$node_id))
})
