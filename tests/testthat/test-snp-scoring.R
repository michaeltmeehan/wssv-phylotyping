source(test_path("../../R/encoding.R"))
source(test_path("../../R/tree_utils.R"))
source(test_path("../../R/snp_scoring.R"))
source(test_path("../../R/target_diagnostics.R"))

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

test_that("stage 03 summaries split all and target clades with transparent ordering", {
  site_map <- data.frame(
    site = c(10L, 20L),
    alignment_position = c(1001L, 1002L)
  )
  node_metadata <- data.frame(
    node_index = c(1L, 2L),
    node_id = c(11L, 12L),
    posterior_support = c(0.95, 0.91),
    clade_size = c(5L, 4L),
    complement_size = c(10L, 11L),
    depth = c(0.4, 0.8),
    balance = c(0.5, 0.4),
    weight = c(0.18, 0.10)
  )
  site_node_scores <- data.frame(
    node_index = c(1L, 2L, 1L, 2L),
    node_id = c(11L, 12L, 11L, 12L),
    site = c(10L, 10L, 20L, 20L),
    gain = c(0.30, 0.25, 0.85, 0.15),
    normalized_gain = c(0.55, 0.60, 0.70, 0.20),
    best_allele = c("A", "A", "G", "G"),
    direction = c("clade", "outside", "clade", "outside"),
    allele_count_inside = c(3L, 2L, 4L, 1L),
    allele_count_outside = c(1L, 1L, 0L, 0L),
    observed_inside = c(4L, 3L, 5L, 2L),
    observed_outside = c(6L, 7L, 5L, 8L),
    total_observed = c(10L, 10L, 10L, 10L),
    weight = c(100, 100, 1, 1)
  )

  site_summary_all <- make_site_summary(
    site_node_scores,
    site_map = site_map,
    n_tip = 12L,
    helped_count_name = "nodes_helped"
  )
  site_summary_targets <- make_site_summary(
    site_node_scores,
    site_map = site_map,
    n_tip = 12L,
    helped_count_name = "target_clades_helped"
  )
  node_summary_targets <- make_node_summary(site_node_scores, node_metadata)

  expect_equal(site_summary_all$site, c(20L, 10L))
  expect_equal(site_summary_all$alignment_position, c(1002L, 1001L))
  expect_true(site_summary_all$weighted_gain_sum[1] < site_summary_all$weighted_gain_sum[2])
  expect_equal(site_summary_all$nodes_helped, c(2L, 2L))
  expect_equal(site_summary_all$best_allele, c("G", "A"))
  expect_equal(site_summary_all$best_direction, c("clade", "outside"))
  expect_equal(site_summary_all$missing_sequences_mean, c(2, 2))
  expect_true("legacy_weighted_gain_sum" %in% names(site_summary_all))
  expect_true("weighted_normalized_gain_sum" %in% names(site_summary_all))
  expect_equal(site_summary_all$weighted_normalized_gain_sum, site_summary_all$weighted_gain_sum)

  expect_equal(site_summary_targets$site, c(20L, 10L))
  expect_equal(site_summary_targets$target_clades_helped, c(2L, 2L))
  expect_true(all(site_summary_targets$best_normalized_gain >= site_summary_targets$normalized_gain_mean))

  expect_equal(nrow(node_summary_targets), 2L)
  expect_true(all(c("posterior_support", "clade_size", "complement_size", "depth") %in% names(node_summary_targets)))
  expect_equal(node_summary_targets$node_id, c(11L, 12L))
  expect_equal(node_summary_targets$informative_sites, c(2L, 2L))
  expect_equal(node_summary_targets$best_site, c(20L, 10L))
  expect_equal(node_summary_targets$best_raw_gain, c(0.85, 0.25))
  expect_equal(node_summary_targets$best_normalized_gain, c(0.70, 0.60))
})

test_that("target-clade checkpoint keeps all configured targets and surfaces strongest SNPs", {
  pre <- list(
    aln_int = matrix(0L, nrow = 5L, ncol = 3L),
    site_map = data.frame(
      site = c(1L, 2L, 3L),
      alignment_position = c(101L, 102L, 103L)
    ),
    polymorphic_sites = c(1L, 2L, 3L),
    target_node_metadata = data.frame(
      node_index = c(1L, 2L, 3L),
      node_id = c(11L, 12L, 13L),
      posterior_support = c(0.95, 0.91, 0.99),
      clade_size = c(3L, 2L, 2L),
      complement_size = c(2L, 3L, 3L),
      depth = c(0.4, 0.8, 0.6),
      balance = c(0.5, 0.4, 0.4),
      weight = c(0.18, 0.10, 0.08)
    ),
    target_clades = list(
      target_clade_spec = list(
        list(target_name = "Alpha"),
        list(node_id = 12L),
        list(target_name = "Gamma")
      )
    )
  )
  site_node_scores_targets <- data.frame(
    node_index = c(1L, 1L, 2L),
    node_id = c(11L, 11L, 12L),
    site = c(1L, 2L, 3L),
    gain = c(0.30, 0.20, 0.40),
    normalized_gain = c(0.60, 0.50, 0.80),
    best_allele = c("A", "G", "T"),
    direction = c("clade", "outside", "clade"),
    allele_count_inside = c(2L, 1L, 1L),
    allele_count_outside = c(0L, 0L, 1L),
    observed_inside = c(3L, 4L, 2L),
    observed_outside = c(2L, 1L, 3L),
    total_observed = c(5L, 5L, 5L),
    stringsAsFactors = FALSE
  )

  checkpoint <- make_target_clade_checkpoint(pre, site_node_scores_targets, top_n = 2L)

  expect_equal(checkpoint$summary$target, c("Alpha", "node 12", "Gamma"))
  expect_equal(checkpoint$summary$node_id, c(11L, 12L, 13L))
  expect_equal(checkpoint$summary$informative_sites_passing_threshold, c(2L, 1L, 0L))
  expect_equal(checkpoint$summary$polymorphic_sites_assessed, c(3L, 3L, 3L))
  expect_equal(checkpoint$summary$best_alignment_position, c(101, 103, NA_real_))
  expect_equal(checkpoint$summary$best_site_missing_sequences, c(0L, 0L, NA_integer_))
  expect_equal(checkpoint$summary$max_normalized_gain, c(0.60, 0.80, NA_real_))
  expect_equal(checkpoint$summary$mean_normalized_gain, c(0.55, 0.80, NA_real_))

  expect_equal(nrow(checkpoint$strongest_snps), 3L)
  expect_equal(checkpoint$strongest_snps$target, c("Alpha", "Alpha", "node 12"))
  expect_equal(checkpoint$strongest_snps$rank_within_target, c(1L, 2L, 1L))
  expect_equal(checkpoint$strongest_snps$best_rule, c("A / clade", "G / outside", "T / clade"))
  expect_equal(round(checkpoint$strongest_snps$accuracy, 2), c(0.90, 1.00, 1.00))
  expect_equal(round(checkpoint$strongest_snps$baseline_accuracy, 2), c(0.60, 0.80, 0.60))
  expect_true(all(c("raw_gain", "normalized_gain", "missing_fraction") %in% names(checkpoint$strongest_snps)))
})
