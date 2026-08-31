source(test_path("../../R/script_utils.R"))
source(test_path("../../R/window_scoring.R"))

synthetic_amplicon_inputs <- function(weight_scale = 1) {
  aln_int <- matrix(1L, nrow = 4L, ncol = 40L)
  rownames(aln_int) <- paste0("t", seq_len(nrow(aln_int)))
  colnames(aln_int) <- as.character(seq_len(ncol(aln_int)))
  aln_int[, 1:7] <- 2L
  aln_int[, 14:19] <- 2L
  aln_int[, 30:35] <- 3L

  config <- list(
    analysis = list(
      amplicons = list(
        primary_target_nodes = c(45L, 66L, 72L, 80L),
        diagnostic_target_nodes = c(67L),
        min_length = 6L,
        max_length = 6L,
        flank_search_width = 4L,
        strong_gain_threshold = 0.8,
        minimum_gain_threshold = 0.5,
        flank_min_non_missing_fraction = 0.9,
        flank_max_variable_fraction = 0.25,
        flank_min_gc_fraction = 0.35,
        flank_min_conserved_fraction = 0.9,
        flank_complete_fraction = 0.95
      )
    )
  )

  site_node_scores <- data.frame(
    node_index = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 5L, 5L),
    node_id = c(45L, 66L, 72L, 80L, 45L, 66L, 72L, 80L, 67L, 67L),
    site = c(10L, 11L, 12L, 13L, 20L, 21L, 22L, 23L, 11L, 21L),
    gain = c(0.5, 0.4, 0.3, 0.6, 0.2, 0.5, 0.45, 0.35, 0.9, 0.8),
    normalized_gain = c(0.9, 0.8, 0.7, 0.95, 0.6, 0.9, 0.85, 0.75, 1.4, 1.3),
    best_allele = c("A", "C", "G", "T", "A", "C", "G", "T", "G", "T"),
    direction = rep("clade", 10L),
    allele_count_inside = c(8L, 7L, 6L, 9L, 8L, 7L, 6L, 9L, 5L, 5L),
    allele_count_outside = c(1L, 2L, 3L, 1L, 1L, 2L, 3L, 1L, 0L, 0L),
    observed_inside = rep(20L, 10L),
    observed_outside = rep(20L, 10L),
    total_observed = rep(40L, 10L),
    weight = c(1, 1, 1, 1, 1, 1, 1, 1, 9, 9) * weight_scale,
    stringsAsFactors = FALSE
  )

  list(
    site_node_scores = site_node_scores,
    polymorphic_sites = c(10L, 11L, 12L, 13L, 20L, 21L, 22L, 23L),
    aln_int = aln_int,
    site_map = data.frame(site = seq_len(ncol(aln_int)), alignment_position = seq_len(ncol(aln_int))),
    config = config
  )
}

test_that("candidate intervals respect min and max length and stay deterministic", {
  input <- synthetic_amplicon_inputs()

  scored_one <- score_candidate_amplicons(
    site_node_scores = input$site_node_scores,
    polymorphic_sites = input$polymorphic_sites,
    aln_int = input$aln_int,
    site_map = input$site_map,
    config = input$config
  )
  scored_two <- score_candidate_amplicons(
    site_node_scores = input$site_node_scores,
    polymorphic_sites = input$polymorphic_sites,
    aln_int = input$aln_int,
    site_map = input$site_map,
    config = input$config
  )

  expect_equal(scored_one$candidate_amplicons, scored_two$candidate_amplicons)
  expect_equal(scored_one$candidate_amplicon_target_scores, scored_two$candidate_amplicon_target_scores)
  expect_equal(scored_one$candidate_amplicon_snps, scored_two$candidate_amplicon_snps)
  expect_true(all(scored_one$candidate_amplicons$length >= 6L))
  expect_true(all(scored_one$candidate_amplicons$length <= 6L))
  expect_true(all(scored_one$candidate_amplicons$end - scored_one$candidate_amplicons$start + 1L == scored_one$candidate_amplicons$length))
  expect_equal(anyDuplicated(paste(scored_one$candidate_amplicons$start, scored_one$candidate_amplicons$end)), 0L)
})

test_that("primary target filtering excludes diagnostic child nodes", {
  input <- synthetic_amplicon_inputs()

  scored <- score_candidate_amplicons(
    site_node_scores = input$site_node_scores,
    polymorphic_sites = input$polymorphic_sites,
    aln_int = input$aln_int,
    site_map = input$site_map,
    config = input$config
  )

  expect_true(all(scored$candidate_amplicon_target_scores$target_node %in% c(45L, 66L, 72L, 80L)))
  expect_false(any(scored$candidate_amplicon_target_scores$target_node == 67L))
})

test_that("legacy weight does not affect the new stage-04 outputs", {
  input_a <- synthetic_amplicon_inputs(weight_scale = 1)
  input_b <- synthetic_amplicon_inputs(weight_scale = 100)

  scored_a <- score_candidate_amplicons(
    site_node_scores = input_a$site_node_scores,
    polymorphic_sites = input_a$polymorphic_sites,
    aln_int = input_a$aln_int,
    site_map = input_a$site_map,
    config = input_a$config
  )
  scored_b <- score_candidate_amplicons(
    site_node_scores = input_b$site_node_scores,
    polymorphic_sites = input_b$polymorphic_sites,
    aln_int = input_b$aln_int,
    site_map = input_b$site_map,
    config = input_b$config
  )

  expect_equal(scored_a$candidate_amplicons, scored_b$candidate_amplicons)
  expect_equal(scored_a$candidate_amplicon_target_scores, scored_b$candidate_amplicon_target_scores)
  expect_equal(scored_a$candidate_amplicon_snps, scored_b$candidate_amplicon_snps)
})

test_that("candidate-target summaries agree with underlying SNP records", {
  input <- synthetic_amplicon_inputs()

  scored <- score_candidate_amplicons(
    site_node_scores = input$site_node_scores,
    polymorphic_sites = input$polymorphic_sites,
    aln_int = input$aln_int,
    site_map = input$site_map,
    config = input$config
  )

  snps <- scored$candidate_amplicon_snps
  targets <- scored$candidate_amplicon_target_scores

  for (i in seq_len(nrow(targets))) {
    row <- targets[i, , drop = FALSE]
    subset <- snps[snps$candidate_id == row$candidate_id & snps$target_node == row$target_node, , drop = FALSE]
    expect_equal(length(unique(subset$site)), row$n_informative_snps)
    expect_equal(max(subset$raw_gain), row$max_raw_gain)
    expect_equal(max(subset$normalized_gain), row$max_normalized_gain)
    expect_equal(mean(subset$normalized_gain), row$mean_normalized_gain)
    expect_equal(stats::median(subset$normalized_gain), row$median_normalized_gain)
    best <- subset[order(-subset$normalized_gain, -subset$raw_gain, subset$site), , drop = FALSE][1L, ]
    expect_equal(best$site, row$best_snp_position)
    expect_equal(best$best_snp_rule, row$best_snp_rule)
  }
})

test_that("candidate summaries preserve the evidence vector by target", {
  input <- synthetic_amplicon_inputs()

  scored <- score_candidate_amplicons(
    site_node_scores = input$site_node_scores,
    polymorphic_sites = input$polymorphic_sites,
    aln_int = input$aln_int,
    site_map = input$site_map,
    config = input$config
  )

  expect_true(all(c("targets_with_any_signal", "targets_with_strong_signal", "target_nodes_represented") %in% names(scored$candidate_amplicons)))
  expect_true(all(scored$candidate_amplicons$targets_with_any_signal <= scored$candidate_amplicons$primary_targets_represented))
  expect_true(all(scored$candidate_amplicons$targets_with_strong_signal <= scored$candidate_amplicons$primary_targets_represented))
})
