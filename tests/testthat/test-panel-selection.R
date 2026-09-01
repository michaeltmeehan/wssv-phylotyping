source(test_path("../../R/script_utils.R"))
source(test_path("../../R/panel_selection.R"))

synthetic_panel_inputs <- function() {
  config <- list(
    analysis = list(
      amplicons = list(
        primary_target_nodes = c(45L, 66L, 72L, 80L),
        minimum_gain_threshold = 0.5,
        strong_gain_threshold = 0.8
      ),
      panel_selection = list(
        min_panel_size = 2L,
        max_panel_size = 3L,
        top_n_panels = 10L,
        exact_enumeration_limit = 10000L,
        approximate_search_budget = 10L,
        max_reduced_candidates = 10L,
        evidence_threshold = 0.5,
        strong_evidence_threshold = 0.8,
        material_gain_improvement = 0.01,
        material_relative_gain_improvement = 0.01,
        redundancy_cap_per_target = 2L,
        max_allowed_overlap = 0.5,
        min_distance = 20L,
        local_redundancy_overlap_threshold = 0.8,
        local_redundancy_min_distance = 20L
      )
    )
  )

  candidate_amplicons <- data.frame(
    candidate_id = c("amp_A", "amp_B", "amp_C", "amp_D", "amp_E"),
    start = c(1L, 121L, 241L, 361L, 5L),
    end = c(100L, 220L, 340L, 460L, 105L),
    length = c(100L, 100L, 100L, 100L, 101L),
    complete_fraction = c(0.98, 0.97, 0.96, 0.95, 0.90),
    missing_fraction = c(0.02, 0.03, 0.04, 0.05, 0.10),
    left_flank_conserved_stretch = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    right_flank_conserved_stretch = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    unique_informative_snp_positions = c(2L, 2L, 2L, 2L, 2L),
    stringsAsFactors = FALSE
  )

  candidate_amplicon_target_scores <- rbind(
    data.frame(candidate_id = "amp_A", target_node = 45L, max_normalized_gain = 0.95, n_informative_snps = 2L, best_snp_position = 10L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_A", target_node = 66L, max_normalized_gain = 0.70, n_informative_snps = 2L, best_snp_position = 12L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_A", target_node = 72L, max_normalized_gain = 0.25, n_informative_snps = 1L, best_snp_position = 14L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_A", target_node = 80L, max_normalized_gain = 0.10, n_informative_snps = 1L, best_snp_position = 16L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 45L, max_normalized_gain = 0.20, n_informative_snps = 1L, best_snp_position = 130L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 66L, max_normalized_gain = 0.60, n_informative_snps = 2L, best_snp_position = 132L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 72L, max_normalized_gain = 0.85, n_informative_snps = 2L, best_snp_position = 134L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 80L, max_normalized_gain = 0.15, n_informative_snps = 1L, best_snp_position = 136L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 45L, max_normalized_gain = 0.25, n_informative_snps = 1L, best_snp_position = 250L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 66L, max_normalized_gain = 0.55, n_informative_snps = 1L, best_snp_position = 252L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 72L, max_normalized_gain = 0.90, n_informative_snps = 2L, best_snp_position = 254L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 80L, max_normalized_gain = 0.20, n_informative_snps = 1L, best_snp_position = 256L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 45L, max_normalized_gain = 0.15, n_informative_snps = 1L, best_snp_position = 370L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 66L, max_normalized_gain = 0.30, n_informative_snps = 1L, best_snp_position = 372L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 72L, max_normalized_gain = 0.40, n_informative_snps = 1L, best_snp_position = 374L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 80L, max_normalized_gain = 0.85, n_informative_snps = 2L, best_snp_position = 376L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 45L, max_normalized_gain = 0.95, n_informative_snps = 2L, best_snp_position = 10L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 66L, max_normalized_gain = 0.70, n_informative_snps = 2L, best_snp_position = 12L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 72L, max_normalized_gain = 0.25, n_informative_snps = 1L, best_snp_position = 14L, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 80L, max_normalized_gain = 0.10, n_informative_snps = 1L, best_snp_position = 16L, stringsAsFactors = FALSE)
  )

  candidate_amplicon_snps <- rbind(
    data.frame(candidate_id = "amp_A", target_node = 45L, site = 10L, normalized_gain = 0.95, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_A", target_node = 45L, site = 11L, normalized_gain = 0.90, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_A", target_node = 66L, site = 12L, normalized_gain = 0.70, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_A", target_node = 72L, site = 14L, normalized_gain = 0.25, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_A", target_node = 80L, site = 16L, normalized_gain = 0.10, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 45L, site = 130L, normalized_gain = 0.20, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 66L, site = 132L, normalized_gain = 0.60, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 72L, site = 134L, normalized_gain = 0.85, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_B", target_node = 80L, site = 136L, normalized_gain = 0.15, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 45L, site = 250L, normalized_gain = 0.25, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 66L, site = 252L, normalized_gain = 0.55, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 72L, site = 254L, normalized_gain = 0.90, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_C", target_node = 80L, site = 256L, normalized_gain = 0.20, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 45L, site = 370L, normalized_gain = 0.15, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 66L, site = 372L, normalized_gain = 0.30, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 72L, site = 374L, normalized_gain = 0.40, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_D", target_node = 80L, site = 376L, normalized_gain = 0.85, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 45L, site = 10L, normalized_gain = 0.95, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 45L, site = 11L, normalized_gain = 0.90, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 66L, site = 12L, normalized_gain = 0.70, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 72L, site = 14L, normalized_gain = 0.25, stringsAsFactors = FALSE),
    data.frame(candidate_id = "amp_E", target_node = 80L, site = 16L, normalized_gain = 0.10, stringsAsFactors = FALSE)
  )

  list(
    config = config,
    candidate_amplicons = candidate_amplicons,
    candidate_amplicon_target_scores = candidate_amplicon_target_scores,
    candidate_amplicon_snps = candidate_amplicon_snps
  )
}

make_panel_metrics <- function(panel_size, min_abs, min_rel, redundancy, total_length, candidate_key, panel_id) {
  data.frame(
    panel_size = panel_size,
    min_absolute_target_gain = min_abs,
    min_target_relative_gain = min_rel,
    capped_target_redundancy = redundancy,
    min_regional_complete_fraction = 0.95,
    both_flanks_acceptable = 2L,
    total_length = total_length,
    candidate_id_key = candidate_key,
    panel_id = panel_id,
    stringsAsFactors = FALSE
  )
}

test_that("candidate universe preserves target evidence vectors and relative gains", {
  input <- synthetic_panel_inputs()
  universe <- build_candidate_universe(
    input$candidate_amplicons,
    input$candidate_amplicon_target_scores,
    input$candidate_amplicon_snps,
    input$config
  )

  expect_true(all(c("gain_45", "gain_66", "gain_72", "gain_80", "relative_gain_45", "relative_gain_66", "relative_gain_72", "relative_gain_80") %in% names(universe)))
  expect_equal(universe$gain_72[universe$candidate_id == "amp_C"], 0.90)
  expect_equal(universe$relative_gain_72[universe$candidate_id == "amp_C"], 1)
  expect_true(universe$gain_66[universe$candidate_id == "amp_B"] < 0.8)
})

test_that("filtering does not discard node 66 merely because it is below a universal 0.8 threshold", {
  input <- synthetic_panel_inputs()
  universe <- build_candidate_universe(
    input$candidate_amplicons,
    input$candidate_amplicon_target_scores,
    input$candidate_amplicon_snps,
    input$config
  )
  filtered <- filter_candidate_universe(universe, resolve_panel_selection_config(input$config))
  expect_true("amp_B" %in% filtered$candidate_id)
  expect_true(any(filtered$gain_66 < 0.8))
})

test_that("exact and approximate search modes switch at the configured limit", {
  input <- synthetic_panel_inputs()
  universe <- build_candidate_universe(
    input$candidate_amplicons,
    input$candidate_amplicon_target_scores,
    input$candidate_amplicon_snps,
    input$config
  )
  cfg <- resolve_panel_selection_config(input$config)
  exact <- evaluate_panel_size(universe, cfg, 2L)
  expect_identical(exact$search_method, "exact")
  expect_equal(exact$panels_evaluated, choose(nrow(universe), 2L))

  cfg$exact_enumeration_limit <- 1L
  approximate <- evaluate_panel_size(universe, cfg, 3L)
  expect_identical(approximate$search_method, "approximate")
  expect_gt(approximate$panels_evaluated, 0L)
  expect_lte(approximate$panels_evaluated, cfg$approximate_search_budget)
})

test_that("within-size Pareto logic respects maximisation and minimisation directions", {
  panels <- rbind(
    make_panel_metrics(2L, 0.8, 0.9, 3L, 200L, "A|B", "panel_k2_001"),
    make_panel_metrics(2L, 0.8, 0.9, 2L, 260L, "C|D", "panel_k2_002")
  )
  within <- identify_pareto_panels(panels, within_size = TRUE)
  expect_true(within[1L])
  expect_false(within[2L])
})

test_that("cross-size Pareto keeps panels that trade size against redundancy", {
  panels <- rbind(
    make_panel_metrics(2L, 0.8, 0.9, 2L, 200L, "A|B", "panel_k2_001"),
    make_panel_metrics(3L, 0.8, 0.9, 4L, 260L, "A|B|C", "panel_k3_001")
  )
  cross <- identify_pareto_panels(panels, within_size = FALSE)
  expect_true(all(cross))
})

test_that("a clear trade-off can yield multiple Pareto panels", {
  panels <- rbind(
    make_panel_metrics(2L, 0.8, 0.9, 2L, 200L, "A|B", "panel_k2_001"),
    make_panel_metrics(2L, 0.8, 0.9, 1L, 150L, "C|D", "panel_k2_002")
  )
  pareto <- identify_pareto_panels(panels, within_size = TRUE)
  expect_true(all(pareto))
})

test_that("dominance does not collapse panel-size trade-offs when size is an objective", {
  panel_a <- make_panel_metrics(2L, 0.8, 0.9, 2L, 200L, "A|B", "panel_k2_001")
  panel_b <- make_panel_metrics(3L, 0.8, 0.9, 4L, 260L, "A|B|C", "panel_k3_001")
  expect_false(dominates_objectives(
    panel_a,
    panel_b,
    c("min_absolute_target_gain", "min_target_relative_gain", "capped_target_redundancy", "min_regional_complete_fraction", "both_flanks_acceptable", "panel_size", "total_length"),
    c("max", "max", "max", "max", "max", "min", "min")
  ))
  expect_false(dominates_objectives(
    panel_b,
    panel_a,
    c("min_absolute_target_gain", "min_target_relative_gain", "capped_target_redundancy", "min_regional_complete_fraction", "both_flanks_acceptable", "panel_size", "total_length"),
    c("max", "max", "max", "max", "max", "min", "min")
  ))
})

test_that("redundancy cap saturates repeated support", {
  input <- synthetic_panel_inputs()
  cfg <- resolve_panel_selection_config(input$config)
  panel <- data.frame(
    candidate_id = c("x1", "x2", "x3", "x4"),
    start = c(1L, 101L, 201L, 301L),
    end = c(80L, 180L, 280L, 380L),
    length = c(80L, 80L, 80L, 80L),
    gain_45 = c(0.9, 0.9, 0.9, 0.9),
    gain_66 = c(0, 0, 0, 0),
    gain_72 = c(0, 0, 0, 0),
    gain_80 = c(0, 0, 0, 0),
    relative_gain_45 = c(1, 1, 1, 1),
    relative_gain_66 = c(0, 0, 0, 0),
    relative_gain_72 = c(0, 0, 0, 0),
    relative_gain_80 = c(0, 0, 0, 0),
    complete_fraction = c(0.95, 0.95, 0.95, 0.95),
    left_flank_conserved_stretch = TRUE,
    right_flank_conserved_stretch = TRUE,
    candidate_snp_positions_key = c("1", "2", "3", "4"),
    stringsAsFactors = FALSE
  )
  metrics <- panel_metrics(panel, cfg)
  expect_equal(metrics$capped_target_redundancy, 2L)
})

test_that("marginal improvement calculations are explicit", {
  input <- synthetic_panel_inputs()
  selection <- select_amplicon_panels(
    input$candidate_amplicons,
    input$candidate_amplicon_target_scores,
    input$candidate_amplicon_snps,
    input$config
  )

  expect_true(all(c(
    "delta_min_absolute_target_gain",
    "delta_min_target_relative_gain",
    "delta_capped_target_redundancy",
    "delta_number_of_targets_with_two_supporting_amplicons",
    "delta_total_length",
    "delta_min_regional_complete_fraction"
  ) %in% names(selection$panel_size_comparison)))
  expect_true(any(selection$panel_size_comparison$panel_size %in% c(2L, 3L)))
})

test_that("recommendation chooses the smallest size after discrimination and redundancy saturate", {
  input <- synthetic_panel_inputs()
  cfg <- resolve_panel_selection_config(input$config)
  comparison <- rbind(
    data.frame(panel_size = 2L, best_panel_id = "panel_k2_001", min_absolute_target_gain = 0.6667, min_target_relative_gain = 0.7778, capped_target_redundancy = 2L, delta_min_absolute_target_gain = NA_real_, delta_min_target_relative_gain = NA_real_, delta_capped_target_redundancy = NA_real_, stringsAsFactors = FALSE),
    data.frame(panel_size = 3L, best_panel_id = "panel_k3_001", min_absolute_target_gain = 0.6667, min_target_relative_gain = 0.7778, capped_target_redundancy = 4L, delta_min_absolute_target_gain = 0, delta_min_target_relative_gain = 0, delta_capped_target_redundancy = 2L, stringsAsFactors = FALSE),
    data.frame(panel_size = 4L, best_panel_id = "panel_k4_001", min_absolute_target_gain = 0.6667, min_target_relative_gain = 0.7778, capped_target_redundancy = 4L, delta_min_absolute_target_gain = 0, delta_min_target_relative_gain = 0, delta_capped_target_redundancy = 0L, stringsAsFactors = FALSE)
  )
  recommendation <- recommend_provisional_panel_size(comparison, cfg)
  expect_equal(recommendation$panel_size, 3L)
})

test_that("enumeration is deterministic and uses the same reduced candidate pool", {
  input <- synthetic_panel_inputs()
  selection_one <- select_amplicon_panels(
    input$candidate_amplicons,
    input$candidate_amplicon_target_scores,
    input$candidate_amplicon_snps,
    input$config
  )
  selection_two <- select_amplicon_panels(
    input$candidate_amplicons,
    input$candidate_amplicon_target_scores,
    input$candidate_amplicon_snps,
    input$config
  )

  expect_equal(selection_one$panel_candidates, selection_two$panel_candidates)
  expect_equal(selection_one$recommended_provisional_panel, selection_two$recommended_provisional_panel)
  expect_true(identical(selection_one$candidate_table$candidate_id, selection_two$candidate_table$candidate_id))
})

test_that("panel outputs do not depend on legacy weighted node score columns", {
  input <- synthetic_panel_inputs()
  selection <- select_amplicon_panels(
    input$candidate_amplicons,
    input$candidate_amplicon_target_scores,
    input$candidate_amplicon_snps,
    input$config
  )

  expect_false("total_weighted_gain" %in% names(selection$panel_candidates))
  expect_false(any(grepl("weighted", names(selection$panel_candidates))))
  expect_false(any(grepl("weighted", names(selection$recommended_provisional_panel))))
})
