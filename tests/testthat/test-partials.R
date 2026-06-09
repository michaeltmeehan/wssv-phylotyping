source(test_path("../../R/encoding.R"))
source(test_path("../../R/preprocess.R"))
source(test_path("../../R/classifier.R"))
source(test_path("../../R/partials.R"))

skip_if_no_pairwise_backend <- function() {
  if (is.null(alignment_pkg())) {
    skip("Neither pwalign nor Biostrings pairwise alignment backend is available")
  }
}

make_partial_test_classifier <- function(...) {
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
  train_classifier(
    aln, target_mask, node_metadata, scores,
    use_selected_panel = FALSE,
    ...
  )
}

make_partial_panel_test_classifier <- function(...) {
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
    node_id = c("10", "11", "12"),
    site = c(2L, 4L, 2L),
    best_allele = c("A", "A", "C"),
    direction = c("clade", "clade", "clade"),
    normalized_gain = c(1, 1, 1)
  )
  selected_panel <- data.frame(
    start = 2L,
    end = 2L
  )
  train_classifier(
    aln, target_mask, node_metadata, scores,
    selected_panel = selected_panel,
    use_selected_panel = TRUE,
    ...
  )
}

test_that("partial FASTA files are read with metadata", {
  dir <- tempfile()
  dir.create(dir)
  fasta <- file.path(dir, "partials.fasta")
  writeLines(c(">NC_001 id text", "ACGT", ">sample two", "TT--"), fasta)

  records <- read_partial_fastas(dir, extensions = "fasta")

  expect_equal(nrow(records), 2L)
  expect_equal(records$sequence_id, c("NC_001", "sample"))
  expect_equal(records$length, c(4L, 4L))
  expect_equal(records$accession[1], "NC_001")
  expect_equal(records$source_file, rep(fasta, 2L))
})

test_that("empty partial input directory returns well-formed records", {
  dir <- tempfile()
  dir.create(dir)

  records <- read_partial_fastas(dir)

  expect_equal(nrow(records), 0L)
  expect_named(records, c("sequence_id", "source_file", "header", "accession", "sequence", "length"))
})

test_that("already aligned full-length partials are accepted", {
  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "-ACG--",
    length = 6L,
    stringsAsFactors = FALSE
  )

  mapped <- map_partial_sequences(records, alignment_length = 6L, mapping_mode = "aligned_full_length")

  expect_true(mapped$mapped)
  expect_equal(mapped$mapped_start, 2L)
  expect_equal(mapped$mapped_end, 4L)
  expect_equal(mapped$mapped_sequence, "-ACG--")
})

test_that("simple exact substrings are mapped against a reference sequence", {
  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "CGT",
    length = 3L,
    stringsAsFactors = FALSE
  )

  mapped <- map_partial_sequences(records, alignment_length = 6L, reference_sequence = "A-CGTA", mapping_mode = "exact_substring")

  expect_true(mapped$mapped)
  expect_equal(mapped$mapped_start, 3L)
  expect_equal(mapped$mapped_end, 5L)
  expect_equal(mapped$mapped_sequence, "--CGT-")
})

test_that("unsupported unaligned partials remain unmapped", {
  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "CCCC",
    length = 4L,
    stringsAsFactors = FALSE
  )

  mapped <- map_partial_sequences(records, alignment_length = 8L, reference_sequence = "AAAA----", mapping_mode = "auto")

  expect_false(mapped$mapped)
  expect_true(is.na(mapped$mapped_sequence))
})

test_that("pairwise local mapping reports orientation and coordinates", {
  skip_if_no_pairwise_backend()

  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "TAAC",
    length = 4L,
    stringsAsFactors = FALSE
  )

  mapped <- map_partial_sequences(
    records,
    alignment_length = 8L,
    reference_sequence = "AA-CGTTA",
    mapping_mode = "pairwise_local",
    pairwise_local_thresholds = list(min_identity = 0.75, min_query_coverage = 1, min_aligned_length = 4L)
  )

  expect_true(mapped$mapped)
  expect_equal(mapped$mapping_strand, "reverse_complement")
  expect_equal(mapped$mapped_start, 5L)
  expect_equal(mapped$mapped_end, 8L)
  expect_equal(mapped$mapping_reference_start, 4L)
  expect_equal(mapped$mapping_reference_end, 7L)
  expect_equal(mapped$mapping_aligned_length, 4L)
  expect_equal(mapped$mapping_query_coverage, 1)
  expect_equal(mapped$mapping_identity, 1)
})

test_that("pairwise local mapping handles forward, mismatched, and gapped local alignments", {
  skip_if_no_pairwise_backend()

  forward <- data.frame(
    sequence_id = "forward",
    source_file = "synthetic.fa",
    header = "forward",
    accession = NA_character_,
    sequence = "CGT",
    length = 3L,
    stringsAsFactors = FALSE
  )
  mismatched <- data.frame(
    sequence_id = "mismatched",
    source_file = "synthetic.fa",
    header = "mismatched",
    accession = NA_character_,
    sequence = "CGCTA",
    length = 5L,
    stringsAsFactors = FALSE
  )
  gapped <- data.frame(
    sequence_id = "gapped",
    source_file = "synthetic.fa",
    header = "gapped",
    accession = NA_character_,
    sequence = "CGTA",
    length = 4L,
    stringsAsFactors = FALSE
  )

  mapped_forward <- map_partial_sequences(
    forward,
    alignment_length = 8L,
    reference_sequence = "AA-CGTTA",
    mapping_mode = "pairwise_local",
    pairwise_local_thresholds = list(min_identity = 1, min_query_coverage = 1, min_aligned_length = 3L)
  )
  mapped_mismatched <- map_partial_sequences(
    mismatched,
    alignment_length = 8L,
    reference_sequence = "AA-CGTTA",
    mapping_mode = "pairwise_local",
    pairwise_local_thresholds = list(min_identity = 0.75, min_query_coverage = 1, min_aligned_length = 5L)
  )
  mapped_gapped <- map_partial_sequences(
    gapped,
    alignment_length = 8L,
    reference_sequence = "AA-CGTTA",
    mapping_mode = "pairwise_local",
    pairwise_local_thresholds = list(min_identity = 0.75, min_query_coverage = 1, min_aligned_length = 5L, gap_opening = 0)
  )

  expect_true(mapped_forward$mapped)
  expect_equal(mapped_forward$mapping_strand, "forward")
  expect_equal(mapped_forward$mapping_reference_start, 3L)
  expect_equal(mapped_forward$mapping_reference_end, 5L)
  expect_equal(mapped_forward$mapped_start, 4L)
  expect_equal(mapped_forward$mapped_end, 6L)

  expect_true(mapped_mismatched$mapped)
  expect_equal(mapped_mismatched$mapping_strand, "forward")
  expect_equal(mapped_mismatched$mapping_reference_start, 3L)
  expect_equal(mapped_mismatched$mapping_reference_end, 7L)
  expect_equal(mapped_mismatched$mapped_start, 4L)
  expect_equal(mapped_mismatched$mapped_end, 8L)
  expect_equal(mapped_mismatched$mapping_identity, 0.8)

  expect_true(mapped_gapped$mapped)
  expect_equal(mapped_gapped$mapping_reference_start, 3L)
  expect_equal(mapped_gapped$mapping_reference_end, 7L)
  expect_equal(mapped_gapped$mapped_start, 4L)
  expect_equal(mapped_gapped$mapped_end, 8L)
  expect_equal(mapped_gapped$mapping_identity, 0.8)
  expect_equal(mapped_gapped$mapping_query_coverage, 1)
})

test_that("pairwise local mapping rejects below identity threshold", {
  skip_if_no_pairwise_backend()

  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "CGCTA",
    length = 5L,
    stringsAsFactors = FALSE
  )

  mapped <- map_partial_sequences(
    records,
    alignment_length = 8L,
    reference_sequence = "AA-CGTTA",
    mapping_mode = "pairwise_local",
    pairwise_local_thresholds = list(min_identity = 0.99, min_query_coverage = 1, min_aligned_length = 5L)
  )

  expect_false(mapped$mapped)
  expect_true(is.na(mapped$mapping_identity))
  expect_match(mapped$mapping_note, "unmapped")
})

test_that("pairwise local mapping rejects below query coverage threshold", {
  skip_if_no_pairwise_backend()

  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "CGAAAA",
    length = 6L,
    stringsAsFactors = FALSE
  )

  mapped <- map_partial_sequences(
    records,
    alignment_length = 8L,
    reference_sequence = "AA-CGTTA",
    mapping_mode = "pairwise_local",
    pairwise_local_thresholds = list(min_identity = 1, min_query_coverage = 0.9, min_aligned_length = 2L)
  )

  expect_false(mapped$mapped)
  expect_true(is.na(mapped$mapping_query_coverage))
  expect_match(mapped$mapping_note, "unmapped")
})

test_that("auto mapping falls back from aligned length to exact substring to pairwise local", {
  skip_if_no_pairwise_backend()

  full <- data.frame(
    sequence_id = "full",
    source_file = "synthetic.fa",
    header = "full",
    accession = NA_character_,
    sequence = "-ACG--",
    length = 6L,
    stringsAsFactors = FALSE
  )
  exact <- data.frame(
    sequence_id = "exact",
    source_file = "synthetic.fa",
    header = "exact",
    accession = NA_character_,
    sequence = "CGT",
    length = 3L,
    stringsAsFactors = FALSE
  )
  local <- data.frame(
    sequence_id = "local",
    source_file = "synthetic.fa",
    header = "local",
    accession = NA_character_,
    sequence = "TAAC",
    length = 4L,
    stringsAsFactors = FALSE
  )

  mapped_full <- map_partial_sequences(full, alignment_length = 6L, reference_sequence = "A-CGTA", mapping_mode = "auto")
  mapped_exact <- map_partial_sequences(exact, alignment_length = 6L, reference_sequence = "A-CGTA", mapping_mode = "auto")
  mapped_local <- map_partial_sequences(local, alignment_length = 8L, reference_sequence = "AA-CGTTA", mapping_mode = "auto", pairwise_local_thresholds = list(min_identity = 0.75, min_query_coverage = 1, min_aligned_length = 4L))

  expect_true(mapped_full$mapped)
  expect_equal(mapped_full$mapping_note, "accepted aligned full-length sequence")
  expect_true(mapped_exact$mapped)
  expect_match(mapped_exact$mapping_note, "exact substring")
  expect_true(mapped_local$mapped)
  expect_match(mapped_local$mapping_note, "pairwise local alignment")
})

test_that("pairwise local mapping no longer calls subjectStart", {
  skip_if_no_pairwise_backend()

  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "TAAC",
    length = 4L,
    stringsAsFactors = FALSE
  )

  expect_silent(
    mapped <- map_partial_sequences(
      records,
      alignment_length = 8L,
      reference_sequence = "AA-CGTTA",
      mapping_mode = "pairwise_local",
      pairwise_local_thresholds = list(min_identity = 0.75, min_query_coverage = 1, min_aligned_length = 4L)
    )
  )
  expect_true(mapped$mapped)
})

test_that("partials implementation does not call unsafe Biostrings alignment accessors", {
  partials_source <- readLines(test_path("../../R/partials.R"), warn = FALSE)
  unsafe <- c(
    "Biostrings::subject(",
    "Biostrings::pattern(",
    "Biostrings::subjectStart(",
    "Biostrings::subjectEnd("
  )

  for (call in unsafe) {
    expect_false(any(grepl(call, partials_source, fixed = TRUE)), info = call)
  }
})

test_that("mapped partials are encoded with missing states", {
  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "ACGTN-",
    length = 6L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 6L, mapping_mode = "aligned_full_length")

  encoded <- encode_mapped_partials(mapped, alignment_length = 6L)

  expect_equal(unname(encoded[1, ]), c(1L, 2L, 3L, 4L, 0L, 0L))
})

test_that("partial with enough informative sites is classified", {
  classifier <- make_partial_test_classifier()
  records <- data.frame(
    sequence_id = "p1",
    source_file = "synthetic.fa",
    header = "p1",
    accession = NA_character_,
    sequence = "AAGA",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(mapped, classifier, min_total_informative_sites = 1L)

  expect_equal(out$classifications$status, "resolved")
  expect_equal(out$classifications$assigned_node, "11")
  expect_true(nrow(out$node_evidence) > 0L)
})

test_that("partials report no informative sites and unmapped statuses", {
  classifier <- make_partial_test_classifier()
  records <- data.frame(
    sequence_id = c("noinfo", "unmapped"),
    source_file = "synthetic.fa",
    header = c("noinfo", "unmapped"),
    accession = NA_character_,
    sequence = c("A---", "CCCC"),
    length = c(4L, 4L),
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, reference_sequence = "AAAA", mapping_mode = "aligned_full_length")
  mapped$mapped[2] <- FALSE
  mapped$mapped_sequence[2] <- NA_character_
  mapped$mapped_start[2] <- NA_integer_
  mapped$mapped_end[2] <- NA_integer_

  out <- classify_mapped_partials(mapped, classifier, min_total_informative_sites = 1L)

  expect_equal(out$classifications$status, c("no_informative_sites", "unmapped"))
  expect_equal(out$classifications$observed_informative_sites, c(0L, 0L))
})

test_that("opportunistic mode classifies a partial outside the selected panel", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = c("10", "11", "12"),
    site = c(2L, 4L, 2L),
    best_allele = c("A", "A", "C"),
    direction = c("clade", "clade", "clade"),
    normalized_gain = c(1, 1, 1)
  )
  records <- data.frame(
    sequence_id = "opp",
    source_file = "synthetic.fa",
    header = "opp",
    accession = NA_character_,
    sequence = "---A",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(
    mapped,
    classifier,
    classification_mode = "opportunistic",
    site_node_scores = full_scores,
    opportunistic_settings = list(min_observed_informative_sites = 1L)
  )

  expect_equal(out$classifications$status, "resolved_opportunistic")
  expect_equal(out$classifications$assigned_node, "11")
  expect_equal(out$classifications$sites_available_in_region, 1L)
  expect_equal(out$classifications$observed_informative_sites, 1L)
  expect_true(nrow(out$opportunistic_node_evidence) > 0L)
})

test_that("toy mapped partial with one scored site produces one site-evidence row", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = "11",
    site = 4L,
    best_allele = "A",
    direction = "clade",
    normalized_gain = 1
  )
  records <- data.frame(
    sequence_id = "one-site",
    source_file = "synthetic.fa",
    header = "one-site",
    accession = NA_character_,
    sequence = "---A",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(
    mapped,
    classifier,
    classification_mode = "opportunistic",
    site_node_scores = full_scores,
    opportunistic_settings = list(min_observed_informative_sites = 1L)
  )

  expect_equal(nrow(out$opportunistic_site_evidence), 1L)
  expect_equal(out$opportunistic_site_evidence$alignment_site, 4L)
  expect_true(out$opportunistic_site_evidence$supports_assigned_node)
  expect_equal(out$opportunistic_site_summary$unique_observed_scored_sites, out$classifications$observed_informative_sites)
})

test_that("toy mapped partial with several scored sites produces expected site-evidence rows", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = c("10", "11", "12"),
    site = c(2L, 4L, 2L),
    best_allele = c("A", "A", "C"),
    direction = c("clade", "clade", "clade"),
    normalized_gain = c(1, 1, 1)
  )
  records <- data.frame(
    sequence_id = "several-sites",
    source_file = "synthetic.fa",
    header = "several-sites",
    accession = NA_character_,
    sequence = "AA-A",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(
    mapped,
    classifier,
    classification_mode = "opportunistic",
    site_node_scores = full_scores,
    opportunistic_settings = list(min_observed_informative_sites = 1L)
  )

  expect_equal(nrow(out$opportunistic_site_evidence), 3L)
  expect_equal(length(unique(out$opportunistic_site_evidence$alignment_site)), 2L)
  expect_equal(out$opportunistic_site_summary$unique_observed_scored_sites, out$classifications$observed_informative_sites)
})

test_that("opportunistic mode reports no informative sites when interval has no scored SNPs", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = "11",
    site = 4L,
    best_allele = "A",
    direction = "clade",
    normalized_gain = 1
  )
  records <- data.frame(
    sequence_id = "empty-region",
    source_file = "synthetic.fa",
    header = "empty-region",
    accession = NA_character_,
    sequence = "A---",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(mapped, classifier, classification_mode = "opportunistic", site_node_scores = full_scores)

  expect_equal(out$classifications$status, "no_informative_sites")
  expect_equal(out$classifications$sites_available_in_region, 0L)
  expect_equal(nrow(out$opportunistic_site_evidence), 0L)
  expect_named(out$opportunistic_site_evidence, names(empty_partial_opportunistic_site_evidence_table()))
})

test_that("ambiguous bases at opportunistic informative sites are ignored", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = "11",
    site = 4L,
    best_allele = "A",
    direction = "clade",
    normalized_gain = 1
  )
  mapped <- data.frame(
    sequence_id = "ambiguous",
    source_file = "synthetic.fa",
    header = "ambiguous",
    accession = NA_character_,
    sequence = "---N",
    length = 4L,
    mapped = TRUE,
    mapped_start = 4L,
    mapped_end = 4L,
    mapped_sequence = "---N",
    mapping_note = "synthetic mapped interval",
    mapping_strand = "forward",
    mapping_identity = 1,
    mapping_query_coverage = 1,
    mapping_aligned_length = 1L,
    mapping_reference_start = 4L,
    mapping_reference_end = 4L,
    stringsAsFactors = FALSE
  )

  out <- classify_mapped_partials(mapped, classifier, classification_mode = "opportunistic", site_node_scores = full_scores)

  expect_equal(out$classifications$status, "no_observed_informative_sites")
  expect_equal(out$classifications$sites_available_in_region, 1L)
  expect_equal(out$classifications$observed_informative_sites, 0L)
  expect_equal(nrow(out$opportunistic_site_evidence), 0L)
  expect_named(out$opportunistic_site_evidence, names(empty_partial_opportunistic_site_evidence_table()))
})

test_that("one site supporting multiple nodes is represented explicitly", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = c("10", "11"),
    site = c(4L, 4L),
    best_allele = c("A", "A"),
    direction = c("clade", "clade"),
    normalized_gain = c(1, 1)
  )
  records <- data.frame(
    sequence_id = "multi-node-site",
    source_file = "synthetic.fa",
    header = "multi-node-site",
    accession = NA_character_,
    sequence = "---A",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(
    mapped,
    classifier,
    classification_mode = "opportunistic",
    site_node_scores = full_scores,
    opportunistic_settings = list(min_observed_informative_sites = 1L)
  )

  expect_equal(nrow(out$opportunistic_site_evidence), 2L)
  expect_equal(length(unique(out$opportunistic_site_evidence$alignment_site)), 1L)
  expect_equal(out$classifications$observed_informative_sites, 1L)
  expect_equal(out$opportunistic_site_summary$unique_observed_scored_sites, 1L)
})

test_that("site-level evidence row counts agree with classification summaries where expected", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = c("10", "11"),
    site = c(2L, 4L),
    best_allele = c("A", "A"),
    direction = c("clade", "clade"),
    normalized_gain = c(1, 1)
  )
  records <- data.frame(
    sequence_id = "summary-match",
    source_file = "synthetic.fa",
    header = "summary-match",
    accession = NA_character_,
    sequence = "AA-A",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(
    mapped,
    classifier,
    classification_mode = "opportunistic",
    site_node_scores = full_scores,
    opportunistic_settings = list(min_observed_informative_sites = 1L)
  )

  expect_equal(length(unique(out$opportunistic_site_evidence$alignment_site)), out$classifications$observed_informative_sites)
  expect_equal(out$opportunistic_site_summary$supporting_sites_for_assigned_node, out$classifications$supporting_sites)
})

test_that("panel classification mode preserves selected-panel behavior", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = "11",
    site = 4L,
    best_allele = "A",
    direction = "clade",
    normalized_gain = 1
  )
  records <- data.frame(
    sequence_id = "panel-only",
    source_file = "synthetic.fa",
    header = "panel-only",
    accession = NA_character_,
    sequence = "---A",
    length = 4L,
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(mapped, classifier, classification_mode = "panel", site_node_scores = full_scores)

  expect_equal(out$classifications$status, "no_informative_sites")
  expect_true(is.na(out$classifications$assigned_node))
  expect_equal(out$classifications$classification_source, "panel")
})

test_that("auto mode falls back only when panel sites are unavailable", {
  classifier <- make_partial_panel_test_classifier()
  full_scores <- data.frame(
    node_id = c("10", "11", "12"),
    site = c(2L, 4L, 2L),
    best_allele = c("A", "A", "C"),
    direction = c("clade", "clade", "clade"),
    normalized_gain = c(1, 1, 1)
  )
  records <- data.frame(
    sequence_id = c("fallback", "panel_signal"),
    source_file = "synthetic.fa",
    header = c("fallback", "panel_signal"),
    accession = NA_character_,
    sequence = c("---A", "AC-A"),
    length = c(4L, 4L),
    stringsAsFactors = FALSE
  )
  mapped <- map_partial_sequences(records, alignment_length = 4L, mapping_mode = "aligned_full_length")

  out <- classify_mapped_partials(mapped, classifier, classification_mode = "auto", site_node_scores = full_scores)

  expect_equal(out$classifications$classification_source, c("panel_then_opportunistic", "panel"))
  expect_equal(out$classifications$status[[1L]], "resolved_opportunistic")
  expect_equal(out$classifications$assigned_node[[1L]], "11")
  expect_equal(out$classifications$assigned_node[[2L]], "12")
})
