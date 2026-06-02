source(test_path("../../R/encoding.R"))
source(test_path("../../R/preprocess.R"))
source(test_path("../../R/classifier.R"))
source(test_path("../../R/partials.R"))

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
