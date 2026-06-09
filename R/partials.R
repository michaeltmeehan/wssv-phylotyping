partial_default_extensions <- c("fa", "fas", "fasta", "fna")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

read_partial_fastas <- function(input_dir, extensions = partial_default_extensions) {
  if (!dir.exists(input_dir)) {
    return(empty_partial_records())
  }
  pattern <- paste0("\\.(", paste(gsub("^\\.", "", extensions), collapse = "|"), ")$")
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0L) {
    return(empty_partial_records())
  }

  records <- lapply(files, read_partial_fasta_file)
  out <- do.call(rbind, records)
  rownames(out) <- NULL
  out
}

read_partial_fasta_file <- function(path) {
  seqs <- read_fasta_character(path)
  if (length(seqs) == 0L) {
    return(empty_partial_records())
  }
  headers <- names(seqs)
  data.frame(
    sequence_id = make.unique(vapply(headers, partial_sequence_id, character(1L))),
    source_file = path,
    header = headers,
    accession = vapply(headers, parse_partial_accession, character(1L)),
    sequence = as.character(seqs),
    length = nchar(as.character(seqs), type = "bytes"),
    stringsAsFactors = FALSE
  )
}

empty_partial_records <- function() {
  data.frame(
    sequence_id = character(),
    source_file = character(),
    header = character(),
    accession = character(),
    sequence = character(),
    length = integer(),
    stringsAsFactors = FALSE
  )
}

partial_sequence_id <- function(header) {
  id <- strsplit(trimws(header), "\\s+")[[1L]][[1L]]
  if (is.na(id) || id == "") header else id
}

parse_partial_accession <- function(header) {
  first <- partial_sequence_id(header)
  patterns <- c(
    "^[A-Z]{1,3}_[0-9]+(\\.[0-9]+)?$",
    "^[A-Z]{1,4}[0-9]{5,}(\\.[0-9]+)?$"
  )
  if (any(vapply(patterns, function(p) grepl(p, first), logical(1L)))) {
    return(first)
  }
  m <- regexpr("([A-Z]{1,3}_[0-9]+|[A-Z]{1,4}[0-9]{5,})(\\.[0-9]+)?", header)
  if (m[[1L]] < 0L) NA_character_ else regmatches(header, m)
}

reference_sequence_from_alignment <- function(aln_int, reference_id = NULL) {
  if (!is.matrix(aln_int) || nrow(aln_int) == 0L) {
    stop("aln_int must be a non-empty integer matrix.", call. = FALSE)
  }
  row_i <- 1L
  if (!is.null(reference_id)) {
    match_i <- match(reference_id, rownames(aln_int))
    if (!is.na(match_i)) {
      row_i <- match_i
    }
  }
  chars <- base_code_to_allele(aln_int[row_i, ])
  chars[is.na(chars)] <- "-"
  paste0(chars, collapse = "")
}

reference_alignment_coordinate_map <- function(reference_sequence) {
  ref_chars <- strsplit(toupper(reference_sequence), "", fixed = TRUE)[[1L]]
  ref_bases <- ref_chars %in% c("A", "C", "G", "T")
  data.frame(
    reference_ungapped_position = seq_len(sum(ref_bases)),
    alignment_position = which(ref_bases),
    stringsAsFactors = FALSE
  )
}

map_partial_sequences <- function(partials, alignment_length, reference_sequence = NULL,
                                  mapping_mode = "auto",
                                  pairwise_local_thresholds = list()) {
  if (nrow(partials) == 0L) {
    partials$mapped <- logical()
    partials$mapped_start <- integer()
    partials$mapped_end <- integer()
    partials$mapped_sequence <- character()
    partials$mapping_mode <- character()
    partials$mapping_note <- character()
    partials$mapping_strand <- character()
    partials$mapping_identity <- numeric()
    partials$mapping_query_coverage <- numeric()
    partials$mapping_aligned_length <- integer()
    partials$mapping_reference_start <- integer()
    partials$mapping_reference_end <- integer()
    return(partials)
  }
  allowed <- c("aligned_full_length", "exact_substring", "pairwise_local", "auto")
  if (!mapping_mode %in% allowed) {
    stop("mapping_mode must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  pairwise_local_thresholds <- normalize_pairwise_local_thresholds(pairwise_local_thresholds)
  mapped <- lapply(seq_len(nrow(partials)), function(i) {
    map_partial_sequence(
      partials[i, , drop = FALSE],
      alignment_length,
      reference_sequence,
      mapping_mode,
      pairwise_local_thresholds = pairwise_local_thresholds
    )
  })
  out <- cbind(partials, do.call(rbind, mapped))
  rownames(out) <- NULL
  out
}

map_partial_sequence <- function(record, alignment_length, reference_sequence = NULL,
                                 mapping_mode = "auto",
                                 pairwise_local_thresholds = list()) {
  seq <- toupper(gsub("\\s+", "", record$sequence[[1L]]))
  if (mapping_mode %in% c("aligned_full_length", "auto") && nchar(seq, type = "bytes") == alignment_length) {
    encoded <- encode_classifier_query(seq)
    observed <- which(encoded != 0L)
    return(data.frame(
      mapped = TRUE,
      mapped_start = if (length(observed) > 0L) min(observed) else NA_integer_,
      mapped_end = if (length(observed) > 0L) max(observed) else NA_integer_,
      mapped_sequence = seq,
      mapping_mode = "aligned_full_length",
      mapping_note = "accepted aligned full-length sequence",
      mapping_strand = "forward",
      mapping_identity = 1,
      mapping_query_coverage = 1,
      mapping_aligned_length = nchar(seq, type = "bytes"),
      mapping_reference_start = 1L,
      mapping_reference_end = alignment_length,
      stringsAsFactors = FALSE
    ))
  }

  if (mapping_mode %in% c("exact_substring", "auto") && !is.null(reference_sequence)) {
    placed <- map_exact_substring(seq, reference_sequence, alignment_length)
    if (!is.null(placed)) {
      return(placed)
    }
  }

  if (mapping_mode %in% c("pairwise_local", "auto") && !is.null(reference_sequence)) {
    placed <- map_pairwise_local(
      seq,
      reference_sequence,
      alignment_length,
      thresholds = pairwise_local_thresholds
    )
    if (!is.null(placed)) {
      return(placed)
    }
  }

  data.frame(
    mapped = FALSE,
    mapped_start = NA_integer_,
    mapped_end = NA_integer_,
    mapped_sequence = NA_character_,
    mapping_mode = NA_character_,
    mapping_strand = NA_character_,
    mapping_identity = NA_real_,
    mapping_query_coverage = NA_real_,
    mapping_aligned_length = NA_integer_,
    mapping_reference_start = NA_integer_,
    mapping_reference_end = NA_integer_,
    mapping_note = "unmapped: unsupported length or no accepted reference mapping",
    stringsAsFactors = FALSE
  )
}

map_exact_substring <- function(sequence, reference_sequence, alignment_length = nchar(reference_sequence, type = "bytes")) {
  query <- gsub("[-.]", "", toupper(sequence))
  if (nchar(query, type = "bytes") == 0L) {
    return(NULL)
  }
  if (grepl("[^ACGT]", query)) {
    return(NULL)
  }

  ref_chars <- strsplit(toupper(reference_sequence), "", fixed = TRUE)[[1L]]
  ref_bases <- ref_chars %in% c("A", "C", "G", "T")
  ref_ungapped <- paste0(ref_chars[ref_bases], collapse = "")
  hit <- regexpr(query, ref_ungapped, fixed = TRUE)
  if (hit[[1L]] < 0L) {
    return(NULL)
  }
  ungapped_start <- as.integer(hit[[1L]])
  ungapped_end <- ungapped_start + nchar(query, type = "bytes") - 1L
  ref_map <- reference_alignment_coordinate_map(reference_sequence)
  ref_cols <- ref_map$alignment_position
  aln_start <- ref_cols[[ungapped_start]]
  aln_end <- ref_cols[[ungapped_end]]

  full <- rep("-", alignment_length)
  full[ref_cols[ungapped_start:ungapped_end]] <- strsplit(query, "", fixed = TRUE)[[1L]]
  data.frame(
    mapped = TRUE,
    mapped_start = aln_start,
    mapped_end = aln_end,
    mapped_sequence = paste0(full, collapse = ""),
    mapping_mode = "exact_substring",
    mapping_strand = "forward",
    mapping_identity = 1,
    mapping_query_coverage = nchar(query, type = "bytes") / nchar(query, type = "bytes"),
    mapping_aligned_length = nchar(query, type = "bytes"),
    mapping_reference_start = ungapped_start,
    mapping_reference_end = ungapped_end,
    mapping_note = "mapped by exact substring against reference sequence",
    stringsAsFactors = FALSE
  )
}

normalize_pairwise_local_thresholds <- function(thresholds) {
  if (is.null(thresholds)) {
    thresholds <- list()
  }
  defaults <- list(
    min_identity = 0.85,
    min_query_coverage = 0.6,
    min_aligned_length = 20L,
    gap_opening = -10,
    gap_extension = -0.5
  )
  out <- modifyList(defaults, thresholds)
  out$min_identity <- as.numeric(out$min_identity)
  out$min_query_coverage <- as.numeric(out$min_query_coverage)
  out$min_aligned_length <- as.integer(out$min_aligned_length)
  out$gap_opening <- as.numeric(out$gap_opening)
  out$gap_extension <- as.numeric(out$gap_extension)
  out
}

alignment_pkg <- function() {
  if (requireNamespace("pwalign", quietly = TRUE)) {
    return("pwalign")
  }
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    biostrings_exports <- c("pairwiseAlignment", "nucleotideSubstitutionMatrix", "alignedPattern", "alignedSubject")
    if (all(vapply(biostrings_exports, function(x) exists(x, envir = asNamespace("Biostrings"), inherits = FALSE), logical(1L)))) {
      return("Biostrings")
    }
  }
  NULL
}

pairwise_alignment <- function(..., pkg = alignment_pkg()) {
  if (is.null(pkg)) {
    stop("Neither pwalign nor the required Biostrings pairwise-alignment functions are available.", call. = FALSE)
  }
  getExportedValue(pkg, "pairwiseAlignment")(...)
}

nucleotide_substitution_matrix <- function(..., pkg = alignment_pkg()) {
  if (is.null(pkg)) {
    stop("Neither pwalign nor the required Biostrings pairwise-alignment functions are available.", call. = FALSE)
  }
  getExportedValue(pkg, "nucleotideSubstitutionMatrix")(...)
}

aligned_pattern <- function(aln, pkg = alignment_pkg()) {
  if (is.null(pkg)) {
    stop("Neither pwalign nor the required Biostrings pairwise-alignment functions are available.", call. = FALSE)
  }
  as.character(getExportedValue(pkg, "alignedPattern")(aln))
}

aligned_subject <- function(aln, pkg = alignment_pkg()) {
  if (is.null(pkg)) {
    stop("Neither pwalign nor the required Biostrings pairwise-alignment functions are available.", call. = FALSE)
  }
  as.character(getExportedValue(pkg, "alignedSubject")(aln))
}

alignment_score <- function(aln, pkg = alignment_pkg()) {
  if (is.null(pkg) || !"score" %in% getNamespaceExports(pkg)) {
    return(NA_real_)
  }
  as.numeric(getExportedValue(pkg, "score")(aln))
}

anchor_aligned_reference <- function(aligned_reference, reference_ungapped) {
  ref_fragment <- gsub("-", "", toupper(aligned_reference), fixed = TRUE)
  if (nchar(ref_fragment, type = "bytes") == 0L) {
    return(list(start = NA_integer_, end = NA_integer_, ambiguous = TRUE))
  }

  hits <- gregexpr(ref_fragment, reference_ungapped, fixed = TRUE)[[1L]]
  if (identical(hits, -1L)) {
    return(list(start = NA_integer_, end = NA_integer_, ambiguous = FALSE))
  }
  if (length(hits) > 1L) {
    return(list(start = NA_integer_, end = NA_integer_, ambiguous = TRUE))
  }
  start <- as.integer(hits[[1L]])
  list(
    start = start,
    end = start + nchar(ref_fragment, type = "bytes") - 1L,
    ambiguous = FALSE
  )
}

pairwise_local_alignment_summary <- function(query, reference_ungapped, strand,
                                             thresholds, pkg = alignment_pkg()) {
  aln <- tryCatch(
    pairwise_alignment(
      Biostrings::DNAString(query),
      Biostrings::DNAString(reference_ungapped),
      type = "local",
      substitutionMatrix = nucleotide_substitution_matrix(match = 1, mismatch = -1, baseOnly = TRUE, pkg = pkg),
      gapOpening = thresholds$gap_opening,
      gapExtension = thresholds$gap_extension,
      scoreOnly = FALSE,
      pkg = pkg
    ),
    error = function(e) NULL
  )
  if (is.null(aln)) {
    return(NULL)
  }

  pattern_aln <- aligned_pattern(aln, pkg = pkg)
  subject_aln <- aligned_subject(aln, pkg = pkg)
  aligned_query <- strsplit(pattern_aln, "", fixed = TRUE)[[1L]]
  aligned_reference <- strsplit(subject_aln, "", fixed = TRUE)[[1L]]
  if (length(aligned_query) == 0L || length(aligned_reference) == 0L) {
    return(NULL)
  }

  query_positions <- which(aligned_query != "-")
  comparable <- aligned_query != "-" & aligned_reference != "-"
  matched <- sum(aligned_query[comparable] == aligned_reference[comparable])
  aligned_length <- length(aligned_reference)
  identity <- if (aligned_length > 0L) matched / aligned_length else NA_real_
  query_coverage <- if (nchar(query, type = "bytes") > 0L) length(query_positions) / nchar(query, type = "bytes") else NA_real_
  anchored <- anchor_aligned_reference(subject_aln, reference_ungapped)
  if (isTRUE(anchored$ambiguous) || is.na(anchored$start) || is.na(anchored$end)) {
    return(NULL)
  }

  list(
    aligned_query = pattern_aln,
    aligned_reference = subject_aln,
    reference_start = anchored$start,
    reference_end = anchored$end,
    strand = strand,
    identity = identity,
    query_coverage = query_coverage,
    aligned_length = as.integer(aligned_length),
    score = alignment_score(aln, pkg = pkg)
  )
}

map_pairwise_local <- function(sequence, reference_sequence, alignment_length,
                               thresholds = list()) {
  pkg <- alignment_pkg()
  if (is.null(pkg)) {
    return(NULL)
  }
  query <- gsub("[-.\\s]", "", toupper(sequence))
  if (nchar(query, type = "bytes") == 0L) {
    return(NULL)
  }

  ref_chars <- strsplit(toupper(reference_sequence), "", fixed = TRUE)[[1L]]
  ref_ungapped <- paste0(ref_chars[ref_chars %in% c("A", "C", "G", "T")], collapse = "")
  if (nchar(ref_ungapped, type = "bytes") == 0L) {
    return(NULL)
  }

  thresholds <- normalize_pairwise_local_thresholds(thresholds)

  ref_map <- reference_alignment_coordinate_map(reference_sequence)
  best <- NULL
  best_score <- NA_real_
  if (as.numeric(nchar(query, type = "bytes")) * as.numeric(nchar(ref_ungapped, type = "bytes")) > .Machine$integer.max) {
    return(NULL)
  }
  for (strand in c("forward", "reverse_complement")) {
    q <- if (strand == "forward") query else as.character(Biostrings::reverseComplement(Biostrings::DNAString(query)))
    summary <- pairwise_local_alignment_summary(q, ref_ungapped, strand, thresholds, pkg = pkg)
    if (is.null(summary)) next
    aligned_pattern <- strsplit(summary$aligned_query, "", fixed = TRUE)[[1L]]
    aligned_subject <- strsplit(summary$aligned_reference, "", fixed = TRUE)[[1L]]
    identity <- summary$identity
    query_coverage <- summary$query_coverage
    aligned_length <- summary$aligned_length
    if (is.na(identity) || is.na(query_coverage) || aligned_length < as.integer(thresholds$min_aligned_length) ||
      identity < thresholds$min_identity || query_coverage < thresholds$min_query_coverage) {
      next
    }

    subject_start <- summary$reference_start
    subject_end <- summary$reference_end
    ref_positions <- seq.int(subject_start, subject_end)
    ref_alignment_cols <- ref_map$alignment_position[match(ref_positions, ref_map$reference_ungapped_position)]
    if (any(is.na(ref_alignment_cols))) next
    mapped_sequence <- rep("-", alignment_length)
    query_used <- if (strand == "forward") query else as.character(Biostrings::reverseComplement(Biostrings::DNAString(query)))
    query_chars <- strsplit(query_used, "", fixed = TRUE)[[1L]]
    subject_non_gap_cols <- which(aligned_subject != "-")
    shared_cols <- intersect(which(aligned_pattern != "-"), subject_non_gap_cols)
    if (length(shared_cols) == 0L) next
    shared_query_idx <- cumsum(aligned_pattern != "-")[shared_cols]
    shared_subject_idx <- cumsum(aligned_subject != "-")[shared_cols]
    ref_shared_cols <- ref_alignment_cols[shared_subject_idx]
    if (any(is.na(ref_shared_cols)) || any(is.na(shared_query_idx))) next
    mapped_sequence[ref_shared_cols] <- query_chars[shared_query_idx]
    candidate <- data.frame(
      mapped = TRUE,
      mapped_start = min(ref_shared_cols, na.rm = TRUE),
      mapped_end = max(ref_shared_cols, na.rm = TRUE),
      mapped_sequence = paste0(mapped_sequence, collapse = ""),
      mapping_mode = "pairwise_local",
      mapping_strand = strand,
      mapping_identity = identity,
      mapping_query_coverage = query_coverage,
      mapping_aligned_length = as.integer(aligned_length),
      mapping_reference_start = as.integer(subject_start),
      mapping_reference_end = as.integer(subject_end),
      mapping_note = paste0("mapped by pairwise local alignment on ", strand, " using ", pkg),
      stringsAsFactors = FALSE
    )
    if (is.null(best) || candidate$mapping_identity[[1L]] > best$mapping_identity[[1L]] ||
      (candidate$mapping_identity[[1L]] == best$mapping_identity[[1L]] &&
        candidate$mapping_query_coverage[[1L]] > best$mapping_query_coverage[[1L]]) ||
      (candidate$mapping_identity[[1L]] == best$mapping_identity[[1L]] &&
        candidate$mapping_query_coverage[[1L]] == best$mapping_query_coverage[[1L]] &&
        isTRUE(summary$score > best_score))) {
      best <- candidate
      best_score <- summary$score
    }
  }
  best
}

encode_mapped_partials <- function(mapped_partials, alignment_length = NULL) {
  if (nrow(mapped_partials) == 0L) {
    ncol <- if (is.null(alignment_length)) 0L else as.integer(alignment_length)
    return(matrix(integer(), nrow = 0L, ncol = ncol))
  }
  if (is.null(alignment_length)) {
    alignment_length <- nchar(mapped_partials$mapped_sequence[which(mapped_partials$mapped)[[1L]]], type = "bytes")
  }
  encoded <- matrix(0L, nrow = nrow(mapped_partials), ncol = alignment_length)
  for (i in seq_len(nrow(mapped_partials))) {
    if (isTRUE(mapped_partials$mapped[[i]])) {
      encoded[i, ] <- encode_classifier_query(mapped_partials$mapped_sequence[[i]])
    }
  }
  rownames(encoded) <- mapped_partials$sequence_id
  colnames(encoded) <- as.character(seq_len(alignment_length))
  encoded
}

classify_mapped_partials <- function(mapped_partials, classifier,
                                     min_total_informative_sites = NULL,
                                     write_unresolved_evidence = FALSE,
                                     classification_mode = "panel",
                                     site_node_scores = NULL,
                                     opportunistic_settings = list()) {
  classification_mode <- normalize_partial_classification_mode(classification_mode)
  empty <- empty_partial_classification_table()
  if (nrow(mapped_partials) == 0L) {
    return(list(
      classifications = empty,
      node_evidence = empty_partial_node_evidence_table(),
      panel_classifications = empty,
      opportunistic_classifications = empty,
      opportunistic_node_evidence = empty_partial_opportunistic_node_evidence_table(),
      opportunistic_site_evidence = empty_partial_opportunistic_site_evidence_table(),
      opportunistic_site_summary = empty_partial_opportunistic_site_summary_table(),
      region_diagnostics = empty_partial_region_diagnostics_table()
    ))
  }

  panel <- classify_mapped_partials_panel(
    mapped_partials,
    classifier,
    min_total_informative_sites = min_total_informative_sites,
    write_unresolved_evidence = write_unresolved_evidence
  )
  if (classification_mode == "panel") {
    return(list(
      classifications = panel$classifications,
      node_evidence = panel$node_evidence,
      panel_classifications = panel$classifications,
      opportunistic_classifications = empty,
      opportunistic_node_evidence = empty_partial_opportunistic_node_evidence_table(),
      opportunistic_site_evidence = empty_partial_opportunistic_site_evidence_table(),
      opportunistic_site_summary = empty_partial_opportunistic_site_summary_table(),
      region_diagnostics = empty_partial_region_diagnostics_table()
    ))
  }

  if (is.null(site_node_scores)) {
    stop("site_node_scores is required for opportunistic partial classification.", call. = FALSE)
  }
  opportunistic <- classify_mapped_partials_opportunistic(
    mapped_partials,
    classifier,
    site_node_scores,
    min_total_informative_sites = min_total_informative_sites,
    settings = opportunistic_settings
  )
  if (classification_mode == "opportunistic") {
    return(list(
      classifications = opportunistic$classifications,
      node_evidence = opportunistic$node_evidence,
      panel_classifications = panel$classifications,
      opportunistic_classifications = opportunistic$classifications,
      opportunistic_node_evidence = opportunistic$node_evidence,
      opportunistic_site_evidence = opportunistic$site_evidence,
      opportunistic_site_summary = opportunistic$site_summary,
      region_diagnostics = opportunistic$region_diagnostics
    ))
  }

  combined <- panel$classifications
  fallback <- combined$mapped & combined$observed_informative_sites == 0L
  if (any(fallback)) {
    combined[fallback, ] <- harmonize_partial_classification_columns(opportunistic$classifications[fallback, , drop = FALSE], combined)
    combined$classification_source[fallback] <- "panel_then_opportunistic"
  }
  list(
    classifications = combined,
    node_evidence = panel$node_evidence,
    panel_classifications = panel$classifications,
    opportunistic_classifications = opportunistic$classifications,
    opportunistic_node_evidence = opportunistic$node_evidence,
    opportunistic_site_evidence = opportunistic$site_evidence,
    opportunistic_site_summary = opportunistic$site_summary,
    region_diagnostics = opportunistic$region_diagnostics
  )
}

classify_mapped_partials_panel <- function(mapped_partials, classifier,
                                           min_total_informative_sites = NULL,
                                           write_unresolved_evidence = FALSE) {
  empty <- empty_partial_classification_table()
  if (nrow(mapped_partials) == 0L) {
    return(list(classifications = empty, node_evidence = empty_partial_node_evidence_table()))
  }

  encoded <- encode_mapped_partials(mapped_partials, classifier$alignment_length)
  rows <- vector("list", nrow(mapped_partials))
  evidence_rows <- list()
  for (i in seq_len(nrow(mapped_partials))) {
    record <- mapped_partials[i, , drop = FALSE]
    observed_non_missing <- sum(encoded[i, ] != 0L)
    if (!isTRUE(record$mapped[[1L]])) {
      rows[[i]] <- partial_classification_row(record, observed_non_missing, 0L, 0L, 0L, 0L, NA_character_, NA_real_, NA_real_, "unmapped", record$mapping_note, "panel", "panel")
      next
    }

    pred <- classify_tree_path(
      encoded[i, ],
      classifier,
      min_total_informative_sites = min_total_informative_sites
    )
    note <- partial_status_note(pred$status)
    rows[[i]] <- partial_classification_row(
      record,
      observed_non_missing,
      length(classifier$informative_sites[classifier$informative_sites >= record$mapped_start[[1L]] & classifier$informative_sites <= record$mapped_end[[1L]]]),
      pred$observed_informative_sites,
      {
        assigned_ev <- if (!is.na(pred$assigned_node)) pred$evidence[pred$evidence$node_id == pred$assigned_node, , drop = FALSE] else pred$evidence[0, , drop = FALSE]
        if (nrow(assigned_ev) > 0L) assigned_ev$support_sites[[1L]] else 0L
      },
      {
        assigned_ev <- if (!is.na(pred$assigned_node)) pred$evidence[pred$evidence$node_id == pred$assigned_node, , drop = FALSE] else pred$evidence[0, , drop = FALSE]
        if (nrow(assigned_ev) > 0L) assigned_ev$conflict_sites[[1L]] else 0L
      },
      pred$assigned_node,
      pred$support_score,
      prediction_support_margin(pred),
      pred$status,
      note,
      "panel",
      "panel"
    )
    if (isTRUE(write_unresolved_evidence) || pred$status %in% c("resolved", "weak_support", "conflicting")) {
      ev <- pred$evidence
      ev <- ev[!is.na(ev$observed_sites) & ev$observed_sites > 0L, , drop = FALSE]
      if (nrow(ev) > 0L) {
        ev$sequence_id <- record$sequence_id[[1L]]
        ev$source_file <- record$source_file[[1L]]
        evidence_rows[[length(evidence_rows) + 1L]] <- ev
      }
    }
  }
  classifications <- do.call(rbind, rows)
  evidence <- if (length(evidence_rows) > 0L) do.call(rbind, evidence_rows) else empty_partial_node_evidence_table()
  if (nrow(evidence) > 0L) {
    evidence <- evidence[, c("sequence_id", "source_file", setdiff(names(evidence), c("sequence_id", "source_file"))), drop = FALSE]
  }
  list(classifications = classifications, node_evidence = evidence)
}

classify_mapped_partials_opportunistic <- function(mapped_partials, classifier,
                                                   site_node_scores,
                                                   min_total_informative_sites = NULL,
                                                   settings = list()) {
  settings <- normalize_opportunistic_settings(settings, classifier, min_total_informative_sites)
  opportunistic_classifier <- make_opportunistic_classifier(classifier, site_node_scores, settings$use_weighted_support)
  encoded <- encode_mapped_partials(mapped_partials, classifier$alignment_length)
  rows <- vector("list", nrow(mapped_partials))
  evidence_rows <- list()
  diagnostic_rows <- vector("list", nrow(mapped_partials))
  site_evidence_rows <- list()

  for (i in seq_len(nrow(mapped_partials))) {
    record <- mapped_partials[i, , drop = FALSE]
    observed_non_missing <- sum(encoded[i, ] != 0L)
    if (!isTRUE(record$mapped[[1L]])) {
      rows[[i]] <- partial_classification_row(record, observed_non_missing, 0L, 0L, 0L, 0L, NA_character_, NA_real_, NA_real_, "unmapped", record$mapping_note, "opportunistic", "opportunistic")
      diagnostic_rows[[i]] <- partial_region_diagnostics_row(record, 0L, 0L, 0L)
      next
    }

    interval_rules <- opportunistic_classifier$rules[
      opportunistic_classifier$rules$site >= record$mapped_start[[1L]] &
        opportunistic_classifier$rules$site <= record$mapped_end[[1L]],
      ,
      drop = FALSE
    ]
    sites_available <- length(unique(interval_rules$site))
    observed_sites <- extract_observed_scored_sites(encoded[i, ], interval_rules)
    observed_count <- length(unique(observed_sites$site))
    diagnostic_rows[[i]] <- partial_region_diagnostics_row(record, sites_available, observed_count, observed_non_missing)

    if (sites_available == 0L) {
      rows[[i]] <- partial_classification_row(record, observed_non_missing, 0L, 0L, 0L, 0L, NA_character_, NA_real_, NA_real_, "no_informative_sites", "mapped interval has no scored informative SNPs", "opportunistic", "opportunistic")
      next
    }
    if (observed_count == 0L) {
      rows[[i]] <- partial_classification_row(record, observed_non_missing, sites_available, 0L, 0L, 0L, NA_character_, NA_real_, NA_real_, "no_observed_informative_sites", "scored informative SNPs were present but missing, ambiguous, or gapped in the query", "opportunistic", "opportunistic")
      next
    }
    if (observed_count < settings$min_observed_informative_sites) {
      rows[[i]] <- partial_classification_row(record, observed_non_missing, sites_available, observed_count, 0L, 0L, NA_character_, NA_real_, NA_real_, "weak_support", "too few observed scored informative SNPs in the mapped interval", "opportunistic", "opportunistic")
      next
    }

    pred <- classify_tree_path(
      encoded[i, ],
      opportunistic_classifier,
      min_total_informative_sites = settings$min_observed_informative_sites,
      min_support = settings$min_support,
      max_conflict = settings$conflict_support_threshold,
      support_margin = settings$min_support_margin
    )
    status <- if (pred$status == "resolved") "resolved_opportunistic" else pred$status
    note <- partial_status_note(status)
    assigned_ev <- if (!is.na(pred$assigned_node)) pred$evidence[pred$evidence$node_id == pred$assigned_node, , drop = FALSE] else pred$evidence[0, , drop = FALSE]
    rows[[i]] <- partial_classification_row(
      record,
      observed_non_missing,
      sites_available,
      observed_count,
      if (nrow(assigned_ev) > 0L) assigned_ev$support_sites[[1L]] else 0L,
      if (nrow(assigned_ev) > 0L) assigned_ev$conflict_sites[[1L]] else 0L,
      pred$assigned_node,
      pred$support_score,
      prediction_support_margin(pred),
      status,
      note,
      "opportunistic",
      "opportunistic"
    )
    ev <- summarise_opportunistic_node_evidence(record, pred$evidence, pred$assigned_node, sites_available, observed_count, opportunistic_classifier$target_mask)
    if (nrow(ev) > 0L) {
      evidence_rows[[length(evidence_rows) + 1L]] <- ev
    }
  }

  classifications <- do.call(rbind, rows)
  for (i in seq_len(nrow(mapped_partials))) {
    ev <- build_partial_opportunistic_site_evidence(
      mapped_partials[i, , drop = FALSE],
      classifications[i, , drop = FALSE],
      opportunistic_classifier
    )
    if (nrow(ev) > 0L) {
      site_evidence_rows[[length(site_evidence_rows) + 1L]] <- ev
    }
  }
  site_evidence <- if (length(site_evidence_rows) > 0L) {
    do.call(rbind, site_evidence_rows)
  } else {
    empty_partial_opportunistic_site_evidence_table()
  }

  list(
    classifications = classifications,
    node_evidence = if (length(evidence_rows) > 0L) do.call(rbind, evidence_rows) else empty_partial_opportunistic_node_evidence_table(),
    site_evidence = site_evidence,
    site_summary = summarise_partial_opportunistic_site_evidence(site_evidence, classifications),
    region_diagnostics = do.call(rbind, diagnostic_rows)
  )
}

build_partial_opportunistic_site_evidence <- function(record, classification,
                                                      opportunistic_classifier,
                                                      reference_sequence = NULL) {
  empty <- empty_partial_opportunistic_site_evidence_table()
  if (nrow(record) == 0L || !isTRUE(record$mapped[[1L]])) {
    return(empty)
  }
  rules <- opportunistic_classifier$rules
  if (nrow(rules) == 0L) {
    return(empty)
  }
  interval_rules <- rules[
    rules$site >= record$mapped_start[[1L]] &
      rules$site <= record$mapped_end[[1L]],
    ,
    drop = FALSE
  ]
  if (nrow(interval_rules) == 0L) {
    return(empty)
  }

  x <- encode_classifier_query(record$mapped_sequence[[1L]])
  interval_rules <- interval_rules[interval_rules$site <= length(x), , drop = FALSE]
  query_code <- x[interval_rules$site]
  observed <- query_code != 0L
  if (!any(observed)) {
    return(empty)
  }
  interval_rules <- interval_rules[observed, , drop = FALSE]
  query_code <- query_code[observed]
  query_base <- base_code_to_allele(query_code)
  supports <- ifelse(
    interval_rules$direction == "clade",
    query_code == interval_rules$allele_code,
    query_code != interval_rules$allele_code
  )

  assigned_node <- as.character(classification$assigned_node[[1L]])
  if (is.na(assigned_node) || assigned_node == "") {
    assigned_node <- NA_character_
  }
  node_id <- as.character(interval_rules$node_id)
  supports_assigned_node <- !is.na(assigned_node) & node_id == assigned_node & supports
  supports_assigned_path <- if (is.na(assigned_node)) {
    rep(NA, length(node_id))
  } else {
    supports & nodes_nested(assigned_node, node_id, opportunistic_classifier$target_mask)
  }
  contradicts_assigned_node <- !is.na(assigned_node) & node_id == assigned_node & !supports
  on_assigned_path <- if (is.na(assigned_node)) {
    rep(NA, length(node_id))
  } else {
    nodes_nested(assigned_node, node_id, opportunistic_classifier$target_mask)
  }
  role <- ifelse(
    supports_assigned_node,
    "supporting",
    ifelse(!is.na(on_assigned_path) & on_assigned_path & supports, "supporting_path",
      ifelse(supports, "off_path", "uninformative_observed")
    )
  )
  site_reason <- ifelse(
    supports,
    "query base matches this node rule",
    "query base does not match this node rule"
  )
  if (any(interval_rules$direction == "outside")) {
    site_reason[interval_rules$direction == "outside" & supports] <- "query base differs from outside-directed informative allele"
    site_reason[interval_rules$direction == "outside" & !supports] <- "query base matches outside-directed informative allele"
  }
  reference_base <- rep(NA_character_, nrow(interval_rules))
  if (!is.null(reference_sequence)) {
    ref <- strsplit(toupper(reference_sequence), "", fixed = TRUE)[[1L]]
    ok <- interval_rules$site <= length(ref)
    reference_base[ok] <- ref[interval_rules$site[ok]]
  }
  assigned_label <- classification$assigned_node_label[[1L]]
  if (length(assigned_label) == 0L || is.na(assigned_label) || assigned_label == "") {
    assigned_label <- assigned_node
  }
  node_label <- if ("node_label" %in% names(interval_rules)) {
    as.character(interval_rules$node_label)
  } else {
    node_id
  }

  out <- data.frame(
    record_id = record$sequence_id[[1L]],
    classification_mode = "opportunistic",
    classification_source = classification$classification_source[[1L]],
    mapped = as.logical(record$mapped[[1L]]),
    mapping_mode = as.character(record$mapping_mode[[1L]]),
    strand = as.character(record$mapping_strand[[1L]]),
    alignment_start = as.integer(record$mapped_start[[1L]]),
    alignment_end = as.integer(record$mapped_end[[1L]]),
    alignment_site = as.integer(interval_rules$site),
    query_base = query_base,
    reference_base = reference_base,
    node = node_id,
    node_label = node_label,
    node_allele = as.character(interval_rules$best_allele),
    site_score = if ("normalized_gain" %in% names(interval_rules)) as.numeric(interval_rules$normalized_gain) else as.numeric(interval_rules$rule_weight),
    site_weight = as.numeric(interval_rules$rule_weight),
    supports_node = as.logical(supports),
    supports_assigned_node = as.logical(supports_assigned_node),
    supports_assigned_path = as.logical(supports_assigned_path),
    contradicts_assigned_node = as.logical(contradicts_assigned_node),
    assigned_node = assigned_node,
    assigned_node_label = assigned_label,
    assigned_status = as.character(classification$status[[1L]]),
    evidence_role = role,
    site_reason = site_reason,
    stringsAsFactors = FALSE
  )
  out[order(out$alignment_site, out$node), names(empty), drop = FALSE]
}

summarise_partial_opportunistic_site_evidence <- function(site_evidence, classifications = NULL) {
  empty <- empty_partial_opportunistic_site_summary_table()
  if (is.null(classifications) || nrow(classifications) == 0L) {
    return(empty)
  }
  rows <- lapply(seq_len(nrow(classifications)), function(i) {
    cls <- classifications[i, , drop = FALSE]
    ev <- site_evidence[site_evidence$record_id == cls$record_id[[1L]], , drop = FALSE]
    supporting_assigned <- ev[!is.na(ev$supports_assigned_node) & ev$supports_assigned_node, , drop = FALSE]
    supporting_path <- ev[!is.na(ev$supports_assigned_path) & ev$supports_assigned_path, , drop = FALSE]
    off_path <- ev[ev$supports_node & (is.na(ev$supports_assigned_path) | !ev$supports_assigned_path), , drop = FALSE]
    supported_ev <- ev[ev$supports_node, , drop = FALSE]
    node_counts <- if (nrow(supported_ev) == 0L) {
      data.frame(node = character(), count = integer())
    } else {
      counts <- aggregate(alignment_site ~ node, supported_ev, function(x) length(unique(x)))
      names(counts) <- c("node", "count")
      counts[order(-counts$count, counts$node), , drop = FALSE]
    }
    top_nodes <- if (nrow(node_counts) == 0L) "" else {
      paste(utils::head(paste0(node_counts$node, ":", node_counts$count), 5L), collapse = ";")
    }
    interpretation <- if (nrow(ev) == 0L) {
      "no observed scored SNP evidence in mapped interval"
    } else if (cls$status[[1L]] == "resolved_opportunistic") {
      "opportunistic assignment supported by observed region-specific scored SNPs"
    } else if (cls$status[[1L]] == "weak_support") {
      "weak opportunistic signal; too few observed scored SNPs or insufficient support"
    } else {
      cls$reason[[1L]]
    }
    data.frame(
      record_id = cls$record_id[[1L]],
      assigned_node = cls$assigned_node[[1L]],
      assigned_status = cls$status[[1L]],
      unique_scored_sites_in_region = as.integer(cls$sites_available_in_region[[1L]]),
      unique_observed_scored_sites = length(unique(ev$alignment_site)),
      supporting_sites_for_assigned_node = length(unique(supporting_assigned$alignment_site)),
      supporting_sites_for_assigned_path = length(unique(supporting_path$alignment_site)),
      off_path_supporting_sites = length(unique(off_path$alignment_site)),
      top_supported_nodes = top_nodes,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[, names(empty), drop = FALSE]
}

extract_observed_scored_sites <- function(query, scored_rules) {
  x <- encode_classifier_query(query)
  if (nrow(scored_rules) == 0L) {
    return(data.frame(site = integer(), allele_code = integer()))
  }
  sites <- sort(unique(as.integer(scored_rules$site)))
  sites <- sites[sites >= 1L & sites <= length(x)]
  observed <- sites[x[sites] != 0L]
  data.frame(site = observed, allele_code = as.integer(x[observed]))
}

summarise_opportunistic_node_evidence <- function(record, evidence, assigned_node,
                                                  sites_available_in_region,
                                                  sites_observed,
                                                  target_mask = NULL) {
  ev <- evidence[!is.na(evidence$observed_sites) & evidence$observed_sites > 0L, , drop = FALSE]
  if (nrow(ev) == 0L) {
    return(empty_partial_opportunistic_node_evidence_table())
  }
  ev$record_id <- record$sequence_id[[1L]]
  ev$sequence_id <- record$sequence_id[[1L]]
  ev$source_file <- record$source_file[[1L]]
  ev$node <- ev$node_id
  ev$classification_mode <- "opportunistic"
  ev$sites_available_in_region <- as.integer(sites_available_in_region)
  ev$sites_observed <- as.integer(sites_observed)
  ev$sites_supporting_node <- as.integer(ev$support_sites)
  ev$sites_contradicting_node <- as.integer(ev$conflict_sites)
  ev$weighted_supporting_score <- as.numeric(ev$support_weight)
  ev$weighted_contradicting_score <- as.numeric(ev$conflict_weight)
  ev$support_margin <- vapply(seq_len(nrow(ev)), function(i) {
    competitors <- ev$support[ev$node_id != ev$node_id[[i]]]
    if (length(competitors) == 0L || all(is.na(competitors))) {
      return(NA_real_)
    }
    ev$support[[i]] - max(competitors, na.rm = TRUE)
  }, numeric(1L))
  ev$is_on_assigned_path <- if (is.na(assigned_node) || is.null(target_mask)) {
    NA
  } else {
    nodes_nested(assigned_node, ev$node_id, target_mask)
  }
  ev$evidence_summary <- paste0(ev$support_sites, "/", ev$observed_sites, " observed sites support node")
  cols <- c(
    "record_id", "sequence_id", "source_file", "node", "node_id", "classification_mode",
    "sites_available_in_region", "sites_observed", "sites_supporting_node",
    "sites_contradicting_node", "weighted_supporting_score", "weighted_contradicting_score",
    "support", "support_margin", "is_on_assigned_path", "evidence_summary"
  )
  ev[, intersect(cols, names(ev)), drop = FALSE]
}

normalize_partial_classification_mode <- function(classification_mode) {
  classification_mode <- as.character(classification_mode %||% "panel")
  allowed <- c("panel", "opportunistic", "auto")
  if (!classification_mode %in% allowed) {
    stop("classification_mode must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  classification_mode
}

normalize_opportunistic_settings <- function(settings, classifier, min_total_informative_sites = NULL) {
  if (is.null(settings)) {
    settings <- list()
  }
  settings$min_observed_informative_sites <- as.integer(settings$min_observed_informative_sites %||%
    min_total_informative_sites %||%
    classifier$settings$min_total_informative_sites %||%
    1L)
  settings$min_support <- as.numeric(settings$min_support %||% classifier$settings$min_support %||% 0.8)
  settings$min_support_margin <- as.numeric(settings$min_support_margin %||% classifier$settings$support_margin %||% 0.05)
  settings$conflict_support_threshold <- as.numeric(settings$conflict_support_threshold %||% classifier$settings$max_conflict %||% 0.2)
  settings$use_weighted_support <- isTRUE(settings$use_weighted_support %||% TRUE)
  settings
}

make_opportunistic_classifier <- function(classifier, site_node_scores, use_weighted_support = TRUE) {
  rules <- make_classifier_rules(site_node_scores, selected_panel = NULL, use_selected_panel = FALSE)
  if (!isTRUE(use_weighted_support)) {
    rules$rule_weight <- 1
  }
  metadata_cols <- intersect(
    c("node_id", "node_index", "clade_size", "complement_size", "depth", "balance", "weight"),
    names(site_node_scores)
  )
  metadata <- unique(site_node_scores[, metadata_cols, drop = FALSE])
  if (!"node_id" %in% names(metadata) || nrow(metadata) == 0L) {
    metadata <- classifier$node_table[, intersect(c("node_id", "node_index", "clade_size", "complement_size", "depth", "balance", "weight"), names(classifier$node_table)), drop = FALSE]
  }
  metadata$node_id <- as.character(metadata$node_id)
  fallback_metadata <- classifier$node_table[, intersect(c("node_id", "node_index", "clade_size", "complement_size", "depth", "balance", "weight"), names(classifier$node_table)), drop = FALSE]
  fallback_metadata$node_id <- as.character(fallback_metadata$node_id)
  for (name in setdiff(names(fallback_metadata), names(metadata))) {
    metadata[[name]] <- fallback_metadata[[name]][match(metadata$node_id, fallback_metadata$node_id)]
  }
  metadata <- metadata[metadata$node_id %in% as.character(rules$node_id), , drop = FALSE]
  node_table <- make_classifier_node_table(classifier$target_mask, metadata, rules)
  classifier$rules <- rules
  classifier$node_table <- node_table
  classifier$target_mask <- classifier$target_mask[match(node_table$node_id, rownames(classifier$target_mask)), , drop = FALSE]
  classifier$informative_sites <- sort(unique(as.integer(rules$site)))
  classifier$settings$use_selected_panel <- FALSE
  classifier
}

prediction_support_margin <- function(prediction) {
  assigned <- prediction$support_score
  competing <- prediction$competing_support
  if (is.na(assigned) || is.na(competing)) {
    return(NA_real_)
  }
  assigned - competing
}

harmonize_partial_classification_columns <- function(rows, target) {
  missing <- setdiff(names(target), names(rows))
  for (name in missing) {
    rows[[name]] <- NA
  }
  rows[, names(target), drop = FALSE]
}

partial_region_diagnostics_row <- function(record, sites_available, sites_observed, observed_non_missing) {
  data.frame(
    record_id = record$sequence_id,
    sequence_id = record$sequence_id,
    source_file = record$source_file,
    mapped = as.logical(record$mapped),
    alignment_start = as.integer(record$mapped_start),
    alignment_end = as.integer(record$mapped_end),
    sites_available_in_region = as.integer(sites_available),
    sites_observed = as.integer(sites_observed),
    observed_non_missing_sites = as.integer(observed_non_missing),
    stringsAsFactors = FALSE
  )
}

partial_classification_row <- function(record, observed_non_missing, sites_available_in_region,
                                       observed_informative, supporting_sites, contradictory_sites,
                                       assigned_node, support_score, support_margin,
                                       status, note, classification_mode = "panel",
                                       classification_source = "panel") {
  mapping_col <- function(name, default) {
    if (name %in% names(record)) record[[name]] else default
  }
  data.frame(
    record_id = record$sequence_id,
    sequence_id = record$sequence_id,
    source_file = record$source_file,
    raw_length = as.integer(record$length),
    mapped = as.logical(record$mapped),
    mapping_mode = as.character(mapping_col("mapping_mode", mapping_col("mapping_note", NA_character_))),
    strand = as.character(mapping_col("mapping_strand", NA_character_)),
    identity = as.numeric(mapping_col("mapping_identity", NA_real_)),
    query_coverage = as.numeric(mapping_col("mapping_query_coverage", NA_real_)),
    aligned_length = as.integer(mapping_col("mapping_aligned_length", NA_integer_)),
    alignment_start = as.integer(record$mapped_start),
    alignment_end = as.integer(record$mapped_end),
    mapped_start = as.integer(record$mapped_start),
    mapped_end = as.integer(record$mapped_end),
    classification_mode = classification_mode,
    classification_source = classification_source,
    sites_available_in_region = as.integer(sites_available_in_region),
    observed_non_missing_sites = as.integer(observed_non_missing),
    observed_informative_sites = as.integer(observed_informative),
    supporting_sites = as.integer(supporting_sites),
    contradictory_sites = as.integer(contradictory_sites),
    assigned_node = assigned_node,
    assigned_node_label = assigned_node,
    assigned_support = as.numeric(support_score),
    support_margin = as.numeric(support_margin),
    support_score = as.numeric(support_score),
    status = status,
    reason = note,
    note = note,
    stringsAsFactors = FALSE
  )
}

empty_partial_classification_table <- function() {
  data.frame(
    record_id = character(),
    sequence_id = character(),
    source_file = character(),
    raw_length = integer(),
    mapped = logical(),
    mapping_mode = character(),
    strand = character(),
    identity = numeric(),
    query_coverage = numeric(),
    aligned_length = integer(),
    alignment_start = integer(),
    alignment_end = integer(),
    mapped_start = integer(),
    mapped_end = integer(),
    classification_mode = character(),
    classification_source = character(),
    sites_available_in_region = integer(),
    observed_non_missing_sites = integer(),
    observed_informative_sites = integer(),
    supporting_sites = integer(),
    contradictory_sites = integer(),
    assigned_node = character(),
    assigned_node_label = character(),
    assigned_support = numeric(),
    support_margin = numeric(),
    support_score = numeric(),
    status = character(),
    reason = character(),
    note = character(),
    stringsAsFactors = FALSE
  )
}

empty_partial_node_evidence_table <- function() {
  data.frame(
    sequence_id = character(),
    source_file = character(),
    node_id = character(),
    observed_sites = integer(),
    support_sites = integer(),
    conflict_sites = integer(),
    support_weight = numeric(),
    conflict_weight = numeric(),
    support = numeric(),
    conflict = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_partial_opportunistic_node_evidence_table <- function() {
  data.frame(
    record_id = character(),
    sequence_id = character(),
    source_file = character(),
    node = character(),
    node_id = character(),
    classification_mode = character(),
    sites_available_in_region = integer(),
    sites_observed = integer(),
    sites_supporting_node = integer(),
    sites_contradicting_node = integer(),
    weighted_supporting_score = numeric(),
    weighted_contradicting_score = numeric(),
    support = numeric(),
    support_margin = numeric(),
    is_on_assigned_path = logical(),
    evidence_summary = character(),
    stringsAsFactors = FALSE
  )
}

empty_partial_opportunistic_site_evidence_table <- function() {
  data.frame(
    record_id = character(),
    classification_mode = character(),
    classification_source = character(),
    mapped = logical(),
    mapping_mode = character(),
    strand = character(),
    alignment_start = integer(),
    alignment_end = integer(),
    alignment_site = integer(),
    query_base = character(),
    reference_base = character(),
    node = character(),
    node_label = character(),
    node_allele = character(),
    site_score = numeric(),
    site_weight = numeric(),
    supports_node = logical(),
    supports_assigned_node = logical(),
    supports_assigned_path = logical(),
    contradicts_assigned_node = logical(),
    assigned_node = character(),
    assigned_node_label = character(),
    assigned_status = character(),
    evidence_role = character(),
    site_reason = character(),
    stringsAsFactors = FALSE
  )
}

empty_partial_opportunistic_site_summary_table <- function() {
  data.frame(
    record_id = character(),
    assigned_node = character(),
    assigned_status = character(),
    unique_scored_sites_in_region = integer(),
    unique_observed_scored_sites = integer(),
    supporting_sites_for_assigned_node = integer(),
    supporting_sites_for_assigned_path = integer(),
    off_path_supporting_sites = integer(),
    top_supported_nodes = character(),
    interpretation = character(),
    stringsAsFactors = FALSE
  )
}

empty_partial_region_diagnostics_table <- function() {
  data.frame(
    record_id = character(),
    sequence_id = character(),
    source_file = character(),
    mapped = logical(),
    alignment_start = integer(),
    alignment_end = integer(),
    sites_available_in_region = integer(),
    sites_observed = integer(),
    observed_non_missing_sites = integer(),
    stringsAsFactors = FALSE
  )
}

partial_status_note <- function(status) {
  switch(
    status,
    resolved = "classification resolved by conservative tree-path classifier",
    resolved_opportunistic = "classification resolved with region-dependent opportunistic SNP evidence",
    weak_support = "best available tree-path evidence is below support threshold",
    no_informative_sites = "mapped sequence has too few observed classifier-informative sites",
    no_observed_informative_sites = "mapped sequence overlaps informative sites but has no observed unambiguous bases at those sites",
    conflicting = "observed informative sites support competing or conflicting nodes",
    unresolved = "mapped sequence has informative signal but no conservative node assignment",
    unmapped = "sequence could not be mapped to alignment coordinates",
    "classification completed"
  )
}
