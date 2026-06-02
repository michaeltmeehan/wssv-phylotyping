partial_default_extensions <- c("fa", "fas", "fasta", "fna")

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

map_partial_sequences <- function(partials, alignment_length, reference_sequence = NULL,
                                  mapping_mode = "auto") {
  if (nrow(partials) == 0L) {
    partials$mapped <- logical()
    partials$mapped_start <- integer()
    partials$mapped_end <- integer()
    partials$mapped_sequence <- character()
    partials$mapping_note <- character()
    return(partials)
  }
  allowed <- c("aligned_full_length", "exact_substring", "auto")
  if (!mapping_mode %in% allowed) {
    stop("mapping_mode must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  mapped <- lapply(seq_len(nrow(partials)), function(i) {
    map_partial_sequence(partials[i, , drop = FALSE], alignment_length, reference_sequence, mapping_mode)
  })
  out <- cbind(partials, do.call(rbind, mapped))
  rownames(out) <- NULL
  out
}

map_partial_sequence <- function(record, alignment_length, reference_sequence = NULL,
                                 mapping_mode = "auto") {
  seq <- toupper(gsub("\\s+", "", record$sequence[[1L]]))
  if (mapping_mode %in% c("aligned_full_length", "auto") && nchar(seq, type = "bytes") == alignment_length) {
    encoded <- encode_classifier_query(seq)
    observed <- which(encoded != 0L)
    return(data.frame(
      mapped = TRUE,
      mapped_start = if (length(observed) > 0L) min(observed) else NA_integer_,
      mapped_end = if (length(observed) > 0L) max(observed) else NA_integer_,
      mapped_sequence = seq,
      mapping_note = "accepted aligned full-length sequence",
      stringsAsFactors = FALSE
    ))
  }

  if (mapping_mode %in% c("exact_substring", "auto") && !is.null(reference_sequence)) {
    placed <- map_exact_substring(seq, reference_sequence, alignment_length)
    if (!is.null(placed)) {
      return(placed)
    }
  }

  data.frame(
    mapped = FALSE,
    mapped_start = NA_integer_,
    mapped_end = NA_integer_,
    mapped_sequence = NA_character_,
    mapping_note = "unmapped: unsupported length or no exact reference substring match",
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
  ref_cols <- which(ref_bases)
  aln_start <- ref_cols[[ungapped_start]]
  aln_end <- ref_cols[[ungapped_end]]

  full <- rep("-", alignment_length)
  full[ref_cols[ungapped_start:ungapped_end]] <- strsplit(query, "", fixed = TRUE)[[1L]]
  data.frame(
    mapped = TRUE,
    mapped_start = aln_start,
    mapped_end = aln_end,
    mapped_sequence = paste0(full, collapse = ""),
    mapping_note = "mapped by exact substring against reference sequence",
    stringsAsFactors = FALSE
  )
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
      rows[[i]] <- partial_classification_row(record, observed_non_missing, 0L, NA_character_, NA_real_, "unmapped", record$mapping_note)
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
      pred$observed_informative_sites,
      pred$assigned_node,
      pred$support_score,
      pred$status,
      note
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

partial_classification_row <- function(record, observed_non_missing, observed_informative,
                                       assigned_node, support_score, status, note) {
  data.frame(
    sequence_id = record$sequence_id,
    source_file = record$source_file,
    raw_length = as.integer(record$length),
    mapped = as.logical(record$mapped),
    mapped_start = as.integer(record$mapped_start),
    mapped_end = as.integer(record$mapped_end),
    observed_non_missing_sites = as.integer(observed_non_missing),
    observed_informative_sites = as.integer(observed_informative),
    assigned_node = assigned_node,
    support_score = as.numeric(support_score),
    status = status,
    note = note,
    stringsAsFactors = FALSE
  )
}

empty_partial_classification_table <- function() {
  data.frame(
    sequence_id = character(),
    source_file = character(),
    raw_length = integer(),
    mapped = logical(),
    mapped_start = integer(),
    mapped_end = integer(),
    observed_non_missing_sites = integer(),
    observed_informative_sites = integer(),
    assigned_node = character(),
    support_score = numeric(),
    status = character(),
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

partial_status_note <- function(status) {
  switch(
    status,
    resolved = "classification resolved by conservative tree-path classifier",
    weak_support = "best available tree-path evidence is below support threshold",
    no_informative_sites = "mapped sequence has too few observed classifier-informative sites",
    conflicting = "observed informative sites support competing or conflicting nodes",
    unresolved = "mapped sequence has informative signal but no conservative node assignment",
    unmapped = "sequence could not be mapped to alignment coordinates",
    "classification completed"
  )
}
