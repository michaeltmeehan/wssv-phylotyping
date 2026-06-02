base_code_map <- c(A = 1L, C = 2L, G = 3L, T = 4L)

base_code_to_allele <- function(code) {
  alleles <- c("A", "C", "G", "T")
  out <- rep(NA_character_, length(code))
  ok <- !is.na(code) & code %in% seq_along(alleles)
  out[ok] <- alleles[code[ok]]
  out
}

alignment_to_character <- function(alignment) {
  if (inherits(alignment, "DNAStringSet")) {
    seqs <- as.character(alignment)
    names(seqs) <- names(alignment)
    return(seqs)
  }

  if (is.character(alignment)) {
    return(alignment)
  }

  stop("Alignment must be a DNAStringSet or named character vector.", call. = FALSE)
}

validate_alignment_lengths <- function(alignment) {
  seqs <- alignment_to_character(alignment)
  if (length(seqs) < 2L) {
    stop("Alignment contains fewer than 2 sequences.", call. = FALSE)
  }
  if (is.null(names(seqs)) || anyNA(names(seqs)) || any(names(seqs) == "")) {
    stop("Alignment sequences must be named.", call. = FALSE)
  }
  if (anyDuplicated(names(seqs))) {
    stop("Duplicate sequence names in alignment.", call. = FALSE)
  }

  lengths <- nchar(seqs, type = "bytes")
  if (length(unique(lengths)) != 1L) {
    stop("Sequences have different lengths; expected a fixed-width alignment.", call. = FALSE)
  }

  invisible(lengths[[1L]])
}

encode_alignment_int <- function(alignment) {
  seqs <- alignment_to_character(alignment)
  aln_len <- validate_alignment_lengths(seqs)

  lut <- integer(256L)
  lut[as.integer(charToRaw("A")) + 1L] <- 1L
  lut[as.integer(charToRaw("C")) + 1L] <- 2L
  lut[as.integer(charToRaw("G")) + 1L] <- 3L
  lut[as.integer(charToRaw("T")) + 1L] <- 4L
  lut[as.integer(charToRaw("a")) + 1L] <- 1L
  lut[as.integer(charToRaw("c")) + 1L] <- 2L
  lut[as.integer(charToRaw("g")) + 1L] <- 3L
  lut[as.integer(charToRaw("t")) + 1L] <- 4L

  encoded <- matrix(0L, nrow = length(seqs), ncol = aln_len)
  for (i in seq_along(seqs)) {
    raw_seq <- charToRaw(seqs[[i]])
    encoded[i, ] <- lut[as.integer(raw_seq) + 1L]
  }

  rownames(encoded) <- names(seqs)
  colnames(encoded) <- as.character(seq_len(aln_len))
  encoded
}

is_polymorphic_site <- function(x) {
  observed <- x[x != 0L]
  length(unique(observed)) > 1L
}

find_polymorphic_sites <- function(aln_int) {
  if (!is.matrix(aln_int)) {
    stop("aln_int must be an integer matrix.", call. = FALSE)
  }
  which(apply(aln_int, 2L, is_polymorphic_site))
}
