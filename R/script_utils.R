#' Return a default value when a config/list value is NULL.
value_or <- function(x, default) {
  if (is.null(x)) default else x
}

#' Read a generated output table by stem, preferring RDS over CSV.
#'
#' Used by scripts so downstream stages report clear "run upstream script first"
#' errors when expected generated outputs are missing. Optional fallback stems
#' preserve compatibility with renamed outputs while later stages migrate to the
#' new filenames.
read_table_output <- function(table_dir, stem, upstream_script = NULL, fallback_stems = NULL) {
  stems <- unique(c(stem, fallback_stems))
  stems <- stems[!is.na(stems) & nzchar(stems)]
  for (current_stem in stems) {
    rds_path <- file.path(table_dir, paste0(current_stem, ".rds"))
    csv_path <- file.path(table_dir, paste0(current_stem, ".csv"))
    if (file.exists(rds_path)) {
      return(readRDS(rds_path))
    }
    if (file.exists(csv_path)) {
      return(read.csv(csv_path, stringsAsFactors = FALSE))
    }
  }

  hint <- if (is.null(upstream_script)) "" else paste0(" Run ", upstream_script, " first.")
  primary_csv <- file.path(table_dir, paste0(stem, ".csv"))
  primary_rds <- file.path(table_dir, paste0(stem, ".rds"))
  fallback_note <- if (length(stems) > 1L) {
    paste0(" Accepted fallback stems: ", paste(stems[-1L], collapse = ", "), ".")
  } else {
    ""
  }
  stop("Missing ", primary_csv, " or ", primary_rds, ".", fallback_note, hint, call. = FALSE)
}

format_count_list <- function(x) {
  if (length(x) == 0L) "none" else paste(x, collapse = ", ")
}
