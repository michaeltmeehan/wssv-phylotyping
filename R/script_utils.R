#' Return a default value when a config/list value is NULL.
value_or <- function(x, default) {
  if (is.null(x)) default else x
}

#' Read a generated output table by stem, preferring RDS over CSV.
#'
#' Used by scripts so downstream stages report clear "run upstream script first"
#' errors when expected generated outputs are missing.
read_table_output <- function(table_dir, stem, upstream_script = NULL) {
  rds_path <- file.path(table_dir, paste0(stem, ".rds"))
  csv_path <- file.path(table_dir, paste0(stem, ".csv"))
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  if (file.exists(csv_path)) {
    return(read.csv(csv_path, stringsAsFactors = FALSE))
  }

  hint <- if (is.null(upstream_script)) "" else paste0(" Run ", upstream_script, " first.")
  stop("Missing ", csv_path, " or ", rds_path, ".", hint, call. = FALSE)
}

format_count_list <- function(x) {
  if (length(x) == 0L) "none" else paste(x, collapse = ", ")
}
