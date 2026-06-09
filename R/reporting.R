read_optional_output <- function(table_dir, stem) {
  rds_path <- file.path(table_dir, paste0(stem, ".rds"))
  csv_path <- file.path(table_dir, paste0(stem, ".csv"))
  if (file.exists(rds_path)) {
    return(list(data = readRDS(rds_path), path = rds_path, missing = FALSE))
  }
  if (file.exists(csv_path)) {
    return(list(data = read.csv(csv_path, stringsAsFactors = FALSE), path = csv_path, missing = FALSE))
  }
  list(data = NULL, path = csv_path, missing = TRUE)
}

table_or_note <- function(x, columns = NULL, n = 10L, missing_note = "Missing output.") {
  if (is.null(x)) {
    return(missing_note)
  }
  if (nrow(x) == 0L) {
    return("No rows were available.")
  }
  if (!is.null(columns)) {
    columns <- intersect(columns, names(x))
    x <- x[, columns, drop = FALSE]
  }
  x <- utils::head(x, n)
  markdown_table(x)
}

markdown_table <- function(x) {
  if (is.null(x) || nrow(x) == 0L || ncol(x) == 0L) {
    return("No rows were available.")
  }
  format_cell <- function(value) {
    if (is.na(value)) {
      return("")
    }
    if (is.numeric(value)) {
      return(format(round(value, 4L), trim = TRUE, scientific = FALSE))
    }
    value <- trimws(gsub("\\|", "/", as.character(value)))
    if (nchar(value) > 80L) {
      value <- paste0(substr(value, 1L, 77L), "...")
    }
    value
  }
  rows <- apply(x, 1L, function(row) paste(vapply(row, format_cell, character(1L)), collapse = " | "))
  paste(
    paste(names(x), collapse = " | "),
    paste(rep("---", ncol(x)), collapse = " | "),
    paste(rows, collapse = "\n"),
    sep = "\n"
  )
}

count_values <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(data.frame(value = character(), count = integer()))
  }
  out <- as.data.frame(table(x, useNA = "ifany"), stringsAsFactors = FALSE)
  names(out) <- c("value", "count")
  out[order(-out$count, out$value), , drop = FALSE]
}

numeric_summary <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(data.frame(n = 0L, min = NA_real_, median = NA_real_, mean = NA_real_, max = NA_real_))
  }
  data.frame(
    n = length(x),
    min = min(x),
    median = stats::median(x),
    mean = mean(x),
    max = max(x)
  )
}

interpret_signal_concentration <- function(site_summary, top_n = 10L) {
  if (is.null(site_summary) || nrow(site_summary) == 0L || !"weighted_gain_sum" %in% names(site_summary)) {
    return("SNP concentration could not be assessed because site summary output was unavailable.")
  }
  total <- sum(site_summary$weighted_gain_sum, na.rm = TRUE)
  if (total <= 0) {
    return("Weighted SNP signal is zero or unavailable, so concentration could not be assessed.")
  }
  ordered <- site_summary[order(-site_summary$weighted_gain_sum), , drop = FALSE]
  share <- sum(utils::head(ordered$weighted_gain_sum, top_n), na.rm = TRUE) / total
  if (share >= 0.5) {
    paste0("Signal is relatively concentrated: the top ", min(top_n, nrow(ordered)), " sites account for about ", round(100 * share, 1), "% of total weighted gain.")
  } else {
    paste0("Signal is spread across many sites: the top ", min(top_n, nrow(ordered)), " sites account for about ", round(100 * share, 1), "% of total weighted gain.")
  }
}

panel_size_subset <- function(panel_summary, sizes = c(1L, 2L, 3L, 5L, 10L)) {
  if (is.null(panel_summary) || nrow(panel_summary) == 0L || !"panel_size" %in% names(panel_summary)) {
    return(data.frame())
  }
  panel_summary[panel_summary$panel_size %in% sizes, , drop = FALSE]
}

opportunistic_assigned_site_summary <- function(site_summary) {
  if (is.null(site_summary) || nrow(site_summary) == 0L) {
    return(site_summary)
  }
  if (!"assigned_status" %in% names(site_summary)) {
    return(site_summary)
  }
  assigned <- site_summary[site_summary$assigned_status == "resolved_opportunistic", , drop = FALSE]
  if (nrow(assigned) > 0L) {
    return(assigned)
  }
  site_summary
}

write_report_extract <- function(x, path, columns = NULL, n = NULL) {
  if (is.null(x)) {
    return(FALSE)
  }
  if (!is.null(columns)) {
    columns <- intersect(columns, names(x))
    x <- x[, columns, drop = FALSE]
  }
  if (!is.null(n)) {
    x <- utils::head(x, n)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
  TRUE
}
