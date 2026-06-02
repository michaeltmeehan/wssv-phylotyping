source(test_path("../../R/script_utils.R"))

test_that("read_table_output prefers RDS over CSV", {
  dir <- tempfile()
  dir.create(dir)
  saveRDS(data.frame(x = 1L), file.path(dir, "example.rds"))
  write.csv(data.frame(x = 2L), file.path(dir, "example.csv"), row.names = FALSE)

  out <- read_table_output(dir, "example")

  expect_equal(out$x, 1L)
})

test_that("read_table_output falls back to CSV", {
  dir <- tempfile()
  dir.create(dir)
  write.csv(data.frame(x = 2L), file.path(dir, "example.csv"), row.names = FALSE)

  out <- read_table_output(dir, "example")

  expect_equal(out$x, 2L)
})

test_that("read_table_output reports missing upstream output clearly", {
  dir <- tempfile()
  dir.create(dir)

  expect_error(
    read_table_output(dir, "example", "scripts/04_score_windows.R"),
    "Missing .*example.csv.*example.rds.*Run scripts/04_score_windows.R first"
  )
})

test_that("value_or only replaces NULL values", {
  expect_equal(value_or(NULL, "fallback"), "fallback")
  expect_equal(value_or(FALSE, TRUE), FALSE)
  expect_true(is.na(value_or(NA, "fallback")))
})

test_that("format_count_list names empty vectors clearly", {
  expect_equal(format_count_list(character()), "none")
  expect_equal(format_count_list(c(1L, 3L)), "1, 3")
})
