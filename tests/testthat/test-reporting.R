source(test_path("../../R/reporting.R"))

test_that("numeric summaries handle values and missing-only vectors", {
  summary <- numeric_summary(c(1, 2, 3, NA))

  expect_equal(summary$n, 3L)
  expect_equal(summary$median, 2)
  expect_true(is.na(numeric_summary(c(NA_real_))$mean))
})

test_that("signal concentration interpretation distinguishes concentrated signal", {
  concentrated <- data.frame(weighted_normalized_gain_sum = c(9, 1, 0))
  spread <- data.frame(weighted_normalized_gain_sum = rep(1, 20))

  expect_match(interpret_signal_concentration(concentrated, top_n = 1L), "concentrated")
  expect_match(interpret_signal_concentration(spread, top_n = 2L), "spread")
})

test_that("markdown tables escape pipe characters", {
  out <- markdown_table(data.frame(name = "a|b", value = 1))

  expect_match(out, "a/b", fixed = TRUE)
})
