test_that("test results match pre-calculated results", {
  results <- calcClassScot(data = fcs2::demo_data)

  expect_equal(results, demo_results)
})