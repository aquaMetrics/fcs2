test_that("test results match pre-calculated results", {
  results <- calcClassScot(data = fcs2::demo_data, seed = 42)
  expect_equal(results, fcs2::demo_results)
})
