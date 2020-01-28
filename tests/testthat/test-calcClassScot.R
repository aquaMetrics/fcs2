test_that("test results match pre-calculated results", {
  # this test was set in R < 3.6.0 which changed the default RNG kind method
  # for sampling: https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17494
  if (getRversion() >= '3.6.0') RNGkind(sample.kind = "Rounding")
  results <- calcClassScot(data = fcs2::demo_data, seed = 42)
  expect_equal(results, fcs2::demo_results)
})
