test_that("test works", {
# simulate random dataset
data <- data.frame(SurveyArea = rlnorm(100, 4.6, 0.5)) # random survey area
data$Catch <- rzinbinom(100,
  size = 1.1,
  zeroprob = 0.3,
  nbmean = 0.3 * data$SurveyArea
) # single catch per survey

# fit approximate model with INLA - model contains no regression terms
fit1 <- fcs2FitModel("Catch", dataFrame = data, surveyAreaVar = "SurveyArea")

# fit full model with OpenBUGS
# (more iterations may be required to achieve convergence)
fit1 <- fcs2FitModel(
  fit = fit1, runBUGS = TRUE, n.iter = 1000,
  bugsProgram = "OpenBUGS"
)

# summarise fit
summary(fit1)
})
