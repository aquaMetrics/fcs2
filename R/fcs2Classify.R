#' Probabilistic WFD Classifications from EQR Samples
#'
#' Produces probabilistic Water Framework Directive classifications of
#' \emph{Bad}, \emph{Poor}, \emph{Good}, \emph{Moderate} or \emph{High} from
#' Monte Carlo \acronym{EQR} samples.
#'
#' @param eqr an object of class \code{"fcs2EQR"} containing Monte Carlo
#'   \acronym{EQR} samples, as calculated from \code{\link{fcs2SingleEQR}},
#'   \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}.
#' @param survey index specifying which surveys (or water bodies etc if surveys
#'   joined) to classify.
#' @param species index specifying which species to classify.
#' @param boundaries a vector of length 4 giving the \acronym{EQR} boundaries
#'   separating the classes \emph{Bad}, \emph{Poor}, \emph{Good},
#'   \emph{Moderate} and \emph{High}. If missing, regularly spaced boundaries of
#'   \code{c(0.2, 0.4, 0.6, 0.8)} are used with a warning.
#' @return A matrix or array containing the probabilities of each \acronym{WFD}
#'   class, as estimated from the proportion of \acronym{EQR} samples between
#'   the class boundaries. The classes \emph{Bad}, \emph{Poor}, \emph{Good},
#'   \emph{Moderate} and \emph{High} are given in rows 1 to 5 respectively,
#'   while the columns indicate the selected surveys (or sites/water bodies/etc
#'   if surveys were joined when calculating \acronym{EQR}s). If \code{eqr}
#'   contains \acronym{EQR} samples for multiple fits/species, an array is
#'   returned with species as the third dimension.
#' @seealso \code{\link{plot.fcs2EQR}} with \code{type="class"} for plotting the
#'   probabilistic classifications.#' \code{\link{fcs2SingleEQR}},
#'   \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}} for
#'   producing \code{"fcs2EQR"} objects.
#' @examples
#'     \dontrun{
#' ### Very simple example using a single species EQR
#'
#' # simulate random dataset
#' Data <- data.frame(SurveyArea=rlnorm(100, 4.6, 0.5))
#'
#' # random survey area
#' Data$Catch <- rzinbinom(100, size=1.1, zeroprob=0.3,
#'                   nbmean=0.3 * Data$SurveyArea)  # single catch per survey
#'
#' # define a simple model with no covariates
#' # and fit full model with OpenBUGS
#' fit <- fcs2FitModel("Catch", dataFrame=Data, surveyAreaVar="SurveyArea",
#'                            runBUGS=TRUE, bugsProgram="OpenBUGS", n.iter=1000)
#'
#' # calculate samples of single EQR, using same dataset
#' eqr <- fcs2SingleEQR(fit, Data)
#'
#' # define WFD class boundaries
#' boundaries <- c(0.001, 0.01, 0.25, 0.5)
#'
#' # plot EQR variables for first 9 surveys,
#' # shading the probabilities of each class
#' plot(eqr, 1:9, type="density", boundaries=boundaries)
#'
#' # calculate probabilistic classifications
#' fcs2Classify(eqr, boundaries=boundaries)
#'
#' # plot these probabilities for the first 9 surveys
#' plot(eqr, 1:9, type="class", boundaries=boundaries)
#' }
#' @export


fcs2Classify <-
function(eqr, survey = 1:ncol(eqr), species, boundaries)
{
    # if boundaries missing, use regularly spaced boundaries
    if (missing(boundaries)) {
        boundaries <- seq(0.2, 0.8, by=0.2)
        warning("using regularly spaced class boundaries as 'boundaries' argument not supplied")
    }
    boundaries <- c(-0.001, boundaries, 1)

    # check survey not logical
    if (inherits(survey, "logical"))
        survey <- which(survey)

    # calculate probability of each class
    if (length(dim(eqr)) == 2) {
        # eqr by survey

        # calculate probs
        prob <- array(dim=c(length(boundaries) - 1, length(survey)),
                      dimnames=list(c("Bad", "Poor", "Moderate", "Good", "High"), colnames(eqr[, survey, drop=FALSE])))
        for (i in 1:(length(boundaries) - 1)) {
             for (j in 1:length(survey))
                prob[i, j] <- mean(eqr[, survey[j]] > boundaries[i] & eqr[, survey[j]] <= boundaries[i + 1])
        }

    } else {
        # eqr by survey and fit

        # set species if necessary
        if (missing(species))
            species <- 1:(dim(eqr)[3])

        # check species not logical
        if (inherits(species, "logical"))
            species <- which(species)

        # calculate probs
        prob <- array(dim=c(length(boundaries) - 1, length(survey), length(species)),
                      dimnames=c(list(c("Bad", "Poor", "Moderate", "Good", "High")),
                                 dimnames(eqr[, survey, species, drop=FALSE])[2:3]))
        for (i in 1:(length(boundaries) - 1)) {
            for (j in 1:length(survey))
                for (k in 1:length(species))
                    prob[i, j, k] <- mean(eqr[, survey[j], species[k]] > boundaries[i] &
                                          eqr[, survey[j], species[k]] <= boundaries[i + 1])
        }
    }

    prob
}

