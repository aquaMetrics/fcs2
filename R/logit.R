#' The Logit Function
#'
#' The \code{logit} function and its inverse \code{expit}.  This is used as the
#' link function for the prevalence regression in the \acronym{FCS2} model to
#' transform the interval \eqn{(0, 1)} to the real line.
#'
#' The logit transformation is defined by \eqn{logit(x) = log(x) - log(1 - x)}
#' so that its inverse expit is given by \eqn{expit(x) = 1 / (1 + exp(-x))}.
#'
#' @name logit
#' @aliases logit expit
#' @usage logit(x)
#' @usage expit(x)
#' @param x a numeric vector or matrix
#' @return A vector or matrix the same size as \code{x} containing the
#' transformed values.
#' @seealso \code{\link{log}} for the natural logarithm and
#' \code{\link{prevalence}} for producing samples of the prevalence model
#' component, using the \code{logit} transformation.
#' @export
logit <-
function(x)
    log(x) - log(1 - x)


## expit
#' @rdname logit
#' @export
expit <-
function(x) {
    1 / (1 + exp(-x))
}
