#' Quantiles of an EQR Sample
#'
#' Estimates quantiles of an \acronym{EQR} variable from a collection of Monte
#' Carlo \acronym{EQR} samples.
#'
#'
#' @param x an object of class \code{"fcs2EQR"} containing Monte Carlo
#' \acronym{EQR} samples, as calculated from \code{\link{fcs2SingleEQR}},
#' \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}.
#' @param probs numeric vector of probabilities with values in [0, 1].
#' @param \dots further arguments passed to \code{\link{quantile}}.
#' @return Returns a matrix or array of quantile estimates for each of the
#' surveys (or sites/etc if surveys were joined) and species contained in
#' \code{x}.
#' @seealso \code{\link{mean.fcs2EQR}}
#' @export
quantile.fcs2EQR <-
function(x, probs = seq(0, 1, 0.25), ...)
{
    if (length(dim(x)) == 2)
        apply(x, 2, quantile, probs=probs, ...)
    else
        apply(x, 2:3, quantile, probs=probs, ...)
}

