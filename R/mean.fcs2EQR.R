#' Mean EQR
#'
#' @encoding UTF-8
#' @title Mean EQR
#'
#' @description
#' Estimates the expected \acronym{EQR} from a collection of Monte Carlo
#' \acronym{EQR} samples.
#'
#' @name mean
#' @usage mean(x, ...)
#' @param x an object of class \code{"fcs2EQR"} containing Monte Carlo
#' \acronym{EQR} samples, as calculated from \code{\link{fcs2SingleEQR}},
#' \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}.
#' @param \dots Not currently used.
#' @return Returns a vector or matrix containing the expected \acronym{EQR} for
#' each of the surveys (or sites/etc if surveys were joined) and species
#' contained in \code{x}.
#' @seealso \code{\link{quantile.fcs2EQR}}
#' @export
mean <- function(x, ...) UseMethod("mean")

#' @rdname mean
#' @export
mean.fcs2EQR <-
function(x, ...)
{
    if (length(dim(x)) == 2)
        colMeans(x)
    else
        apply(x, 2:3, mean)
}

