#' Convert EQR object to array
#'
#' @encoding UTF-8
#' @title Convert EQR object to array
#'
#' @description
#' Converts an \code{"fcs2EQR"} object to a standard \R \code{\link{array}} of
#' Monte Carlo \acronym{EQR} samples.
#'
#' @name as.array
#' @usage as.array(x, ...)
#' @param x an object of class \code{"fcs2EQR"} containing Monte Carlo
#'   \acronym{EQR} samples, as calculated from \code{\link{fcs2SingleEQR}},
#'   \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}.
#' @param \dots Not currently used.
#' @return Returns the Monte Carlo \acronym{EQR} samples as an
#'   \code{\link{array}} object.
#' @seealso \code{\link{fcs2SingleEQR}}, \code{\link{fcs2JointEQR}} or
#'   \code{\link{fcs2JointAndSingleEQR}} for producing \code{"fcs2EQR"}
#'   objects.\cr \code{[.fcs2EQR} for extracting a subset of
#'   \acronym{EQR} values.
#' @export
as.array <- function(x, ...) UseMethod("as.array")

#' @export
as.array.fcs2EQR <- function(x, ...)
{
    class(x) <- "array"
    x
}
