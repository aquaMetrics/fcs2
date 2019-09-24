#' Missing Values in EQR Samples
#'
#' Extract information about surveys (or water bodies, sites, etc) that were
#' removed when calculating the \acronym{EQR} due to missing values.
#'
#'
#' @param object an object of class \code{"fcs2EQR"} containing Monte Carlo
#' \acronym{EQR} samples, as calculated from \code{\link{fcs2SingleEQR}},
#' \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}.
#' @return If the \acronym{EQR} object contains \acronym{EQR} samples for
#' multiple surveys, the returned vector indicates surveys in the data frame
#' that were removed due to missing observations or covariates.
#'
#' Alternatively, if the \acronym{EQR} object contains \acronym{EQR} samples
#' corresponding to another variable (the \code{joinByVar} in the call to
#' \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}, e.g.
#' water body ID), the returned vector indicates values of the variable that
#' were removed from the \acronym{EQR} object due to missing values.
#' @seealso \code{\link{fcs2SingleEQR}}, \code{\link{fcs2JointEQR}} and
#' \code{\link{fcs2JointAndSingleEQR}} for generating \acronym{EQR} samples.
#' @export
na.action.fcs2EQR <-
function(object)
    attr(object, "na.action")

