#' Missing Values in FCS2 Model Fit
#'
#' @encoding UTF-8
#' @title Missing Values in FCS2 Model Fit
#'
#' @description
#' Extract information about surveys that were removed when fitting the
#' \acronym{FCS2} statistical model due to missing values.
#'
#'
#' @param object an object of class \code{"fcs2Fit"}, usually returned by
#' \code{\link{fcs2FitModel}}.
#' @param \dots Not currently used.
#' @return the returned vector indicates surveys in the data frame that were
#' removed due to missing observations or covariates.
#' @seealso \code{\link{fcs2FitModel}} for fitting the \acronym{FCS2} model.
#' @export
na.action.fcs2Fit <-
  function(object, ...)
    object$na.action
