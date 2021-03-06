% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.fcs2Fit.R
\name{predict.fcs2Fit}
\alias{predict.fcs2Fit}
\title{Predicted Total Catch}
\usage{
\method{predict}{fcs2Fit}(object, newData, subset = 1:nrow(newData), na.action, mu, rho, ...)
}
\arguments{
\item{object}{an \code{"fcs2Fit"} object containing a full \acronym{FCS2}
model object, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS
= TRUE}.}

\item{newData}{a data frame with surveys as rows and variables as columns.
It should contain all variables required by \code{object}.}

\item{subset}{an optional vector specifying a subset of surveys to calculate
the predicted total catch for.}

\item{na.action}{a function which indicates what should happen when the data
contain missing values (\code{NA}s).  The default is set by the
\code{na.action} setting of \code{\link{options}} and this is usually set to
\code{\link{na.omit}}.  This setting removes surveys that contain missing
data in any required variables.  A vector indicating the rows that were
removed can be extracted from the returned object using
\code{\link{na.action}}.  Alternatively, \code{\link{na.pass}} can be used
to ignore missing values (where possible) or \code{\link{na.fail}} can be
given to signal an error if missing values are found.}

\item{mu}{a matrix of posterior samples of the abundance component \eqn{\mu}
can optionally be given to save recalculation if already available.  This is
assumed to have been calculated from \code{\link{abundance}} using the same
arguments as above.}

\item{rho}{a matrix of posterior samples of the prevalence component
\eqn{\rho} can optionally be given to save recalculation if already
available.  This is assumed to have been calculated from
\code{\link{prevalence}} using the same arguments as above.}

\item{\dots}{Not currently used.}
}
\value{
a vector containing the expected total catch for each survey
selected.
}
\description{
Calculates the expected total catch over all runs for a range of surveys, as
predicted by a fitted \acronym{FCS2} model.
}
\seealso{
\code{\link{fcs2FitModel}} for fitting the \acronym{FCS2} model.\cr
\code{\link{pCatch}} for calculating probabilities relating to the predicted
total catch.
}
