#' Fisheries Classification Scheme 2 For SNIFFER
#'
#' Provides functions that carry out
#' \acronym{SNIFFER}'s implementation of the Environment Agency's
#' \dfn{Fisheries Classification Scheme 2} (\acronym{FCS2}).  This package was
#' developed for use in Scotland, Northern Ireland and the Republic of Ireland
#' as part of \acronym{SNIFFER} project WFD68c: Science Work.
#'
#' The main functions are \code{\link{fcs2FitModel}} which fits the
#' \acronym{FCS2} statistical model, \code{\link{fcs2JointEQR}} which
#' calculates samples of the joint \acronym{EQR} and \code{\link{fcs2Classify}}
#' which uses these to produce probabilistic \acronym{WFD} classifications.
#'
#' The \acronym{FCS2} approach consists primarily of two main tasks:
#' \enumerate{ \item Fit the statistical model for each species.  \item Produce
#' \acronym{EQR}s and classifications.  } This package provides functions for
#' carrying out all stages of this analysis.  Other packages may be of use, for
#' example \pkg{RODBC} can be used to read data into from a database.
#'
#' Fitting the model for a single species consists of the following steps:
#' \enumerate{ \item Select covariate terms for prevalence and abundance
#' regressions.\cr For this, \code{\link{fcs2FitModel}} can be used with
#' \code{runBUGS = FALSE} to find approximate parameter estimates using
#' \acronym{INLA} with which to judge the significance of each suggested term.
#' Alternatively \code{\link{fcs2ModelSelection}} can be used to automatically
#' select a set of significant regression terms.  \item Set prior
#' distributions.\cr Default priors are given but \code{\link{plot.fcs2Fit}}
#' and \code{\link{fcs2Priors}} can be used to check the priors so that they
#' can be modified if necessary before fitting the full Bayesian model.  \item
#' Fit the full model.\cr \code{\link{fcs2FitModel}} used with \code{runBUGS =
#' TRUE} will fit the model using \acronym{MCMC} via either WinBUGS or
#' OpenBUGS.  This can take much time.  \item Check Monte Carlo samples.\cr The
#' convergence of the \acronym{MCMC} chains can be checked with
#' \code{\link{plotBUGSTrace}} and the sample can be thinned if necessary with
#' \code{\link{thinBUGSSamples}}.  }
#'
#' After fitting the model for every species, the \acronym{EQR}s and
#' classifications can be found by the following steps: \enumerate{ \item
#' Select values for pressure variables at reference conditions.\cr Observed
#' values of pressure variables are replaced with reference values for the
#' model to make predictions of reference conditions.  \item Calculate
#' \acronym{EQR} samples.\cr \code{\link{fcs2SingleEQR}} can calculate
#' \acronym{EQR}s for single species and surveys and \code{\link{fcs2JointEQR}}
#' calculates combined \acronym{EQR} values.  \item Select \acronym{EQR} class
#' boundaries.\cr These may be selected by comparing the mean \acronym{EQR}
#' values produced from surveys of each class.  \item Calculate probability of
#' each class.\cr This can be calculated using \code{\link{fcs2Classify}}.\cr
#' Alternatively, the function \code{\link{fcs2EQRSummaryMatrix}} can be used
#' to provide a tabulated summary of joint and single \acronym{EQR}s as well as
#' probabilities of class.  }
#'
#' @name fcs2-package
#' @aliases fcs2-package fcs2
#' @docType package
#' @author David Wyncoll \email{d.wyncoll@@hrwallingford.co.uk} of HR
#' Wallingford \url{http://www.hrwallingford.co.uk} for \acronym{SNIFFER}
#' \url{http://www.sniffer.org.uk}
#' @keywords package
NULL
