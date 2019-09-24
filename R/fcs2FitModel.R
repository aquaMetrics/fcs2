#' Fit the FCS2 Statistical Model
#'
#' Fits the both the original \acronym{EA} \acronym{FCS2} statistical model and
#' also \acronym{SNIFFER}'s multiple-pass extension that allows multiple runs
#' per survey.  \cr The model allows observations to be given either as fish
#' counts in each run of the survey, the total count over all runs or as an
#' upper and lower bound for the all runs total.  The fish counts depend upon
#' abundance \eqn{\mu} and prevalence \eqn{\rho} components that are each
#' modelled by a covariate regression of linear, non-linear and spatial
#' terms.\cr The Bayesian model is fitted by sampling a set of potential
#' parameter values from the posterior distribution using \acronym{MCMC} via
#' WinBUGS or OpenBUGS.  However, this can take much time so an approximate
#' model fit can be quickly obtained using \acronym{INLA}.
#'
#' The full model (described below) is fitted by drawing Monte Carlo samples
#' from the posterior distribution of the unknown paramters.  These represent
#' equally plausable values of the parameters. Properties of the fitted model
#' can be estimated by averaging with these values.
#'
#' The values are sampled using Markov Chain Monte Carlo (\acronym{MCMC}) via
#' WinBUGS or OpenBUGS. This can take much time to run and convergence of the
#' markov chains should always be verified using \code{\link{plotBUGSTrace}}
#' before the posterior samples can be trusted.
#'
#' Approximate parameter estimates can first be produced using Integrated
#' Nested Laplace Approximations (\acronym{INLA}).  This method produces
#' estimates in seconds but cannot be used to fit the full \acronym{FCS2}
#' model. Instead \acronym{INLA} fits two approximate models, one for the
#' prevalence regression and one for the abundance.  The \acronym{INLA} fits
#' may be used to predict the significance of regression terms in the model.
#' See \code{\link{fcs2RunINLA}} for further details of these models.
#'
#' @param runTotalVars a character vector of columns in \code{dataFrame} that
#' give the number of fish caught in each run of a multiple-run survey.  These
#' should be in order from the first run to the last.
#' @param allRunsTotalVar the name of a column in \code{dataFrame} that gives
#' the total number of fish caught over all runs of a multiple-run survey.
#' @param allRunsRangeVars the names of two columns in \code{dataFrame} that
#' give the minimum and the maximum value of a range of possible values for the
#' total number of fish caught over all runs in a multiple-run survey.
#' @param muFormula a \code{\link{formula}} specifying which terms should
#' appear in the abundance regression equation.  This should take the form
#' \code{~ term1 + term2 + } \ldots{} where each term is either a standard
#' linear term, a non-linear first-order random walk term (given by
#' \code{\link{rw1}}), a second-order random walk term (given by
#' \code{\link{rw2}}) or a spatial term (given by \code{\link{spatial}}).  A
#' constant term is always present.
#' @param rhoFormula a \code{\link{formula}} specifying which terms should
#' appear in the prevalence regression equation.  As with \code{muFormula}, the
#' regression can contain linear terms but also non-linear terms and a single
#' spatial term.
#' @param dataFrame a data frame with surveys as rows and variables as columns.
#' It should contain all variables specified through other arguments.
#' @param surveyAreaVar the name of a column in \code{dataFrame} that gives the
#' survey area.  If not specified, the function will search for the default
#' value \code{"SurveyArea"}.
#' @param nRunsVar the name of a column in \code{dataFrame} that gives the
#' number of runs in each survey.  If missing, the number of runs is assumed to
#' be the number of non-missing run total entries, unless \code{runTotalVars}
#' has length 1 or is missing in which case a single-run model is used.  Number
#' of runs values greater than the number of run total columns are clipped with
#' a warning.
#' @param subset an optional vector specifying a subset of surveys to be used
#' in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#' contain missing values (\code{NA}s).  The default is set by the
#' \code{na.action} setting of \code{\link{options}} and this is usually set to
#' \code{\link{na.omit}}.  This setting removes surveys that contain missing
#' data in any required variables.  A vector indicating the rows that were
#' removed can be extracted from the returned object using
#' \code{\link{na.action.fcs2Fit}}.  Alternatively, \code{\link{na.pass}} can
#' be used to ignore missing values (where possible) or \code{\link{na.fail}}
#' can be given to signal an error if missing values are found.
#' @param prior.parameters an optional named list of named vectors giving the
#' parameter values to use for the prior distribution of a variable.  See
#' \code{\link{fcs2Priors}} for further details on the prior distributions and
#' how to specify prior parameters.  Priors that are not given parameter values
#' have defaults given to them by \code{\link{.fcs2SetDefaultPriors}}.  Note
#' that priors for linear variables are ignored by \acronym{INLA} but all
#' priors are used by \acronym{BUGS} to fit the full model.
#' @param initial.values an optional named list of initial values to use for
#' variables when fitting the full model using WinBUGS or OpenBUGS.  Initial
#' values that are not specified are estimated using the median estimate from
#' \acronym{INLA}, unless \code{runINLA} is \code{FALSE} in which case
#' \acronym{BUGS} will simulate an initial value from the prior.
#' @param runINLA whether to run \acronym{INLA} to fit approximate abundance
#' and prevalence models.  By default, \acronym{INLA} is run for both models
#' only if necessary for providing initial values for all variables.
#' \code{runINLA} can be \code{"rho"} or \code{"mu"} to run only the prevalence
#' or abundance models respectively.  However, \code{"mu"} can only be given if
#' the \acronym{INLA} fit for prevalence has already been run and is provided
#' through \code{fit}.  See \code{\link{fcs2RunINLA}} for details.
#' @param runBUGS whether to run WinBUGS or OpenBUGS to fit the full
#' \acronym{FCS2} model by sampling from the posterior distribution of the
#' parameters.  Note that \acronym{BUGS} can often take weeks to run with a
#' large data set, many regression terms, and a large number of iterations
#' \code{n.iter}.  By default, \code{runBUGS} is \code{FALSE} and the model is
#' simply approximated using \acronym{INLA}.  However, posterior samples from
#' \acronym{BUGS} are required for calculating \acronym{EQR}s using
#' \code{\link{fcs2SingleEQR}}, \code{\link{fcs2JointEQR}} or
#' \code{\link{fcs2JointAndSingleEQR}}.
#' @param estAllRunsTotalVar the name of a column in \code{dataFrame} that
#' gives an estimate of the total number of fish caught over all runs for
#' surveys where only range data is available.  This is used in the approximate
#' abundance \acronym{INLA} fit only and is not used in the full model as
#' fitted with \acronym{BUGS}.  If this is not provided and range data are
#' present, the central value of each range is used with a warning.
#' @param n.chains the number of \acronym{MCMC} chains to run using
#' \acronym{BUGS}. The default is \code{2} and at least 2 chains are required
#' if convergence measures are required.
#' @param n.iter the total number of iterations per chain (including burn-in).
#' The default is \code{2000} but it is likely that considerably more will be
#' required if a few regression terms are included.
#' @param n.burnin the length of burn-in, i.e. number of iterations to discard
#' at the beginning.  The default is \code{n.iter / 2}, that is to discarding
#' the first half of the simulations.
#' @param n.thin the thinning rate. Must be a positive integer. Set
#' \code{n.thin > 1} to save memory and computation time if \code{n.iter} is
#' large.  The default is \code{max(1, floor(n.chains * (n.iter - n.burnin) /
#' 1000))} which will only thin if there are at least \code{2000} simulations.
#' @param n.sims the approximate number of simulations to keep after thinning.
#' @param bugsFilename the name to use for the \acronym{BUGS} model file. The
#' default is \code{"model.bug"}.
#' @param bugsProgram the \acronym{BUGS} program, either \code{winbugs},
#' \code{WinBUGS}, \code{openbugs} or \code{OpenBUGS}.  If using OpenBUGS, the
#' package \pkg{BRugs} must be installed.
#' @param verbose whether to print progress to screen.
#' @param fit an \code{"fcs2Fit"} object, as returned by this function, can
#' alternatively be provided instead of the arguments that specify the model
#' and data.  This allows the object returned after setting \code{runINLA} or
#' \code{runBUGS} to \code{FALSE} to be extended, for example to view
#' approximate INLA fits before running \acronym{BUGS}.
#' @param \dots further arguments that are passed to \code{openbugs}.
#' @return Returns an \code{"fcs2Fit"} object that contains the model
#' specification, the required data and any \acronym{INLA} and \acronym{BUGS}
#' model fits.  The \code{"fcs2Fit"} object is essentially a list with the
#' following items: \item{call}{ the matched call.  } \item{modelMatrix}{ a
#' matrix of variables required by the model, extracted from \code{dataFrame}.
#' } \item{muLinearVars}{ a character vector of linear terms in the abundance
#' regression, each of which appears as a column in \code{modelMatrix} }
#' \item{muRW1Vars, muRW2Vars}{ a character vector of variables with first and
#' second-order random walk terms in the abundance regression.  }
#' \item{muSpatialVar}{ the name of the abundance spatial variable, if present.
#' } \item{muAdjacency}{ the adjacency information for the abundance spatial
#' term, if present.  } \item{rhoLinearVars, rhoRW1Vars, rhoRW2Vars,
#' rhoSpatialVar, rhoAdjacency}{ as their \code{mu} counterparts but applying
#' to the prevalence regression rather than abundance } \item{dataType}{ a
#' \code{\link{factor}} with levels \code{run}, \code{total} and \code{range}
#' indicating whether each survey contains the number of fish caught in each
#' \emph{run}, the \emph{total} number over all runs, or a minimum and maximum
#' describing a \emph{range} of possible values of the total catch.  }
#' \item{rwNoLevels}{ a named list giving the number of discrete levels to use
#' to represent each random walk term.  } \item{rwBoundaries}{ a named list
#' giving the locations of the discrete points that are used to represent each
#' non-linear covariate term.  } \item{covariateMin, covariateMax}{ the minimum
#' and maximum value for each covariate in the model.  } \item{N}{ the number
#' of surveys.  } \item{multiRun}{ whether the multiple-run extension to the
#' \acronym{FCS2} model is used.  This is necessary if more than one run total
#' variable is specified by \code{runTotalVars}, or if the number of runs
#' variable \code{nRunsVar} exceeds 1.  } \item{inlaFits}{ if \acronym{INLA}
#' has been run, a list with two components: \tabular{ll}{ \code{rhoFit} \tab
#' an \code{"inla"} object containing the approximate prevalence model fit \cr
#' \code{muFit} \tab an \code{"inla"} object containing the approximate
#' abundance model fit \cr } See \code{inla} from the package \pkg{INLA}
#' for a description of \code{"inla"} objects.  } \item{bugsFit}{ if
#' \acronym{BUGS} has been run, a \code{"bugs"} object containing the full
#' \acronym{FCS2} model fit.  See \code{bugs} from the package
#' \pkg{R2WinBUGS} for a description of this object.  } In addition, the
#' following arguments are stored: \code{runTotalVars}, \code{allRunsTotalVar},
#' \code{allRunsRangeVars}, \code{muFormula}, \code{rhoFormula},
#' \code{surveyAreaVar} and \code{nRunsVar}.  Also, \code{prior.parameters} and
#' \code{initial.values} are corrected and completed.
#' @note WinBUGS
#' (\url{http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/contents.shtml}) or the
#' newer OpenBUGS (\url{http://www.openbugs.info/w/}) must be separately
#' installed in order to fit the full \acronym{FCS2} model.  If OpenBUGS is
#' used, the package \pkg{BRugs} must also be installed.  The latest version of
#' \pkg{BRugs} comes with OpenBUGS bundled within it so that an external
#' installation is no longer necessary.
#'
#' Note that \acronym{INLA} only uses the priors for the shape parameter
#' \eqn{r} and the hyperparameter for each non-linear term.
#' @section FCS2 model: The \acronym{FCS2} model is a Bayesian hierarchical
#' model of fish counts over a number of surveys.  The original \acronym{EA}
#' \acronym{FCS2} model assumes a single fish count \eqn{C} per survey that is
#' modelled by a zero-inflated Negative Binomial (\acronym{ZINB}) distribution
#' (see \code{\link{rzinbinom}}).  Specifically, \deqn{C \sim ZINB(r, m=a \mu, z=1 -
#' \rho)}{C ~ ZINB(r, m=a \mu, z=1 - \rho)} where \eqn{r} is a shape parameter,
#' \eqn{a} is the survey area, \eqn{\mu} is the \dfn{abundance} (mean density)
#' and \eqn{\rho} is the \dfn{prevalence} (probability present).
#'
#' The \acronym{SNIFFER} multiple-pass extension to the \acronym{FCS2} model
#' instead assumes that joint fish count \eqn{C_1, \ldots, C_d} over \eqn{d}
#' runs follows a zero-inflated Negative Multinomial (\acronym{ZINM})
#' distribution (see \code{fcs2:::zinmultinom}).  Specifically, \deqn{(C_1, \ldots,
#' C_d) \sim ZINM(r, m_j = a \mu q (1 - q)^{(j - 1)}, z=1 - \rho)}{(C_1,
#' \ldots, C_d) ~ ZINM(r, m_j = a \mu q (1 - q)^(j - 1), z=1 - \rho)} where the
#' additional parameter \eqn{q} represents the catch probability.
#'
#' The remaining components of the model are common to both the original and
#' multiple-pass \acronym{FCS2}.  The abundance and prevalence components are
#' both modelled using regressions in terms of other covariates
#' \eqn{x_{i,j}}{x_i,j}.  Specifically, for survey \eqn{i}, \deqn{log(\mu_i) =
#' \beta_0 + \ldots + f_j(x_{i,j}; \beta_j) + \ldots}{log(\mu_i) = \beta_0 +
#' \ldots + f_j(x_i,j; \beta_j) + \ldots} \deqn{logit(\rho_i) = \gamma_0 +
#' \ldots + g_j(x_{i,j}; \gamma_j) + \ldots}{log(\mu_i) = \gamma_0 + \ldots +
#' g_j(x_i,j; \gamma_j) + \ldots} where \eqn{f_j} and \eqn{g_j} are regression
#' terms that depend upon parameters \eqn{\beta_j} and \eqn{\gamma_j}
#' respectively.  The functions \link{log} and \link{logit} are used to
#' transform the components to the real line.  The regression equations are
#' specified by setting \code{muFormula} and \code{rhoFormula} using symbolic
#' formulae.
#'
#' The terms can be either linear, non-linear or spatial.  Linear terms such as
#' \eqn{\beta_1 x + \beta_2 x^2} are specified as standard
#' \code{\link{formula}}.  Non-linear relationships may be constructed using
#' random walk terms using \code{\link{rw1}} for first-order or
#' \code{\link{rw2}} for second-order random walks.  A single spatial term may
#' be added using \code{\link{spatial}}.  This relies upon a geographical
#' partition of the land into a number of spatial regions.  The spatial
#' relationship is modelled by assuming correlations between neighbouring
#' regions.
#'
#' Being a Bayesian model, prior distributions are required for every unknown
#' parameter. These should represent your beliefs about their values prior to
#' incorporating the information from the dataset.  It is common to specify
#' wide uninformative prior distributions to represent prior ignorance.  The
#' default priors attempt to do this but try to balance against too wide priors
#' which can cause WinBUGS to fail when fitting the model.  The prior
#' distributions should always be checked before fitting the full model.  See
#' \code{\link{fcs2Priors}} for further details on the prior distributions and
#' how to specify prior parameters via \code{prior.parameters}.
#' @seealso \code{inla} from the package \pkg{INLA} which is used to
#' provide the approximate model fits using Integrated Nested Laplace
#' Approximation;\cr \code{bugs} from the package \pkg{BRUGS} which
#' is used to provide the full model fit using Markov Chain Monte Carlo
#' (\acronym{MCMC});\cr
#'
#' \code{\link{print.fcs2Fit}} and \code{\link{summary.fcs2Fit}} for
#' summarising \code{"fcs2Fit"} objects;\cr \code{\link{plot.fcs2Fit}},
#' \code{\link{plotSpatialTerm}} and \code{\link{plotCatchPMF}} for plotting
#' the \acronym{FCS2} model fit;\cr \code{\link{plotBUGSTrace}} for assessing
#' the MCMC convergence;\cr \code{\link{ppplot}} for assessing the model fit
#' with a \acronym{P-P} plot;\cr
#'
#' \code{\link{fcs2SingleEQR}}, \code{\link{fcs2JointEQR}} or
#' \code{\link{fcs2JointAndSingleEQR}} for producing single or joint
#' \acronym{EQR} samples.
#' @keywords models regression nonlinear multivariate spatial
#' @export
#' @examples
#'
#' ### Example 1: Very simple example with no covariates
#' ###
#'
#' # simulate random dataset
#' Data <- data.frame(SurveyArea=rlnorm(100, 4.6, 0.5))   # random survey area
#' Data$Catch <- rzinbinom(100, size=1.1, zeroprob=0.3,
#'                         nbmean=0.3 * Data$SurveyArea)  # single catch per survey
#'
#' # fit approximate model with INLA - model contains no regression terms
#' fit1 <- fcs2FitModel("Catch", dataFrame=Data, surveyAreaVar="SurveyArea")
#'
#' \dontrun{
#' # fit full model with OpenBUGS
#' # (more iterations may be required to achieve convergence)
#' fit1 <- fcs2FitModel(fit=fit1, runBUGS=TRUE, n.iter=1000,
#'                      bugsProgram="OpenBUGS")  }
#'
#' # summarise fit
#' summary(fit1)
#'
#'
#'
#'
#' ### Example 2: Multiple-run data with prevalence affected by a barrier
#' ###
#'
#' # add to simulated dataset
#' Data$NumberOfRuns <- 3   # 3 runs
#' Data$Barrier <- runif(100) > 0.8   # 20% chance that Barrier = TRUE
#'
#' # 3 run catch with barrier effecting chance of 0
#' Catch <- rzinmultinom(100, size=1.1,
#'                       zeroprob=expit(-1 + 2 * Data$Barrier),
#'                       nmmean=0.3 * Data$SurveyArea %o% 0.4^(0:2))
#' Data$Run1 <- Catch[,1]
#' Data$Run2 <- Catch[,2]
#' Data$Run3 <- Catch[,3]
#'
#' # define model with a single barrier term in prevalence regression
#' # and run INLA to provide an approximate fit
#' fit2 <- fcs2FitModel(c("Run1", "Run2", "Run3"), dataFrame=Data,
#'                      surveyAreaVar="SurveyArea", nRunsVar="NumberOfRuns",
#'                      rhoFormula= ~ Barrier)
#'
#' # summarise fit to see whether barrier term is significant
#' summary(fit2)
#'
#' \dontrun{
#' # fit full model with OpenBUGS
#' fit2 <- fcs2FitModel(fit=fit2, runBUGS=TRUE, n.iter=10000,
#'                      bugsProgram="OpenBUGS")  }
#'
#' # plot the parameter estimates
#' plot(fit2, group=FALSE)
#'
#'
#'
#'
#' ### Example 3: Detailed example with multiple terms
#' ###
#'
#' # define adjacency information for 3 spatial regions appearing in a line
#' adjacency <- list(num=c(1, 2, 1),  # number of regions adjacent to each region
#'                   adj=c(2,
#'                         1, 3,
#'                         2),        # indices of these regions
#'                   sumNumNeigh=4)   # total number of adjacencies
#'
#' # add to simulated dataset
#' # randomly place each survey in one of 3 regions
#' Data$RegionNumber <- sample(1:3, 100, replace=TRUE)
#' Data$Altitude <- rlnorm(100, 4.4, 1.1)   # random altitude
#' Data$WetWidth <- rlnorm(100, 1.5, 0.6)   # random wet width
#'
#' # simulate catch that varies with barrier, region, altitude but not wet width
#' Data$Catch <- rzinbinom(100, size=1.1,
#'                         zeroprob=expit(-1.3 + 2.9 * Data$Barrier +
#'                                        2.55 * (Data$RegionNumber == 1) -
#'                                        3.7 * (Data$RegionNumber == 3)),
#'                         nbmean=0.3 * Data$SurveyArea *
#'                                exp(-9.5 + 4 * log(Data$Altitude) -
#'                                    0.5 * log(Data$Altitude)^2))
#'
#' # define model with barrier and spatial terms in prevalence and
#' #   with log(altitude) and log(wet width) in abundance
#' # run INLA to provide an approximate fit
#' fit3 <- fcs2FitModel("Catch", dataFrame=Data, surveyAreaVar="SurveyArea",
#'                      rhoFormula= ~ Barrier + spatial(RegionNumber, adjacency),
#'                      muFormula= ~ log(Altitude) + I(log(Altitude)^2) +
#'                                 rw2(log(WetWidth)))
#'
#' # summarise fit
#' summary(fit3)
#'
#' # plot fit
#' plot(fit3)
#'
#' # fit again without wet width term and restrict the variability
#' #   between regions in the spatial term
#' fit3 <- fcs2FitModel("Catch", dataFrame=Data, surveyAreaVar="SurveyArea",
#'                      rhoFormula= ~ Barrier +
#'                                    spatial(RegionNumber, adjacency,
#'                                            scale.parameters=c(mean=1, sd=3)),
#'                      muFormula= ~ log(Altitude) + I(log(Altitude)^2))
#'
#' # summarise fit, showing each of the variables in the spatial term
#' summary(fit3, allVars=TRUE)
#'
#' # plot fit
#' plot(fit3)
#'
#' # adjust the priors to make sure non-informative and fit again
#' fit3 <- fcs2FitModel("Catch", dataFrame=Data, surveyAreaVar="SurveyArea",
#'                      rhoFormula= ~ Barrier +
#'                                    spatial(RegionNumber, adjacency,
#'                                            scale.parameters=c(mean=1, sd=3)),
#'                      muFormula= ~ log(Altitude) + I(log(Altitude)^2),
#'                      prior.parameters=list(beta.const=c(mean=0, sd=30),
#'                                            "beta.log(Altitude)"=c(mean=0, sd=10),
#'                                            gamma.BarrierTRUE=c(mean=0, sd=5)))
#'
#' # plot fit
#' plot(fit3)
#'
#' \dontrun{
#'
#' # fit full model with OpenBUGS using a large number of iterations
#' #   - this may take some time or
#' #     even crash if some terms are not significant
#' fit3 <- fcs2FitModel(fit=fit3, runBUGS=TRUE, n.iter=10000,
#'                      bugsProgram="OpenBUGS")
#'
#' # check MCMC chains converged and mixed well
#' plotBUGSTrace(fit3)
#'
#' # plot fit
#' plot(fit3)
#'
#' # plot predicted catch distribution, for first 9 surveys
#' plotCatchPMF(fit3, Data, 1:9, boundaries=NULL)
#' }
#'
fcs2FitModel <-
function(runTotalVars = NULL, allRunsTotalVar = NULL, allRunsRangeVars = NULL, muFormula = ~1, rhoFormula = ~1, dataFrame,
         surveyAreaVar = "SurveyArea", nRunsVar = NULL, subset = 1:nrow(dataFrame), na.action,
         prior.parameters = list(), initial.values = list(), runINLA, runBUGS = FALSE, estAllRunsTotalVar = NULL,
         n.chains = 2, n.iter = 2000, n.burnin = floor(n.iter/2), n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)),
         n.sims = 1000, bugsFilename = "model.bug", bugsProgram = "OpenBUGS", verbose = TRUE, fit, ...)  # further arguments to BUGS
{
    # get default na.action if missing
    if (missing(na.action)) {
        na.action <- getOption("na.action")
        if (is.null(na.action))
            na.action <- na.omit  # use 'na.omit' if option not set
        else
            na.action <- eval(parse(text=na.action))
    }

    # if fit object provided, call .fcs2FitModelFromComponent
    if (!missing(fit)) {
        if (missing(runINLA))
            return(.fcs2FitModelFromComponent(call=match.call(), n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, n.sims=n.sims,
                                              bugsFilename=bugsFilename, bugsProgram=bugsProgram, verbose=verbose, fit=fit, runBUGS=runBUGS, ...))
        else
            return(.fcs2FitModelFromComponent(call=match.call(), n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, n.sims=n.sims,
                                              bugsFilename=bugsFilename, bugsProgram=bugsProgram, verbose=verbose, fit=fit, runBUGS=runBUGS, runINLA=runINLA, ...))
    }

    # set subset if missing
    if (missing(subset))
        subset <- 1:nrow(dataFrame)

    # make sure subset is not logical
    ## NOTE: subset may still be names of rows but this is ok
    if (inherits(subset, "logical") || min(subset) == 0)
        subset <- which(as.logical(subset))

    # check catch variables are in data frame
    if (length(var <- which(!(runTotalVars %in% names(dataFrame)))) > 0)
        stop(paste("run total variable '", runTotalVars[var[1]], "' not found in data frame", sep=""))
    if (length(var <- which(!(allRunsTotalVar %in% names(dataFrame)))) > 0)
        stop(paste("all runs total variable '", allRunsTotalVar[var[1]], "' not found in data frame", sep=""))
    if (length(var <- which(!(allRunsRangeVars %in% names(dataFrame)))) > 0)
        stop(paste("all runs range variable '", allRunsRangeVars[var[1]], "' not found in data frame", sep=""))
    if (length(var <- which(!(surveyAreaVar %in% names(dataFrame)))) > 0)
        stop(paste("survey area variable '", surveyAreaVar[var[1]], "' not found in data frame", sep=""))
    if (length(var <- which(!(nRunsVar %in% names(dataFrame)))) > 0)
        stop(paste("number of runs variable '", nRunsVar[var[1]], "' not found in data frame", sep=""))
    if (length(var <- which(!(estAllRunsTotalVar %in% names(dataFrame)))) > 0)
        stop(paste("all runs total estimate variable '", estAllRunsTotalVar[var[1]], "' not found in data frame", sep=""))


    ## remove 'rw1', 'rw2' and 'spatial' from formulae

    # mu
    muTerms <- terms(muFormula, specials=c("rw1", "rw2", "spatial"))
    muTermLabels <- attr(muTerms, "term.labels")
    specials <- attr(muTerms, "specials")
    iMuSpecials <- c(specials$rw1, specials$rw2, specials$spatial)
    if (length(iMuSpecials) > 0) {
        muSpecialTerms <- muTermLabels[iMuSpecials]
        muTermLabels <- muTermLabels[-iMuSpecials]
    }

    # rho
    rhoTerms <- terms(rhoFormula, specials=c("rw1", "rw2", "spatial"))
    rhoTermLabels <- attr(rhoTerms, "term.labels")
    specials <- attr(rhoTerms, "specials")
    iRhoSpecials <- c(specials$rw1, specials$rw2, specials$spatial)
    if (length(iRhoSpecials) > 0) {
        rhoSpecialTerms <- rhoTermLabels[iRhoSpecials]
        rhoTermLabels <- rhoTermLabels[-iRhoSpecials]
    }


    ## extract linear variables

    # mu
    simpleFormula <- formula(paste("~", paste(c("1", muTermLabels), collapse="+")))
    modelMatrix <- model.matrix(simpleFormula, dataFrame)[, -1, drop=FALSE]
    muLinearVars <- colnames(modelMatrix)

    # rho
    simpleFormula <- formula(paste("~", paste(c("1", rhoTermLabels), collapse="+")))
    modelMatrix <- model.matrix(simpleFormula, dataFrame)[, -1, drop=FALSE]
    rhoLinearVars <- colnames(modelMatrix)


    ## add 'rw1', 'rw2' and 'spatial' information

    # mu
    muRW1Vars <- muRW2Vars <- c()
    muSpatialVar <- muAdjacency <- NULL
    rwNoLevels <- rwBoundaries <- list()
    if (length(iMuSpecials) > 0) {
        for (text in muSpecialTerms) {
            # extract term info
            term <- eval(parse(text=text), dataFrame)
            muTermLabels <- c(muTermLabels, term$name)

            # set variable name
            varName <- paste("beta", term$name, sep=".")

            # add variable to data frame (if not already there)
            if (!(term$name %in% colnames(dataFrame))) {
                dataFrame <- cbind(dataFrame, term$val)
                colnames(dataFrame)[ncol(dataFrame)] <- term$name
            }

            # add variable to correct list
            if (term$type == "rw1")
                muRW1Vars <- c(muRW1Vars, term$name)
            else if (term$type == "rw2")
                muRW2Vars <- c(muRW2Vars, term$name)
            else if (term$type == "spatial") {
                muSpatialVar <- term$name
                muAdjacency <- term$adjacency
            }

            # add boundaries or number of boundary levels
            if (term$type != "spatial") {
                if (is.null(term$boundaries)) {
                    rwNoLevels <- c(rwNoLevels, list(term$noLevels))
                    names(rwNoLevels)[length(rwNoLevels)] <- varName

                } else {
                    rwBoundaries <- c(rwBoundaries, list(term$boundaries))
                    names(rwBoundaries)[length(rwBoundaries)] <- varName
                }
            }

            # add parameter info to prior.parameters
            if (!is.null(term$scale.parameters)) {
                prior.parameters <- c(prior.parameters, list(term$scale.parameters))
                names(prior.parameters)[length(prior.parameters)] <- paste("sigma", term$name, sep=".")

            } else if (!is.null(term$precision.parameters)) {
                prior.parameters <- c(prior.parameters, list(term$precision.parameters))
                names(prior.parameters)[length(prior.parameters)] <- paste("tau", term$name, sep=".")
            }

            # add initial values
            if (!is.null(term$initial.values)) {
                initial.values <- c(initial.values, list(term$initial.values))
                names(initial.values)[length(initial.values)] <- varName
            }
            if (!is.null(term$precision.initial.value)) {
                initial.values <- c(initial.values, list(term$precision.initial.value))
                names(initial.values)[length(initial.values)] <- paste("tau", term$name, sep=".")
            }
        }
    }

    # rho
    rhoRW1Vars <- rhoRW2Vars <- c()
    rhoSpatialVar <- rhoAdjacency <- NULL
    if (length(iRhoSpecials) > 0) {
        for (text in rhoSpecialTerms) {
            # extract term info
            term <- eval(parse(text=text), dataFrame)
            rhoTermLabels <- c(rhoTermLabels, term$name)

            # set variable name
            varName <- paste("gamma", term$name, sep=".")
            # add variable to model matrix (if not already there)
            if (!(term$name %in% colnames(dataFrame))) {
                dataFrame <- cbind(dataFrame, term$val)
                colnames(dataFrame)[ncol(dataFrame)] <- term$name
            }

            # add variable to correct list
            if (term$type == "rw1")
                rhoRW1Vars <- c(rhoRW1Vars, term$name)
            else if (term$type == "rw2")
                rhoRW2Vars <- c(rhoRW2Vars, term$name)
            else if (term$type == "spatial") {
                rhoSpatialVar <- term$name
                rhoAdjacency <- term$adjacency
            }

            # add boundaries or number of boundary levels
            if (term$type != "spatial") {
                if (is.null(term$boundaries)) {
                    rwNoLevels <- c(rwNoLevels, list(term$noLevels))
                    names(rwNoLevels)[length(rwNoLevels)] <- varName

                } else {
                    rwBoundaries <- c(rwBoundaries, list(term$boundaries))
                    names(rwBoundaries)[length(rwBoundaries)] <- varName
                }
            }

            # add parameter info to prior.parameters
            if (!is.null(term$scale.parameters)) {
                prior.parameters <- c(prior.parameters, list(term$scale.parameters))
                names(prior.parameters)[length(prior.parameters)] <- paste("nu", term$name, sep=".")

            } else if (!is.null(term$precision.parameters)) {
                prior.parameters <- c(prior.parameters, list(term$precision.parameters))
                names(prior.parameters)[length(prior.parameters)] <- paste("phi", term$name, sep=".")
            }

            # add initial values
            if (!is.null(term$initial.values)) {
                initial.values <- c(initial.values, list(term$initial.values))
                names(initial.values)[length(initial.values)] <- varName
            }
            if (!is.null(term$precision.initial.value)) {
                initial.values <- c(initial.values, list(term$precision.initial.value))
                names(initial.values)[length(initial.values)] <- paste("phi", term$name, sep=".")
            }
        }
    }


    # if multiple run data ...
    if (!is.null(runTotalVars) && length(runTotalVars) > 1) {

        # create nRunsVar if missing (but only bother for subset)
        if (is.null(nRunsVar)) {
            nRunsVar <- "NumberOfRuns"
            nRuns <- rep(NA, nrow(dataFrame))
            for (i in (1:nrow(dataFrame))[subset])
                nRuns[i] <- sum(!is.na(dataFrame[i, runTotalVars]))
            dataFrame <- cbind(dataFrame, NumberOfRuns=nRuns)
        }

        # check number of runs variable isn't larger than number of catch variables
        if (max(dataFrame[subset, nRunsVar], na.rm=TRUE) > length(runTotalVars)) {
            warning(paste("Number of runs variable '", nRunsVar, "' clipped to not exceed number of catch variables (", length(runTotalVars), ")", sep=""))
            dataFrame[subset, nRunsVar][!is.na(dataFrame[subset, nRunsVar]) & dataFrame[subset, nRunsVar] > length(runTotalVars)] <- length(runTotalVars)
        }

        # fill in runs that didn't take place with 0s (but only bother for subset)
#         for (i in (1:nrow(dataFrame))[subset])
#             if (!is.na(dataFrame[i, nRunsVar]) && dataFrame[i, nRunsVar] < length(runTotalVars))
#                 dataFrame[i, runTotalVars[(dataFrame[i, nRunsVar] + 1):length(runTotalVars)]] <- 0
        zero <- array(FALSE, dim=c(length(subset), length(runTotalVars)))
        for (j in 1:ncol(zero))
            zero[, j] <- j > dataFrame[subset, nRunsVar]
        dataFrame[subset, runTotalVars][zero] <- 0


        # check nRunsVar present if needed for all run total or range
    } else if (is.null(nRunsVar) && (!is.null(allRunsTotalVar) || !is.null(allRunsRangeVars))) {
        # create nRuns = 1 with warning
        nRunsVar <- "NumberOfRuns"
        dataFrame <- cbind(dataFrame, NumberOfRuns=1)
        warning("Number of runs variable 'nRunsVar' is missing so a single run is assumed for all surveys")
    }


    # create variable that stores the best type of catch data available for each survey (as "run", "total", "range" factor)
    dataType <- factor(rep(NA, length(subset)), c("run", "total", "range"))
    if (!is.null(runTotalVars))
        dataType[!is.na(dataFrame[subset, runTotalVars[1]])] <- "run"
    if (!is.null(allRunsTotalVar))
        dataType[!is.na(dataFrame[subset, allRunsTotalVar]) & is.na(dataType)] <- "total"
    if (!is.null(allRunsRangeVars))
        dataType[!is.na(dataFrame[subset, allRunsRangeVars[1]]) & !is.na(dataFrame[subset, allRunsRangeVars[2]]) & is.na(dataType)] <- "range"


    # if multi-run data or all run info given ...
    if ((!is.null(runTotalVars) && length(runTotalVars) > 1) || !is.null(nRunsVar)) {

        # create allRunsTotalVar if missing
        if (is.null(allRunsTotalVar)) {
            allRunsTotalVar <- "AllRunsTotal"
            dataFrame <- cbind(dataFrame, AllRunsTotal=NA)
        }

        # fill in all run total blanks with catch run data, if present
        if (!is.null(runTotalVars)) {
            allRunTotal <- apply(dataFrame[subset, runTotalVars, drop=FALSE], 1, sum, na.rm=TRUE)  # NOTE: this will give 0 when all missing
            dataFrame[subset[is.na(dataFrame[subset, allRunsTotalVar]) & !is.na(dataFrame[subset, runTotalVars[1]])], allRunsTotalVar] <-
                        allRunTotal[is.na(dataFrame[subset, allRunsTotalVar]) & !is.na(dataFrame[subset, runTotalVars[1]])]
        }

        # fill in all run total blanks with estimated total or with central value of all runs range data (rounded down)
        if (!is.null(allRunsRangeVars)) {
            if (!is.null(estAllRunsTotalVar))
                allRunTotal <- dataFrame[subset, estAllRunsTotalVar]
            else
                allRunTotal <- floor(apply(dataFrame[subset, allRunsRangeVars, drop=FALSE], 1, mean, na.rm=TRUE))
            if (sum(is.na(dataFrame[subset, allRunsTotalVar]) & !is.nan(allRunTotal)) > 0) {
                # issue warning if taking central value of range
                if (is.null(estAllRunsTotalVar))
                    warning("using centre of range as all runs total estimate for INLA abundance fit")

                dataFrame[subset[is.na(dataFrame[subset, allRunsTotalVar]) & !is.nan(allRunTotal)], allRunsTotalVar] <-
                        allRunTotal[is.na(dataFrame[subset, allRunsTotalVar]) & !is.nan(allRunTotal)]
            }
        }

    } else {
        # if allRunsTotalVar still missing, set to runTotalVars[1]
        ## NOTE: can only get here if there is a single run total
        if (is.null(allRunsTotalVar))
            allRunsTotalVar <- runTotalVars[1]
    }


    # fill missing catch data with 0s when other type is available
    if (!is.null(runTotalVars) && sum(!is.na(dataType) & dataType != "run") > 0)
        dataFrame[subset[!is.na(dataType) & dataType != "run"], runTotalVars] <- 0
    if (!is.null(allRunsRangeVars) && sum(!is.na(dataType) & dataType != "range") > 0)
        dataFrame[subset[!is.na(dataType) & dataType != "range"], allRunsRangeVars] <- 0


    # create model matrix containing required variables only
    simpleFormula <- formula(paste("~", paste(c(runTotalVars, allRunsTotalVar, allRunsRangeVars,
                                                surveyAreaVar, nRunsVar, muTermLabels, rhoTermLabels), collapse="+")))
    modelFrame <- model.frame(simpleFormula, dataFrame, subset=subset, na.action=na.action)
    na.action <- attr(modelFrame, "na.action")
    modelMatrix <- model.matrix(simpleFormula, modelFrame)[, -1, drop=FALSE]

    # remove entries or set NA values in dataType
    if (!is.null(na.action))
        dataType <- dataType[-na.action]
    else if (!is.null(na.action.omit <- attr(model.frame(simpleFormula, dataFrame, subset=subset, na.action=na.omit), "na.action")))
        dataType[na.action.omit] <- NA

    # call .fcs2FitModelFromComponent
    .fcs2FitModelFromComponent(runTotalVars, allRunsTotalVar, allRunsRangeVars, modelMatrix, surveyAreaVar, nRunsVar,
                               muLinearVars, muRW1Vars, muRW2Vars, muSpatialVar, muAdjacency,
                               rhoLinearVars, rhoRW1Vars, rhoRW2Vars, rhoSpatialVar, rhoAdjacency,
                               1:nrow(modelMatrix), na.action, dataType,
                               rwNoLevels, rwBoundaries, prior.parameters, initial.values,
                               runINLA, runBUGS, n.chains, n.iter, n.burnin, n.thin, n.sims,
                               bugsFilename, bugsProgram, verbose, muFormula, rhoFormula, match.call(), ...)
}

