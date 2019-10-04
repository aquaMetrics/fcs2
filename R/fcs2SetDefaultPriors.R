#' Set Default Priors
#'
#' Sets default prior parameters for variables missing in
#' \code{fit$prior.parameters}.  %$
#'
#' The default priors are as follows:
#'
#' The shape parameter \eqn{r} follows a log-Normal distribution with
#' \code{mean=1} and \code{sd=10}.\cr The catch probability \eqn{q} follows a
#' Uniform distribution between 0 and 1.\cr The constant regression terms
#' \code{beta.const} and \code{gamma.const} follow a Normal distribution with
#' \code{mean=0} and \code{sd=5}.\cr All other linear terms \code{beta.}* or
#' \code{gamma.}* follow a Normal distribution with \code{mean=0} and
#' \code{sd=1}.\cr The hyperparameters for random walk terms \code{\link{rw1}}
#' or \code{\link{rw2}} are specified for the scale parameters \code{sigma.}*
#' or \code{nu.}* which follow the \code{\link{hyperprior}} distribution with
#' \code{mean=0.5} and \code{sd=3}.\cr The hyperparameters for
#' \code{\link{spatial}} terms are specified for the precision parameters
#' \code{tau.}* or \code{phi.}* and follow a Gamma distribution with
#' \code{shape=1} and \code{rate=0.001}.
#'
#' The priors should always be checked (with \code{\link{plot.fcs2Fit}} or
#' \code{\link{fcs2Priors}}) before fitting the full \acronym{FCS2} model with
#' \acronym{BUGS} and adjusted if necessary through the \code{prior.parameters}
#' argument of \code{\link{fcs2FitModel}}.
#'
#' @param fit an object of class \code{"fcs2Fit"}, usually returned by
#' \code{\link{fcs2FitModel}}.
#' @return the list \code{fit$prior.parameters} modified to contain parameters
#' for every variable in the model. %$
#' @note The user should not call this function directly but it is called by
#' \code{\link{fcs2FitModel}} when creating a complete \code{"fcs2Fit"} object.
#' @seealso \code{\link{fcs2FitModel}}
#' @export
.fcs2SetDefaultPriors <-
function(fit)
{
  with(fit, {

    # catch probability q
    if (fit$multiRun) {
        if (!("q" %in% names(prior.parameters)))
            prior.parameters <- c(prior.parameters, list(q=c(a=1, b=1)))  # q ~ Beta(a=1, b=1) = Uniform(0, 1)
    }

    # shape parameter r
    if (!("r" %in% names(prior.parameters)))
        prior.parameters <- c(prior.parameters, list(r=c(mean=1, sd=10)))  # r ~ Gamma(mean=1, sd=10)

    # constant terms
    for (var in c("beta.const", "gamma.const")) {
        if (!(var %in% names(prior.parameters))) {
            prior.parameters <- c(prior.parameters, list(c(mean=0, sd=5)))  # *.const ~ Normal(mean=0, sd=5)
            names(prior.parameters)[length(prior.parameters)] <- var
        }
    }

    # linear terms
    if (length(muLinearVars) > 0) {
        for (var in paste("beta.", muLinearVars, sep='')) {
            if (!(var %in% names(prior.parameters))) {
                prior.parameters <- c(prior.parameters, list(c(mean=0, sd=1)))  # beta.* ~ Normal(mean=0, sd=1)
                names(prior.parameters)[length(prior.parameters)] <- var
            }
        }
    }
    if (length(rhoLinearVars) > 0) {
        for (var in paste("gamma.", rhoLinearVars, sep='')) {
            if (!(var %in% names(prior.parameters))) {
                prior.parameters <- c(prior.parameters, list(c(mean=0, sd=1)))  # gamma.* ~ Normal(mean=0, sd=1)
                names(prior.parameters)[length(prior.parameters)] <- var
            }
        }
    }

    # hyperparameters for RW terms
    if (length(c(muRW1Vars, muRW2Vars)) > 0) {
        for (var in paste("sigma.", c(muRW1Vars, muRW2Vars), sep='')) {
            precVar <- sub("sigma", "tau", var)
            if (!(var %in% names(prior.parameters)) && !(precVar %in% names(prior.parameters))) {
                prior.parameters <- c(prior.parameters, list(c(mean=0.5, sd=3)))  # sigma.* ~ hyperprior(mean=0.5, sd=3)
                names(prior.parameters)[length(prior.parameters)] <- var
            }
        }
    }
    if (length(c(rhoRW1Vars, rhoRW2Vars)) > 0) {
        for (var in paste("nu.", c(rhoRW1Vars, rhoRW2Vars), sep='')) {
            precVar <- sub("nu", "phi", var)
            if (!(var %in% names(prior.parameters)) && !(precVar %in% names(prior.parameters))) {
                prior.parameters <- c(prior.parameters, list(c(mean=0.5, sd=3)))  # nu.* ~ hyperprior(mean=0.5, sd=3)
                names(prior.parameters)[length(prior.parameters)] <- var
            }
        }
    }

    # hyperparameters for spatial terms
    if (length(muSpatialVar) > 0) {
        for (var in paste("tau.", muSpatialVar, sep='')) {
            scaleVar <- sub("tau", "sigma", var)
            if (!(var %in% names(prior.parameters)) && !(scaleVar %in% names(prior.parameters))) {
                prior.parameters <- c(prior.parameters, list(c(a=1, b=0.001)))  # tau.* ~ Gamma(a=1, b=0.001)
                names(prior.parameters)[length(prior.parameters)] <- var
            }
        }
    }
    if (length(rhoSpatialVar) > 0) {
        for (var in paste("phi.", rhoSpatialVar, sep='')) {
            scaleVar <- sub("phi", "nu", var)
            if (!(var %in% names(prior.parameters)) && !(scaleVar %in% names(prior.parameters))) {
                prior.parameters <- c(prior.parameters, list(c(a=1, b=0.001)))  # phi.* ~ Gamma(a=1, b=0.001)
                names(prior.parameters)[length(prior.parameters)] <- var
            }
        }
    }

    prior.parameters
  })
}
