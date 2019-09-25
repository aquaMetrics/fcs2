#' @encoding UTF-8
#' @title Variable Names in FCS2 Model
#'
#' @description
#' Returns a vector of selected variables in an \acronym{FCS2} model.
#' @name variable.names
#' @param object an object of class \code{"fcs2Fit"}, usually returned by \code{\link{fcs2FitModel}}.
#' @param q whether to include the catch probability variable \eqn{q}.
#' @param r whether to include the shape parameter \eqn{r} (given by \code{size}
#' in \code{\link{dzinbinom}} and \code{\link{dzinmultinom}}).
#' @param abundance whether to include any variables related to the abundance regression.
#' @param prevalence whether to include any variables related to the prevalence regression.
#' @param linear whether to include linear regression terms.
#' @param rw whether and how to include random walk terms. If \code{TRUE} or \code{"group"},
#' a single variable name is returned for each random walk term. If \code{"singluar"},
#' each individual variable is included.
#' @param spatial whether and how to include spatial terms. If \code{TRUE} or \code{"group"},
#' a single variable name is returned for each spatial term. If \code{"singluar"},
#' each individual variable is included.
#' @param hyperparams whether and how to include hyperparameters. If \code{"scale"},
#' scale hyperparameters (\eqn{\sigma} or \eqn{\nu}) are included. If \code{"precision"},
#' hyperparameters are specified on the precision scale (\eqn{\tau} or \eqn{\phi}).
#' @param \dots Not currently used.
#' @return A character vector of selected variables in the \acronym{FCS2} model specified by \code{object}.
#' @seealso \code{\link{fcs2FitModel}} for fitting the \acronym{FCS2} model.

#' @export
variable.names.fcs2Fit <-
  function(object, q = object$multiRun, r = TRUE, abundance = TRUE, prevalence = TRUE, linear = TRUE,
           rw = "group", spatial = "group", hyperparams = "scale", ...)
{
    ##attach(object, warn.conflicts=FALSE)
    ##on.exit(detach(object))
  with(object, {
    params <- c()

    # Catch probability q
    if (q)
        params <- c(params, "q")

    # Shape parameter r
    if (r)
        params <- c(params, "r")

    # Abundance terms:
    if (abundance) {
        # Constant and linear terms
        if (linear)
            params <- c(params, paste("beta", c("const", muLinearVars), sep="."))

        # Random walk terms
        if (length(c(muRW1Vars, muRW2Vars)) > 0) {
            if (!is.na(pmatch(rw, c("group", TRUE))))
                params <- c(params, paste("beta", c(muRW1Vars, muRW2Vars), sep="."))
            else if (!is.na(pmatch(rw, "singular"))) {
                for (var in c(muRW1Vars, muRW2Vars))
                    params <- c(params, paste("beta.", var, "[", 1:length(object$rwBoundaries[[paste("beta", var, sep=".")]]), "]", sep=""))
            }
        }

        # Spatial term
        if (!is.null(muSpatialVar)) {
            if (!is.na(pmatch(spatial, c("group", TRUE))))
                params <- c(params, paste("beta", muSpatialVar, sep="."))
            else if (!is.na(pmatch(spatial, "singular")))
                params <- c(params, paste("beta.", muSpatialVar, "[", 1:length(object$muAdjacency$num), "]", sep=""))
        }

        # Random walk hyperparameters
        if (length(c(muRW1Vars, muRW2Vars)) > 0) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("sigma", c(muRW1Vars, muRW2Vars), sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("tau", c(muRW1Vars, muRW2Vars), sep="."))
        }

        # Spatial hyperparameter
        if (!is.null(muSpatialVar)) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("sigma", muSpatialVar, sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("tau", muSpatialVar, sep="."))
        }
    }

    # Prevalence terms:
    if (prevalence) {
        # Constant and linear terms
        if (linear)
            params <- c(params, paste("gamma", c("const", rhoLinearVars), sep="."))

        # Random walk terms
        if (length(c(rhoRW1Vars, rhoRW2Vars)) > 0) {
            if (!is.na(pmatch(rw, c("group", TRUE))))
                params <- c(params, paste("gamma", c(rhoRW1Vars, rhoRW2Vars), sep="."))
            else if (!is.na(pmatch(rw, "singular"))) {
                for (var in c(rhoRW1Vars, rhoRW2Vars))
                    params <- c(params, paste("gamma.", var, "[", 1:length(object$rwBoundaries[[paste("gamma", var, sep=".")]]), "]", sep=""))
            }
        }

        # Spatial terms
        if (!is.null(rhoSpatialVar)) {
            if (!is.na(pmatch(spatial, c("group", TRUE))))
                params <- c(params, paste("gamma", rhoSpatialVar, sep="."))
            else if (!is.na(pmatch(spatial, "singular")))
                params <- c(params, paste("gamma.", rhoSpatialVar, "[", 1:length(object$rhoAdjacency$num), "]", sep=""))
        }

        # Random walk hyperparameterw
        if (length(c(rhoRW1Vars, rhoRW2Vars)) > 0) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("nu", c(rhoRW1Vars, rhoRW2Vars), sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("phi", c(rhoRW1Vars, rhoRW2Vars), sep="."))
        }

        # Spatial hyperparameter
        if (!is.null(rhoSpatialVar)) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("nu", rhoSpatialVar, sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("phi", rhoSpatialVar, sep="."))
        }
    }

    params
  })
}

