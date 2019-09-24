#' Fit Approximate FCS2 Models with INLA
#'
#' Fits two approximate \acronym{FCS2} sub-models using \dfn{Integrated Nested
#' Laplace Approximations} (\acronym{INLA}).
#'
#' \acronym{INLA} is used to fit two models that together approximate the full
#' \acronym{FCS2} statistical model (see \code{\link{fcs2FitModel}}.
#'
#' Firstly, an approximate presence/absence model is fitted to estimate the
#' prevalence regression.  This assumes the observed presence follows a
#' Bernoulli distribution with probability of presence \eqn{\rho} (the
#' prevalence).  The \code{\link{logit}} regression for \eqn{\rho} takes the
#' same form as in the full \acronym{FCS2} model.
#'
#' After fitting the prevalence sub-model, the presence of the species can be
#' estimated.  Data for surveys either observed or predicted to be present are
#' then modelled by the second model that attempts to estimate the abundance
#' regression.  This models the total catch \eqn{T} over all passes by the
#' Negative Binomial distribution with shape \eqn{r} and mean given by \eqn{a
#' \mu}, where \eqn{a} is the survey area and \eqn{\mu} is the abundance.  The
#' regression for \eqn{\mu} takes the same form as in the full \acronym{FCS2}
#' model.
#'
#' This model is exact for the original \acronym{FCS2} model but only an
#' approximation of the total catch when present in the multiple-pass model
#' since it ignores the dependence upon the catch probability \eqn{q} and the
#' number of runs \eqn{d}. However, since \eqn{q} is constant and if the number
#' of runs does not vary significantly, this discrepancy should be mostly
#' absorbed into an error in the constant parameter \eqn{\beta_0}.  If this is
#' true, the model estimates can still be used to estimate the regression terms
#' so that their significance can be judged.
#'
#' @param fit an object of class \code{"fcs2Fit"}, usually returned by
#' \code{\link{fcs2FitModel}}.
#' @param run which approximate \acronym{FCS2} models to fit with
#' \acronym{INLA}.  If \code{TRUE} (the default) or \code{"both"},
#' \acronym{INLA} is run for both models.  Alternatively, if \code{"rho"} or
#' \code{"prevalence"} the prevalence model only is run and if \code{"mu"} or
#' \code{"abundance"} the approximate abundance model only is run.  The
#' abundance model can only be run if the prevalence model is either run first
#' or is provided through \code{fit}.
#' @param verbose whether to print progress to screen.
#' @return a list with two components: \item{rhoFit}{ an \code{"inla"} object
#' containing the approximate prevalence model fit } \item{muFit}{ an
#' \code{"inla"} object containing the approximate abundance model fit, or
#' \code{NULL} if this was not calculated or provided } See \code{inla}
#' from the package \pkg{INLA} for a description of the \code{"inla"} object.
#' @note The user would not usually call this function directly as it is called
#' if requested by \code{\link{fcs2FitModel}}.
#' @seealso \code{\link{fcs2FitModel}} for the preferred way of calling this
#' function. \cr \code{inla} from the package \pkg{INLA} for details of
#' the \acronym{INLA} estimation method and the \code{"inla"} objects it
#' returns.
#' @keywords models regression nonlinear spatial
#' @export
fcs2RunINLA <-
function(fit, run = TRUE, verbose = TRUE)
{
    # attach required data (with corrected names)
    if (length(fit$rhoRW1Vars) > 0 || length(fit$rhoRW2Vars) > 0 || length(fit$muRW1Vars) > 0 || length(fit$muRW2Vars) > 0) {
        inlaRwBoundaries <- fit$rwBoundaries
        names(inlaRwBoundaries) <- make.names(names(inlaRwBoundaries))
        attach(inlaRwBoundaries, warn.conflicts=FALSE)
        on.exit(detach(inlaRwBoundaries), add=TRUE)
    }
    if (length(fit$rhoRW1Vars) > 0 || length(fit$rhoRW2Vars) > 0 || !is.null(fit$rhoSpatialVar) ||
        length(fit$muRW1Vars) > 0 || length(fit$muRW2Vars) > 0 || !is.null(fit$muSpatialVar)) {
        inlaPriorParameters <- fit$prior.parameters
        names(inlaPriorParameters) <- make.names(names(inlaPriorParameters))
        attach(inlaPriorParameters, warn.conflicts=FALSE)
        on.exit(detach(inlaPriorParameters), add=TRUE)
    }

    ## Fit presence/absence using prevalence (rho) terms:
    if (!is.na(pmatch(run, c(TRUE, "rho", "prevalence", "both")))) {

        # if using spatial term, create INLA graph file
        if(!is.null(fit$rhoSpatialVar)) {
            filename <- "inla_graph_file.txt"
            geobugs2inla(fit$rhoAdjacency$adj, fit$rhoAdjacency$num, filename)
            on.exit(unlink(filename))
        }

        # create model formula
        formulaString <- paste("Present ~ ",
                paste(collapse=' + ',
                      c("1",
                        if(length(fit$rhoLinearVars) > 0)
                            paste(make.names(fit$rhoLinearVars), collapse=" + ", sep=''),

                        if(length(fit$rhoRW1Vars) > 0)
                            paste("f(", make.names(fit$rhoRW1Vars), ", model='rw1', ",
                                  "values=", names(inlaRwBoundaries)[match(paste("gamma", fit$rhoRW1Vars, sep="."), names(fit$rwBoundaries))], ", ",
                                  "param=", names(inlaPriorParameters)[match(paste("phi", fit$rhoRW1Vars, sep="."), names(fit$prior.parameters))], ")",
                                  collapse=" + ", sep=''),

                        if(length(fit$rhoRW2Vars) > 0)
                            paste("f(", make.names(fit$rhoRW2Vars), ", model='rw2', ",
                                  "values=", names(inlaRwBoundaries)[match(paste("gamma", fit$rhoRW2Vars, sep="."), names(fit$rwBoundaries))], ", ",
                                  "param=", names(inlaPriorParameters)[match(paste("phi", fit$rhoRW2Vars, sep="."), names(fit$prior.parameters))], ")",
                                  collapse=" + ", sep=''),

                        if(!is.null(fit$rhoSpatialVar))
                            paste("f(", make.names(fit$rhoSpatialVar), ", model='besag', graph.file='", filename, "', ",
                                  "param=", names(inlaPriorParameters)[match(paste("phi", fit$rhoSpatialVar, sep="."), names(fit$prior.parameters))], ")",
                                  collapse=" + ", sep='')
                        )), sep='')
        fm <- formula(formulaString)


        # set data for INLA (with corrected names)
        inlaData <- data.frame(Present=as.numeric(fit$modelMatrix[, fit$allRunsTotalVar] > 0),
                               fit$modelMatrix[, c(fit$rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars, fit$rhoSpatialVar), drop=FALSE])
        if (!is.null(fit$rhoSpatialVar))
            inlaData[, fit$rhoSpatialVar] <- inlaData[, fit$rhoSpatialVar] - 1  # reduce spatial no by 1 as INLA is zero-based
        colnames(inlaData) <- make.names(colnames(inlaData))
        Ntrials <- rep(1, fit$N)

        # call INLA
        ### NOTE: does not specify prior parameters for linear variables
        if (verbose)
            cat("Fitting prevalence with INLA ...")
        flush.console()
        suppressWarnings({
            prevalenceFit <- inla(fm, "binomial", inlaData,
                                  Ntrials=Ntrials, control.compute=list(dic=TRUE))
        })

        # check that INLA ran ok
        if (is.null(prevalenceFit$summary.fitted.values)) {
            if (sum(inlaData$Present, na.rm=TRUE) == 0)
                warning("All observations are 0")
            stop("INLA failed to fit prevalence")
        }

        if (verbose)
            cat(" done\n")

    } else {
        # INLA not run for prevalence, but old run can be used for abundance, if present
        if (!is.null(fit$inlaFits) && !is.null(fit$inlaFits$rhoFit)) {
            prevalenceFit <- fit$inlaFits$rhoFit
            inlaData <- data.frame(Present=as.numeric(fit$modelMatrix[, fit$allRunsTotalVar] > 0))

        } else {
            if (!is.na(pmatch(run, c(TRUE, "mu", "abundance", "both"))))
                stop("cannot fit abundance with INLA without prevalence fit")
            else
                return(fit$inlaFits)
        }
    }


    ## Fit total catches estimated as present using abundance (mu) terms:
    if (!is.na(pmatch(run, c(TRUE, "mu", "abundance", "both")))) {
        subset <- !is.na(inlaData$Present) & (inlaData$Present | (prevalenceFit$summary.fitted.values[, 4] > 0.5))  # observed present or predicted prevalence > 0.5

        # check there are observations to fit
        if (sum(subset) == 0) {
            warning("Unable to fit abundance with INLA as no fish were observed or predicted to be present")
            return(list(rhoFit=prevalenceFit, muFit=fit$inlaFits$muFit))

        } else {

            # if using spatial term, create INLA graph file, if not made above
            if(!is.null(fit$muSpatialVar) && (is.null(fit$rhoSpatialVar) || fit$muSpatialVar != fit$rhoSpatialVar)) {
                filename <- "inla_graph_file.txt"  # Note: tempfile() didn't seem to work
                geobugs2inla(fit$muAdjacency$adj, fit$muAdjacency$num, filename)
                on.exit(unlink(filename), add=TRUE)
            }

            # create model formula
            formulaString <- paste("CatchTotal ~ ",
                    paste(collapse=' + ',
                          c("1",
                            if(length(fit$muLinearVars) > 0)
                                paste(make.names(fit$muLinearVars), collapse=" + ", sep=''),

                            if(length(fit$muRW1Vars) > 0)
                                paste("f(", make.names(fit$muRW1Vars), ", model='rw1', ",
                                      "values=", names(inlaRwBoundaries)[match(paste("beta", fit$muRW1Vars, sep="."), names(fit$rwBoundaries))], ", ",
                                      "param=", names(inlaPriorParameters)[match(paste("tau", fit$muRW1Vars, sep="."), names(fit$prior.parameters))], ")",
                                      collapse=" + ", sep=''),

                            if(length(fit$muRW2Vars) > 0)
                                paste("f(", make.names(fit$muRW2Vars), ", model='rw2', ",
                                      "values=", names(inlaRwBoundaries)[match(paste("beta", fit$muRW2Vars, sep="."), names(fit$rwBoundaries))], ", ",
                                      "param=", names(inlaPriorParameters)[match(paste("tau", fit$muRW2Vars, sep="."), names(fit$prior.parameters))], ")",
                                      collapse=" + ", sep=''),

                            if(!is.null(fit$muSpatialVar))
                                paste("f(", make.names(fit$muSpatialVar), ", model='besag', graph.file='", filename, "', ",
                                      "param=", names(inlaPriorParameters)[match(paste("tau", fit$muSpatialVar, sep="."), names(fit$prior.parameters))], ")",
                                      collapse=" + ", sep='')
                            )), sep='')
            fm <- formula(formulaString)

            # set data for INLA (with corrected names)
            inlaData <- data.frame(CatchTotal=fit$modelMatrix[subset, fit$allRunsTotalVar],
                                   fit$modelMatrix[subset, c(fit$muLinearVars, fit$muRW1Vars, fit$muRW2Vars, fit$muSpatialVar), drop=FALSE])
            if (!is.null(fit$muSpatialVar))
                inlaData[, fit$muSpatialVar] <- inlaData[, fit$muSpatialVar] - 1  # reduce spatial no by 1 as INLA is zero-based
            colnames(inlaData) <- make.names(colnames(inlaData))
            E <- fit$modelMatrix[subset, fit$surveyAreaVar]

            # call INLA
            ## NOTE: Priors not set for linear variables
            if (verbose)
                cat("Fitting abundance with INLA ...")
            flush.console()
            suppressWarnings({
                abundanceFit <- inla(fm, "nbinomial", inlaData, E=E, control.compute=list(dic=TRUE) ,
                                     control.data=list(prior="gaussian", param=fit$prior.parameters[["r"]]))
            })

            # check that INLA ran ok
            if (is.null(abundanceFit$summary.fitted.values))
                stop("INLA failed to fit abundance")

            if (verbose)
                cat(" done\n")

            # return both model fits
            list(rhoFit=prevalenceFit, muFit=abundanceFit)
        }

    } else {
        # INLA not run for abundance so return
        list(rhoFit=prevalenceFit, muFit=fit$inlaFits$muFit)
    }
}

