#' Combine BUGS fits
#'
#' Combines multiple BUGS fits to the same model into a single fit object. This
#' may be useful for combining multiple \acronym{MCMC} chains that were split
#' between computers to reduce computational time. \code{\link{thinBUGSSamples}}
#' is used when necessary with a warning to thin samples so that the number per
#' chain is consistent across fits.
#'
#' @param fit1 an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#'   fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS =
#'   TRUE}.
#' @param ... further \code{"fcs2Fit"} objects, each resulting from a full
#'   \acronym{FCS2} model fit with the same model and data but perhaps different
#'   initial values or random seed. The modified \code{"fcs2Fit"} object
#'   \code{fit} containing the combined sample.
#' @seealso \code{\link{thinBUGSSamples}} for thinning samples to reduce
#'   autocorrelation
#' @export

combineBUGSFits <-
function(fit1, ...)
{
    # create list of fits
    if (missing(fit1))
        fits <- list(...)
    else
        fits <- list(fit1, ...)

    # check each object is a bugs fit and extract number of samples per chain
    k <- length(fits)
    n.keep <- numeric(k)
    for (i in 1:k) {
        # check fcs2Fit object
        if (!inherits(fits[[i]], "fcs2Fit"))
            stop(paste("unrecognised argument", names(fits)[i]))

        # check fit contains bugs fit
        if (is.null(fits[[i]]$bugsFit))
            stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

        # store number of samples per chain
        n.keep[i] <- fits[[i]]$bugsFit$n.keep
    }

    # extract first fit
    fit1 <- fits[[1]]

    # return fit if only one provided
    if (k == 1)
        return(fit1)

    # check fits all have identical model setup and data (by testing against fit1)
    elements <- c("runTotalVars", "allRunsTotalVar", "allRunsRangeVars", "modelMatrix", "surveyAreaVar", "nRunsVar",
                  "rhoLinearVars", "rhoRW1Vars", "rhoRW2Vars", "rhoSpatialVar", "rhoAdjacency",
                  "rhoLinearVars", "rhoRW1Vars", "rhoRW2Vars", "rhoSpatialVar", "rhoAdjacency",
                  "dataType", "rwBoundaries", "prior.parameters", "multiRun")
    for (i in 2:k) {
        if (!identical(fit1[elements], fits[[i]][elements]))
            stop("fits do not have identical models and data")
    }

    # thin samples to match smallest no samples if necessary
    if (sd(n.keep) > 0) {
        # give warning
        warning("thinning samples so that the number per chain is consistent across fits")

        for (i in 1:k) {
            if (n.keep[i] > min(n.keep))
                fits[[i]] <- thinBUGSSamples(fits[[i]], which=round(seq(1, n.keep[i], length=min(n.keep))))
        }
    }

    # update ref to first fit
    fit1 <- fits[[1]]

    # combine sims arrays
    sims.array <- fit1$bugsFit$sims.array
    for (i in 2:k) {
        # extend current array
        sims.array <- sims.array[, c(1:ncol(sims.array), rep(NA, fits[[i]]$bugsFit$n.chains)), , drop=FALSE]

        # copy samples
        sims.array[, ncol(sims.array) - (fits[[i]]$bugsFit$n.chains - 1):0, ] <- fits[[i]]$bugsFit$sims.array
    }

    # calculate mean DIC
    DIC <- NULL
    for (i in 1:k)
        DIC <- c(DIC, fits[[i]]$bugsFit$DIC)
    if (is.null(DIC))
        DICOutput <- NULL
    else
        DICOutput <- rbind(c(fit1$bugsFit$isDIC, fit1$bugsFit$DICbyR, mean(DIC), fit1$bugsFit$pD))

    # calculate min n.iter and n.burnin
    n.iter <- n.burnin <- Inf
    for (i in 1:k) {
        n.iter <- min(n.iter, fits[[i]]$bugsFit$n.iter)
        n.burnin <- min(n.burnin, fits[[i]]$bugsFit$n.burnin)
    }

    # calculate n.thin
    n.thin <- (n.iter - n.burnin) / min(n.keep)  # not sure about this but doesn't really matter

    # recreate BUGS fit
    fit1$bugsFit <- as.bugs.array(sims.array, fit1$bugsFit$model.file, fit1$bugsFit$program, TRUE, DICOutput,
                                  n.iter, n.burnin, n.thin)

    # correct call
    fit1$call <- match.call()

    # return
    fit1
}

