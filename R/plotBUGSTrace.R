#' Trace Plots of MCMC Chains from BUGS
#'
#' Produces trace plots of the \acronym{MCMC} samples for each parameter in the
#' fitted \acronym{FCS2} model.  These can be used to visualise the convergence
#' of the \acronym{MCMC} chains and the autocorrelation.
#'
#' If the \acronym{MCMC} chains have converged, the multiple chains (in
#' different colours) should agree; Their values should not be identical but
#' they should cover the same range of values as each other.
#'
#' The \acronym{MCMC} chains should also mix well meaning that they should move
#' about the parameter space quickly.  A chain that moves slowly across the
#' parameter space has high autocorrelation.  This can also be verified by
#' plotting the autocorrelation of a parameter using \code{\link{acf}}.  The
#' autocorrelation can be reduced by thinning the sample using
#' \code{\link{thinBUGSSamples}} or by specifying a larger number for
#' \code{n.thin} when fitting.
#'
#' @param fit an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#' fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS = TRUE}.
#' @param variables an optional character vector giving the names of the model
#' variables to plot.  \code{fcs2:::variable.names.fcs2Fit} can be used for
#' this, but only singular variable names can be given.
#' @param chains which \acronym{MCMC} chains to plot. Default to all.
#' @param nPerPage the number of trace plots to create on each page.
#' @param col a vector giving the colour to use for each chain in
#' \code{chains}.
#' @param ylab,ylim,main,...  further arguments passed to \code{\link{plot}}
#' when creating each plot.
#' @seealso \code{\link{thinBUGSSamples}}
#' @keywords hplot
#' @export
plotBUGSTrace <-
function(fit, variables = variable.names(fit, rw="singular", spatial="singular", hyperpar="prec"),
         chains = 1:fit$bugsFit$n.chains, nPerPage = 12, col = chains + 1, ylab, ylim, main, ...)
{
    # check fit contains bugs fit
    if (is.null(fit$bugsFit))
        stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

    # set up multiple plots
    if (length(variables) > 1) {
        if (length(variables) <= nPerPage) {
            w <- ceiling(sqrt(length(variables)))
            par.old <- par(mfrow=c(w, w - (length(variables) <= w*(w-1))))

        } else {
            w <- ceiling(sqrt(nPerPage))
            par.old <- par(mfrow=c(w, w - (nPerPage <= w*(w-1))), ask=TRUE)
        }
        on.exit(par(par.old))
    }

    # set up colours
    if (length(col) < length(chains))
        col <- rep(col, length(chains))

    # plot trace for each variable
    for (var in variables) {
        # set plot options with defaults or given values
        ylab.var <- ifelse(missing(ylab), var, ylab)
        main.var <- ifelse(missing(main), var, main)
        if (missing(ylim))
            ylim.var <- range(fit$bugsFit$sims.array[, chains, var])
        else
            ylim.var <- ylim

        plot(fit$bugsFit$sims.array[, chains[1], var], ylim=ylim.var, type='l', col=col[1], ylab=ylab.var, main=main.var, ...)
        if (length(chains) > 1) {
            for (i in 2:length(chains))
                lines(fit$bugsFit$sims.array[, chains[i], var], col=col[i])
        }
    }
}

