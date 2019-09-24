#' Thin BUGS Samples
#'
#' Used when necessary with a warning to thin samples so that the number per
#' chain is consistent across fits.
#'
#' @param fit fit
#' @param n.thin number to thin
#' @param first first
#' @param last last
#' @param which which
#' @param chains chains
#'
#' @return thinner samples
#' @export
#'
#' @examples
#' \dontrun{
#' thinBUGSSamples(fit,
#'                 n.thin = 1,
#'                 first = 1,
#'                 last = fit$bugsFit$n.keep,
#'                 which,
#'                 chains = 1:fit$bugsFit$n.chains)
#' }
thinBUGSSamples <-
function(fit, n.thin = 1, first = 1, last = fit$bugsFit$n.keep, which, chains = 1:fit$bugsFit$n.chains)
{
    # check fit contains bugs fit
    if (is.null(fit$bugsFit))
        stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

    # extract sims array
    sims.array <- fit$bugsFit$sims.array

    # thin sample
    if (missing(which))
        which <- seq(first, last, by=n.thin)
    else
        n.thin <- nrow(sims.array) / length(which)
    sims.array <- sims.array[which, chains, , drop=FALSE]  # keeping all chains for now
    n.thin <- fit$bugsFit$n.thin * n.thin

    # recreate BUGS fit
    if (is.null(fit$bugsFit$DIC))
        DICOutput <- NULL
    else
        DICOutput <- rbind(c(fit$bugsFit$isDIC, fit$bugsFit$DICbyR, fit$bugsFit$DIC, fit$bugsFit$pD))
    fit$bugsFit <- as.bugs.array(sims.array, fit$bugsFit$model.file, fit$bugsFit$program, TRUE, DICOutput,
                                 fit$bugsFit$n.iter, fit$bugsFit$n.burnin, n.thin)

    # correct call
    fit$call <- match.call()

    # return
    fit
}

