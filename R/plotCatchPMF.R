#' Plot Probability Mass Function of Predicted Total Catch
#'
#' Plots a probability mass function representing the predicted probability of
#' catching a given number of fish in total over all runs.  The plot is
#' coloured to indicate the mean \acronym{EQR} achieved if each number was
#' observed in terms of the \acronym{WFD} classification this would indicate.
#'
#'
#' @param fit an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#' fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS = TRUE}.
#' @param newData a data frame with surveys as rows and variables as columns.
#' It should contain all variables required by \code{fit}.
#' @param subset an vector specifying which surveys to plot the predicted
#' \acronym{PMF} for.
#' @param na.action a function which indicates what should happen when the data
#' contain missing values (\code{NA}s).  The default is set by the
#' \code{na.action} setting of \code{\link{options}} and this is usually set to
#' \code{\link{na.omit}}.  This setting removes surveys that contain missing
#' data in any required variables.  Alternatively, \code{\link{na.pass}} can be
#' used to ignore missing values (where possible) or \code{\link{na.fail}} can
#' be given to signal an error if missing values are found.
#' @param boundaries a vector of length 4 giving the \acronym{EQR} boundaries
#' separating the classes \emph{Bad}, \emph{Poor}, \emph{Good}, \emph{Moderate}
#' and \emph{High}.  These are used only to colour the plot with \emph{Bad} red
#' and \emph{High} blue.  If missing, regularly spaced boundaries of
#' \code{c(0.2, 0.4, 0.6, 0.8)} are used with a warning.  If \code{NULL}, the
#' probability that defines the single \acronym{EQR} is coloured blue.
#' @param mu a matrix of posterior samples of the abundance component
#' \eqn{\mu}{mu} can optionally be given to save recalculation if already
#' available.  This is assumed to have been calculated from
#' \code{\link{abundance}} using the same arguments as above.
#' @param rho a matrix of posterior samples of the prevalence component
#' \eqn{\rho}{rho} can optionally be given to save recalculation if already
#' available.  This is assumed to have been calculated from
#' \code{\link{prevalence}} using the same arguments as above.
#' @param xmax an optional x-axis upper limit giving the largest value of the
#' total catch to plot the probability of.
#' @param xmaxscale if \code{xmax} is not specified, the x-axis upper limit is
#' set to \code{xmaxscale *} the expected total catch.
#' @param maxpts the maximum number of points to plot on the x-axis. If
#' \code{xmax >= maxpts} a simpler plot with lines rather than bars is drawn
#' with \code{maxpts} points only.
#' @param title an optional character vector giving the title for each plot.
#' Should have one entry per survey selected by \code{subset}.
#' @param \dots further arguments passed to \code{\link{barplot}} if \code{xmax
#' < maxpts} of \code{\link{plot}} otherwise.
#' @seealso \code{\link{fcs2FitModel}} for fitting the model.\cr
#' \code{\link{pCatch}} for obtaining the probabilities plotted by this
#' function.\cr
#'
#' \code{\link{fcs2InteractivePrediction}} can be used to interactively view
#' how the predicted catch distribution changes with covariates.
#' @keywords hplot
#' @export
plotCatchPMF <-
function(fit, newData, subset = 1, na.action, boundaries, mu, rho, xmax, xmaxscale = 2, maxpts = 200, title, ...)
{
    # check 'newData' provided
    if (missing(newData))
        stop("a data frame 'newData' must be provided")

    if (inherits(subset, "logical") || min(subset) == 0)
        subset <- which(as.logical(subset))

    # get default na.action if missing
    if (missing(na.action)) {
        na.action <- getOption("na.action")
        if (is.null(na.action))
            na.action <- na.omit  # use 'na.omit' if option not set
        else
            na.action <- eval(parse(text=na.action))
    }

    # calculate or extract posterior samples of shape r, prevalence, abundance and q
    r <- fit$bugsFit$sims.list$r
    if (missing(mu))
        mu <- abundance(fit, newData, subset, na.action)
    if (missing(rho))
        rho <- prevalence(fit, newData, subset, na.action)
    if (fit$multiRun)
        q <- fit$bugsFit$sims.list$q

    # if multiple runs, calculate or extract number of runs
    if (fit$multiRun) {
        if (fit$nRunsVar %in% colnames(newData))
            nRuns <- newData[, fit$nRunsVar]
        else {
            # check run total vars available
            if (sum(fit$runTotalVars %in% colnames(newData)) < length(fit$runTotalVars))
                stop("Unable to calculate number of runs as nRunsVar and runTotalVars unavailable")

            nRuns <- rep(NA, nrow(newData))
            for (i in 1:nrow(newData))
                nRuns[i] <- sum(!is.na(newData[i, fit$runTotalVars]))
        }

        # check number of runs variable isn't larger than number of catch variables
        if (max(nRuns, na.rm=TRUE) > length(fit$runTotalVars)) {
            warning(paste("Number of runs ", if (fit$nRunsVar %in% colnames(newData))
                                                 paste("variable '", fit$nRunsVar, "' ", sep=""),
                                             "clipped to not exceed number of catch variables (", length(fit$runTotalVars), ")", sep=""))
            nRuns[!is.na(nRuns) & nRuns > length(fit$runTotalVars)] <- length(fit$runTotalVars)
        }
    }

    # extract (total) catch
    if (!fit$multiRun) {
        # use allRunsTotalVar if available
        if (fit$allRunsTotalVar %in% colnames(newData))
            catchMin <- newData[, fit$allRunsTotalVar]
        else
            catchMin <- rep(NA, nrow(newData))

    } else {
        catchMin <- rep(NA, nrow(newData))

        # from runs
        ## NOTE: can probably rewrite this to speed up by removing for loop
        if (!is.null(fit$runTotalVars) && sum(fit$runTotalVars %in% colnames(newData)) == length(fit$runTotalVars)) {
            for (i in 1:nrow(newData))
                if (!is.na(nRuns[i]))
                    catchMin[i] <- sum(as.numeric(newData[i, fit$runTotalVars[1:nRuns[i]]]))
        }

        # from all runs total
        if (!is.null(fit$allRunsTotalVar) && fit$allRunsTotalVar %in% colnames(newData) && sum(is.na(catchMin)) > 0)
            catchMin[is.na(catchMin)] <- newData[is.na(catchMin), fit$allRunsTotalVar]
    }
    catchMax <- catchMin
    # from all runs range
    if (!is.null(fit$allRunsRangeVars) && sum(fit$allRunsRangeVars %in% colnames(newData)) == 2 && sum(is.na(catchMin)) > 0) {
        catchMin[is.na(catchMin)] <- newData[is.na(catchMin), fit$allRunsRangeVars[1]]
        catchMax[is.na(catchMax)] <- newData[is.na(catchMax), fit$allRunsRangeVars[2]]
    }

    # find joint na.action attribute
    na.action.mu <- attr(mu, "na.action")
    na.action.rho <- attr(rho, "na.action")
    na.action.surveyArea <- attr(na.action(newData[subset, fit$surveyAreaVar]), "na.action")
    na.action <- union(union(na.action.mu, na.action.rho), na.action.surveyArea)
    if (!is.null(na.action))
        na.action <- sort(na.action)

    # set names and na.action class
    if (!is.null(na.action)) {
        names(na.action) <- rownames(newData)[subset][na.action]

        if (!is.null(na.action.mu))
            class(na.action) <- class(na.action.mu)
        else if (!is.null(na.action.rho))
            class(na.action) <- class(na.action.rho)
        else
            class(na.action) <- class(na.action.surveyArea)
    }

    # remove entries missing in other matrix
    if (length(setdiff(na.action, na.action.rho)) > 0) {
        which <- rep(TRUE, length(subset))
        which[na.action] <- FALSE
        if (length(na.action.rho) > 0)
            which <- which[-na.action.rho]
        rho <- rho[, which, drop=FALSE]
    }
    if (length(setdiff(na.action, na.action.mu)) > 0) {
        which <- rep(TRUE, length(subset))
        which[na.action] <- FALSE
        if (length(na.action.mu) > 0)
            which <- which[-na.action.mu]
        mu <- mu[, which, drop=FALSE]
    }

    # check that there are values to plot
    if (ncol(mu) == 0)
        stop("cannot calculate predictive distribution for selected surveys")

    # multiply abundance mu by survey area
    isubset <- subset[setdiff(1:length(subset), na.action)]
    mu <- mu * matrix(newData[isubset, fit$surveyAreaVar], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu))

    # if multiple runs, multiply again by (1 - (1 - q)^nRuns)
    if (fit$multiRun)
        mu <- mu * (1 - (1 - matrix(q, nrow=length(q), ncol=ncol(mu))) ^ matrix(nRuns[isubset], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu)))

    # calculate mean catch
    meanC <- apply(mu * rho, 2, mean)

    # calculate zeroprob = 1 - rho
    rho <- 1 - rho


    # if boundaries missing, use regularly spaced boundaries
    if (missing(boundaries)) {
        boundaries <- seq(0.2, 0.8, by=0.2)
        warning("using regularly spaced class boundaries to colour the plot as 'boundaries' argument not supplied")
    }
    if (missing(title))
        title <- paste("Catch PMF for", rownames(newData)[subset])
    if (length(title) == 1)
        title <- rep(title, length(subset))
    title <- title[setdiff(1:length(subset), na.action)]
    names(title) <- rownames(newData)[isubset]

    # set up multiple plots
    if (length(isubset) > 1) {
        w <- ceiling(sqrt(length(isubset)))
        par.old <- graphics::par(mfrow=c(w - (length(isubset) <= w*(w-1)), w))
        on.exit(graphics::par(par.old))
    }

    # plot pmf
    for (i in 1:length(isubset)) {
        # set x axis limit, xmaxi
        if (missing(xmax))
            xmaxi <- ceiling(xmaxscale * meanC[i])
        else
            xmaxi <- xmax

        # if catchMax missing, set to xmaxi
        if (is.na(catchMax[isubset[i]]))
            catchMax[isubset[i]] <- xmaxi

        # if catchMin missing, set to 0
        if (is.na(catchMin[isubset[i]]))
            catchMin[isubset[i]] <- 0

        # check xmaxi <= catchMax
        if (xmaxi < catchMax[isubset[i]])
            xmaxi <- catchMax[isubset[i]] + 1

        # set range of catch values
        if (xmaxi >= maxpts)
            xpts <- round(seq(0, xmaxi, length=maxpts))
        else
            xpts <- 0:xmaxi

        # calculate prob in each class
        mcProb <- array(dim=c(nrow(mu), length(xpts)))
        for (j in 1:length(xpts))
            mcProb[, j] <- dzinbinom(xpts[j], size=r, zeroprob=as.vector(rho[, i]), nbmean=as.vector(mu[, i]))

        # find mean prob and cumulative prob
        prob <- apply(mcProb, 2, mean)
        cumProb <- cumsum(prob)

        # calculate plot colour
        if (is.null(boundaries)) {
            # if boundaries = NULL, colour blue and grey
            col <- grDevices::hsv(h=0.66, s=(xpts <= catchMax[isubset[i]]) - 0.1 * (xpts < catchMin[isubset[i]]),
                       v=0.8 + 0.2 * (xpts <= catchMax[isubset[i]]) - 0.1 * (xpts < catchMin[isubset[i]]))

        } else {
            # colour by WFD classification
            hue <- rep(0, length(cumProb))
            hue[cumProb > boundaries[1] & cumProb <= boundaries[2]] <- 1/10
            hue[cumProb > boundaries[2] & cumProb <= boundaries[3]] <- 2/10
            hue[cumProb > boundaries[3] & cumProb <= boundaries[4]] <- 3/10
            hue[cumProb > boundaries[4]] <- 5/10
            col <- grDevices::hsv(h=hue, s=0.3 + 0.4 * (xpts <= catchMax[isubset[i]]) - 0.2 * (xpts < catchMin[isubset[i]]),
                       v=0.9 + 0.1 * (xpts <= catchMax[isubset[i]]) - 0.05 * (xpts < catchMin[isubset[i]]))
        }

        # plot
        if (xmaxi >= maxpts) {
            # line plot
            graphics::plot(xpts, prob, t='h', col=col, yaxs="i", bty="n", ylim=c(0, max(prob)), lwd=3,
                 xlab="Total catch", ylab="Probability", main=title[i], ...)
            graphics::abline(h=0, col="grey50")

        } else {
            # bar plot
            graphics::barplot(prob, names.arg=xpts, space=0, col=col,
                    border=grDevices::hsv(s=0, v=0.5 + 0.2 * (xpts > catchMax[isubset[i]])),
                    xlab="Total catch", ylab="Probability", main=title[i], ...)
            if (catchMax[isubset[i]] < length(prob))
                graphics::lines(rep(catchMax[isubset[i]] + 1, 2), c(0, prob[catchMax[isubset[i]] + 1]), col=grDevices::hsv(s=0, v=0.5)) # draw over boarder of data prob
            graphics::abline(h=0, col="grey50")
        }
    }
}

