#' Probability-Probability Plot of FCS2 Model Fit
#'
#' Produces a \dfn{Probability-Probability} (\acronym{P-P}) plot of the fitted
#' total catch that can be used to assess the model fit.
#'
#' If \code{addCI} is \code{TRUE}, a confidence interval is found by Monte
#' Carlo simulation.  \code{n.sims} datasets of fish catches are simulated from
#' the fitted model, each using the same observed covariates, and the
#' \acronym{P-P} plot line is found for each.  The confidence interval then
#' represents the extent covered by the central \code{ciprob} proportion of
#' these.
#'
#' If the observed data is typical of data simulated from the model, we would
#' expect the observed \acronym{P-P} plot line to fall within the confidence
#' interval \code{ciprob * 100}\% of the time.
#'
#' @param fit an \code{"fcs2Fit"} object, usually returned by
#' \code{\link{fcs2FitModel}}.
#' @param dataFrame a data frame with surveys as rows and variables as columns.
#' It should contain all variables required by \code{fit}.  This is usually the
#' same data frame used to create \code{fit}.
#' @param subset an optional vector specifying a subset of surveys to be used
#' to create the \acronym{P-P} plot.
#' @param na.action a function which indicates what should happen when the data
#' contain missing values (\code{NA}s).  The default is set by the
#' \code{na.action} setting of \code{\link{options}} and this is usually set to
#' \code{\link{na.omit}}.  This setting removes surveys that contain missing
#' data in any required variables.  A vector indicating the rows that were
#' removed can be extracted from the returned object using
#' \code{\link{na.action.fcs2Fit}}.  Alternatively, \code{\link{na.pass}} can
#' be used to ignore missing values (where possible) or \code{\link{na.fail}}
#' can be given to signal an error if missing values are found.
#' @param n.sims the number of datasets to simulate from the fitted model to
#' generate the confidence interval if \code{addCI} is \code{TRUE}.
#' @param ciprob the desired limits for the confidence interval added if
#' \code{addCI} is \code{TRUE}.  The default is \code{0.95} representing a 95\%
#' interval.
#' @param title an optional title for the plot.
#' @param addCI whether to add a confidence interval to the plot (default is
#' \code{TRUE}).
#' @param progressBar whether to show a progress bar when calculating the
#' confidence interval, since this can take some time.
#' @param mu a matrix of abundance samples can optionally be provided, as
#' generated from \code{\link{abundance}} using the same \code{fit},
#' \code{dataFrame} and \code{subset}.
#' @param rho a matrix of prevalence samples can optionally be provided, as
#' generated from \code{\link{prevalence}} using the same \code{fit},
#' \code{dataFrame} and \code{subset}.
#' @return If \code{addCI} is \code{TRUE}, the proportion of observed
#' \acronym{P-P} plot points within the confidence limits is printed to screen
#' and invisibly returned.
#' @section Warning: Generating the confidence interval can take some time,
#' especially if \code{n.sims} is large.
#' @seealso \code{\link{plot.fcs2Fit}}
#' @keywords hplot
#' @export
ppplot <-
function(fit, dataFrame, subset = 1:nrow(dataFrame), na.action, n.sims = 100, ciprob = 0.95, title = "", addCI = TRUE, progressBar = TRUE, mu, rho)
{
    if (class(subset) == "logical" || min(subset) == 0)
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
        mu <- abundance(fit, dataFrame, subset, na.action)
    if (missing(rho))
        rho <- prevalence(fit, dataFrame, subset, na.action)
    if (fit$multiRun)
        q <- fit$bugsFit$sims.list$q

    # if multiple runs, calculate or extract number of runs
    if (fit$multiRun) {
        if (fit$nRunsVar %in% colnames(dataFrame))
            nRuns <- dataFrame[, fit$nRunsVar]
        else {
            nRuns <- rep(NA, nrow(dataFrame))
            for (i in 1:nrow(dataFrame))
                nRuns[i] <- sum(!is.na(dataFrame[i, fit$runTotalVars]))
        }

        # check number of runs variable isn't larger than number of catch variables
        if (max(nRuns, na.rm=TRUE) > length(fit$runTotalVars)) {
            warning(paste("Number of runs ", if (fit$nRunsVar %in% colnames(dataFrame))
                                                 paste("variable '", fit$nRunsVar, "' ", sep=""),
                                             "clipped to not exceed number of catch variables (", length(fit$runTotalVars), ")", sep=""))
            nRuns[!is.na(nRuns) & nRuns > length(fit$runTotalVars)] <- length(fit$runTotalVars)
        }
    }

    ## calculate mean eqr
    eqr <- fcs2SingleEQR(fit, dataFrame, subset, na.action, mu=mu, rho=rho)
    na.action <- attr(eqr, "na.action")
    eqr <- mean(eqr)

    # remove entries missing in other matrix
    na.action.mu <- attr(mu, "na.action")
    na.action.rho <- attr(rho, "na.action")
    if (length(setdiff(na.action, na.action.rho)) > 0) {
        which <- rep(TRUE, length(subset))
        which[na.action] <- FALSE
        if (length(na.action.rho) > 0)
            which <- which[-na.action.rho]
        rho <- rho[, which]
    }
    if (length(setdiff(na.action, na.action.mu)) > 0) {
        which <- rep(TRUE, length(subset))
        which[na.action] <- FALSE
        if (length(na.action.mu) > 0)
            which <- which[-na.action.mu]
        mu <- mu[, which]
    }

    # multiply abundance mu by survey area
    isubset <- subset[setdiff(1:length(subset), na.action)]
    mu <- mu * matrix(dataFrame[isubset, fit$surveyAreaVar], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu))

    # if multiple runs, multiply again by (1 - (1 - q)^nRuns)
    if (fit$multiRun)
        mu <- mu * (1 - (1 - matrix(q, nrow=length(q), ncol=ncol(mu))) ^ matrix(nRuns[isubset], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu)))

    # calculate zeroprob = 1 - rho
    rho <- 1 - rho

    # sort eqr
    eqr <- sort(eqr)
    n <- length(eqr)
    prob <- ppoints(n)
    plot(0, 0, col="white", xlim=c(0, 1), ylim=c(0, 1),
         xlab="Theoretical Probabilities", ylab="Sample Probabilities", main=title)  #paste("P-P plot:", title))
    abline(0, 1, col="grey50")
    abline(v=c(0, 1), h=c(0, 1), col="grey80")
    lines(prob, eqr, col="blue")


    ## add CI
    if (addCI) {
        if (progressBar)
            bar <- txtProgressBar(style=3)

        # sample indicies corresponding to Monte Carlo sample
        index <- sample(1:nrow(mu), n.sims, replace=TRUE)

        # sample from fitted distribution
        sample <- matrix(rzinbinom(n.sims*n,
                                   r[rep(index, n)],
                                   zeroprob=as.vector(rho[cbind(rep(index, n), rep(1:n, rep(n.sims,n)))]),
                                   nbmean=as.vector(mu[cbind(rep(index, n), rep(1:n, rep(n.sims,n)))])), nrow = n.sims)

        # for each sample, calculate mean eqr values
        sampleEQRs <- array(NA, dim=c(n.sims, n))
        for (i in 1:n.sims) {
            sampleEQR <- matrix(pzinbinom(sample[i, rep(1:n, rep(nrow(mu), n))],
                                          r[rep(1:nrow(mu), n)],
                                          zeroprob=as.vector(rho[1:nrow(mu), ]),
                                          nbmean=as.vector(mu[1:nrow(mu), ])), nrow=nrow(mu))
            sampleEQRs[i, ] <- colMeans(sampleEQR)
            if (progressBar)
                setTxtProgressBar(bar, i / n.sims)
        }
        if (progressBar)
            close(bar)

        # sort each sample
        apply(sampleEQRs, 1, sort)
        for(i in 1:n.sims)
            sampleEQRs[i,] <- sort(sampleEQRs[i,])

        # find CI
        ci <- array(dim=c(2, n))
        cirange <- (1 - ciprob) / 2
        cirange <- c(cirange, 1 - cirange)
        for(j in 1:n)
            ci[, j] <- quantile(sampleEQRs[,j], cirange, name=FALSE)

        # add to plot
        lines(prob, ci[1,], lty=2)
        lines(prob, ci[2,], lty=2)
        lines(prob, eqr, col="blue") # make sure blue line on top

        # count proportion of data points outside CI
        prop <- sum( eqr >= ci[1,] & eqr <= ci[2,] ) / n
        cat("proportion within ", 100*ciprob, "% CI = ", prop, "\n", sep='')

        return(invisible(prop))
    }
}

