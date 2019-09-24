#' Probabilities of the total catch distribution
#'
#' Density and distribution function of the predicted total catch over all
#' passes.
#'
#' @aliases dCatch pCatch Catch catch
#' @param fit an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#'   fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS =
#'   TRUE}.
#' @param x vector of (non-negative integer) quantiles.
#' @param newData a data frame with surveys as rows and variables as columns. It
#'   should contain all variables required by \code{fit}.
#' @param subset an optional vector specifying a subset of surveys to calculate
#'   the predicted total catch for.
#' @param na.action a function which indicates what should happen when the data
#'   contain missing values (\code{NA}s). The default is set by the
#'   \code{na.action} setting of \code{\link{options}} and this is usually set
#'   to \code{\link{na.omit}}. This setting removes surveys that contain missing
#'   data in any required variables. A vector indicating the rows that were
#'   removed can be extracted from the returned object using
#'   \code{\link{na.action}}. Alternatively, \code{\link{na.pass}} can be used
#'   to ignore missing values (where possible) or \code{\link{na.fail}} can be
#'   given to signal an error if missing values are found.
#' @param mu a matrix of posterior samples of the abundance component
#'   \eqn{\mu}{mu} can optionally be given to save recalculation if already
#'   available. This is assumed to have been calculated from
#'   \code{\link{abundance}} using the same arguments as above.
#' @param rho a matrix of posterior samples of the prevalence component
#'   \eqn{\rho}{rho} can optionally be given to save recalculation if already
#'   available. This is assumed to have been calculated from
#'   \code{\link{prevalence}} using the same arguments as above.
#' @return \code{dCatch} gives the probability mass function and \code{pCatch}
#'   the probability distribution function evaluated at values \code{x} for the
#'   selected surveys. In each case, a matrix is returned with rows for \code{x}
#'   values and columns for selected surveys.
#' @seealso \code{\link{fcs2FitModel}} for fitting the \acronym{FCS2} model.\cr
#'   \code{\link{predict.fcs2Fit}} for calculating the expected total catch.\cr
#'   \code{\link{plotCatchPMF}} for plotting these probabilities.\cr
#'   \code{\link{dzinbinom}} for the zero-inflated negative binomial
#'   distribution, which these functions use.
#' @keywords distribution
#' @export

dCatch <-
function(fit, x, newData, subset = 1:nrow(newData), na.action, mu, rho)
{
    # check 'newData' provided
    if (missing(newData))
        stop("a data frame 'newData' must be provided")

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

    # find joint na.action attribute
    na.action.mu <- attr(mu, "na.action")
    na.action.rho <- attr(rho, "na.action")
    na.action.surveyArea <- attr(na.action(newData[subset, fit$surveyAreaVar]), "na.action")
    na.action <- union(union(na.action.mu, na.action.rho), na.action.surveyArea)
    if (!is.null(na.action))
        na.action <- sort(na.action)

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

    # multiply abundance mu by survey area
    isubset <- subset[setdiff(1:length(subset), na.action)]
    mu <- mu * matrix(newData[isubset, fit$surveyAreaVar], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu))

    # if multiple runs, multiply again by (1 - (1 - q)^nRuns)
    if (fit$multiRun)
        mu <- mu * (1 - (1 - matrix(q, nrow=length(q), ncol=ncol(mu))) ^ matrix(nRuns[isubset], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu)))

    # calculate zeroprob = 1 - rho
    rho <- 1 - rho

    # calculate probabilities
    ret <- array(NA, dim=c(length(x), ncol(mu)), dimnames=list(x, colnames(mu)))
    for (i in 1:nrow(ret)) {
        for (j in 1:ncol(ret))
            ret[i, j] <- mean(dzinbinom(rep(x[i], nrow(mu)), size=r, zeroprob=as.vector(rho[, j]), nbmean=as.vector(mu[, j])))
    }

    # add na.action
    if (!is.null(na.action))
        attr(ret, "na.action") <- na.action

    # return
    ret
}


## pCatch
pCatch <-
function(fit, x, newData, subset = 1:nrow(newData), na.action, mu, rho)
{
    # check 'newData' provided
    if (missing(newData))
        stop("a data frame 'newData' must be provided")

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

    # find joint na.action attribute
    na.action.mu <- attr(mu, "na.action")
    na.action.rho <- attr(rho, "na.action")
    na.action.surveyArea <- attr(na.action(newData[subset, fit$surveyAreaVar]), "na.action")
    na.action <- union(union(na.action.mu, na.action.rho), na.action.surveyArea)
    if (!is.null(na.action))
        na.action <- sort(na.action)

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

    # multiply abundance mu by survey area
    isubset <- subset[setdiff(1:length(subset), na.action)]
    mu <- mu * matrix(newData[isubset, fit$surveyAreaVar], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu))

    # if multiple runs, multiply again by (1 - (1 - q)^nRuns)
    if (fit$multiRun)
        mu <- mu * (1 - (1 - matrix(q, nrow=length(q), ncol=ncol(mu))) ^ matrix(nRuns[isubset], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu)))

    # calculate zeroprob = 1 - rho
    rho <- 1 - rho

    # calculate probabilities
    ret <- array(NA, dim=c(length(x), ncol(mu)), dimnames=list(x, colnames(mu)))
    for (i in 1:nrow(ret)) {
        for (j in 1:ncol(ret))
            ret[i, j] <- mean(pzinbinom(rep(x[i], nrow(mu)), size=r, zeroprob=as.vector(rho[, j]), nbmean=as.vector(mu[, j])))
    }
    # add na.action
    if (!is.null(na.action))
        attr(ret, "na.action") <- na.action

    # return
    ret
}




