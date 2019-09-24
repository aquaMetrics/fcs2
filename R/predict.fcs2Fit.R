#' Predicted Total Catch
#'
#' Calculates the expected total catch over all runs for a range of surveys, as
#' predicted by a fitted \acronym{FCS2} model.
#'
#'
#' @param object an \code{"fcs2Fit"} object containing a full \acronym{FCS2}
#' model object, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS
#' = TRUE}.
#' @param newData a data frame with surveys as rows and variables as columns.
#' It should contain all variables required by \code{object}.
#' @param subset an optional vector specifying a subset of surveys to calculate
#' the predicted total catch for.
#' @param na.action a function which indicates what should happen when the data
#' contain missing values (\code{NA}s).  The default is set by the
#' \code{na.action} setting of \code{\link{options}} and this is usually set to
#' \code{\link{na.omit}}.  This setting removes surveys that contain missing
#' data in any required variables.  A vector indicating the rows that were
#' removed can be extracted from the returned object using
#' \code{\link{na.action}}.  Alternatively, \code{\link{na.pass}} can be used
#' to ignore missing values (where possible) or \code{\link{na.fail}} can be
#' given to signal an error if missing values are found.
#' @param mu a matrix of posterior samples of the abundance component \eqn{\mu}
#' can optionally be given to save recalculation if already available.  This is
#' assumed to have been calculated from \code{\link{abundance}} using the same
#' arguments as above.
#' @param rho a matrix of posterior samples of the prevalence component
#' \eqn{\rho} can optionally be given to save recalculation if already
#' available.  This is assumed to have been calculated from
#' \code{\link{prevalence}} using the same arguments as above.
#' @return a vector containing the expected total catch for each survey
#' selected.
#' @seealso \code{\link{fcs2FitModel}} for fitting the \acronym{FCS2} model.\cr
#' \code{\link{pCatch}} for calculating probabilities relating to the predicted
#' total catch.
#' @export
predict.fcs2Fit <-
function(object, newData, subset = 1:nrow(newData), na.action, mu, rho)
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

    # calculate or extract posterior samples of shape r, prevalence and abundance
    r <- object$bugsFit$sims.list$r
    if (missing(mu))
        mu <- abundance(object, newData, subset, na.action)
    if (missing(rho))
        rho <- prevalence(object, newData, subset, na.action)
    if (object$multiRun)
        q <- object$bugsFit$sims.list$q

    # if multiple runs, calculate or extract number of runs
    if (object$multiRun) {
        if (object$nRunsVar %in% colnames(newData))
            nRuns <- newData[, object$nRunsVar]
        else {
            nRuns <- rep(NA, nrow(newData))
            for (i in 1:nrow(newData))
                nRuns[i] <- sum(!is.na(newData[i, object$runTotalVars]))
        }

        # check number of runs variable isn't larger than number of catch variables
        if (max(nRuns, na.rm=TRUE) > length(object$runTotalVars)) {
            warning(paste("Number of runs ", if (object$nRunsVar %in% colnames(newData))
                                                 paste("variable '", object$nRunsVar, "' ", sep=""),
                                             "clipped to not exceed number of catch variables (", length(object$runTotalVars), ")", sep=""))
            nRuns[!is.na(nRuns) & nRuns > length(object$runTotalVars)] <- length(object$runTotalVars)
        }
    }

    # find joint na.action attribute
    na.action.mu <- attr(mu, "na.action")
    na.action.rho <- attr(rho, "na.action")
    na.action.surveyArea <- attr(na.action(newData[subset, object$surveyAreaVar]), "na.action")
    na.action <- union(union(na.action.mu, na.action.rho), na.action.surveyArea)

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

    # multiply abundance mu by survey area
    isubset <- subset[setdiff(1:length(subset), na.action)]
    mu <- mu * matrix(newData[isubset, object$surveyAreaVar], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu))

    # if multiple runs, multiply again by (1 - (1 - q)^nRuns)
    if (object$multiRun)
        mu <- mu * (1 - (1 - matrix(q, nrow=length(q), ncol=ncol(mu))) ^ matrix(nRuns[isubset], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu)))

    # find expected catch by multiplying mu by rho
    ret <- apply(mu * rho, 2, mean)
    names(ret) <- colnames(mu)

    # add na.action attribute
    attr(ret, "na.action") <- na.action

    ret
}

