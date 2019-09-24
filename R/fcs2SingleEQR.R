#' Single EQR
#'
#' Produces Monte Carlo samples of the \acronym{FCS2} \dfn{Ecological Quality
#' Ratio} (\acronym{EQR}) for a single species and for each survey.  The
#' \acronym{EQR} is a number between 0 and 1 which is found by comparing the
#' observed catch with the model's prediction at reference conditions.  Higher
#' \acronym{EQR} values are better and these are caused by catching larger
#' numbers of fish.  The variability in the Monte Carlo samples of the
#' \acronym{EQR} is due to the uncertainty in the parameters of the
#' \acronym{FCS2} model.  This uncertainty can be used by
#' \code{\link{fcs2Classify}} to produce probabilistic \acronym{WFD}
#' classifications from each \acronym{EQR}.
#'
#'
#' @param fit an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#' fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS = TRUE}.
#' @param newData a data frame with surveys as rows and variables as columns.
#' It should contain all variables required by \code{fit}.  Covariates which
#' are related to human disturbance (\dfn{pressure variables}) should have
#' their values set to the value expected at the site if it were undisturbed
#' (\dfn{reference conditions}) rather than the observed value for each of
#' these variables.
#' @param subset an optional vector specifying a subset of surveys to calculate
#' \acronym{EQR} samples for.
#' @param na.action a function which indicates what should happen when the data
#' contain missing values (\code{NA}s).  The default is set by the
#' \code{na.action} setting of \code{\link{options}} and this is usually set to
#' \code{\link{na.omit}}.  This setting removes surveys that contain missing
#' data in any required variables.  A vector indicating the rows that were
#' removed can be extracted from the returned object using
#' \code{\link{na.action.fcs2EQR}}.  Alternatively, \code{\link{na.pass}} can
#' be used to ignore missing values (where possible) or \code{\link{na.fail}}
#' can be given to signal an error if missing values are found.
#' @param mu a matrix of posterior samples of the abundance component \eqn{\mu}
#' can optionally be given to save recalculation if already available.  This is
#' assumed to have been calculated from \code{\link{abundance}} using the same
#' arguments as above.
#' @param rho a matrix of posterior samples of the prevalence component
#' \eqn{\rho} can optionally be given to save recalculation if already
#' available.  This is assumed to have been calculated from
#' \code{\link{prevalence}} using the same arguments as above.
#' @return Returns an \code{"fcs2EQR"} object that contains the \acronym{EQR}
#' samples.  The \code{"fcs2EQR"} object is essentially a matrix with Monte
#' Carlo samples as rows and surveys as columns.
#' @seealso \code{\link{print.fcs2EQR}}, \code{\link{summary.fcs2EQR}} and
#' \code{\link{fcs2EQRSummaryMatrix}} for summarising \code{"fcs2EQR"}
#' objects;\cr \code{\link{plot.fcs2EQR}} for plotting \acronym{EQR}
#' variables;\cr \code{\link{mean.fcs2EQR}} and \code{\link{quantile.fcs2EQR}}
#' for calculating means and quantiles of \acronym{EQR} variables
#' respectively;\cr
#'
#' \code{\link{fcs2FitModel}} for producing the required \acronym{FCS2} model
#' fit;\cr \code{\link{fcs2Classify}} for using \acronym{EQR} samples to
#' produce probabilistic \acronym{WFD} classifications.
#'
#' \code{\link{fcs2JointEQR}} for producing \acronym{EQR} samples that combine
#' multiple species and/or surveys;\cr \code{\link{fcs2JointAndSingleEQR}} for
#' calculating joint and single \acronym{EQR} samples simultaneously.
#' @export
#' @examples
#'
#' \dontrun{
#'
#' ### Very simple example with no covariates
#' ###
#'
#' # simulate random dataset
#' Data <- data.frame(SurveyArea=rlnorm(100, 4.6, 0.5))   # random survey area
#' Data$Catch <- rzinbinom(100, size=1.1, zeroprob=0.3,
#'                         nbmean=0.3 * Data$SurveyArea)  # single catch per survey
#'
#' # fit full model with OpenBUGS
#' fit <- fcs2FitModel("Catch", dataFrame=Data, surveyAreaVar="SurveyArea",
#'                     runBUGS=TRUE, n.iter=1000, bugsProgram="OpenBUGS")
#'
#' # calculate samples of single EQR, using same dataset
#' eqr <- fcs2SingleEQR(fit, Data)
#'
#' # plot EQR variables for first 9 surveys
#' plot(eqr, 1:9, boundaries=NULL)
#'
#' # calculate mean EQR values
#' mean(eqr)
#' }
#'
fcs2SingleEQR <-
function(fit, newData, subset = 1:nrow(newData), na.action, mu, rho)
{
    # check fit contains bugs fit
    if (is.null(fit$bugsFit))
        stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

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

    # extract (total) catch
    if (!fit$multiRun) {
        # use allRunsTotalVar if available
        if (fit$allRunsTotalVar %in% colnames(newData))
            catch <- newData[, fit$allRunsTotalVar]
        else
            catch <- rep(NA, nrow(newData))

    } else {
        catch <- rep(NA, nrow(newData))

        # from runs
        if (!is.null(fit$runTotalVars) && sum(fit$runTotalVars %in% colnames(newData)) == length(fit$runTotalVars)) {
            for (i in 1:nrow(newData))
                if (!is.na(nRuns[i]))
                    catch[i] <- sum(as.numeric(newData[i, fit$runTotalVars[1:nRuns[i]]]))
        }

        # from all runs total
        if (!is.null(fit$allRunsTotalVar) && fit$allRunsTotalVar %in% names(newData) && sum(is.na(catch)) > 0)
            catch[is.na(catch)] <- newData[is.na(catch), fit$allRunsTotalVar]
    }
    # from all runs range
    if (!is.null(fit$allRunsRangeVars) && sum(fit$allRunsRangeVars %in% colnames(newData)) == 2 && sum(is.na(catch)) > 0) {
        # calculate whether min < max
        irange <- is.na(catch) & !is.na(newData[, fit$allRunsRangeVars[1]]) &
                    newData[, fit$allRunsRangeVars[1]] < newData[, fit$allRunsRangeVars[2]]

        # add min to catch
        catch[is.na(catch)] <- newData[is.na(catch), fit$allRunsRangeVars[1]]

    } else
        irange <- NULL

    # find joint na.action attribute
    na.action.mu <- attr(mu, "na.action")
    na.action.rho <- attr(rho, "na.action")
    na.action.catch <- attr(na.action(catch[subset]), "na.action")
    na.action.surveyArea <- attr(na.action(newData[subset, fit$surveyAreaVar]), "na.action")
    na.action <- union(union(union(na.action.mu, na.action.rho), na.action.catch), na.action.surveyArea)
    if (!is.null(na.action))
        na.action <- sort(na.action)

    # set names and na.action class
    if (!is.null(na.action)) {
        names(na.action) <- rownames(newData)[subset][na.action]

        if (!is.null(na.action.mu))
            class(na.action) <- class(na.action.mu)
        else if (!is.null(na.action.rho))
            class(na.action) <- class(na.action.rho)
        else if (!is.null(na.action.catch))
            class(na.action) <- class(na.action.catch)
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
    mu <- mu * matrix(newData[isubset, fit$surveyAreaVar], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu))

    # if multiple runs, multiply again by (1 - (1 - q)^nRuns)
    if (fit$multiRun)
        mu <- mu * (1 - (1 - matrix(q, nrow=length(q), ncol=ncol(mu))) ^ matrix(nRuns[isubset], byrow=TRUE, nrow=nrow(mu), ncol=ncol(mu)))

    # calculate zeroprob = 1 - rho
    rho <- 1 - rho

    # create indices for MC samples and catch
    ir <- rep(1:fit$bugsFit$n.sims, ncol(mu))
    im <- 1:fit$bugsFit$n.sims
    ic <- rep(isubset, rep(fit$bugsFit$n.sims, ncol(mu)))

    # if range data, sample all runs total from posterior restricted to range and modify catch index
    if (sum(irange[isubset]) > 0) {
        # create index of indices with range data
        iix <- which(ic %in% isubset[irange[isubset]])

        # sample n.sims all runs totals for each survey in irange
        catch <- c(catch, rzinbinom_constrained(sum(irange[isubset]) * fit$bugsFit$n.sims,
                                                r[ir[iix]], zeroprob=as.vector(rho[im, irange[isubset]]), nbmean=as.vector(mu[im, irange[isubset]]),
                                                min=catch[ic[iix]], max=newData[ic[iix], fit$allRunsRangeVars[2]]))

        # replace NaNs by Inf (so that gives EQR 1, though I think these are caused by zeroprob=1 so this always happens)
        catch[is.nan(catch)] <- Inf

        # update catch index
        ic[iix] <- nrow(newData) + 1:(sum(irange[isubset]) * fit$bugsFit$n.sims)
    }

    # calculate probability of obtaining less or equal the observed catch
    ret <- pzinbinom(catch[ic], r[ir], zeroprob=as.vector(rho[im, ]), nbmean=as.vector(mu[im, ]))
    dim(ret) <- c(length(im), ncol(mu))
    dimnames(ret) <- list(NULL, survey=colnames(mu))

    # add na.action attribute
    attr(ret, "na.action") <- na.action

    # set class to fcs2EQR
    class(ret) <- "fcs2EQR"

    ret
}

