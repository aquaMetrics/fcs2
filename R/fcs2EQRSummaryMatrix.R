#' EQR Summary Matrix
#'
#' Produces a matrix summarising the single and joint \acronym{EQR}s and
#' optionally comparing these to the observed catches and the \acronym{FCS2}
#' model's predictions.
#'
#' @param fit1 an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#'   fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS =
#'   TRUE}.
#' @param ... further \code{"fcs2Fit"} objects, each resulting from a full
#'   \acronym{FCS2} model fit with a different species.
#' @param newData a data frame with surveys as rows and variables as columns. It
#'   should contain all variables required by each of the model fits. Covariates
#'   which are related to human disturbance (\dfn{pressure variables}) should
#'   have their values set to the value expected at the site if it were
#'   undisturbed (\dfn{reference conditions}) rather than the observed value for
#'   each of these variables.
#' @param joinByVar the name of a column in \code{dataFrame} to be used for
#'   joining multiple surveys together when calculating the \acronym{EQR}. For
#'   example, \code{joinByVar = "WBId"} would produce \acronym{EQR} samples for
#'   each water body, with every survey in a water body contributing towards its
#'   \acronym{EQR}.\cr If \code{NULL} (the default), \acronym{EQR}s are
#'   calculated for each survey only.
#' @param subset an optional vector specifying a subset of surveys to calculate
#'   \acronym{EQR} samples for.
#' @param na.action a function which indicates what should happen when the data
#'   contain missing values (\code{NA}s). The default is set by the
#'   \code{na.action} setting of \code{\link{options}} and this is usually set
#'   to \code{\link{na.omit}}. This setting removes surveys that contain missing
#'   data in any required variables. Alternatively, \code{\link{na.pass}} can be
#'   used to ignore missing values (where possible) or \code{\link{na.fail}} can
#'   be given to signal an error if missing values are found.
#' @param classify whether to produce probabilistic classifications of each
#'   \acronym{WFD} class.  If \code{joinByVar = NULL} (the default),
#'   classifications are made using the joint \acronym{EQR} for each survey.  If
#'   \code{joinByVar} is provided, classifications are made using the combined
#'   \acronym{EQR}s across multiple surveys.
#' @param boundaries a vector of length 4 giving the \acronym{EQR} boundaries
#'   separating the classes \emph{Bad}, \emph{Poor}, \emph{Good},
#'   \emph{Moderate} and \emph{High}. This is used if \code{classify = TRUE} to
#'   produce probabilistic classifications of each class. If missing, regularly
#'   spaced boundaries of \code{c(0.2, 0.4, 0.6, 0.8)} are used with a warning.
#' @param summary a character vector listing how to summarise the \acronym{EQR}
#'   variables. Can contain \code{"mean"} for the mean and/or \code{"sd"} to
#'   calculate the standard deviation.
#' @param observations whether to give the observed total catch for each survey
#'   and species.
#' @param predictions whether to summarise the model's predicted total catch for
#'   each survey and species. If \code{TRUE} (the default), the expected total
#'   catch is given but if \code{"detail"}, the probability of presence
#'   (prevalence) and the expected total catch if present are alternatively
#'   given. If \code{"all"}, all three columns are shown for each species.
#' @param n.samples the number of Monte Carlo \acronym{EQR} samples to produce
#'   for each survey (or joining variable).
#' @param n.sims the number of Monte Carlo simulations to make for each
#'   \acronym{EQR} sample. These internal samples are used for approximating the
#'   probability that defines the joint \acronym{EQR}.
#' @param eqrs the result of \code{\link{fcs2JointAndSingleEQR}} with the above
#'   arguments can alternatively be given, to save recalculation if available.
#' @return a matrix with surveys as rows and columns containing the selected
#'   information.
#' @seealso \code{\link{fcs2JointAndSingleEQR}}
#' @export


fcs2EQRSummaryMatrix <-
function(fit1, ..., newData, joinByVar = NULL, subset = 1:nrow(newData), na.action,
         classify = !missing(boundaries), boundaries, summary = "mean", observations = TRUE,
         predictions = TRUE, n.samples = 1000, n.sims = 1000, eqrs)
{
    # get default na.action if missing
    if (missing(na.action)) {
        na.action <- getOption("na.action")
        if (is.null(na.action))
            na.action <- na.omit  # use 'na.omit' if option not set
        else
            na.action <- eval(parse(text=na.action))
    }

    # calculate joint and single EQRs, unless provided
    if (missing(eqrs))
        eqrs <- fcs2JointAndSingleEQR(fit1, ..., newData=newData, joinByVar=joinByVar, subset=subset, na.action=na.action,
                                      n.samples=n.samples, n.sims=n.sims, both=TRUE)

    # check whether two sets of eqrs
    if (class(eqrs) == "list") {
        # EQRs by survey and joined, so summarise and combine

        # extract names of fits from EQR
        fitNames <- dimnames(eqrs[[1]])[[3]]

        # check we have joinByVar (as eqrs may have been provided)
        if (missing(joinByVar))
            joinByVar <- sub("EQRBy", "", names(eqrs)[2])

        # create initial survey and join data frames with joinByVar
        ret <- newData[colnames(eqrs[[1]]), joinByVar, drop=FALSE]  # by survey
        joinMx <- data.frame(colnames(eqrs[[2]]))  # by join variable
        colnames(joinMx) <- joinByVar

        # calculate classifications from combined EQRs for all species
        if (classify)
            joinMx <- cbind(joinMx, t(fcs2Classify(eqrs[[2]][,,1], 1:ncol(eqrs[[2]]), boundaries=boundaries)))

        # calculate summaries of eqrs
        if ("mean" %in% summary) {
            ret <- cbind(ret, apply(eqrs[[1]], 2:3, mean))
            joinMx <- cbind(joinMx, apply(eqrs[[2]], 2:3, mean))
            colnames(ret)[ncol(ret) - (dim(eqrs[[1]])[3] - 1):0] <- paste(colnames(ret)[ncol(ret) - (dim(eqrs[[1]])[3] - 1):0], "survey EQR mean")
            colnames(joinMx)[ncol(joinMx) - (dim(eqrs[[2]])[3] - 1):0] <- paste(colnames(joinMx)[ncol(joinMx) - (dim(eqrs[[2]])[3] - 1):0], joinByVar, "EQR mean")
        }
        if ("sd" %in% summary) {
            ret <- cbind(ret, apply(eqrs[[1]], 2:3, sd))
            joinMx <- cbind(joinMx, apply(eqrs[[2]], 2:3, sd))
            colnames(ret)[ncol(ret) - (dim(eqrs[[1]])[3] - 1):0] <- paste(colnames(ret)[ncol(ret) - (dim(eqrs[[1]])[3] - 1):0], "survey EQR sd")
            colnames(joinMx)[ncol(joinMx) - (dim(eqrs[[2]])[3] - 1):0] <- paste(colnames(joinMx)[ncol(joinMx) - (dim(eqrs[[2]])[3] - 1):0], joinByVar, "EQR sd")
        }

        # merge survey and join matrices, but keep rownames of ret (so that can match with rows in newData)
        ret <- cbind(ret, rowNames=rownames(ret))  # add rownames as a column
        ret <- merge(joinMx, ret)  # join with joinByVar first
        rownames(ret) <- ret$rowNames
        ret <- ret[, -match("rowNames", colnames(ret))]  # remove rownames column

    } else {
        # summarise EQRs by survey only

        # check size of eqr
        if (length(dim(eqrs)) == 2)
            dim(eqrs) <- c(dim(eqrs), 1)

        # extract names of fits from EQR
        fitNames <- dimnames(eqrs)[[3]]
        if (is.null(fitNames))
            fitNames <- ""

        # create initial survey data frame (with no columns)
        ret <- as.data.frame(array(dim=c(ncol(eqrs), 0), dimnames=list(colnames(eqrs), NULL)))

        # calculate classifications from all species EQRs
        if (classify)
            ret <- cbind(ret, t(fcs2Classify(eqrs[,,1], 1:ncol(eqrs), boundaries=boundaries)))

        # calculate summaries of eqrs
        if ("mean" %in% summary) {
            ret <- cbind(ret, apply(eqrs, 2:3, mean))
            if (dim(eqrs)[3] == 1)
                colnames(ret)[ncol(ret)] <- paste("Survey EQR mean")
            else
                colnames(ret)[ncol(ret) - (dim(eqrs)[3] - 1):0] <- paste(colnames(ret)[ncol(ret) - (dim(eqrs)[3] - 1):0], "survey EQR mean")
        }
        if ("sd" %in% summary) {
            ret <- cbind(ret, apply(eqrs, 2:3, sd))
            if (dim(eqrs)[3] == 1)
                colnames(ret)[ncol(ret)] <- paste("Survey EQR sd")
            else
                colnames(ret)[ncol(ret) - (dim(eqrs)[3] - 1):0] <- paste(colnames(ret)[ncol(ret) - (dim(eqrs)[3] - 1):0], "survey EQR sd")
        }
    }

    # add observations for each survey / species
    if (observations) {
        ## create list of fits, named by their (first) catch variable
        if (missing(fit1))
            fits <- list(...)
        else
            fits <- list(fit1, ...)

        # check each non-matched object is a fit and extract fit names
        k <- length(fits)
        if (k == 0)
            stop("fits must be provided to summarise observations")
        for (i in 1:k) {
            if (class(fits[[i]]) != "fcs2Fit")
                stop(paste("unrecognised argument", names(fits)[i]))

            if (is.null(names(fits)[i]) || is.na(names(fits)[i]) || names(fits)[i] == "") {
                if (!is.null(fits[[i]]$runTotalVars))
                    names(fits)[i] <- sub(".Run1Total", "", fits[[i]]$runTotalVars[1], fixed=TRUE)
                else if (!is.null(fits[[i]]$allRunsRangeVars))
                    names(fits)[i] <- sub(".AllRunsTotalMin", "", fits[[i]]$allRunsRangeVars[1], fixed=TRUE)
                else
                    names(fits)[i] <- sub(".AllRunsTotal", "", fits[[i]]$allRunsTotalVar[1], fixed=TRUE)
            }
        }

        # add total catch for each survey / species
        for (i in 1:k) {
            ret <- cbind(ret, .allRunsTotal(newData, fits[[i]]$runTotalVars, fits[[i]]$allRunsTotalVar, fits[[i]]$allRunsRangeVars,
                                            fits[[i]]$nRunsVar, subset=rownames(ret), allowRange=TRUE))
            colnames(ret)[ncol(ret)] <- paste(names(fits)[i], "observed total catch")
        }
    }

    # add predictions for each survey / species
    if (predictions != FALSE) {
        if (!observations) {
            ## create list of fits, named by their (first) catch variable
            if (missing(fit1))
                fits <- list(...)
            else
                fits <- list(fit1, ...)

            # check each non-matched object is a fit and extract fit names
            k <- length(fits)
            if (k == 0)
                stop("fits must be provided to summarise predictions")
            for (i in 1:k) {
                if (class(fits[[i]]) != "fcs2Fit") {
                    obName <- names(fits)[i]
                    stop(paste("unrecognised argument", obName))
                }

                if (is.null(names(fits)[i]) || names(fits)[i] == "") {
                    if (!is.null(fits[[i]]$runTotalVars))
                        names(fits)[i] <- sub(".Run1Total", "", fits[[i]]$runTotalVars[1], fixed=TRUE)
                    else if (!is.null(fits[[i]]$allRunsRangeVars))
                        names(fits)[i] <- sub(".AllRunsTotalMin", "", fits[[i]]$allRunsRangeVars[1], fixed=TRUE)
                    else
                        names(fits)[i] <- sub(".AllRunsTotal", "", fits[[i]]$allRunsTotalVar[1], fixed=TRUE)
                }
            }
        }

        # calculate expected catch
        for (i in 1:k) {
            # calculate abundance and prevalence
            mean <- abundance(fits[[i]], newData, rownames(ret), na.pass)
            probPresent <- prevalence(fits[[i]], newData, rownames(ret), na.pass)

            # multiply by survey area
            mean <- mean * matrix(newData[rownames(ret), fits[[i]]$surveyAreaVar], byrow=TRUE, nrow=nrow(mean), ncol=ncol(mean))

            # if multiple runs, multiply again by (1 - (1 - q)^nRuns)
            if (fits[[i]]$multiRun) {
                # extract number of runs
                if (fits[[i]]$nRunsVar %in% colnames(newData))
                    nRuns <- newData[rownames(ret), fits[[i]]$nRunsVar]
                else {
                    nRuns <- rep(NA, nrow(ret))
                    for (j in 1:nrow(ret))
                        nRuns[j] <- sum(!is.na(newData[rownames(ret)[j], fits[[i]]$runTotalVars]))
                }
                mean <- mean * (1 - (1 - matrix(fits[[i]]$bugsFit$sims.list$q, nrow=nrow(mean), ncol=ncol(mean))) ^
                                         matrix(nRuns, byrow=TRUE, nrow=nrow(mean), ncol=ncol(mean)))
            }

            # calculate mean of mean samples
            expProbPresent <- apply(probPresent, 2, mean)
            expTotalIfPresent <- apply(mean, 2, mean)
            expTotal <- apply(mean * probPresent, 2, mean)

            # add probPresent, mean if present and mean to ret
            if (sum(!is.na(pmatch(predictions, c("detail", "all")))) > 0) {
                ret <- cbind(ret, expProbPresent, expTotalIfPresent)
                colnames(ret)[ncol(ret) - (1:0)] <- paste(names(fits)[i], c("probability present", "expected total catch if present"))
            }
            if (sum(!is.na(pmatch(predictions, c("TRUE", "all")))) > 0) {
                ret <- cbind(ret, expTotal)
                colnames(ret)[ncol(ret)] <- paste(names(fits)[i], "expected total catch")
            }
        }
    }

    # reorder columns by grouping together species by searching for each in turn
    ix <- numeric(0)
    if (class(eqrs) == "list")
        ix <- 1
    if (classify)
        ix <- c(ix, grep("Bad", colnames(ret)) + 0:4)
    for (i in 1:length(fitNames))
        ix <- c(ix, grep(fitNames[i], colnames(ret)))
    ix <- c(ix, setdiff(1:ncol(ret), ix))  # put any missed columns at end (though shouldn't be any)
    ret <- ret[, ix, drop=FALSE]

    # return
    ret
}

