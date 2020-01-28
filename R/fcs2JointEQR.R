#' Joint EQR
#'
#' Produces Monte Carlo samples of the joint \dfn{Ecological Quality Ratio}
#' (\acronym{EQR}) that combines multiple species and optionally multiple
#' surveys within a group (such as a water body).  This provides a single
#' \acronym{EQR} variable that is found by comparing the observed fish catches
#' with the model's predictions at reference conditions for a number of species
#' and optionally a number of surveys.\cr This function can be used to provide
#' a single \acronym{EQR} variable for each water body that can be passed to
#' \code{\link{fcs2Classify}} to produce probabilistic \acronym{WFD}
#' classifications.
#'
#'
#' @param fit1 an \code{"fcs2Fit"} object containing a full \acronym{FCS2}
#' model fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS =
#' TRUE}.
#' @param \dots further \code{"fcs2Fit"} objects, each resulting from a full
#' \acronym{FCS2} model fit with a different species.
#' @param newData a data frame with surveys as rows and variables as columns.
#' It should contain all variables required by each of the model fits.
#' Covariates which are related to human disturbance (\dfn{pressure variables})
#' should have their values set to the value expected at the site if it were
#' undisturbed (\dfn{reference conditions}) rather than the observed value for
#' each of these variables.
#' @param joinByVar the name of a column in \code{dataFrame} to be used for
#' joining multiple surveys together when calculating the \acronym{EQR}.  For
#' example, \code{joinByVar = "WBId"} would produce \acronym{EQR} samples for
#' each water body, with every survey in a water body contributing towards its
#' \acronym{EQR}.\cr If \code{NULL} (the default), \acronym{EQR} samples are
#' calculated for each survey.
#' @param subset an optional vector specifying a subset of surveys to calculate
#' \acronym{EQR} samples for.
#' @param na.action a function which indicates what should happen when the data
#' contain missing values (\code{NA}s).  The default is set by the
#' \code{na.action} setting of \code{\link{options}} and this is usually set to
#' \code{\link{na.omit}}.  This setting removes surveys that contain missing
#' data in any required variables.  A vector indicating the rows that were
#' removed can be extracted from the returned object using
#' \code{\link{na.action.fcs2EQR}} if surveys are not joined.  Alternatively,
#' \code{\link{na.pass}} can be used to ignore missing values (where possible)
#' or \code{\link{na.fail}} can be given to signal an error if missing values
#' are found.
#' @param n.samples the number of Monte Carlo \acronym{EQR} samples to produce
#' for each survey (or joining variable).
#' @param n.sims the number of Monte Carlo simulations to make for each
#' \acronym{EQR} sample.  These internal samples are used for approximating the
#' probability that defines the joint \acronym{EQR}.
#' @param both if \code{TRUE} and \code{joinByVar} is provided, \acronym{EQR}
#' samples are produced both for each survey and also for each level of the
#' joining variable.  This allows for efficient calculation of the joint
#' \acronym{EQR} for each survey and also the joint \acronym{EQR} that also
#' joins surveys.
#' @param showProgress whether to display the current progress while generating
#' \acronym{EQR} samples.
#' @return If \code{both = FALSE} (the default), a single \code{"fcs2EQR"}
#' object containing Monte Carlo \acronym{EQR} samples is returned.  The
#' \code{"fcs2EQR"} object is essentially a matrix with Monte Carlo samples as
#' rows and columns corresponding to either surveys, if \code{joinByVar} is
#' missing, or to different levels of the joining variable, e.g. different
#' water bodies.
#'
#' If \code{both = TRUE} and \code{joinByVar} is provided, a list containing
#' both of these \code{"fcs2EQR"} objects is returned.
#' @section Warning: Although this function has been written in C for speed, it
#' can still take a long time to produce a large number \code{n.samples} of
#' \acronym{EQR} samples when the number \code{n.sims} of simulations is also
#' high.
#' @seealso \code{\link{print.fcs2EQR}}, \code{\link{summary.fcs2EQR}} and
#' \code{\link{fcs2EQRSummaryMatrix}} for summarising \code{"fcs2EQR"}
#' objects;\cr \code{\link{plot.fcs2EQR}} for plotting \acronym{EQR}
#' variables;\cr \code{\link{mean.fcs2EQR}} and \code{\link{quantile.fcs2EQR}}
#' for calculating means and quantiles of \acronym{EQR} variables
#' respectively;\cr
#'
#' \code{\link{fcs2FitModel}} for producing the required \acronym{FCS2} model
#' fits;\cr \code{\link{fcs2Classify}} for using \acronym{EQR} samples to
#' produce probabilistic \acronym{WFD} classifications.
#'
#' \code{\link{fcs2SingleEQR}} for producing \acronym{EQR} samples for a single
#' species and survey;\cr \code{\link{fcs2JointAndSingleEQR}} for calculating
#' joint and single \acronym{EQR} samples simultaneously.
#' @export
#' @examples
#'
#' \dontrun{
#'
#' ### Example 1: Simple example with two model fits and no covariates
#' ###
#'
#' # simulate random dataset for example
#' Data <- data.frame(SurveyArea=rlnorm(100, 4.6, 0.5))   # random survey area
#' # a single salmon catch per survey
#' Data$Salmon <- rzinbinom(100, size=1.1,
#'                          zeroprob=0.3, nbmean=0.3 * Data$SurveyArea)
#' # a single trout catch per survey
#' Data$Trout <- rzinbinom(100, size=0.87,
#'                         zeroprob=0.19, nbmean=0.42 * Data$SurveyArea)
#'
#' # fit full model for salmon with OpenBUGS
#' salmonFit <- fcs2FitModel("Salmon", dataFrame=Data, surveyAreaVar="SurveyArea",
#'                           runBUGS=TRUE, n.iter=1000, bugsProgram="OpenBUGS")
#'
#' # fit full model for trout with OpenBUGS
#' troutFit <- fcs2FitModel("Trout", dataFrame=Data, surveyAreaVar="SurveyArea",
#'                          runBUGS=TRUE, n.iter=1000, bugsProgram="OpenBUGS")
#'
#' # calculate samples of EQR that combines salmon and trout for each survey,
#' #   using same dataset as no pressure variables to adjust to reference values
#' eqr <- fcs2JointEQR(salmonFit, troutFit, newData=Data,
#'                     n.samples=100, n.sims=100)
#'
#' # plot EQR variables for first 9 surveys
#' plot(eqr, 1:9, boundaries=NULL)
#'
#' # calculate mean EQR values
#' mean(eqr)
#'
#'
#'
#' ### Example 2: Joining surveys as well as species to create water body EQRs
#' ###
#'
#' # extend dataset to include water body indicator
#' # randomly assign water body A, B or C to each survey
#' Data$WaterBody <- sample(c("A", "B", "C"), 100, replace=TRUE)
#'
#' # calculate samples of EQR that combines salmon and trout
#' #   as well as combining all surveys within each water body
#' eqr <- fcs2JointEQR(salmonFit, troutFit, newData=Data,
#'                     n.samples=100, n.sims=100, joinByVar="WaterBody")
#'
#' # plot EQR variables for all 3 water bodies
#' plot(eqr, 1:3, boundaries=NULL)
#'
#' # calculate mean EQR values
#' mean(eqr)
#' }
#'
fcs2JointEQR <-
function(fit1, ..., newData, joinByVar = NULL, subset = 1:nrow(newData), na.action,
         n.samples = 1000, n.sims = 1000, both = FALSE, showProgress = TRUE)
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

    ## check whether joinByVar in newData
    if (!is.null(joinByVar) && !(joinByVar %in% names(newData)))
        stop("join by variable '", joinByVar, "' not found in data frame 'newData'")

    ## create list of fits, named by their (first) catch variable
    if (missing(fit1))
        fits <- list(...)
    else
        fits <- list(fit1, ...)

    # check each non-matched object is a fit and extract fit names
    k <- length(fits)
    for (i in 1:k) {
        # check fcs2Fit object
        if (!inherits(fits[[i]], "fcs2Fit"))
            stop(paste("unrecognised argument", names(fits)[i]))

        # check fit contains bugs fit
        if (is.null(fits[[i]]$bugsFit))
            stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

        if (is.null(names(fits)[i]) || is.na(names(fits)[i]) || names(fits)[i] == "") {
            if (!is.null(fits[[i]]$runTotalVars))
                names(fits)[i] <- sub(".Run1Total", "", fits[[i]]$runTotalVars[1], fixed=TRUE)
            else if (!is.null(fits[[i]]$allRunsRangeVars))
                names(fits)[i] <- sub(".AllRunsTotalMin", "", fits[[i]]$allRunsRangeVars[1], fixed=TRUE)
            else
                names(fits)[i] <- sub(".AllRunsTotal", "", fits[[i]]$allRunsTotalVar[1], fixed=TRUE)
        }
    }

    ## Sample sets of indices corresponding to parameter samples
    ix <- array(dim=c(n.samples, k))
    for (i in 1:k) ix[, i] <- sample(1:fits[[i]]$bugsFit$n.sims, n.samples, replace=TRUE)

    # for each fit, calculate abundance and prevalence, storing only values needed
    nbmean <- zeroprob <- array(NA, dim=c(n.samples, length(subset), k))
    na.action.part <- c()
    for (i in 1:k) {
        # extract abundance mu
        singleMx <- abundance(fits[[i]], newData, subset, na.action)
        na.action.single <- attr(singleMx, "na.action")
        na.action.part <- union(na.action.part, na.action.single)

        # multiply mu by survey area
        if (!is.null(na.action.single))
            nbmean[, -na.action.single, i] <- singleMx[ix[, i], ] *
                        matrix(newData[subset[-na.action.single], fits[[i]]$surveyAreaVar], byrow=TRUE, nrow=n.samples, ncol=ncol(singleMx))
        else
            nbmean[, , i] <- singleMx[ix[, i], ] * matrix(newData[subset, fits[[i]]$surveyAreaVar], byrow=TRUE, nrow=n.samples, ncol=ncol(singleMx))

        # extract prevalence rho
        singleMx <- prevalence(fits[[i]], newData, subset, na.action)
        na.action.single <- attr(singleMx, "na.action")
        na.action.part <- union(na.action.part, na.action.single)

        # calculate zeroprob = 1 - rho
        if (!is.null(na.action.single))
            zeroprob[, -na.action.single, i] <- 1 - singleMx[ix[, i], ]
        else
            zeroprob[, , i] <- 1 - singleMx[ix[, i], ]
    }
    rm(singleMx)

    # if multiple runs, calculate or extract number of runs
    nRuns <- array(NA, dim=c(nrow(newData), k), dimnames=list(rownames(newData), names(fits)))
    for (i in 1:k) {
        if (fits[[i]]$multiRun) {
            if (fits[[i]]$nRunsVar %in% colnames(newData))
                nRuns[, i] <- newData[, fits[[i]]$nRunsVar]
            else {
                for (j in 1:nrow(newData))
                    nRuns[j, i] <- sum(!is.na(newData[j, fits[[i]]$runTotalVars]))
            }

            # check number of runs variable isn't larger than number of catch variables
            if (max(nRuns[, i], na.rm=TRUE) > length(fits[[i]]$runTotalVars)) {
                warning(paste("Number of runs ", if (fits[[i]]$nRunsVar %in% colnames(newData))
                                                     paste("variable '", fits[[i]]$nRunsVar, "' ", sep=""),
                                                 "clipped to not exceed number of catch variables (", length(fits[[i]]$runTotalVars), ") for fit ", i, sep=""))
                nRuns[!is.na(nRuns[, i]) & nRuns[, i] > length(fits[[i]]$runTotalVars), i] <- length(fits[[i]]$runTotalVars)
            }
        }
    }

    # extract (total) catch min and max
    catchMin <- catchMax <- array(dim=c(nrow(newData), k), dimnames=list(rownames(newData), names(fits)))
    for (i in 1:k) {
        if (!fits[[i]]$multiRun)
            catchMin[, i] <- catchMax[, i] <- newData[, fits[[i]]$allRunsTotalVar]
        else {
            # from runs
            if (!is.null(fits[[i]]$runTotalVars)) {
                for (j in 1:nrow(newData))
                    if (!is.na(nRuns[j, i]))
                        catchMin[j, i] <- catchMax[j, i] <- sum(as.numeric(newData[j, fits[[i]]$runTotalVars[1:nRuns[j, i]]]))
            }

            # from all runs total
            if (!is.null(fits[[i]]$allRunsTotalVar) && fits[[i]]$allRunsTotalVar %in% names(newData) && sum(is.na(catchMin[, i])) > 0)
                catchMin[is.na(catchMin[, i]), i] <- catchMax[is.na(catchMin[, i]), i] <- newData[is.na(catchMin[, i]), fits[[i]]$allRunsTotalVar]
        }
        # from all runs range
        if (!is.null(fits[[i]]$allRunsRangeVars) && sum(is.na(catchMin[, i])) > 0) {
            catchMin[is.na(catchMin[, i]), i] <- newData[is.na(catchMin[, i]), fits[[i]]$allRunsRangeVars[1]]
            catchMax[is.na(catchMax[, i]), i] <- newData[is.na(catchMax[, i]), fits[[i]]$allRunsRangeVars[2]]
        }
    }

    # extract survey area variables
    surveyAreaVars <- c()
    for (i in 1:k)
        surveyAreaVars <- c(surveyAreaVars, fits[[i]]$surveyAreaVar)

    # calculate overall list of removed sites/surveys
    na.action <- union(na.action.part,
                       union(attr(na.action(catchMin[subset, ]), "na.action"),
                             union(attr(na.action(newData[subset, unique(surveyAreaVars)]), "na.action"),
                                   attr(na.action(newData[subset, joinByVar]), "na.action"))))
    if (!is.null(na.action)) {
        na.action <- sort(na.action)
        names(na.action) <- rownames(newData)[subset][na.action]
        class(na.action) <- "omit"   # assumes na.action was na.omit (as too much effort to search where it came from)
    }

    # remove these from each abundance and zeroprob
    if (!is.null(na.action)) {
        nbmean <- nbmean[, -na.action, , drop=FALSE]
        zeroprob <- zeroprob[, -na.action, , drop=FALSE]
    }

    # calculate overall subset
    isubset <- subset[setdiff(1:length(subset), na.action)]
    catchMin <- catchMin[isubset, , drop=FALSE]
    catchMax <- catchMax[isubset, , drop=FALSE]

    # if multiple runs, multiply nbmean by (1 - (1 - q)^nRuns)
    for (i in 1:k) {
        if (fits[[i]]$multiRun)
            nbmean[, , i] <- nbmean[, , i] *
                    (1 - (1 - matrix(fits[[i]]$bugsFit$sims.list$q[ix[, i]], nrow=n.samples, ncol=length(isubset))) ^
                                  matrix(nRuns[isubset, i], byrow=TRUE, nrow=n.samples, ncol=length(isubset)))
    }

    # if multiple runs, extract or set index of values to join (eg join sites within same water body)
    if (is.null(joinByVar))
        iJoin <- 1:length(isubset)
    else
        iJoin <- as.numeric(factor(newData[isubset, joinByVar]))
    nJoinLevels <- max(iJoin)

    # extract r samples
    r <- array(NA, dim=c(n.samples, k))
    for (i in 1:k)
        r[, i] <- fits[[i]]$bugsFit$sims.list$r[ix[, i]]

    # create further subset of non-NAs
    subsubset <- rep(TRUE, length(isubset))
    for (i in 1:k)
        subsubset <- subsubset & !is.na(catchMin[, i]) & !is.na(nbmean[1, , i]) & !is.na(zeroprob[1, , i])
    nSurveys <- sum(subsubset)
    catchMin <- catchMin[subsubset, , drop=FALSE]
    catchMax <- catchMax[subsubset, , drop=FALSE]
    nbmean <- nbmean[, subsubset, , drop=FALSE]
    zeroprob <- zeroprob[, subsubset, , drop=FALSE]
    if (!is.null(joinByVar))
        iJoin <- iJoin[subsubset]

    # calculate both by survey and joining surveys
    if (both && !is.null(joinByVar)) {

        # check number of non-NA surveys
        if (nSurveys > 0) {
            ret <- .C("jointEQRJoiningSurveysAndNot", as.integer(catchMin), as.integer(catchMax), as.double(r), as.double(nbmean), as.double(zeroprob),
                      as.integer(iJoin), as.integer(nSurveys), as.integer(k), as.integer(n.samples), as.integer(n.sims), as.integer(nJoinLevels), as.integer(showProgress),
                      EQRBySurvey=double(n.samples * nSurveys), EQRJoiningSurveys=double(n.samples * nJoinLevels),
                      PACKAGE="fcs2")[c("EQRBySurvey", "EQRJoiningSurveys")]

            dim(ret$EQRBySurvey) <- c(n.samples, nSurveys)
            dim(ret$EQRJoiningSurveys) <- c(n.samples, nJoinLevels)

            # extend EQRBySurvey with NAs
            cols <- rep(NA, length(isubset))
            cols[subsubset] <- 1:nSurveys
            ret$EQRBySurvey <- ret$EQRBySurvey[, cols]

        } else {
            ret <- list(EQRBySurvey=array(NA, dim=c(n.samples, length(isubset))),
                        EQRJoiningSurveys=array(NA, dim=c(n.samples, nJoinLevels)))
        }

        # add column names
        dimnames(ret$EQRBySurvey) <- list(NULL, survey=rownames(newData)[isubset])
        dimnames(ret$EQRJoiningSurveys) <- list(NULL, levels(factor(newData[isubset, joinByVar])))
        names(dimnames(ret$EQRJoiningSurveys))[2] <- joinByVar

        # add na.action to EQRBySurvey
        attr(ret$EQRBySurvey, "na.action") <- na.action

        # create na.action to indicate join levels that are excluded for EQRJoiningSurveys
        na.action.names <- setdiff(levels(factor(newData[subset, joinByVar])), levels(factor(newData[isubset, joinByVar])))
        if (length(na.action.names) != 0) {
            na.action <- match(na.action.names, levels(factor(newData[subset, joinByVar])))
            names(na.action) <- na.action.names
            attr(ret$EQRJoiningSurveys, "na.action") <- na.action
        }

        # change name of joined component
        names(ret)[2] <- paste("EQRBy", joinByVar, sep='')

        # set classes to fcs2EQR
        class(ret[[1]]) <- class(ret[[2]]) <- "fcs2EQR"

        # return
        return(ret)
    }

    # check number of non-NA surveys
    if (nSurveys > 0) {

        ## Calculate joint EQR for each index set

        if (is.null(joinByVar)) {
            # by survey
            ret <- .C("jointEQR", as.integer(catchMin), as.integer(catchMax), as.double(r), as.double(nbmean), as.double(zeroprob),
                      as.integer(nSurveys), as.integer(k), as.integer(n.samples), as.integer(n.sims), as.integer(showProgress),
                      PACKAGE="fcs2", ret=double(n.samples * nSurveys))$ret
            dim(ret) <- c(n.samples, nSurveys)

            # extend ret with NAs
            cols <- rep(NA, length(isubset))
            cols[subsubset] <- 1:nSurveys
            ret <- ret[, cols]

        } else {
            # joining surveys
            ret <- .C("jointEQRJoiningSurveys", as.integer(catchMin), as.integer(catchMax), as.double(r), as.double(nbmean), as.double(zeroprob),
                      as.integer(iJoin), as.integer(nSurveys), as.integer(k), as.integer(n.samples), as.integer(n.sims), as.integer(nJoinLevels),
                      as.integer(showProgress), PACKAGE="fcs2", ret=double(n.samples * nJoinLevels))$ret
            dim(ret) <- c(n.samples, nJoinLevels)

            # change na.action to indicate join levels that are excluded
            na.action.names <- setdiff(levels(factor(newData[subset, joinByVar])), levels(factor(newData[isubset, joinByVar])))
            if (length(na.action.names) == 0)
                na.action <- NULL
            else {
                na.action <- match(na.action.names, levels(factor(newData[subset, joinByVar])))
                names(na.action) <- na.action.names
                class(na.action) <- "omit"
            }
        }

    } else
        ret <- array(NA, dim=c(n.samples, nJoinLevels))

    # set dimnames
    if (is.null(joinByVar))
        dimnames(ret) <- list(NULL, survey=rownames(newData)[isubset])
    else {
        dimnames(ret) <- list(NULL, levels(factor(newData[isubset, joinByVar])))
        names(dimnames(ret))[2] <- joinByVar
    }

    # add na.action attribute
    attr(ret, "na.action") <- na.action

    # set class to fcs2EQR
    class(ret) <- "fcs2EQR"

    ret
}

