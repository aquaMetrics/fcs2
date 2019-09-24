#' Automated Model Selection
#'
#' Automatically selects an optimal set of covariate terms for the abundance
#' and prevalence regression equations of the \acronym{FCS2} model.  Terms are
#' attempted sequentially and the approximate abundance and prevalence
#' \acronym{INLA} fits are used to test the significance of each term.
#'
#'
#' @param runTotalVars a character vector of columns in \code{dataFrame} that
#' give the number of fish caught in each run of a multiple-run survey.  These
#' should be in order from the first run to the last.
#' @param allRunsTotalVar the name of a column in \code{dataFrame} that gives
#' the total number of fish caught over all runs of a multiple-run survey.
#' @param allRunsRangeVars the names of two columns in \code{dataFrame} that
#' give the minimum and the maximum value of a range of possible values for the
#' total number of fish caught over all runs in a multiple-run survey.
#' @param dataFrame a data frame with surveys as rows and variables as columns.
#' It should contain all variables specified by other arguments.
#' @param surveyAreaVar the name of a column in \code{dataFrame} that gives the
#' survey area.  If not specified, the function will search for the default
#' value \code{"SurveyArea"}.
#' @param nRunsVar the name of a column in \code{dataFrame} that gives the
#' number of runs in each survey.  If missing, the number of runs is assumed to
#' be the number of non-missing run total entries, unless \code{runTotalVars}
#' has length 1 or is missing in which case a single-run model is used.  Number
#' of runs values greater than the number of run total columns are clipped with
#' a warning.
#' @param muVars a character vector naming variables to use for terms in the
#' abundance regression equation.  Variables are attempted one-by-one in the
#' order given.
#' @param muVarType a character vector of the same length as \code{muVars}
#' indicating the type of abundance term to attempt for each variable.  Each
#' element should be one of \code{"asis"}, \code{"factor"}, \code{"linear"},
#' \code{"continuous"} or \code{"spatial"}.  See \sQuote{Details} for a
#' description of each.
#' @param rhoVars a character vector naming variables to use for terms in the
#' prevalence regression equation.  Variables are attempted one-by-one in the
#' order given.  Defaults to \code{muVars} to use the same variables as
#' specified for abundance.
#' @param rhoVarType a character vector of the same length as \code{muVars}
#' indicating the type of prevalence term to attempt for each variable.  Each
#' element should be one of \code{"asis"}, \code{"factor"}, \code{"linear"},
#' \code{"continuous"} or \code{"spatial"}.  See \sQuote{Details} for a
#' description of each.  Defaults to \code{muVarType} to use the same variables
#' as specified for abundance.
#' @param rhoFormula an optional \code{\link{formula}} specifying which terms
#' should appear in the prevalence regression equation when selecting terms for
#' abundance.  If specified, the prevalence model selection is skipped.
#' @param subset an optional vector specifying a subset of surveys to be used
#' in the fitting process.
#' @param tolerance the threshold for each term's significance probability,
#' below which a term is considered to be significant and is retained.  The
#' default value is \code{0.01} but a smaller value will cause fewer terms to
#' be retained and vice-versa.\cr See \code{\link{summary.fcs2Fit}} for a
#' definition of the significance probability for each variable.  The
#' probabilities corresponding to each variable within a \code{\link{rw2}} or
#' \code{\link{spatial}} term are combined (by rescaling the minimum under the
#' assumption of independence) before comparison with \code{tolerance}.
#' @param maxNoOrders the maximum number of polynomial orders to try for each
#' linear component.  Defaults to 3.
#' @param nSweeps the number of sweeps to make through each list of potential
#' variables.  The default is \code{2} so that a second test is made to each
#' term in case the presence of later variables will make an earlier term
#' significant.
#' @param prior.parameters an optional named list of named vectors giving the
#' parameter values to use for the prior distribution of a variable.  See
#' \code{\link{fcs2Priors}} for further details on the prior distributions and
#' how to specify prior parameters.  Priors that are not given parameter values
#' have defaults given to them by \code{\link{.fcs2SetDefaultPriors}}.  Note
#' that priors for linear variables are ignored by \acronym{INLA}.
#' @param estAllRunsTotalVar the name of a column in \code{dataFrame} that
#' gives an estimate of the total number of fish caught over all runs for
#' surveys where only range data is available.  This is used in the approximate
#' abundance \acronym{INLA} fit only and is not used in the full model as
#' fitted with \acronym{BUGS}.  If this is not provided and range data are
#' present, the central value of each range is used with a warning.
#' @param verbose whether to print progress to screen.
#' @return a list of two matricies, each containing a summary of the terms
#' attempted in each iteration.  The first matrix gives the term history for
#' the prevalence (\eqn{\rho}{rho}) and the second gives the history for the
#' abundance (\eqn{\mu}{mu}).
#'
#' Each row of a matrix summarises the regression formula used for that
#' component, with covariates appearing as columns.  Linear terms are
#' represented as a number giving the order of the term, random walk terms are
#' given by \code{"rw2"} and spatial terms by \code{"spatial"}.
#' \code{fcs2:::termSummary2Formula} can be used to convert a row to an formula
#' for simplier input into \code{\link{fcs2FitModel}}.
#' @note Since the prevalence terms are selected first and then the abundance
#' terms second, it is possible that some prevalence terms are no longer
#' significant in the final fit as these were not attempted with the selected
#' abundance formula.
#' @seealso \code{\link{fcs2FitModel}}, \code{fcs2:::termSummary2Formula}
#' @export
fcs2ModelSelection <-
function(runTotalVars = NULL, allRunsTotalVar = NULL, allRunsRangeVars = NULL, dataFrame, surveyAreaVar = "SurveyArea", nRunsVar = NULL,
         muVars, muVarType, rhoVars = muVars, rhoVarType = muVarType, rhoFormula, subset = 1:nrow(dataFrame),
         tolerance = 0.01, maxNoOrders = 3, nSweeps = 1,
         prior.parameters = list(), estAllRunsTotalVar = NULL, verbose = FALSE)
{
    # check arguments
    if (length(muVars) != length(muVarType))
        stop("'muVars' and 'muVarType' must be the same length")
    if (length(rhoVars) != length(rhoVarType))
        stop("'rhoVars' and 'rhoVarType' must be the same length")

    # convert factors to individual terms
    rhoAllVars <- rhoAllVarType <- character(0)
    for (i in 1:length(rhoVars)) {
        if (rhoVarType[i] == "factor") {
            # add individual factor terms to rhoAllVars
            terms <- levels(factor(dataFrame[subset, rhoVars[i]]))
            rhoAllVars <- c(rhoAllVars, paste('I(', rhoVars[i], ' == "', terms, '")', sep=''))
            rhoAllVarType <- c(rhoAllVarType, rep("asis", length(terms)))

        } else {
            # add term as is to rhoAllVars
            rhoAllVars <- c(rhoAllVars, rhoVars[i])
            rhoAllVarType <- c(rhoAllVarType, rhoVarType[i])
        }
    }
    rhoVars <- rhoAllVars
    rhoVarType <- rhoAllVarType

    muAllVars <- muAllVarType <- character(0)
    for (i in 1:length(muVars)) {
        if (muVarType[i] == "factor") {
            # add individual factor terms to muAllVars
            terms <- levels(factor(dataFrame[subset, muVars[i]]))
            muAllVars <- c(muAllVars, paste('I(', muVars[i], ' == "', terms, '")', sep=''))
            muAllVarType <- c(muAllVarType, rep("asis", length(terms)))

        } else {
            # add term as is to muAllVars
            muAllVars <- c(muAllVars, muVars[i])
            muAllVarType <- c(muAllVarType, muVarType[i])
        }
    }
    muVars <- muAllVars
    muVarType <- muAllVarType

    ## Prevalence rho
    if (missing(rhoFormula)) {
        cat("\nPrevalence terms:\n")

        # create array to store results of attempted fits
        rhoFits <- array("", dim=c(0, length(rhoVars)), dimnames=list(NULL, rhoVars))

        # initial rho terms
        currentTerms <- character(0)

        # try adding rho terms one by one
        for (i in rep(1:length(rhoVars), nSweeps)) {
            term <- rhoVars[i]
            termType <- rhoVarType[i]
            orders <- NULL

            # set possible orders of term to try
            if (termType == "asis")
                orders <- 1
            else if (termType == "spatial")
                orders <- -1
            else if (termType == "linear") {
                # count number of terms in there already
                termCount <- length(grep(term, currentTerms, fix=TRUE))

                # continue from current place
                if (termCount < maxNoOrders)
                    orders <- (termCount + 1):maxNoOrders

            } else if (termType == "continuous") {
                # if rw2 term there, break
                if (length(grep(paste("rw2(", term, ")", sep=""), currentTerms, fix=TRUE)) > 0)
                    next

                # count number of terms in there already
                termCount <- length(grep(term, currentTerms, fix=TRUE))

                # continue from current place
                if (termCount < maxNoOrders)
                    orders <- (termCount + 1):maxNoOrders

                # if no terms there, try rw2 first
                if (termCount == 0)
                    orders <- c(0, orders)

            }

            for (j in orders) {

                # add new term
                if (j == -1)
                    term <- paste("spatial(", rhoVars[i], ", adjacency)", sep="")  ### NOTE: assumes adjacency variable named 'adjacency'
                else if (j == 0)
                    term <- paste("rw2(", rhoVars[i], ")", sep="")
                else if (j == 1)
                    term <- rhoVars[i]
                else
                    term <- paste("I(", rhoVars[i], "^", j, ")", sep="")

                # check that term isn't already there (possible if multiple sweeps)
                if (length(grep(term, currentTerms, fix=TRUE)) > 0)
                    break

                # add term
                cat("\n+ Adding term", term, "\n")
                currentTermsBeforeAdding <- currentTerms
                currentTerms <- c(currentTerms, term)

                # add entry in rhoFits and extract prior parameters
                rhoFits <- rbind(rhoFits, "")
                rhoPriors <- list()
                for (k in 1:length(rhoVars)) {
                    # search for var in current terms
                    pos <- grep(rhoVars[k], currentTerms, fixed=TRUE)
                    no <- length(pos)

                    if (no > 0) {
                        # add rhoFits entry
                        if (no == 1 && length(grep("spatial", currentTerms[pos], fixed=TRUE)) > 0)
                            rhoFits[nrow(rhoFits), k] <- "spatial"
                        else if (no == 1 && length(grep("rw2", currentTerms[pos], fixed=TRUE)) > 0)
                            rhoFits[nrow(rhoFits), k] <- "rw2"
                        else
                            rhoFits[nrow(rhoFits), k] <- no

                        # search for prior parameters
                        if (rhoFits[nrow(rhoFits), k] == "spatial" || rhoFits[nrow(rhoFits), k] == "rw2") {
                            # non-linear term, so search for nu.var or phi.var
                            varName <- paste("nu", rhoVars[k], sep=".")
                            if (varName %in% names(prior.parameters))
                                rhoPriors <- c(rhoPriors, prior.parameters[varName])
                            varName <- paste("phi", rhoVars[k], sep=".")
                            if (varName %in% names(prior.parameters))
                                rhoPriors <- c(rhoPriors, prior.parameters[varName])

                        } else {
                            # linear term(s), so search for gamma.term for each term in currentTerms[pos]
                            for (t in currentTerms[pos]) {
                                varName <- paste("gamma", t, sep=".")
                                if (varName %in% names(prior.parameters))
                                    rhoPriors <- c(rhoPriors, prior.parameters[varName])
                            }
                        }
                    }
                }

                # fit model
                rhoFormula <- formula(paste("~", paste(currentTerms, collapse=" + ")))
                cat("Fitting model with\n")
                print(formula(paste("logit(rho)", paste(deparse(rhoFormula), collapse=" "))), showEnv=FALSE)
                fit <- fcs2FitModel(runTotalVars, allRunsTotalVar, allRunsRangeVars, ~1, rhoFormula, dataFrame, surveyAreaVar, nRunsVar, subset,
                                    na.omit, rhoPriors, verbose=verbose, estAllRunsTotalVar=estAllRunsTotalVar, runINLA=FALSE, runBUGS=FALSE)
                if (class(try({
                    fit <- fcs2FitModel(fit=fit, verbose=verbose, runINLA="rho", runBUGS=FALSE)

                }, silent=TRUE)) == "try-error") {
                    # error fitting model, so reset terms as before added this term
                    cat("\n- Removing term ", term, " as error fitting model with INLA\n", sep="")
                    currentTerms <- currentTermsBeforeAdding
                    if (j > 0)
                        break
                    next
                }

                # calculate summary (with 10 attempts as INLA can be tempermental)
                rhoSum <- NULL
                for (k in 1:10) {
                    try({
                        rhoSum <- summary(fit, prior=FALSE, allVars=TRUE)$inla$rho
                        break
                    }, silent=TRUE)
                }
                if (is.null(rhoSum)) {
                    # error calculating summary, so reset terms as before added this term
                    cat("\n- Removing term ", term, " as error calculating INLA summary (after 10 attempts)\n", sep="")
                    currentTerms <- currentTermsBeforeAdding
                    if (j > 0)
                        break
                    next
                }

                # extract significance prob of term
                rhoSignifPr <- rhoSum[, "Signif. Pr."]
                names(rhoSignifPr) <- rownames(rhoSum)
                rhoSignifPr <- rhoSignifPr[!is.na(rhoSignifPr)]  # keep only vars with prob

                # search for multiple vars with each name
                for (k in 1:length(rhoVars)) {
                    pos <- grep(rhoVars[k], names(rhoSignifPr), fixed=TRUE)

                    # if multiple random, use transformed min
                    if (length(pos2 <- grep("[", names(rhoSignifPr)[pos], fixed=TRUE)) > 0) {
                        pos <- pos[pos2]
                        minPr <- min(rhoSignifPr[pos])
                        rhoSignifPr <- rhoSignifPr[-pos[-1]]
                        rhoSignifPr[pos[1]] <- 1 - ((1 - minPr) ^ length(pos))  # transform min prob so on uniform scale (assuming independence!)
                        if (rhoVarType[k] == "continuous")
                            names(rhoSignifPr)[pos[1]] <- paste("rw2(", rhoVars[k], ")", sep="")
                        else if (rhoVarType[k] == "spatial")
                            names(rhoSignifPr)[pos[1]] <- paste("spatial(", rhoVars[k], ", adjacency)", sep="")

                        # remove lower orders of terms so don't remove these first
                    } else if (length(pos) > 1)
                        rhoSignifPr <- rhoSignifPr[-pos[-length(pos)]]
                }

                # remove insignificant terms, in order of least significance
                while (length(rhoSignifPr) > 0 && max(rhoSignifPr) > tolerance) {
                    # remove term with largest signif prob (which must be > tolerance)
                    worstTerm <- sub("gamma.", "", names(which.max(rhoSignifPr)), fixed=TRUE)
                    worstTerm <- sub("TRUE", "", worstTerm)  # remove TRUE if present
                    worstTerm <- sub("FALSE", "", worstTerm)  # remove TRUE if present
                    cat("\n- Removing term ", worstTerm, " as significance probability = ", max(rhoSignifPr), "\n", sep="")
                    currentTerms <- setdiff(currentTerms, worstTerm)

                    # fit model (unless first term to remove is term just added, in which case break)
                    if (length(currentTerms) > 0 && !(term == worstTerm && (rhoFits[nrow(rhoFits) - 1, i] == "" ||
                                                                            (nchar(rhoFits[nrow(rhoFits), i]) == 1 && as.numeric(rhoFits[nrow(rhoFits), i]) > 1 &&
                                                                             as.numeric(rhoFits[nrow(rhoFits), i]) - 1 == as.numeric(rhoFits[nrow(rhoFits) - 1, i]))))) {
                        # add entry in rhoFits and extract prior parameters
                        rhoFits <- rbind(rhoFits, "")
                        rhoPriors <- list()
                        for (k in 1:length(rhoVars)) {
                            # search for var in current terms
                            pos <- grep(rhoVars[k], currentTerms, fixed=TRUE)
                            no <- length(pos)

                            if (no > 0) {
                                # add rhoFits entry
                                if (no == 1 && length(grep("spatial", currentTerms[pos], fixed=TRUE)) > 0)
                                    rhoFits[nrow(rhoFits), k] <- "spatial"
                                else if (no == 1 && length(grep("rw2", currentTerms[pos], fixed=TRUE)) > 0)
                                    rhoFits[nrow(rhoFits), k] <- "rw2"
                                else
                                    rhoFits[nrow(rhoFits), k] <- no

                                # search for prior parameters
                                if (rhoFits[nrow(rhoFits), k] == "spatial" || rhoFits[nrow(rhoFits), k] == "rw2") {
                                    # non-linear term, so search for nu.var or phi.var
                                    varName <- paste("nu", rhoVars[k], sep=".")
                                    if (varName %in% names(prior.parameters))
                                        rhoPriors <- c(rhoPriors, prior.parameters[varName])
                                    varName <- paste("phi", rhoVars[k], sep=".")
                                    if (varName %in% names(prior.parameters))
                                        rhoPriors <- c(rhoPriors, prior.parameters[varName])

                                } else {
                                    # linear term(s), so search for gamma.term for each term in currentTerms[pos]
                                    for (t in currentTerms[pos]) {
                                        varName <- paste("gamma", t, sep=".")
                                        if (varName %in% names(prior.parameters))
                                            rhoPriors <- c(rhoPriors, prior.parameters[varName])
                                    }
                                }
                            }
                        }

                        # fit model
                        rhoFormula <- formula(paste("~", paste(currentTerms, collapse=" + ")))
                        cat("Fitting model with\n")
                        print(formula(paste("logit(rho)", paste(deparse(rhoFormula), collapse=" "))), showEnv=FALSE)
                        fit <- fcs2FitModel(runTotalVars, allRunsTotalVar, allRunsRangeVars, ~1, rhoFormula, dataFrame, surveyAreaVar, nRunsVar, subset,
                                            na.omit, rhoPriors, verbose=verbose, estAllRunsTotalVar=estAllRunsTotalVar, runINLA=FALSE, runBUGS=FALSE)
                        if (class(try({
                            fit <- fcs2FitModel(fit=fit, verbose=verbose, runINLA="rho", runBUGS=FALSE)

                        }, silent=TRUE)) == "try-error") {
                            # error fitting model, so reset terms as before added this term
                            cat("\n- Removing term ", term, " and its consequences as error fitting model with INLA\n", sep="")
                            currentTerms <- currentTermsBeforeAdding
                            rhoSignifPr <- NULL
                            next  # note: only breaks out of while loop
                        }

                        # calculate summary (with 10 attempts as INLA can be tempermental)
                        rhoSum <- NULL
                        for (k in 1:10) {
                            try({
                                rhoSum <- summary(fit, prior=FALSE, allVars=TRUE)$inla$rho
                                break
                            }, silent=TRUE)
                        }
                        if (is.null(rhoSum)) {
                            # error calculating summary, so reset terms as before added this term
                            cat("\n- Removing term ", term, " as error calculating INLA summary (after 10 attempts)\n", sep="")
                            currentTerms <- currentTermsBeforeAdding
                            rhoSignifPr <- NULL
                            next  # note: only breaks out of while loop
                        }

                        # extract significance prob of term
                        rhoSignifPr <- rhoSum[, "Signif. Pr."]
                        names(rhoSignifPr) <- rownames(rhoSum)
                        rhoSignifPr <- rhoSignifPr[!is.na(rhoSignifPr)]  # keep only vars with prob

                        # search for multiple vars with each name
                        for (k in 1:length(rhoVars)) {
                            pos <- grep(rhoVars[k], names(rhoSignifPr), fixed=TRUE)

                            # if multiple random, keep only mean
                            if (length(pos2 <- grep("[", names(rhoSignifPr)[pos], fixed=TRUE)) > 0) {
                                pos <- pos[pos2]
                                minPr <- min(rhoSignifPr[pos])
                                rhoSignifPr <- rhoSignifPr[-pos[-1]]
                                rhoSignifPr[pos[1]] <- 1 - ((1 - minPr) ^ length(pos))  # transform min prob so on uniform scale (assuming independence!)
                                if (rhoVarType[k] == "continuous")
                                    names(rhoSignifPr)[pos[1]] <- paste("rw2(", rhoVars[k], ")", sep="")
                                else if (rhoVarType[k] == "spatial")
                                    names(rhoSignifPr)[pos[1]] <- paste("spatial(", rhoVars[k], ", adjacency)", sep="")

                                # remove lower orders of terms so don't remove these first
                            } else if (length(pos) > 1)
                                rhoSignifPr <- rhoSignifPr[-pos[-length(pos)]]
                        }

                    } else
                        rhoSignifPr <- NULL
                }

                # if rw2 and still there, break
                if (j == 0 && term %in% currentTerms)
                    break

                # if more orders to try, check this term is still in there and break if not
                if (j > 0 && j < max(orders) && !(term %in% currentTerms))
                    break
             }
        }

        # add entry in rhoFits (unless this is already final entry ie no terms were removed at end)
        finalRhoFit <- rhoFits[nrow(rhoFits), ]
        finalRhoFit[finalRhoFit == "spatial" | finalRhoFit == "rw2"] <- 1
        finalRhoFit <- as.numeric(finalRhoFit)
        if (sum(finalRhoFit, na.rm=TRUE) != length(currentTerms)) {
            rhoFits <- rbind(rhoFits, "")
            rhoPriors <- list()
            for (k in 1:length(rhoVars)) {
                # search for var in current terms
                pos <- grep(rhoVars[k], currentTerms, fixed=TRUE)
                no <- length(pos)

                if (no > 0) {
                    # add rhoFits entry
                    if (no == 1 && length(grep("spatial", currentTerms[pos], fixed=TRUE)) > 0)
                        rhoFits[nrow(rhoFits), k] <- "spatial"
                    else if (no == 1 && length(grep("rw2", currentTerms[pos], fixed=TRUE)) > 0)
                        rhoFits[nrow(rhoFits), k] <- "rw2"
                    else
                        rhoFits[nrow(rhoFits), k] <- no

                    # search for prior parameters
                    if (rhoFits[nrow(rhoFits), k] == "spatial" || rhoFits[nrow(rhoFits), k] == "rw2") {
                        # non-linear term, so search for nu.var or phi.var
                        varName <- paste("nu", rhoVars[k], sep=".")
                        if (varName %in% names(prior.parameters))
                            rhoPriors <- c(rhoPriors, prior.parameters[varName])
                        varName <- paste("phi", rhoVars[k], sep=".")
                        if (varName %in% names(prior.parameters))
                            rhoPriors <- c(rhoPriors, prior.parameters[varName])

                    } else {
                        # linear term(s), so search for gamma.term for each term in currentTerms[pos]
                        for (t in currentTerms[pos]) {
                            varName <- paste("gamma", t, sep=".")
                            if (varName %in% names(prior.parameters))
                                rhoPriors <- c(rhoPriors, prior.parameters[varName])
                        }
                    }
                }
            }
        }

        # print best rho fit
        if (length(currentTerms) == 0)
            currentTerms <- "1"
        rhoFormula <- formula(paste("~", paste(currentTerms, collapse=" + ")))
        cat("\nBest prevalence fit is\n")
        print(formula(paste("logit(rho)", paste(deparse(rhoFormula), collapse=" "))), showEnv=FALSE)
        fitRho <- fcs2FitModel(runTotalVars, allRunsTotalVar, allRunsRangeVars, ~1, rhoFormula, dataFrame, surveyAreaVar, nRunsVar, subset,
                               na.omit, rhoPriors, verbose=verbose, estAllRunsTotalVar=estAllRunsTotalVar, runINLA="rho", runBUGS=FALSE)
        cat("\n")

        # make each rhoFits a dataframe (so that each column is turned to a factor)
        rhoFits <- as.data.frame(rhoFits)

    } else {
        fitRho <- NULL
        rhoFits <- NULL
    }


    ## Abundance mu
    cat("\nAbundance terms:\n")

    # create array to store results of attempted fits
    muFits <- array("", dim=c(0, length(muVars)), dimnames=list(NULL, muVars))

    # initial mu terms
    currentTerms <- character(0)

    # try adding rho terms one by one
    for (i in rep(1:length(muVars), nSweeps)) {
        term <- muVars[i]
        termType <- muVarType[i]
        orders <- NULL

        # set possible orders of term to try
        if (termType == "asis")
            orders <- 1
        else if (termType == "spatial")
            orders <- -1
        else if (termType == "linear") {
            # count number of terms in there already
            termCount <- length(grep(term, currentTerms, fix=TRUE))

            # continue from current place
            if (termCount < maxNoOrders)
                orders <- (termCount + 1):maxNoOrders

        } else if (termType == "continuous") {
            # if rw2 term there, break
            if (length(grep(paste("rw2(", term, ")", sep=""), currentTerms, fix=TRUE)) > 0)
                next

            # count number of terms in there already
            termCount <- length(grep(term, currentTerms, fix=TRUE))

            # continue from current place
            if (termCount < maxNoOrders)
                orders <- (termCount + 1):maxNoOrders

            # if no terms there, try rw2 first
            if (termCount == 0)
                orders <- c(0, orders)

        }

        for (j in orders) {

            # add new term
            if (j == -1)
                term <- paste("spatial(", muVars[i], ", adjacency)", sep="")  ### NOTE: assumes adjacency variable named 'adjacency'
            else if (j == 0)
                term <- paste("rw2(", muVars[i], ")", sep="")
            else if (j == 1)
                term <- muVars[i]
            else
                term <- paste("I(", muVars[i], "^", j, ")", sep="")

            # check that term isn't already there (possible if multiple sweeps)
            if (length(grep(term, currentTerms, fix=TRUE)) > 0)
                break

            # add term
            cat("\n+ Adding term", term, "\n")
            currentTermsBeforeAdding <- currentTerms
            currentTerms <- c(currentTerms, term)

            # add entry in muFits and extract prior parameters
            muFits <- rbind(muFits, "")
            priors <- rhoPriors
            if ("r" %in% names(prior.parameters))
                priors <- c(priors, prior.parameters["r"])
            for (k in 1:length(muVars)) {
                # search for var in current terms
                pos <- grep(muVars[k], currentTerms, fixed=TRUE)
                no <- length(pos)

                if (no > 0) {
                    # add muFits entry
                    if (no == 1 && length(grep("spatial", currentTerms[pos], fixed=TRUE)) > 0)
                        muFits[nrow(muFits), k] <- "spatial"
                    else if (no == 1 && length(grep("rw2", currentTerms[pos], fixed=TRUE)) > 0)
                        muFits[nrow(muFits), k] <- "rw2"
                    else
                        muFits[nrow(muFits), k] <- no

                    # search for prior parameters
                    if (muFits[nrow(muFits), k] == "spatial" || muFits[nrow(muFits), k] == "rw2") {
                        # non-linear term, so search for nu.var or phi.var
                        varName <- paste("sigma", muVars[k], sep=".")
                        if (varName %in% names(prior.parameters))
                            priors <- c(priors, prior.parameters[varName])
                        varName <- paste("tau", muVars[k], sep=".")
                        if (varName %in% names(prior.parameters))
                            priors <- c(priors, prior.parameters[varName])

                    } else {
                        # linear term(s), so search for gamma.term for each term in currentTerms[pos]
                        ## NOTE: this is unnecessary since INLA currently ignores priors for linear variables
                        for (t in currentTerms[pos]) {
                            varName <- paste("beta", t, sep=".")
                            if (varName %in% names(prior.parameters))
                                priors <- c(priors, prior.parameters[varName])
                        }
                    }
                }
            }

            # fit model
            muFormula <- formula(paste("~", paste(currentTerms, collapse=" + ")))
            cat("Fitting model with\n")
            print(formula(paste("log(mu)", paste(deparse(muFormula), collapse=" "))), showEnv=FALSE)
            fit <- fcs2FitModel(runTotalVars, allRunsTotalVar, allRunsRangeVars, muFormula, rhoFormula, dataFrame, surveyAreaVar, nRunsVar, subset,
                                na.omit, priors, verbose=verbose, estAllRunsTotalVar=estAllRunsTotalVar, runINLA=FALSE, runBUGS=FALSE)
            if (class(try({
                if (is.null(fitRho)) {
                    fit <- fcs2FitModel(fit=fit, verbose=verbose, runINLA=TRUE, runBUGS=FALSE)
                    fitRho <- fit
                    fit$Rho$inla$mu <- NULL

                } else if (identical(fit$modelMatrix, fitRho$modelMatrix)) {
                    fit$inlaFits <- fitRho$inlaFits
                    fit <- fcs2FitModel(fit=fit, verbose=verbose, runINLA="mu", runBUGS=FALSE)

                } else
                    fit <- fcs2FitModel(fit=fit, verbose=verbose, runINLA=TRUE, runBUGS=FALSE)

            }, silent=TRUE)) == "try-error") {
                # error fitting model, so reset terms as before added this term
                cat("\n- Removing term ", term, " as error fitting model with INLA\n", sep="")
                currentTerms <- currentTermsBeforeAdding
                if (j > 0)
                    break
                next
            }

            # calculate summary (with 10 attempts as INLA can be tempermental)
            muSum <- NULL
            for (k in 1:10) {
                try({
                    muSum <- summary(fit, prior=FALSE, allVars=TRUE)$inla$mu
                    break
                }, silent=TRUE)
            }
            if (is.null(muSum)) {
                # error calculating summary, so reset terms as before added this term
                cat("\n- Removing term ", term, " as error calculating INLA summary (after 10 attempts)\n", sep="")
                currentTerms <- currentTermsBeforeAdding
                if (j > 0)
                    break
                next
            }

            # extract significance prob of term
            muSignifPr <- muSum[, "Signif. Pr."]
            names(muSignifPr) <- rownames(muSum)
            muSignifPr <- muSignifPr[!is.na(muSignifPr)]  # keep only vars with prob

            # search for multiple vars with each name
            for (k in 1:length(muVars)) {
                pos <- grep(muVars[k], names(muSignifPr), fixed=TRUE)

                # if multiple random, use transformed min
                if (length(pos2 <- grep("[", names(muSignifPr)[pos], fixed=TRUE)) > 0) {
                    pos <- pos[pos2]
                    minPr <- min(muSignifPr[pos])
                    muSignifPr <- muSignifPr[-pos[-1]]
                    muSignifPr[pos[1]] <- 1 - ((1 - minPr) ^ length(pos))  # transform min prob so on uniform scale (assuming independence!)
                    if (muVarType[k] == "continuous")
                        names(muSignifPr)[pos[1]] <- paste("rw2(", muVars[k], ")", sep="")
                    else if (muVarType[k] == "spatial")
                        names(muSignifPr)[pos[1]] <- paste("spatial(", muVars[k], ", adjacency)", sep="")

                    # remove lower orders of terms so don't remove these first
                } else if (length(pos) > 1)
                    muSignifPr <- muSignifPr[-pos[-length(pos)]]
            }

            # remove insignificant terms, in order of least significance
            while (length(muSignifPr) > 0 && max(muSignifPr) > tolerance) {
                # remove term with largest signif prob (which must be > tolerance)
                worstTerm <- sub("beta.", "", names(which.max(muSignifPr)), fixed=TRUE)
                worstTerm <- sub("TRUE", "", worstTerm)  # remove TRUE if present
                worstTerm <- sub("FALSE", "", worstTerm)  # remove TRUE if present
                cat("\n- Removing term ", worstTerm, " as significance probability = ", max(muSignifPr), "\n", sep="")
                currentTerms <- setdiff(currentTerms, worstTerm)

                # fit model (unless first term to remove is term just added, in which case break)
                if (length(currentTerms) > 0 && !(term == worstTerm && (muFits[nrow(muFits) - 1, i] == "" ||
                                                                        (nchar(muFits[nrow(muFits), i]) == 1 && as.numeric(muFits[nrow(muFits), i]) > 1 &&
                                                                         as.numeric(muFits[nrow(muFits), i]) - 1 == as.numeric(muFits[nrow(muFits) - 1, i]))))) {
                    # add entry in muFits and extract prior parameters
                    muFits <- rbind(muFits, "")
                    priors <- rhoPriors
                    if ("r" %in% names(prior.parameters))
                        priors <- c(priors, prior.parameters["r"])
                    for (k in 1:length(muVars)) {
                        # search for var in current terms
                        pos <- grep(muVars[k], currentTerms, fixed=TRUE)
                        no <- length(pos)

                        if (no > 0) {
                            # add muFits entry
                            if (no == 1 && length(grep("spatial", currentTerms[pos], fixed=TRUE)) > 0)
                                muFits[nrow(muFits), k] <- "spatial"
                            else if (no == 1 && length(grep("rw2", currentTerms[pos], fixed=TRUE)) > 0)
                                muFits[nrow(muFits), k] <- "rw2"
                            else
                                muFits[nrow(muFits), k] <- no

                            # search for prior parameters
                            if (muFits[nrow(muFits), k] == "spatial" || muFits[nrow(muFits), k] == "rw2") {
                                # non-linear term, so search for nu.var or phi.var
                                varName <- paste("sigma", muVars[k], sep=".")
                                if (varName %in% names(prior.parameters))
                                    priors <- c(priors, prior.parameters[varName])
                                varName <- paste("tau", muVars[k], sep=".")
                                if (varName %in% names(prior.parameters))
                                    priors <- c(priors, prior.parameters[varName])

                            } else {
                                # linear term(s), so search for gamma.term for each term in currentTerms[pos]
                                ## NOTE: this is unnecessary since INLA currently ignores priors for linear variables
                                for (t in currentTerms[pos]) {
                                    varName <- paste("beta", t, sep=".")
                                    if (varName %in% names(prior.parameters))
                                        priors <- c(priors, prior.parameters[varName])
                                }
                            }
                        }
                    }

                    # fit model
                    muFormula <- formula(paste("~", paste(currentTerms, collapse=" + ")))
                    cat("Fitting model with\n")
                    print(formula(paste("log(mu)", paste(deparse(muFormula), collapse=" "))), showEnv=FALSE)
                    fit <- fcs2FitModel(runTotalVars, allRunsTotalVar, allRunsRangeVars, muFormula, rhoFormula, dataFrame, surveyAreaVar, nRunsVar, subset,
                                        na.omit, priors, verbose=verbose, estAllRunsTotalVar=estAllRunsTotalVar, runINLA=FALSE, runBUGS=FALSE)
                    if (class(try({
                        if (identical(fit$modelMatrix, fitRho$modelMatrix)) {
                            fit$inlaFits <- fitRho$inlaFits
                            fit <- fcs2FitModel(fit=fit, verbose=verbose, runINLA="mu", runBUGS=FALSE)

                        } else
                            fit <- fcs2FitModel(fit=fit, verbose=verbose, runINLA=TRUE, runBUGS=FALSE)

                    }, silent=TRUE)) == "try-error") {
                        # error fitting model, so reset terms as before added this term
                        cat("\n- Removing term ", term, " and its consequences as error fitting model with INLA\n", sep="")
                        currentTerms <- currentTermsBeforeAdding
                        muSignifPr <- NULL
                        next  # note: only breaks out of while loop
                    }

                    # calculate summary (with 10 attempts as INLA can be tempermental)
                    muSum <- NULL
                    for (k in 1:10) {
                        try({
                            muSum <- summary(fit, prior=FALSE, allVars=TRUE)$inla$mu
                            break
                        }, silent=TRUE)
                    }
                    if (is.null(muSum)) {
                        # error calculating summary, so reset terms as before added this term
                        cat("\n- Removing term ", term, " as error calculating INLA summary (after 10 attempts)\n", sep="")
                        currentTerms <- currentTermsBeforeAdding
                        muSignifPr <- NULL
                        next
                    }

                    # extract significance prob of term
                    muSignifPr <- muSum[, "Signif. Pr."]
                    names(muSignifPr) <- rownames(muSum)
                    muSignifPr <- muSignifPr[!is.na(muSignifPr)]  # keep only vars with prob

                    # search for multiple vars with each name
                    for (k in 1:length(muVars)) {
                        pos <- grep(muVars[k], names(muSignifPr), fixed=TRUE)

                        # if multiple random, use transformed min
                        if (length(pos2 <- grep("[", names(muSignifPr)[pos], fixed=TRUE)) > 0) {
                            pos <- pos[pos2]
                            minPr <- min(muSignifPr[pos])
                            muSignifPr <- muSignifPr[-pos[-1]]
                            muSignifPr[pos[1]] <- 1 - ((1 - minPr) ^ length(pos))  # transform min prob so on uniform scale (assuming independence!)
                            if (muVarType[k] == "continuous")
                                names(muSignifPr)[pos[1]] <- paste("rw2(", muVars[k], ")", sep="")
                            else if (muVarType[k] == "spatial")
                                names(muSignifPr)[pos[1]] <- paste("spatial(", muVars[k], ", adjacency)", sep="")

                            # remove lower orders of terms so don't remove these first
                        } else if (length(pos) > 1)
                            muSignifPr <- muSignifPr[-pos[-length(pos)]]
                    }

                } else
                    muSignifPr <- NULL
            }

            # if rw2 and still there, break
            if (j == 0 && term %in% currentTerms)
                break

            # if more orders to try, check this term is still in there and break if not
            if (j > 0 && j < max(orders) && !(term %in% currentTerms))
                break
        }
    }

    # add entry in muFits (unless this is already final entry ie no terms were removed at end)
    finalMuFit <- muFits[nrow(muFits), ]
    finalMuFit[finalMuFit == "spatial" | finalMuFit == "rw2"] <- 1
    finalMuFit <- as.numeric(finalMuFit)
    if (sum(finalMuFit, na.rm=TRUE) != length(currentTerms)) {
        muFits <- rbind(muFits, "")
        for (k in 1:length(muVars)) {
            pos <- grep(muVars[k], currentTerms, fixed=TRUE)
            no <- length(pos)
            if (no == 1 && length(grep("spatial", currentTerms[pos], fixed=TRUE)) > 0)
                muFits[nrow(muFits), k] <- "spatial"
            else if (no == 1 && length(grep("rw2", currentTerms[pos], fixed=TRUE)) > 0)
                muFits[nrow(muFits), k] <- "rw2"
            else if (no > 0)
                muFits[nrow(muFits), k] <- no
        }
    }

    # print best fit
    if (length(currentTerms) == 0)
        currentTerms <- "1"
    muFormula <- formula(paste("~", paste(currentTerms, collapse=" + ")))
    cat("\nBest abundance fit is\n")
    print(formula(paste("log(mu)", paste(deparse(muFormula), collapse=" "))), showEnv=FALSE)
    cat("\n")

    # make each rhoFits a dataframe (so that each column is turned to a factor)
    muFits <- as.data.frame(muFits)

    # invisibly return table of models fitted
    return(list(rhoFits=rhoFits, muFits=muFits))
}

