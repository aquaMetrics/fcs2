#' @import stats
#' @import utils
#' @import INLA
#' @import R2OpenBUGS

# Internal (undocumented) functions that should not typically be called by the user
# Extracts or calculates the total number of fish caught over all runs
.allRunsTotal <-
function(dataFrame, runTotalVars = NULL, allRunsTotalVar = NULL, allRunsRangeVars = NULL, nRunsVar = NULL, subset = 1:nrow(dataFrame), allowRange = FALSE)
{
    if (class(subset) == "logical" || min(subset) == 0)
        subset <- which(as.logical(subset))
    ret <- rep(NA, length(subset))
    names(ret) <- rownames(dataFrame)[subset]

    # run total vars
    if (!is.null(runTotalVars) && runTotalVars[1] %in% names(dataFrame)) {
        # if nRunsVar present, clear missing entries in data frame
        if (!is.null(nRunsVar) && nRunsVar %in% names(dataFrame)) {
            # crop nRuns var if necessary
            if (sum(dataFrame[subset, nRunsVar] > length(runTotalVars)) > 0)
                dataFrame[subset[dataFrame[subset, nRunsVar] > length(runTotalVars)], nRunsVar] <- length(runTotalVars)

            # zero entries
            zero <- array(FALSE, dim=c(length(subset), length(runTotalVars)))
            for (j in 1:ncol(zero))
                zero[, j] <- j > dataFrame[subset, nRunsVar]
            dataFrame[subset, runTotalVars][zero] <- 0

            # sum over runs but don't remove NAs
            ret <- apply(dataFrame[subset, runTotalVars, drop=FALSE], 1, sum, na.rm=FALSE)

        } else {
            # nRunsVar not present so assume NAs were not runs
            ret <- apply(dataFrame[subset, runTotalVars, drop=FALSE], 1, sum, na.rm=TRUE)

            # make sure NA if all runs missing
            ret[is.na(dataFrame[subset, runTotalVars[1]])] <- NA
        }
    }

    # all runs total var
    if (!is.null(allRunsTotalVar) && allRunsTotalVar %in% names(dataFrame))
        ret[is.na(ret)] <- dataFrame[subset[is.na(ret)], allRunsTotalVar]

    # all runs total range
    if (!is.null(allRunsRangeVars) && allRunsRangeVars[1] %in% names(dataFrame)) {
        # fill in entries where min = max
        i <- which(is.na(ret) & !is.na(dataFrame[subset, allRunsRangeVars[1]]) &
                   dataFrame[subset, allRunsRangeVars[1]] == dataFrame[subset, allRunsRangeVars[2]])
        if (length(i) > 0)
            ret[i] <- dataFrame[subset[i], allRunsRangeVars[1]]

        # if other range observations, give range as min-max
        i <- which(is.na(ret) & !is.na(dataFrame[subset, allRunsRangeVars[1]]))
        if (length(i) > 0 && allowRange) {
            ret <- as.character(ret) # convert ret to character
            ret[i] <- paste(dataFrame[subset[i], allRunsRangeVars[1]], "to", dataFrame[subset[i], allRunsRangeVars[2]], sep=" ")
        }
    }

    ret
}


## Corrects BUGS output by renaming elements of 'bugs' object
.fcs2CorrectBUGSOutput <-
function(fit)
{
    # get variable names
    correctNames <- c(variable.names(fit, hyperparams="precision"), "deviance")
    fixedNames <- make.names(correctNames)
    correctNamesAll <- c(variable.names(fit, rw="singular", spatial="singular", hyperparams="precision"), "deviance")
    fixedNamesAll <- make.names(correctNamesAll)

    # correct variable names
    dimnames(fit$bugsFit$sims.array)[[3]] <- correctNamesAll[match(make.names(dimnames(fit$bugsFit$sims.array)[[3]]), fixedNamesAll)]
    colnames(fit$bugsFit$sims.matrix) <- correctNamesAll[match(make.names(colnames(fit$bugsFit$sims.matrix)), fixedNamesAll)]
    rownames(fit$bugsFit$summary) <- correctNamesAll[match(make.names(rownames(fit$bugsFit$summary)), fixedNamesAll)]
    names(fit$bugsFit$sims.list) <- names(fit$bugsFit$mean) <- names(fit$bugsFit$median) <- names(fit$bugsFit$sd) <- fit$bugsFit$root.short <-
            correctNames[match(make.names(names(fit$bugsFit$mean)), fixedNames)]

    # correct variable names in list of last values
    for (i in 1:length(fit$bugsFit$last.values))
        names(fit$bugsFit$last.values[[i]]) <- correctNames[match(make.names(names(fit$bugsFit$last.values[[i]])), fixedNames)]

    fit$bugsFit
}


## Creates data required by BUGS
.fcs2CreateBUGSData <-
function(fit)
{
    # remove NA in model matrix (which may not have happened yet if set na.action=na.pass)
    modelMatrix <- na.omit(fit$modelMatrix)
    na.action <- na.action(modelMatrix)
    if (is.null(na.action))
        dataType <- fit$dataType
    else
        dataType <- fit$dataType[-na.action]

    # enter allRunTotal data in allRunMin/Max
    if (is.null(allRunsRangeVars <- fit$allRunsRangeVars)) {
        allRunsRangeVars <- paste("AllRunsTotal", c("Min", "Max"), sep="")
        modelMatrix <- cbind(modelMatrix, AllRunsTotalMin=NA, AllRunsTotalMax=NA)
    }
    if (sum(dataType == "total") > 0)
        modelMatrix[dataType == "total", allRunsRangeVars] <- modelMatrix[dataType == "total", fit$allRunsTotalVar]

    # reorder model matrix so that runs data first
    o <- order(dataType)
    modelMatrix <- modelMatrix[o, ]

    # start with required variables
    bugsData <- list(N=nrow(modelMatrix),
                     area=modelMatrix[, fit$surveyAreaVar])

    # if both total and range data, add 'Nruns'
    if (sum(dataType == "run") > 0 && sum(dataType == "run") < nrow(modelMatrix))
        bugsData$Nruns <- sum(dataType == "run")

    # add catch total data
    if (sum(dataType == "run") > 0) {
        if (!fit$multiRun)
            bugsData$catch <- modelMatrix[, fit$runTotalVars[1]]
        else
            bugsData$catch <- as.matrix(modelMatrix[, fit$runTotalVars])

        # add additional variables for multiRun
        if (fit$multiRun) {
            bugsData$maxNRuns <- max(modelMatrix[, fit$nRuns])
            bugsData$catchSum <- modelMatrix[, fit$allRunsTotalVar]
        }
    }

    # add ranged catch data
    if (sum(dataType != "run") > 0) {
        bugsData$totalCatchMin <- modelMatrix[, allRunsRangeVars[1]]
        bugsData$totalCatchMax <- modelMatrix[, allRunsRangeVars[2]]
    }

    # add no runs if multiRun
    if (fit$multiRun)
        bugsData$nRuns <- modelMatrix[, fit$nRuns]


    # add covariates
    for (var in c(fit$muLinearVars, fit$muRW1Vars, fit$muRW2Vars, fit$muSpatialVar,
                  fit$rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars, fit$rhoSpatialVar)) {
        if (!(var %in% names(bugsData))) {
            bugsData <- c(bugsData, list(modelMatrix[, var]))
            names(bugsData)[length(bugsData)] <- make.names(var)
        }
    }


    ## add boundaries and weights etc for mu RW terms

    # Weight and adjacency matrix corresponding to RW(1) prior
    for (var in fit$muRW1Vars) {
        d <- length(fit$rwBoundaries[[paste("beta.", var, sep='')]])
        weights <- adj <- rep(NA, (d-2)*2 + 2)
        num <- rep(NA, d)

        weights[1] <- 1; adj[1] <- 1+1; num[1] <- 1

        for(t in 2:(d-1)) {
            weights[2+(t-2)*2] <- 1; adj[2+(t-2)*2] <- t-1
            weights[3+(t-2)*2] <- 1; adj[3+(t-2)*2] <- t+1; num[t] <- 2
        }

        weights[(d-2)*2 + 2] <- 1; adj[(d-2)*2 + 2] <- d-1; num[d] <- 1

        bugsData <- c(bugsData, list(fit$rwBoundaries[[paste("beta.", var, sep='')]], adj, weights, num, d))
        names(bugsData)[length(bugsData) - 4:0] <- paste(c("boundaries", "adj", "weights", "num", "nLevels"), ".beta.", make.names(var), sep='')
    }

    # Weight and adjacency matrix corresponding to RW(2) prior
    for (var in fit$muRW2Vars) {
        d <- length(fit$rwBoundaries[[paste("beta.", var, sep='')]])
        weights <- adj <- rep(NA, (d-4) * 4 + 10)
        num <- rep(NA, d)

        weights[1] <- 2;  adj[1] <- 2
        weights[2] <- -1; adj[2] <- 3; num[1] <- 2
        weights[3] <- 2;  adj[3] <- 1
        weights[4] <- 4;  adj[4] <- 3
        weights[5] <- -1; adj[5] <- 4; num[2] <- 3

        for (t in 3:(d-2)) {
            weights[6+(t-3)*4] <- -1; adj[6+(t-3)*4] <- t-2
            weights[7+(t-3)*4] <- 4;  adj[7+(t-3)*4] <- t-1
            weights[8+(t-3)*4] <- 4;  adj[8+(t-3)*4] <- t+1
            weights[9+(t-3)*4] <- -1; adj[9+(t-3)*4] <- t+2; num[t] <- 4
        }

        weights[(d-4)*4 + 6] <- 2;   adj[(d-4)*4 + 6] <- d
        weights[(d-4)*4 + 7] <- 4;   adj[(d-4)*4 + 7] <- d-2
        weights[(d-4)*4 + 8] <- -1;  adj[(d-4)*4 + 8] <- d-3;  num[d-1] <- 3
        weights[(d-4)*4 + 9] <- 2;   adj[(d-4)*4 + 9] <- d-1
        weights[(d-4)*4 + 10] <- -1; adj[(d-4)*4 + 10] <- d-2; num[d] <- 2

        bugsData <- c(bugsData, list(fit$rwBoundaries[[paste("beta.", var, sep='')]], adj, weights, num, d))
        names(bugsData)[length(bugsData) - 4:0] <- paste(c("boundaries", "adj", "weights", "num", "nLevels"), ".beta.", make.names(var), sep='')
    }


    ## add boundaries and weights etc for rho RW terms

    # Weight and adjacency matrix corresponding to RW(1) prior
    for (var in fit$rhoRW1Vars) {
        d <- length(fit$rwBoundaries[[paste("gamma.", var, sep='')]])
        weights <- adj <- rep(NA, (d-2)*2 + 2)
        num <- rep(NA, d)

        weights[1] <- 1; adj[1] <- 1+1; num[1] <- 1

        for(t in 2:(d-1)) {
           weights[2+(t-2)*2] <- 1; adj[2+(t-2)*2] <- t-1
           weights[3+(t-2)*2] <- 1; adj[3+(t-2)*2] <- t+1; num[t] <- 2
        }

        weights[(d-2)*2 + 2] <- 1; adj[(d-2)*2 + 2] <- d-1; num[d] <- 1

        bugsData <- c(bugsData, list(fit$rwBoundaries[[paste("gamma.", var, sep='')]], adj, weights, num, d))
        names(bugsData)[length(bugsData) - 4:0] <- paste(c("boundaries", "adj", "weights", "num", "nLevels"), ".gamma.", make.names(var), sep='')
    }

    # Weight and adjacency matrix corresponding to RW(2) prior
    for (var in fit$rhoRW2Vars) {
        d <- length(fit$rwBoundaries[[paste("gamma.", var, sep='')]])
        weights <- adj <- rep(NA, (d-4) * 4 + 10)
        num <- rep(NA, d)

        weights[1] <- 2;  adj[1] <- 2
        weights[2] <- -1; adj[2] <- 3; num[1] <- 2
        weights[3] <- 2;  adj[3] <- 1
        weights[4] <- 4;  adj[4] <- 3
        weights[5] <- -1; adj[5] <- 4; num[2] <- 3
        for (t in 3:(d-2)) {
            weights[6+(t-3)*4] <- -1; adj[6+(t-3)*4] <- t-2
            weights[7+(t-3)*4] <- 4;  adj[7+(t-3)*4] <- t-1
            weights[8+(t-3)*4] <- 4;  adj[8+(t-3)*4] <- t+1
            weights[9+(t-3)*4] <- -1; adj[9+(t-3)*4] <- t+2; num[t] <- 4
        }
        weights[(d-4)*4 + 6] <- 2;   adj[(d-4)*4 + 6] <- d
        weights[(d-4)*4 + 7] <- 4;   adj[(d-4)*4 + 7] <- d-2
        weights[(d-4)*4 + 8] <- -1;  adj[(d-4)*4 + 8] <- d-3;  num[d-1] <- 3
        weights[(d-4)*4 + 9] <- 2;   adj[(d-4)*4 + 9] <- d-1
        weights[(d-4)*4 + 10] <- -1; adj[(d-4)*4 + 10] <- d-2; num[d] <- 2

        bugsData <- c(bugsData, list(fit$rwBoundaries[[paste("gamma.", var, sep='')]], adj, weights, num, d))
        names(bugsData)[length(bugsData) - 4:0] <- paste(c("boundaries", "adj", "weights", "num", "nLevels"), ".gamma.", make.names(var), sep='')
    }

    # add adjacency info for spatial terms
    if (!is.null(fit$muSpatialVar)) {
        bugsData <- c(bugsData, list(fit$muAdjacency$adj, rep(1, length(fit$muAdjacency$adj)), fit$muAdjacency$num, length(fit$muAdjacency$num)))
        names(bugsData)[length(bugsData) - 3:0] <- paste(c("adj", "weights", "num", "nSpatials"), ".", make.names(fit$muSpatialVar), sep='')
    }
    if (!is.null(fit$rhoSpatialVar) && (is.null(fit$muSpatialVar) || fit$muSpatialVar != fit$rhoSpatialVar)) {
        bugsData <- c(bugsData, list(fit$rhoAdjacency$adj, rep(1, length(fit$rhoAdjacency$adj)), fit$rhoAdjacency$num, length(fit$rhoAdjacency$num)))
        names(bugsData)[length(bugsData) - 3:0] <- paste(c("adj", "weights", "num", "nSpatials"), ".", make.names(fit$rhoSpatialVar), sep='')
    }


    ## add priors

    # only add q.a and q.b if not both 1 (giving uniform)
    if (fit$multiRun && !(fit$prior.parameters[["q"]][1] == 1 && fit$prior.parameters[["q"]][2] == 1))
        bugsData <- c(bugsData, list(q.a=fit$prior.parameters[["q"]][1], q.b=fit$prior.parameters[["q"]][2]))

    bugsData <- c(bugsData, list(mu.r=fit$prior.parameters[["r"]][1], tau.r=fit$prior.parameters[["r"]][2]))

    for (var in c(grep("beta", names(fit$prior.parameters), value=TRUE), grep("gamma", names(fit$prior.parameters), value=TRUE))) {
        bugsData <- c(bugsData, list(fit$prior.parameters[[var]][1], fit$prior.parameters[[var]][2]))
        names(bugsData)[length(bugsData) - 1:0] <- paste(c("mean", "prec"), ".", make.names(var), sep='')
    }

    for (var in c(grep("tau", names(fit$prior.parameters), value=TRUE), grep("phi", names(fit$prior.parameters), value=TRUE))) {
        bugsData <- c(bugsData, list(fit$prior.parameters[[var]][1], fit$prior.parameters[[var]][2]))
        names(bugsData)[length(bugsData) - 1:0] <- paste(c("a", "b"), ".", make.names(var), sep='')
    }

    # return
    bugsData
}


## Fills in missing initial values with median estimates from INLA fits
.fcs2SetInitialValues <-
function(fit)
{
    inits <- fit$initial.values

    if (!is.null(fit$inlaFits)) {
        ## Set initial values using INLA fits to prevalence and abundance components
        ### NOTE: extraction of medians relies on default settings for INLA

        # catch probability q
        if (fit$multiRun)
            inits$q <- fit$prior.parameters[["q"]][1] / (fit$prior.parameters[["q"]][1] + fit$prior.parameters[["q"]][2])


        if (!is.null(fit$inlaFits$muFit)) {
            # shape parameter r
            if (!("r" %in% names(inits)))
                inits <- c(inits, list(r=fit$inlaFits$muFit$summary.hyperpar[pmatch("Size", rownames(fit$inlaFits$muFit$summary.hyperpar)), 4]))

            # mu constant term
            if (!("beta.const" %in% names(inits)))
                inits <- c(inits, list(beta.const=fit$inlaFits$muFit$summary.fixed["(Intercept)", 4]))

            # mu linear terms
            for (var in fit$muLinearVars) {
                varLong <- paste("beta", var, sep='.')
                if (!(varLong %in% names(inits))) {
                    inits <- c(inits, list(fit$inlaFits$muFit$summary.fixed[make.names(var), 4]))
                    names(inits)[length(inits)] <- varLong
                }
            }

            # mu RW and spatial terms
            for (var in c(fit$muRW1Vars, fit$muRW2Vars, fit$muSpatialVar)) {
                varLong <- paste(c("beta", "tau"), var, sep='.')
                if (!(varLong[1] %in% names(inits))) {
                    vals <- fit$inlaFits$muFit$summary.random[[make.names(var)]][, 5]
                    vals <- vals - sum(vals) / length(vals)  # shift so that sum to 0
                    vals[1] <- -sum(vals[-1])  # tiny tweak to be sure sums to 0
                    inits <- c(inits, list(vals))
                    names(inits)[length(inits)] <- varLong[1]
                }
                if (!(varLong[2] %in% names(inits))) {
                    inits <- c(inits, list(fit$inlaFits$muFit$summary.hyperpar[paste("Precision for ", make.names(var), sep=''), 4]))
                    names(inits)[length(inits)] <- varLong[2]
                }
            }
        }

        if (!is.null(fit$inlaFits$rhoFit)) {
            # rho constant term
            if (!("gamma.const" %in% names(inits)))
                inits <- c(inits, list(gamma.const=fit$inlaFits$rhoFit$summary.fixed["(Intercept)", 4]))

            # rho linear terms
            for (var in fit$rhoLinearVars) {
                varLong <- paste("gamma", var, sep='.')
                if (!(varLong %in% names(inits))) {
                    inits <- c(inits, list(fit$inlaFits$rhoFit$summary.fixed[make.names(var), 4]))
                    names(inits)[length(inits)] <- varLong
                }
            }

            # rho RW and spatial terms
            for (var in c(fit$rhoRW1Vars, fit$rhoRW2Vars, fit$rhoSpatialVar)) {
                varLong <- paste(c("gamma", "phi"), var, sep='.')
                if (!(varLong[1] %in% names(inits))) {
                    vals <- fit$inlaFits$rhoFit$summary.random[[make.names(var)]][, 5]
                    vals <- vals - sum(vals) / length(vals)  # shift so that sum to 0
                    vals[1] <- -sum(vals[-1])  # tiny tweak to be sure sums to 0
                    inits <- c(inits, list(vals))
                    names(inits)[length(inits)] <- varLong[1]
                }
                if (!(varLong[2] %in% names(inits))) {
                    inits <- c(inits, list(fit$inlaFits$rhoFit$summary.hyperpar[paste("Precision for ", make.names(var), sep=''), 4]))
                    names(inits)[length(inits)] <- varLong[2]
                }
            }
        }

    } else {
        warning("Initial values not set")
    }

    inits
}


## Creates missing random walk boundaries from specifications of the number of levels
.fcs2SetRWBoundaries <-
function(fit)
{
  with(fit, {

    if (length(rwNoLevels) > 0) {
        for (i in 1:length(rwNoLevels)) {
            if (!(names(rwNoLevels)[i] %in% names(rwBoundaries))) {
                varName <- sub("beta.", "", sub("gamma.", "", names(rwNoLevels[i])))
                rwBoundaries <- c(rwBoundaries, list(seq(min(modelMatrix[, varName], na.rm=TRUE),
                                                         max(modelMatrix[, varName], na.rm=TRUE), length=rwNoLevels[[i]])))
                names(rwBoundaries)[length(rwBoundaries)] <- names(rwNoLevels[i])
            }
        }
    } else {
        for (var in c(muRW1Vars, muRW2Vars)) {
            varName <- paste("beta.", var, sep='')
            if (!(varName %in% names(rwBoundaries))) {
                rwBoundaries <- c(rwBoundaries, list(seq(min(modelMatrix[, var], na.rm=TRUE),
                                                         max(modelMatrix[, var], na.rm=TRUE), length=rwNoLevels)))
                names(rwBoundaries)[length(rwBoundaries)] <- varName
            }
        }
        for (var in c(rhoRW1Vars, rhoRW2Vars)) {
            varName <- paste("gamma.", var, sep='')
            if (!(varName %in% names(rwBoundaries))) {
                rwBoundaries <- c(rwBoundaries, list(seq(min(modelMatrix[, var], na.rm=TRUE),
                                                         max(modelMatrix[, var], na.rm=TRUE), length=rwNoLevels)))
                names(rwBoundaries)[length(rwBoundaries)] <- varName
            }
        }
    }

    rwBoundaries
  })
}


# ## Function called when 'fcs2' package is loaded.
# ## Attaches compiled code for calculating joint EQRs
# .First.lib <-
# function(libname, pkgname)
# {
#     library.dynam("fcs2", package=pkgname, lib.loc=libname)
# }


## Function called when unloading 'fcs2' package.
## Detaches compiled code
#.onUnload <- function(libpath)
#{
#    library.dynam.unload("fcs2", libpath)
#}


## Calculates the significance probability from BUGS samples
.significanceProbability.bugs <-
function(samples)
{
    prob <- mean(samples < 0)
    2 * min(prob, 1 - prob)
}


## Calculates the significance probability from INLA marginal estimate
.significanceProbability.inla <-
function(marginal)
{
    prob <- inla.pmarginal(0, marginal)
    2 * min(prob, 1 - prob)
}


## Produces significance stars from probabilities
.significanceStars <-
function(prob)
{
    sigStars <- character(length(prob))

    for (i in 1:length(prob)) {
        if (is.na(prob[i]))
            sigStars[i] <- "   "
        else if (prob[i] <= 0.001)
            sigStars[i] <- "***"
        else if (prob[i] <= 0.01)
            sigStars[i] <- "** "
        else if (prob[i] <= 0.05)
            sigStars[i] <- "*  "
        else if (prob[i] <= 0.1)
            sigStars[i] <- ".  "
        else
            sigStars[i] <- "   "
    }

    sigStars
}





