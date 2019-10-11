#' Posterior Samples of Abundance and Prevalence Model Components
#'
#' Calculates posterior samples of the abundance or prevalence model components
#'
#' @name abundance
#' @aliases abundance prevalence
#' @usage abundance(fit, newData, subset = 1:nrow(newData), na.action, breakdown = FALSE)
#' @usage prevalence(fit, newData, subset = 1:nrow(newData), na.action, breakdown = FALSE)
#' @param fit an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#'   fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS =
#'   TRUE}.
#' @param newData a data frame with surveys as rows and variables as columns. It
#'   should contain all variables required by \code{fit}. If missing, the model
#'   matrix contained within \code{fit} is used.
#' @param subset an optional vector specifying a subset of surveys to calculate
#'   samples for.
#' @param na.action a function which indicates what should happen when the data
#'   contain missing values (\code{NA}s). The default is set by the
#'   \code{na.action} setting of \code{\link{options}} and this is usually set
#'   to \code{\link{na.omit}}. This setting removes surveys that contain missing
#'   data in any required variables. A vector indicating the rows that were
#'   removed can be extracted from the returned object using
#'   \code{\link{na.action}}. Alternatively, \code{\link{na.pass}} can be used
#'   to ignore missing values (where possible) or \code{\link{na.fail}} can be
#'   given to signal an error if missing values are found.
#' @param breakdown logical; If \code{FALSE} (the default), samples of the
#'   abundance or prevalence are returned for each survey. If \code{TRUE}, an
#'   array is returned contining each additive term in the abundance or
#'   prevalence regression equation. If \code{breakdown = FALSE} (the default),
#'   a matrix of posterior samples of the abundance \eqn{\mu} or prevalence
#'   \eqn{\rho} model component is returned for the selected surveys. Monte
#'   Carlo samples appear as rows and surveys as columns.
#'
#'   If \code{breakdown = TRUE}, a three-dimensional array of regression
#'   component samples is returned. Samples are given in the first dimension,
#'   surveys in the second and the third dimension corresponds to each additive
#'   term in the abundance or prevalence regression. The output from
#'   \code{breakdown = FALSE} can be calculated by summing over the third
#'   dimension and transforming with \code{\link{exp}} for abundance or
#'   \code{\link{expit}} for prevalence.
#'
#'   If surveys were removed by \code{na.action = na.omit}, these can be
#'   recovered by applying \code{\link{na.action}} to the matrix or array.
#' @seealso \code{\link{fcs2FitModel}} for producing the required \acronym{FCS2}
#'   model fit.\cr
#' @export

abundance <-
function(fit, newData, subset = 1:nrow(newData), na.action, breakdown = FALSE)
{
    # check fit contains bugs fit
    if (is.null(fit$bugsFit))
        stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

    # get default na.action if missing
    if (missing(na.action)) {
        na.action <- getOption("na.action")
        if (is.null(na.action))
            na.action <- na.omit  # use 'na.omit' if option not set
        else
            na.action <- eval(parse(text=na.action))
    }

    if (missing(newData)) {
        ## use model matrix from fit
        modelMatrix <- fit$modelMatrix
        na.action.data <- fit$na.action

    } else {
        ## create model matrix from newData, subset and formulae

        # correct factors in linear terms
        muLinearVars <- fit$muLinearVars
        if (length(muLinearVars) > 0) {
            for (i in 1:length(muLinearVars))
                if (length(grep("factor", muLinearVars[i])) > 0) {
                    # temporarily swap 'factor(var)level' with 'I(var == "level")' and swap back later
                    ## NOTE: will not work if factor contained additional arguments! - could be fixed by parsing expression and splitting terms but hard!
                    var <- substr(muLinearVars[i], regexpr("(", muLinearVars[i], fixed=TRUE) + 1, regexpr(")", muLinearVars[i], fixed=TRUE) - 1)
                    level <- substr(muLinearVars[i], regexpr(")", muLinearVars[i], fixed=TRUE) + 1, nchar(muLinearVars[i]))
                    muLinearVars[i] <- paste('I(', var, ' == "', level, '")', sep='')
                }
        } else
            muLinearVars <- "1"

        # remove TRUE or FALSE
        muLinearVars <- sub("TRUE", "", muLinearVars)
        muLinearVars <- sub("FALSE", "", muLinearVars)

        # calculate new model matrix
        simpleFormula <- formula(paste("~", paste(c(muLinearVars, fit$muRW1Vars, fit$muRW2Vars, fit$muSpatialVar), collapse="+")))
        modelFrame <- model.frame(simpleFormula, newData, subset=subset, na.action=na.action)
        modelMatrix <- model.matrix(simpleFormula, modelFrame)[, -1, drop=FALSE]
        if (is.null(rownames(modelMatrix)))
            rownames(modelMatrix) <- rownames(modelFrame)
        na.action.data <- attr(modelFrame, "na.action")

        # replace factor variable names
        if (length(muLinearVars) > 0) {
            for (i in 1:length(muLinearVars))
                if (length(grep("factor", fit$muLinearVars[i])) > 0)
                    colnames(modelMatrix)[grep(muLinearVars[i], colnames(modelMatrix), fixed=TRUE)] <- fit$muLinearVars[i]
        }
    }

    if (breakdown) {
        # store each term separately

        # create empty array to store samples at each survey of each term
        ret <- array(NA, dim=c(fit$bugsFit$n.sims, nrow(modelMatrix),
                               1 + length(fit$muLinearVars) + length(fit$muRW1Vars) + length(fit$muRW2Vars) + length(fit$muSpatialVar)),
                     dimnames=list(NULL, rownames(modelMatrix), make.unique(c("const", fit$muLinearVars, fit$muRW1Vars, fit$muRW2Vars, fit$muSpatialVar))))
        k <- 1

        # check size of model matrix
        if (nrow(modelMatrix) > 0) {

            ## calculate each abundance term
            l <- fit$bugsFit$sims.list

            # constant term
            ret[, , k] <- matrix(l$beta.const, fit$bugsFit$n.sims, nrow(modelMatrix), byrow=FALSE, dimnames=list(NULL, rownames(modelMatrix)))
            k <- k + 1

            # linear terms
            for (var in fit$muLinearVars) {
                varLong <- paste("beta", var, sep=".")
                ret[, , k] <- tcrossprod(l[[varLong]], modelMatrix[, var])
                k <- k + 1
            }

            # random walk terms
            for (var in c(fit$muRW1Vars, fit$muRW2Vars)) {
                varLong <- paste("beta", var, sep=".")
                for (i in 1:nrow(ret))
                    ret[i, , k] <- approx(fit$rwBoundaries[[varLong]], l[[varLong]][i, ], modelMatrix[, var])$y
                k <- k + 1
            }

            # spatial term
            for (var in fit$muSpatialVar) {
                varLong <- paste("beta", var, sep=".")
                ret[, , k] <- l[[varLong]][, modelMatrix[, var]]
                k <- k + 1
            }

        }

        # remove missing values and add to na.action
        ret[!is.finite(ret)] <- NA  # make infinite values NA
        retMean <- apply(ret, 2:3, mean)  # find mean value over samples
        retMean <- na.action(retMean)
        na.action.ret <- attr(retMean, "na.action")
        if (length(na.action.ret) > 0)
            ret <- ret[-na.action.ret, , ]  # apply na.action to ret
        if (is.null(na.action.data))
            modelMatrixIndices <- 1:nrow(modelMatrix)
        else
            modelMatrixIndices <- (1:(length(na.action.data) + nrow(modelMatrix)))[-na.action.data]
        na.action.ret[] <- modelMatrixIndices[na.action.ret]
        na.action <- c(na.action.ret, na.action.data)
        if (!is.null(na.action)) {
            na.action <- sort(na.action)
            class(na.action) <- ifelse(is.null(na.action.ret), class(na.action.data), class(na.action.ret))
            attr(ret, "na.action") <- na.action
        }

    } else {
        # store summed terms only

        # check size of model matrix
        if (nrow(modelMatrix) > 0) {

            ## calculate abundance
            l <- fit$bugsFit$sims.list

            # constant term
            ret <- matrix(l$beta.const, fit$bugsFit$n.sims, nrow(modelMatrix), byrow=FALSE, dimnames=list(NULL, rownames(modelMatrix)))

            # linear terms
            for (var in fit$muLinearVars) {
                varLong <- paste("beta", var, sep=".")
                ret <- ret + tcrossprod(l[[varLong]], modelMatrix[, var])
            }

            # random walk terms
            for (var in c(fit$muRW1Vars, fit$muRW2Vars)) {
                varLong <- paste("beta", var, sep=".")
                for (i in 1:nrow(ret))
                    ret[i, ] <- ret[i, ] + approx(fit$rwBoundaries[[varLong]], l[[varLong]][i, ], modelMatrix[, var])$y
            }

            # spatial term
            for (var in fit$muSpatialVar) {
                varLong <- paste("beta", var, sep=".")
                ret <- ret + l[[varLong]][, modelMatrix[, var]]
            }

            # transform
            ret <- exp(ret)

        } else
            ret <- array(dim=c(fit$bugsFit$n.sims, 0))

        # remove missing values and add to na.action
        ret[!is.finite(ret)] <- NA  # make infinite values NA
        ret <- t(na.action(t(ret)))
        na.action.ret <- attr(ret, "na.action")
        if (is.null(na.action.data))
            modelMatrixIndices <- 1:nrow(modelMatrix)
        else
            modelMatrixIndices <- (1:(length(na.action.data) + nrow(modelMatrix)))[-na.action.data]
        na.action.ret[] <- modelMatrixIndices[na.action.ret]
        na.action <- c(na.action.ret, na.action.data)
        if (!is.null(na.action)) {
            na.action <- sort(na.action)
            class(na.action) <- ifelse(is.null(na.action.ret), class(na.action.data), class(na.action.ret))
            attr(ret, "na.action") <- na.action
        }

    }

    # return
    ret
}


#' @rdname abundance
#' @export
prevalence <-
function(fit, newData, subset = 1:nrow(newData), na.action, breakdown = FALSE)
{
    # check fit contains bugs fit
    if (is.null(fit$bugsFit))
        stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

    # get default na.action if missing
    if (missing(na.action)) {
        na.action <- getOption("na.action")
        if (is.null(na.action))
            na.action <- na.omit  # use 'na.omit' if option not set
        else
            na.action <- eval(parse(text=na.action))
    }

    if (missing(newData)) {
        ## use model matrix from fit
        modelMatrix <- fit$modelMatrix
        na.action.data <- fit$na.action

    } else {
        ## create model matrix from newData, subset and formulae

        # correct factors in linear terms
        rhoLinearVars <- fit$rhoLinearVars
        if (length(rhoLinearVars) > 0) {
            for (i in 1:length(rhoLinearVars))
                if (length(grep("factor", rhoLinearVars[i])) > 0) {
                    # temporarily swap 'factor(var)level' with 'I(var == "level")' and swap back later
                    ## NOTE: will not work if factor contained additional arguments! - could be fixed by parsing expression and splitting terms but hard!
                    var <- substr(rhoLinearVars[i], regexpr("(", rhoLinearVars[i], fixed=TRUE) + 1, regexpr(")", rhoLinearVars[i], fixed=TRUE) - 1)
                    level <- substr(rhoLinearVars[i], regexpr(")", rhoLinearVars[i], fixed=TRUE) + 1, nchar(rhoLinearVars[i]))
                    rhoLinearVars[i] <- paste('I(', var, ' == "', level, '")', sep='')
                }
        } else
            rhoLinearVars <- "1"

        # remove TRUE or FALSE
        rhoLinearVars <- sub("TRUE", "", rhoLinearVars)
        rhoLinearVars <- sub("FALSE", "", rhoLinearVars)

        # calculate new model matrix
        simpleFormula <- formula(paste("~", paste(c(rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars, fit$rhoSpatialVar), collapse="+")))
        modelFrame <- model.frame(simpleFormula, newData, subset=subset, na.action=na.action)
        modelMatrix <- model.matrix(simpleFormula, modelFrame)[, -1, drop=FALSE]
        if (is.null(rownames(modelMatrix)))
            rownames(modelMatrix) <- rownames(modelFrame)
        na.action.data <- attr(modelFrame, "na.action")

        # replace factor variable names
        if (length(rhoLinearVars) > 0) {
            for (i in 1:length(rhoLinearVars))
                if (length(grep("factor", fit$rhoLinearVars[i])) > 0)
                    colnames(modelMatrix)[grep(rhoLinearVars[i], colnames(modelMatrix), fixed=TRUE)] <- fit$rhoLinearVars[i]
        }
    }

    if (breakdown) {
        # store each term separately

        # create empty array to store samples at each survey of each term
        ret <- array(NA, dim=c(fit$bugsFit$n.sims, nrow(modelMatrix),
                               1 + length(fit$rhoLinearVars) + length(fit$rhoRW1Vars) + length(fit$rhoRW2Vars) + length(fit$rhoSpatialVar)),
                     dimnames=list(NULL, rownames(modelMatrix), make.unique(c("const", fit$rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars, fit$rhoSpatialVar))))
        k <- 1

        # check size of model matrix
        if (nrow(modelMatrix) > 0) {

            ## calculate abundance
            l <- fit$bugsFit$sims.list

            # constant term
            ret[, , k] <- matrix(l$gamma.const, fit$bugsFit$n.sims, nrow(modelMatrix), byrow=FALSE, dimnames=list(NULL, rownames(modelMatrix)))
            k <- k + 1

            # linear terms
            for (var in fit$rhoLinearVars) {
                varLong <- paste("gamma", var, sep=".")
                ret[, , k] <- tcrossprod(l[[varLong]], modelMatrix[, var])
                k <- k + 1
            }

            # random walk terms
            for (var in c(fit$rhoRW1Vars, fit$rhoRW2Vars)) {
                varLong <- paste("gamma", var, sep=".")
                for (i in 1:nrow(ret))
                    ret[i, , k] <- approx(fit$rwBoundaries[[varLong]], l[[varLong]][i, ], modelMatrix[, var])$y
                k <- k + 1
            }

            # spatial term
            for (var in fit$rhoSpatialVar) {
                varLong <- paste("gamma", var, sep=".")
                ret[, , k] <- l[[varLong]][, modelMatrix[, var]]
                k <- k + 1
            }

        }

        # remove missing values and add to na.action
        ret[!is.finite(ret)] <- NA  # make infinite values NA
        retMean <- apply(ret, 2:3, mean)  # find mean value over samples
        retMean <- na.action(retMean)
        na.action.ret <- attr(retMean, "na.action")
        if (length(na.action.ret) > 0)
            ret <- ret[-na.action.ret, , ]  # apply na.action to ret
        if (is.null(na.action.data))
            modelMatrixIndices <- 1:nrow(modelMatrix)
        else
            modelMatrixIndices <- (1:(length(na.action.data) + nrow(modelMatrix)))[-na.action.data]
        na.action.ret[] <- modelMatrixIndices[na.action.ret]
        na.action <- c(na.action.ret, na.action.data)
        if (!is.null(na.action)) {
            na.action <- sort(na.action)
            class(na.action) <- ifelse(is.null(na.action.ret), class(na.action.data), class(na.action.ret))
            attr(ret, "na.action") <- na.action
        }

     } else {
        # store summed terms only

        # check size of model matrix
        if (nrow(modelMatrix) > 0) {

            ## calculate abundance
            l <- fit$bugsFit$sims.list

            # constant term
            ret <- matrix(l$gamma.const, fit$bugsFit$n.sims, nrow(modelMatrix), byrow=FALSE, dimnames=list(NULL, rownames(modelMatrix)))

            # linear terms
            for (var in fit$rhoLinearVars) {
                varLong <- paste("gamma", var, sep=".")
                ret <- ret + tcrossprod(l[[varLong]], modelMatrix[, var])
            }

            # random walk terms
            for (var in c(fit$rhoRW1Vars, fit$rhoRW2Vars)) {
                varLong <- paste("gamma", var, sep=".")
                for (i in 1:nrow(ret))
                    ret[i, ] <- ret[i, ] + approx(fit$rwBoundaries[[varLong]], l[[varLong]][i, ], modelMatrix[, var])$y
            }

            # spatial term
            for (var in fit$rhoSpatialVar) {
                varLong <- paste("gamma", var, sep=".")
                ret <- ret + l[[varLong]][, modelMatrix[, var]]
            }

            # transform
            ret <- expit(ret)

        } else
            ret <- array(dim=c(fit$bugsFit$n.sims, 0))

        # remove missing values and add to na.action
        ret[!is.finite(ret)] <- NA  # make infinite values NA
        ret <- t(na.action(t(ret)))
        na.action.ret <- attr(ret, "na.action")
        if (is.null(na.action.data))
            modelMatrixIndices <- 1:nrow(modelMatrix)
        else
            modelMatrixIndices <- (1:(length(na.action.data) + nrow(modelMatrix)))[-na.action.data]
        na.action.ret[] <- modelMatrixIndices[na.action.ret]
        na.action <- c(na.action.ret, na.action.data)
        if (!is.null(na.action)) {
            na.action <- sort(na.action)
            class(na.action) <- ifelse(is.null(na.action.ret), class(na.action.data), class(na.action.ret))
            attr(ret, "na.action") <- na.action
        }

    }

    ret
}

