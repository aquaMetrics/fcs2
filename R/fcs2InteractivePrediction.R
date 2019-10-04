#' Interactive FCS2 model prediction
#'
#' Plots the posterior distribution of the abundance and prevalence components
#' as well as the predictive distribution of a new fish catch at a survey
#' characterised by model covariates that can be modified interactively.
#' Controls are provided for each covariate in the model so that the change in
#' predictive probabilities can be visualised as covariates are adjusted.  The
#' single \acronym{EQR} variable can additionally be displayed and class
#' boundaries can be provided to colour the probabilities of each class.
#'
#'
#' @param fit an \code{"fcs2Fit"} object containing a full \acronym{FCS2} model
#' fit, as returned from \code{\link{fcs2FitModel}} with \code{runBUGS = TRUE}.
#' @param data a data frame with surveys as rows and variables as columns. It
#' should contain all variables required by \code{fit} and is used primarily to
#' specify the upper and lower limits for the interactive controls for each
#' covariate.
#' @param init.row optionally, an integer corresponding to a row in \code{data}
#' that should be used as the initial state for each covariate.  If omitted or
#' the selected value is missing, the initial state for a covariate is selected
#' as the average value in \code{data}.
#' @param eqr whether to give the distribution of the \acronym{EQR} variable.
#' Defaults to \code{TRUE}.
#' @param boundaries a vector of length 4 giving the \acronym{EQR} boundaries
#' separating the classes \emph{Bad}, \emph{Poor}, \emph{Good}, \emph{Moderate}
#' and \emph{High}.  These are used only to colour the plots with \emph{Bad}
#' red and \emph{High} blue.  If \code{NULL} (default), the probability that
#' defines the single \acronym{EQR} is coloured blue.
#' @section Warning: This function requires the additional package \pkg{rpanel}
#' for producing interactive controls.
#' @seealso \code{\link{fcs2InteractiveLikelihood}} which produces a simpler
#' interactive demonstration of how the \acronym{FCS2} \acronym{ZINB}
#' likelihood varies with model parameters.
#' @keywords hplot
#' @export
fcs2InteractivePrediction <- function(fit, data, init.row, eqr = TRUE, boundaries = NULL)
{
    # load package 'rpanel'
    if (!requireNamespace("rpanel", quietly = TRUE))
        stop("`fcs2InteractivePrediction' requires R package `rpanel' - please install using `install.packages'")

    # check fit contains bugs fit
    if (is.null(fit$bugsFit))
        stop("need posterior samples from BUGS model fit - use 'fcs2FitModel(fit=fit, runBUGS=TRUE, ...)'")

    # check 'newData' provided
    if (missing(data))
        stop("a data frame 'data' must be provided to set ranges of variables for controls")

    # convert data to data frame
    data <- as.data.frame(data)


    ## extract names of covariates via get_all_vars and simple formula

    # begin with survey area
    covariates <- data[, fit$surveyAreaVar, drop=FALSE]

    # add number of runs if multiple pass survey
    if (fit$multiRun) {
        if (fit$nRunsVar %in% colnames(data)) {
            # use nRunsVar
            nRuns <- data[, fit$nRunsVar]

        } else if (sum(fit$runTotalVars %in% colnames(data)) == length(fit$runTotalVars)) {
            # calculate from runTotalVars
            nRuns <- rep(NA, nrow(data))
            for (i in 1:nrow(data))
                nRuns[i] <- sum(!is.na(data[i, fit$runTotalVars]))

        } else {
            # set nRuns to range from 1 to 7
            nRuns <- round(seq(1, 7, length=nrow(data)))
        }

        # add column and name
        covariates <- cbind(covariates, nRuns)
        colnames(covariates)[ncol(covariates)] <- fit$nRunsVar

        # make sure class of column is integer
        class(covariates[[ncol(covariates)]]) <- "integer"
    }

    # if showing EQR, add total catch
    if (eqr) {
        if (!fit$multiRun) {
            # use allRunsTotalVar if available
            if (fit$allRunsTotalVar %in% colnames(data))
                catch <- data[, fit$allRunsTotalVar]
            else
                catch <- rep(NA, nrow(data))

        } else {
            catch <- rep(NA, nrow(data))

            # from runs
            if (!is.null(fit$runTotalVars) && sum(fit$runTotalVars %in% colnames(data)) == length(fit$runTotalVars)) {
                for (i in 1:nrow(data))
                    if (!is.na(nRuns[i]))
                        catch[i] <- sum(as.numeric(data[i, fit$runTotalVars[1:nRuns[i]]]))
            }

            # from all runs total
            if (!is.null(fit$allRunsTotalVar) && fit$allRunsTotalVar %in% colnames(data) && sum(is.na(catch)) > 0)
                catch[is.na(catch)] <- data[is.na(catch), fit$allRunsTotalVar]
        }

        # from all runs range
        if (!is.null(fit$allRunsRangeVars) && sum(fit$allRunsRangeVars %in% colnames(data)) == 2 && sum(is.na(catch)) > 0)
            catch[is.na(catch)] <- data[is.na(catch), fit$allRunsRangeVars[1]]

        # if all missing, use range from 0 to 20
        if (sum(is.na(catch)) == nrow(data))
            catch <- round(seq(0, 20, length=nrow(data)))

        # add column and name
        covariates <- cbind(covariates, catch)
        if (is.null(fit$allRunsTotalVar))
            fit$allRunsTotalVar <- "AllRunsTotal"
        colnames(covariates)[ncol(covariates)] <- fit$allRunsTotalVar

        # make sure class of column is integer
        class(covariates[[ncol(covariates)]]) <- "integer"
    }

    # correct factors in linear mu terms
    muLinearVars <- fit$muLinearVars
    if (length(muLinearVars) > 0) {
        for (i in 1:length(muLinearVars))
            if (length(grep("factor", muLinearVars[i])) > 0)
                muLinearVars[i] <- substr(muLinearVars[i], 0, regexpr(")", muLinearVars[i]))
    } else
        muLinearVars <- "1"

    # remove TRUE or FALSE
    muLinearVars <- sub("TRUE", "", muLinearVars)
    muLinearVars <- sub("FALSE", "", muLinearVars)

    # calculate new model matrix
    simpleFormula <- formula(paste("~", paste(c(muLinearVars, fit$muRW1Vars, fit$muRW2Vars, fit$muSpatialVar), collapse="+")))
    covariates <- cbind(covariates, get_all_vars(simpleFormula, data))

    # correct factors in linear rho terms
    rhoLinearVars <- fit$rhoLinearVars
    if (length(rhoLinearVars) > 0) {
        for (i in 1:length(rhoLinearVars))
            if (length(grep("factor", rhoLinearVars[i])) > 0)
                rhoLinearVars[i] <- substr(rhoLinearVars[i], 0, regexpr(")", rhoLinearVars[i]))
    } else
        rhoLinearVars <- "1"

    # remove TRUE or FALSE
    rhoLinearVars <- sub("TRUE", "", rhoLinearVars)
    rhoLinearVars <- sub("FALSE", "", rhoLinearVars)

    # calculate new model matrix
    simpleFormula <- formula(paste("~", paste(c(rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars, fit$rhoSpatialVar), collapse="+")))
    covariates <- cbind(covariates, get_all_vars(simpleFormula, data))

    # remove duplicates
    covariates <- covariates[, unique(colnames(covariates))]
    covnames <- colnames(covariates)


    # calculate range and best values of each parameter
    limit <- best <- list()
    for (cov in covnames)
        if (class(covariates[[cov]]) == "numeric") {
            limit[[cov]] <- range(covariates[[cov]], na.rm=TRUE)
            best[[cov]] <- mean(covariates[[cov]])

        } else if (class(covariates[[cov]]) == "integer") {
            limit[[cov]] <- range(covariates[[cov]], na.rm=TRUE)
            best[[cov]] <- round(mean(covariates[[cov]]))

        } else if (class(covariates[[cov]]) == "factor") {
            best[[cov]] <- names(which.max(table(covariates[[cov]])))
        } else if (class(covariates[[cov]]) == "logical") {
            best[[cov]] <- as.logical(names(which.max(table(covariates[[cov]]))))
        }

    # select initial values
    if (missing(init.row))
        init <- best
    else {
        init <- covariates[init.row,]
        for (cov in covnames)
            if (is.na(init[[cov]]))
                init[[cov]] <- best[[cov]]
    }

    # create panel
    text <- "rp.control("
    for (cov in covnames)
        if (class(covariates[[cov]]) == "factor")
            text <- c(text, paste(cov, ' = "', init[[cov]], '", ', sep=''))
        else
            text <- c(text, paste(cov, "=", init[[cov]], ","))
    text <- c(text, "title='FCS2 prediction')")
    panel <- eval(parse(text=text))

    # create function to redraw plots (pdf of abundance and prevalence for now)
    draw <- function(panel) {
        # create data frame
        data <- as.data.frame(panel[covnames])

        # set up plot
        par.old <- graphics::par(mfrow=c(2,2), ask=FALSE)

        # mu
        mu <- abundance(fit, data)
        if (length(mu) > 0) {
            if (sd(mu) == 0)
                dmu <- density(jitter(mu), from=0)  # prevents bw being too large when all identical
            else
                dmu <- density(mu, from=0)
            graphics::plot(dmu, ylim=c(0, max(dmu$y)), main=expression(paste("Abundance ", mu)), xlab="Abundance")

        } else
            graphics::frame()

        # rho
        rho <- prevalence(fit, data)
        if (length(rho) > 0) {
            if (sd(rho) == 0)
                drho <- density(jitter(rho), from=0, to=1)  # prevents bw being too large when all identical
            else
                drho <- density(rho, from=0, to=1)
            graphics::plot(drho, ylim=c(0, max(drho$y)), main=expression(paste("Prevalence ", rho)), xlab="Prevalence")

            if (length(mu) > 0) {
                # catch
                meanCatch <- predict(fit, data, mu=mu, rho=rho)
                if (fit$multiRun)
                    title <- paste("Total catch: mean =", signif(meanCatch, 3))
                else
                    title <- paste("Catch: mean =", signif(meanCatch, 3))
                plotCatchPMF(fit, data, boundaries=boundaries, mu=mu, rho=rho, title=title, maxpts=100)

                # calculate eqr
                if (eqr) {
                    eqrSamples <- fcs2SingleEQR(fit, data, mu=mu, rho=rho)
                    graphics::plot(eqrSamples, boundaries=boundaries, title=paste("EQR: mean =", signif(mean(eqrSamples), 3)))
                }
            }
        }

        graphics::par(par.old)

        panel
    }

    # add sliders to control parameters
    nrows <- ceiling(sqrt(length(covnames)))
    i <- j <- 1
    for (cov in covnames) {
        if (class(covariates[[cov]]) == "numeric") {
            text <- paste("rp.slider(panel,", cov, ",", limit[[cov]][1], ",", limit[[cov]][2],
                          ", draw, showvalue=TRUE, resolution=", signif(diff(limit[[cov]])/500, 3),
                          ", row=", i, ", column=", j, ")")
            i <- i + 1

        } else if (class(covariates[[cov]]) == "integer") {
            text <- paste("rp.slider(panel,", cov, ",", limit[[cov]][1], ",", limit[[cov]][2],
                          ", draw, showvalue=TRUE, resolution=1, row=",
                          i, ", column=", j, ")")
            i <- i + 1

        } else if (class(covariates[[cov]]) == "factor") {
            text <- paste("rp.radiogroup(panel,", cov, ",levels(covariates[[cov]]), action=draw, row=",
                          i, ", column=", j, ", rowspan=", ceiling((nlevels(covariates[[cov]])+1)/3), ")")
            i <- i + ceiling((nlevels(covariates[[cov]])+1)/3)

        } else if (class(covariates[[cov]]) == "logical") {
            text <- paste("rp.checkbox(panel,", cov, ", action=draw, row=",
                          i, ", column=", j, ")")
            i <- i + 1
        }
        if (i > nrows) {
            i <- 1
            j <- j + 1
        }

        eval(parse(text=text))
    }
}

