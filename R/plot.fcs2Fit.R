#' Plot FCS2 Model Fit
#'
#' @encoding UTF-8
#' @title Plot FCS2 Model Fit
#'
#' @description
#' Produces density plots of the prior and posterior estimates of each variable
#' in the \acronym{FCS2} model.  Can also produce series plots that display the
#' relationship between each covariate and the corresponding abundance or
#' prevalence component.
#'
#'
#' @param x an \code{"fcs2Fit"} object, as returned from
#' \code{\link{fcs2FitModel}}.
#' @param variables an optional character vector giving the names of the model
#' variables to plot.  \code{fcs2:::variable.names.fcs2Fit} can be used for
#' this.
#' @param prior whether to plot the prior distribution. These are only
#' available for plots of the density of a single variable.
#' @param posterior whether to plot the posterior distribution, as estimated
#' from \acronym{BUGS} samples (if available).
#' @param inla whether to plot the approximate posterior estimates produced by
#' \acronym{INLA} (if available).
#' @param legend whether to add a legend to each plot.
#' @param samples whether to display the individual posterior samples rather
#' than summarising each posterior term.
#' @param groupLinearVars whether to group multiple linear terms of the same
#' covariate together to display the overall covariate relationship rather than
#' each density estimate.  This is not available for the prior and
#' \acronym{INLA} only produces a mean estimate of the relationship based on
#' the false assumption of independence between the linear parameters.
#' Therefore the \acronym{INLA} relationship line should only be used as a
#' rough guide when multiple terms exist.
#' @param \dots Not used in this function
#' @seealso \code{\link{fcs2FitModel}}
#' @keywords hplot

#' @export
plot.fcs2Fit <-
function(x, variables = variable.names(x), prior = TRUE, posterior = !is.null(x$bugsFit),
         inla = !posterior, legend = TRUE, samples = FALSE, groupLinearVars = posterior, ...)
{
    # check inputs
    if (posterior && is.null(x$bugsFit))
        posterior <- FALSE
    if (inla && is.null(x$inlaFit))
        inla <- FALSE
    muINLA <- inla && !is.null(x$inlaFits$muFit)
    rhoINLA <- inla && !is.null(x$inlaFits$rhoFit)

    # create variable vectors
    rwGroups <- variable.names(x, q=FALSE, r=FALSE, linear=FALSE, hyperparam=FALSE, rw="group", spatial=FALSE)
    spatialGroups <- variable.names(x, q=FALSE, r=FALSE, linear=FALSE, hyperparam=FALSE, rw=FALSE, spatial="group")

    # if grouping linear terms (and can), remove variables with '^' in name
    #   and create vector of all linear variables
    ## NOTE: this assumes that linear term collections all have first order term written without '^'
    if (groupLinearVars && ((!posterior || is.null(x$bugsFit)) && (!inla || is.null(x$inlaFits)))) {
        groupLinearVars <- FALSE
        warning("series plot of grouped linear variables only available for INLA estimates or posterior samples")
    }
    if (groupLinearVars) {
        if (posterior && !is.null(x$bugsFit)) {
            toRemove <- grep("^", variables, fixed=TRUE)
            if (length(toRemove) > 0)
                variables <- variables[-toRemove]
            linearVars <- variable.names(x, q=FALSE, r=FALSE, linear=TRUE, hyperparam=FALSE, rw=FALSE, spatial=FALSE)

        } else {
            linearVars <- NULL
            if (muINLA) {
                muVars <- grep("beta", variables, fixed=TRUE, value=TRUE)
                toRemove <- grep("^", muVars, fixed=TRUE, value=TRUE)
                toRemove <- match(toRemove, variables)
                if (length(toRemove) > 0)
                    variables <- variables[-toRemove]
                linearVars <- variable.names(x, q=FALSE, r=FALSE, linear=TRUE, hyperparam=FALSE, rw=FALSE, spatial=FALSE, prevalence=FALSE)
            }
            if (rhoINLA) {
                rhoVars <- grep("gamma", variables, fixed=TRUE, value=TRUE)
                toRemove <- grep("^", rhoVars, fixed=TRUE, value=TRUE)
                toRemove <- match(toRemove, variables)
                if (length(toRemove) > 0)
                    variables <- variables[-toRemove]
                linearVars <- c(linearVars, variable.names(x, q=FALSE, r=FALSE, linear=TRUE, hyperparam=FALSE, rw=FALSE, spatial=FALSE, abundance=FALSE))
            }
        }
    }

    # set up multiple plots
    if (length(variables) > 1) {
        w <- ceiling(sqrt(length(variables)))
        par.old <- graphics::par(mfrow=c(w - (length(variables) <= w*(w-1)), w))
        on.exit(graphics::par(par.old))
    }

    # set up plot (colours, lines, etc)
    muEstimates <- c(prior, muINLA, posterior)
    rhoEstimates <- c(prior, rhoINLA, posterior)
    estimateNames <- c("Prior", "INLA", "Posterior")
    lty <- c(2, 4, 1)
    col <- c("black", "red", "blue")  #c("grey30", "grey15", "black")
    prob <- 0.95  ## NOTE: this value is hard-coded in INLA

    for (var in variables) {
        if (!(var %in% c(rwGroups, spatialGroups))) {
            ## PDF plot:

            if (var == "q") {
                ## Catch probability q:

                # set plot limits and find values
                xlim <- c(0, 1)
                ylim <- c(Inf, -Inf)
                if (posterior) {
                    bugsPosterior <- density(x$bugsFit$sims.matrix[, var], from=xlim[1], to=xlim[2])
                    ylim <- c(0, max(bugsPosterior$y))
                }
                if (prior) {
                    priorParameters <- x$prior.parameters[[var]]
                    xplt <- seq(xlim[1], xlim[2], l=200)
                    ypltPrior <- dbeta(xplt, priorParameters["a"], priorParameters["b"])
                    ylim <- c(0, max(ylim[2], ypltPrior[ypltPrior < Inf]))
                }

                # plot
                graphics::plot(0, 0, col="white", xlim=xlim, ylim=ylim, xlab="", ylab="Density", main=var)
                graphics::abline(h=0, v=c(0, 1), col="grey80")
                if (prior)
                    graphics::lines(xplt, ypltPrior, lty=lty[1], col=col[1])
                if (posterior) {
                    graphics::lines(bugsPosterior, lty=lty[3], col=col[3])

                    if (samples) {
                        samplesCol <- grDevices::rgb2hsv(grDevices::col2rgb(col[3]))
                        samplesCol <- grDevices::hsv(h=samplesCol[1], s=samplesCol[2], v=samplesCol[3], a=min(1, 100 / nrow(x$bugsFit$sims.matrix[, var])))  #  OR: 1 / log10(nrow(x$bugsFit$sims.matrix[, var]))))
                        graphics::rug(x$bugsFit$sims.matrix[, var], col=samplesCol)
                    }
                }
                est <- muEstimates
                est[2] <- FALSE  # no INLA
                if (legend)
                    legend("topright", legend=estimateNames[est], lty=lty[est], col=col[est], bg=grDevices::hsv(s=0, a=0.5))


            } else if (var == "r") {
                ## Shape parameter r:

                # set plot limits and find values
                xlim <- ylim <- c(Inf, -Inf)
                if (prior) {
                    priorParameters <- x$prior.parameters[[var]]
                    xlim <- qlnorm(c(0.025, 0.975), priorParameters["mu"], 1 / sqrt(priorParameters["tau"]))
                    if (xlim[1] < 1)
                        xlim[1] <- 0
                }
                if (muINLA) {
                    inlaMarginal <- x$inlaFits$muFit$marginals.hyperpar[[1]]
                    xlim <- c(min(xlim[1], inla.qmarginal(0.01, inlaMarginal), na.rm=TRUE),
                              max(xlim[2], inla.qmarginal(0.99, inlaMarginal), na.rm=TRUE))
                    if (xlim[1] < 1)
                        xlim[1] <- 0
                }
                if (posterior) {
                    xlim <- c(min(xlim[1], x$bugsFit$sims.matrix[, var]), max(xlim[2], x$bugsFit$sims.matrix[, var]))
                    if (xlim[1] < 1)
                        xlim[1] <- 0
                    bugsPosterior <- density(x$bugsFit$sims.matrix[, var], from=xlim[1], to=xlim[2])
                    ylim <- c(0, max(bugsPosterior$y))
                }
                if (prior) {
                    xplt <- seq(xlim[1], xlim[2], l=200)
                    ypltPrior <- dlnorm(xplt, priorParameters["mu"], 1 / sqrt(priorParameters["tau"]))
                    ylim <- c(0, max(ylim[2], ypltPrior[ypltPrior < Inf]))
                }
                if (muINLA) {
                    xplt <- seq(xlim[1], xlim[2], l=200)
                    ypltINLA <- inla.dmarginal(xplt, inlaMarginal)
                    ylim <- c(0, max(ylim[2], ypltINLA))
                }

                # plot
                graphics::plot(0, 0, col="white", xlim=xlim, ylim=ylim, xlab="", ylab="Density", main=var)
                graphics::abline(h=0, v=0, col="grey80")
                if (prior)
                    graphics::lines(xplt, ypltPrior, lty=lty[1], col=col[1])
                if (muINLA)
                    graphics::lines(xplt, ypltINLA, lty=lty[2], col=col[2])
                if (posterior) {
                    graphics::lines(bugsPosterior, lty=lty[3], col=col[3])

                    if (samples) {
                        samplesCol <- grDevices::rgb2hsv(grDevices::col2rgb(col[3]))
                        samplesCol <- grDevices::hsv(h=samplesCol[1], s=samplesCol[2], v=samplesCol[3], a=min(1, 100 / nrow(x$bugsFit$sims.matrix[, var])))  #  OR: 1 / log10(nrow(x$bugsFit$sims.matrix[, var]))))
                        graphics::rug(x$bugsFit$sims.matrix[, var], col=samplesCol)
                    }
                }
                if (legend)
                    legend("topright", legend=estimateNames[muEstimates], lty=lty[muEstimates], col=col[muEstimates], bg=grDevices::hsv(s=0, a=0.5))


            } else if (max(pmatch(c("beta", "gamma"), var, nomatch=0)) > 0) {
                ## Linear terms or individual CAR parameters:

                # check whether variable should be grouped
                if (groupLinearVars && var %in% linearVars &&
                    length(grep(".const", var, fixed=TRUE)) == 0 && length(grep("I(", var, fixed=TRUE)) == 0) {
                    ## Series plot for group of linear variables

                    # extract covariate name (by removing "beta." or "gamma.")
                    if (pmatch("beta", var, nomatch=0))
                        covariate <- sub("beta.", "", var)
                    else
                        covariate <- sub("gamma.", "", var)

                    # search for all linear variables with this covariate (making sure also has beta/gamma as appropriate)
                    vars <- grep(covariate, linearVars, fixed=TRUE, value=TRUE)
                    if (pmatch("beta", var, nomatch=0))
                        vars <- grep("beta.", vars, fixed=TRUE, value=TRUE)
                    else
                        vars <- grep("gamma.", vars, fixed=TRUE, value=TRUE)

                    # for each var, extract order
                    orders <- rep(0, length(vars))
                    for (i in 1:length(vars)) {
                        if (length(grep("^", vars[i], fixed=TRUE)) == 0)
                            orders[i] <- 1
                        else {
                            varSplit <- strsplit(vars[i], "^", fixed=TRUE)
                            orders[i] <- as.numeric(substr(varSplit[[1]][2], 1, 1))
                        }
                    }

                    # extract x-axis range
                    xlim <- range(x$modelMatrix[, covariate])
                    xpoints <- seq(xlim[1], xlim[2], l=30)

                    # calculate term for each posterior sample
                    if (posterior) {
                        bugsSamples <- array(0, dim=c(nrow(x$bugsFit$sims.matrix), length(xpoints)))
                        for (i in 1:length(vars))
                            for (j in 1:length(xpoints))
                                bugsSamples[, j] <- bugsSamples[, j] + (xpoints[j] ^ orders[i]) * x$bugsFit$sims.matrix[, vars[i]]
                    }

                    # calculate term using average marginal values
                    groupINLA <- FALSE
                    if (inla) {
                        inlaEst <- rep(0, length(xpoints))
                        if (pmatch("beta", var, nomatch=0) && muINLA) {
                            for (i in 1:length(vars)) {
                                inlaMarginal <- x$inlaFits$muFit$marginals.fixed[[sub("beta.", "", make.names(vars[i]))]]
                                inlaEst <- inlaEst + (xpoints ^ orders[i]) * inla.emarginal(identity, inlaMarginal)
                            }
                            groupINLA <- TRUE

                        } else if (rhoINLA) {
                            for (i in 1:length(vars)) {
                                inlaMarginal <- x$inlaFits$rhoFit$marginals.fixed[[sub("gamma.", "", make.names(vars[i]))]]
                                inlaEst <- inlaEst + (xpoints ^ orders[i]) * inla.emarginal(identity, inlaMarginal)
                            }
                            groupINLA <- TRUE
                        }
                    }

                    # calculate y-axis limits
                    ylim <- c(Inf, -Inf)
                    if (posterior) {
                        if (samples) {
                            ylim <- range(bugsSamples)

                        } else {
                            bugsMeans <- colMeans(bugsSamples)
                            bugsCI <- cbind(bottom=apply(bugsSamples, 2, quantile, prob=((1 - prob) / 2)),
                                            top=apply(bugsSamples, 2, quantile, prob=1 - ((1 - prob) / 2)))
                            ylim <- range(bugsCI)
                        }
                    }
                    if (groupINLA)
                        ylim <- c(min(ylim[1], inlaEst), max(ylim[2], inlaEst))

                    # plot either samples or summaries
                    graphics::plot(0, 0, col="white", xlim=xlim, ylim=ylim, xlab=covariate, ylab="", main=var)
                    graphics::abline(h=0, col="grey80")
                    if (groupINLA)
                        graphics::lines(xpoints, inlaEst, col=col[2], lty=1)
                    if (posterior) {
                        if (samples) {
                            samplesCol <- grDevices::rgb2hsv(grDevices::col2rgb(col[3]))
                            samplesCol <- grDevices::hsv(h=samplesCol[1], s=samplesCol[2], v=samplesCol[3], a=min(1, 100 / nrow(bugsSamples)))  #  OR: 1 / log10(nrow(bugsSamples))))
                            graphics::matplot(add=TRUE, xpoints, t(bugsSamples), col=samplesCol, t='l', pch=20, lty=1)

                        } else {
                            graphics::lines(xpoints, bugsCI[, 1], col=col[3], lty=2)
                            graphics::lines(xpoints, bugsMeans, col=col[3], lty=1)
                            graphics::lines(xpoints, bugsCI[, 2], col=col[3], lty=2)
                        }
                    }
                    graphics::rug(x$modelMatrix[, covariate])
                    est <- c(FALSE, groupINLA, posterior)
                    if (legend)
                        legend("topright", legend=estimateNames[est], lty=lty[est], col=col[est], bg=grDevices::hsv(s=0, a=0.5))

                } else {
                    ## PDF of single linear variable

                    # set inla
                    if (pmatch("beta", var, nomatch=0))
                        inla <- muINLA
                    else
                        inla <- rhoINLA

                    # set plot limits
                    xlim <- ylim <- c(Inf, -Inf)
                    if (prior && substr(var, nchar(var), nchar(var)) != "]") {
                        priorParameters <- x$prior.parameters[[var]]
                        xlim <- qnorm(c(0.01, 0.99), priorParameters["mean"], 1 / sqrt(priorParameters["precision"]))
                    }
                    if (inla) {
                        # find marginal estimate
                        if (pmatch("beta", var, nomatch=0)) {
                            if (var == "beta.const")
                                inlaMarginal <- x$inlaFits$muFit$marginals.fixed[["(Intercept)"]]
                            else
                                inlaMarginal <- x$inlaFits$muFit$marginals.fixed[[sub("beta.", "", make.names(var))]]
                        } else {
                            if (var == "gamma.const")
                                inlaMarginal <- x$inlaFits$rhoFit$marginals.fixed[["(Intercept)"]]
                            else
                                inlaMarginal <- x$inlaFits$rhoFit$marginals.fixed[[sub("gamma.", "", make.names(var))]]
                        }
                        xlim <- c(min(xlim[1], inla.qmarginal(0.01, inlaMarginal), na.rm=TRUE),
                                  max(xlim[2], inla.qmarginal(0.99, inlaMarginal), na.rm=TRUE))
                    }
                    if (posterior) {
                        xlim <- c(min(xlim[1], x$bugsFit$sims.matrix[, var]), max(xlim[2], x$bugsFit$sims.matrix[, var]))
                        bugsPosterior <- density(x$bugsFit$sims.matrix[, var], from=xlim[1], to=xlim[2])
                        ylim <- c(0, max(bugsPosterior$y))
                    }
                    if (prior && substr(var, nchar(var), nchar(var)) != "]") {
                        xplt <- seq(xlim[1], xlim[2], l=200)
                        ypltPrior <- dnorm(xplt, priorParameters["mean"], 1 / sqrt(priorParameters["precision"]))
                        ylim <- c(0, max(ylim[2], ypltPrior))
                    }
                    if (inla) {
                        xplt <- seq(xlim[1], xlim[2], l=200)
                        ypltINLA <- inla.dmarginal(xplt, inlaMarginal)
                        ylim <- c(0, max(ylim[2], ypltINLA))
                    }

                    # plot
                    graphics::plot(0, 0, col="white", xlim=xlim, ylim=ylim, xlab="", ylab="Density", main=var)
                    graphics::abline(h=0, v=0, col="grey80")
                    if (prior && substr(var, nchar(var), nchar(var)) != "]")
                        graphics::lines(xplt, ypltPrior, lty=lty[1], col=col[1])
                    if (inla)
                        graphics::lines(xplt, ypltINLA, lty=lty[2], col=col[2])
                    if (posterior) {
                        graphics::lines(bugsPosterior, lty=lty[3], col=col[3])

                        if (samples) {
                            samplesCol <- grDevices::rgb2hsv(grDevices::col2rgb(col[3]))
                            samplesCol <- grDevices::hsv(h=samplesCol[1], s=samplesCol[2], v=samplesCol[3], a=min(1, 100 / nrow(x$bugsFit$sims.matrix[, var])))  #  OR: 1 / log10(nrow(x$bugsFit$sims.matrix[, var]))))
                            graphics::rug(x$bugsFit$sims.matrix[, var], col=samplesCol)
                        }
                    }
                    if (legend) {
                        if (pmatch("beta", var, nomatch=0))
                            est <- muEstimates
                        else
                            est <- rhoEstimates
                        est[1] <- (prior && substr(var, nchar(var), nchar(var)) != "]")
                        legend("topright", legend=estimateNames[est], lty=lty[est], col=col[est], bg=grDevices::hsv(s=0, a=0.5))
                    }
                }

            } else if (max(pmatch(c("sigma", "nu"), var, nomatch=0)) > 0) {
                ## Hyperparameters for CAR terms (on sd scale):

                # find name of precision variable and set inla
                if (pmatch("sigma", var, nomatch=0)) {
                    precVarName <- sub("sigma", "tau", var)
                    inla <- muINLA
                    estimates <- muEstimates

                } else {
                    precVarName <- sub("nu", "phi", var)
                    inla <- rhoINLA
                    estimates <- rhoEstimates
                }

                # set plot limits
                xlim <- ylim <- c(Inf, -Inf)
                if (prior) {
                    priorParameters <- x$prior.parameters[[precVarName]]
                    xlim <- qhyperprior(c(0.025, 0.975), priorParameters["a"], priorParameters["b"])
                    if (xlim[1] < 1)
                        xlim[1] <- 0
                }
                if (inla) {
                    # find marginal estimate for precision variable
                    if (pmatch("sigma", var, nomatch=0))
                        inlaMarginal <- x$inlaFits$muFit$marginals.hyperpar[[sub("sigma.", "Precision for ", make.names(var))]]
                    else
                        inlaMarginal <- x$inlaFits$rhoFit$marginals.hyperpar[[sub("nu.", "Precision for ", make.names(var))]]
                    xlim <- c(min(xlim[1], 1 / sqrt(inla.qmarginal(0.99, inlaMarginal)), na.rm=TRUE),
                              max(xlim[2], 1 / sqrt(inla.qmarginal(0.01, inlaMarginal)), na.rm=TRUE))
                }
                if (posterior) {
                    xlim <- c(min(xlim[1], 1 / sqrt(x$bugsFit$sims.matrix[, precVarName]), na.rm=TRUE),
                              max(xlim[2], 1 / sqrt(x$bugsFit$sims.matrix[, precVarName]), na.rm=TRUE))
                    if (xlim[1] < 1)
                        xlim[1] <- 0
                    bugsPosterior <- density(1 / sqrt(x$bugsFit$sims.matrix[, precVarName]), from=xlim[1], to=xlim[2])
                    ylim <- c(0, max(bugsPosterior$y))
                }
                if (prior) {
                    xplt <- seq(xlim[1], xlim[2], l=200)
                    ypltPrior <- dhyperprior(xplt, priorParameters["a"], priorParameters["b"])
                    ylim <- c(0, max(ylim[2], ypltPrior[ypltPrior < Inf]))
                }
                if (inla) {
                    xplt <- seq(xlim[1], xlim[2], l=200)
                    ypltINLA <- inla.dmarginal(1 / (xplt^2), inlaMarginal) * 2 / (xplt^3)
                    ypltINLA[is.nan(ypltINLA)] <- Inf
                    ylim <- c(0, max(ylim[2], ypltINLA[ypltINLA < Inf]))
                }

                # plot
                graphics::plot(0, 0, col="white", xlim=xlim, ylim=ylim, xlab="", ylab="Density", main=var)
                graphics::abline(h=0, v=0, col="grey80")
                if (prior)
                    graphics::lines(xplt, ypltPrior, lty=lty[1], col=col[1])
                if (inla)
                    graphics::lines(xplt, ypltINLA, lty=lty[2], col=col[2])
                if (posterior) {
                    graphics::lines(bugsPosterior, lty=lty[3], col=col[3])

                    if (samples) {
                        samplesCol <- grDevices::rgb2hsv(grDevices::col2rgb(col[3]))
                        samplesCol <- grDevices::hsv(h=samplesCol[1], s=samplesCol[2], v=samplesCol[3], a=min(1, 100 / nrow(x$bugsFit$sims.matrix[, precVarName])))  #  OR: 1 / log10(nrow(x$bugsFit$sims.matrix[, precVarName]))))
                        graphics::rug(1 / sqrt(x$bugsFit$sims.matrix[, precVarName]), col=samplesCol)
                    }
                }
                if (legend)
                    legend("topright", legend=estimateNames[estimates], lty=lty[estimates], col=col[estimates], bg=grDevices::hsv(s=0, a=0.5))

            } else if (max(pmatch(c("tau", "phi"), var, nomatch=0)) > 0) {
                ## Hyperparameters for CAR terms (on precision scale):

                # set inla
                if (pmatch("sigma", var, nomatch=0)) {
                    inla <- rhoINLA
                    estimates <- rhoEstimates

                } else {
                    inla <- muINLA
                    estimates <- muEstimates
                }

                # set plot limits
                xlim <- ylim <- c(Inf, -Inf)
                 if (prior) {
                    priorParameters <- x$prior.parameters[[var]]
                    xlim <- qgamma(c(0.025, 0.975), priorParameters["a"], priorParameters["b"])
                    if (xlim[1] < 1)
                        xlim[1] <- 0
                }
                if (inla) {
                    # find marginal estimate
                    if (pmatch("tau", var, nomatch=0))
                        inlaMarginal <- x$inlaFits$muFit$marginals.hyperpar[[sub("tau.", "Precision for ", make.names(var))]]
                    else
                        inlaMarginal <- x$inlaFits$rhoFit$marginals.hyperpar[[sub("phi.", "Precision for ", make.names(var))]]
                    xlim <- c(min(xlim[1], inla.qmarginal(0.01, inlaMarginal), na.rm=TRUE),
                              max(xlim[2], inla.qmarginal(0.99, inlaMarginal), na.rm=TRUE))
                }
                if (posterior) {
                    xlim <- c(min(xlim[1], x$bugsFit$sims.matrix[, var]), max(xlim[2], x$bugsFit$sims.matrix[, var]))
                    if (xlim[1] < 1)
                        xlim[1] <- 0
                    bugsPosterior <- density(x$bugsFit$sims.matrix[, var], from=xlim[1], to=xlim[2])
                    ylim <- c(0, max(bugsPosterior$y))
                }
                if (prior) {
                    xplt <- seq(xlim[1], xlim[2], l=200)
                    ypltPrior <- dgamma(xplt, priorParameters["a"], priorParameters["b"])
                    ylim <- c(0, max(ylim[2], ypltPrior[ypltPrior < Inf]))
                }
                if (inla) {
                    xplt <- seq(xlim[1], xlim[2], l=200)
                    ypltINLA <- inla.dmarginal(xplt, inlaMarginal)
                    ylim <- c(0, max(ylim[2], ypltINLA[ypltINLA < Inf]))
                }

                # plot
                graphics::plot(0, 0, col="white", xlim=xlim, ylim=ylim, xlab="", ylab="Density", main=var)
                graphics::abline(h=0, v=0, col="grey80")
                if (prior)
                    graphics::lines(xplt, ypltPrior, lty=lty[1], col=col[1])
                if (inla)
                    graphics::lines(xplt, ypltINLA, lty=lty[2], col=col[2])
                if (posterior) {
                    graphics::lines(bugsPosterior, lty=lty[3], col=col[3])

                    if (samples) {
                        samplesCol <- grDevices::rgb2hsv(grDevices::col2rgb(col[3]))
                        samplesCol <- grDevices::hsv(h=samplesCol[1], s=samplesCol[2], v=samplesCol[3], a=min(1, 100 / nrow(x$bugsFit$sims.matrix[, var])))  #  OR: 1 / log10(nrow(x$bugsFit$sims.matrix[, var]))))
                        graphics::rug(x$bugsFit$sims.matrix[, var], col=samplesCol)
                    }
                }
                if (legend)
                    legend("topright", legend=estimateNames[estimates], lty=lty[estimates], col=col[estimates], bg=grDevices::hsv(s=0, a=0.5))

            } else
                warning("Variable '", var, "' not supported for PDF plots", sep="")


        } else {
            ## Series plot (showing non-linear relationship of rw terms and spread of spatial terms):

            # extract covariate name and set inla
            if (pmatch("beta", var, nomatch=0)) {
                covariate <- sub("beta.", "", var)
                fixedCovariate <- sub("beta.", "", make.names(var))
                inla <- muINLA

            } else {
                covariate <- sub("gamma.", "", var)
                fixedCovariate <- sub("gamma.", "", make.names(var))
                inla <- rhoINLA
            }

            # if only prior, skip plot
            if (!posterior && !inla) {
                warning("Cannot plot prior for '", var, "'", sep="")
                next()
            }

            # extract boundaries
            if (var %in% rwGroups)
                boundaries <- x$rwBoundaries[[var]]
            else if (pmatch("beta", var, nomatch=0))
                boundaries <- 1:length(x$muAdjacency$num)
            else
                boundaries <- 1:length(x$rhoAdjacency$num)

            # extract plotting data
            if (inla) {
                if (pmatch(c("beta"), var, nomatch=0)) {
                    inlaMeans <- x$inlaFits$muFit$summary.random[[fixedCovariate]][, 2]
                    inlaCI <- x$inlaFits$muFit$summary.random[[fixedCovariate]][, c(4, 6)]
                } else {
                    inlaMeans <- x$inlaFits$rhoFit$summary.random[[fixedCovariate]][, 2]
                    inlaCI <- x$inlaFits$rhoFit$summary.random[[fixedCovariate]][, c(4, 6)]
                }
            }
            if (posterior) {
                bugsSamples <- x$bugsFit$sims.matrix[, grep(var, colnames(x$bugsFit$sims.matrix), fixed=TRUE)]
                if (!samples) {
                    ## NOTE: could just use BUGS summary for this!
                    bugsMeans <- colMeans(bugsSamples)
                    bugsCI <- cbind(bottom=apply(bugsSamples, 2, quantile, prob=((1 - prob) / 2)),
                                    top=apply(bugsSamples, 2, quantile, prob=1 - ((1 - prob) / 2)))
                } else {
                    samplesCol <- grDevices::rgb2hsv(grDevices::col2rgb(col[3]))
                    samplesCol <- grDevices::hsv(h=samplesCol[1], s=samplesCol[2], v=samplesCol[3], a=min(1, 100 / nrow(bugsSamples)))  #  OR: 1 / log10(nrow(bugsSamples))))
                }
            }

            # set plot limits
            ylim <- c(Inf, -Inf)
            if (inla)
                ylim <- range(c(inlaMeans, inlaCI), na.rm=TRUE)
            if (posterior) {
                if (samples)
                    ylim <- c(min(ylim[1], bugsSamples), max(ylim[2], bugsSamples))
                else
                    ylim <- c(min(ylim[1], bugsCI), max(ylim[2], bugsCI))
            }

            # plot
            graphics::plot(0, 0, col="white", xlim=range(boundaries), ylim=ylim, xlab=covariate, ylab="", main=var)
            graphics::abline(h=0, col="grey80")
            if (posterior && samples)
                graphics::matplot(add=TRUE, boundaries, t(bugsSamples), col=samplesCol, t='l', pch=20, lty=1)
            if (var %in% rwGroups) {
                if (inla) {
                    graphics::lines(boundaries, inlaCI[, 1], col=col[2], lty=2)
                    graphics::lines(boundaries, inlaMeans, col=col[2], lty=1)
                    graphics::lines(boundaries, inlaCI[, 2], col=col[2], lty=2)
                }
                if (posterior && !samples) {
                    graphics::lines(boundaries, bugsCI[, 1], col=col[3], lty=2)
                    graphics::lines(boundaries, bugsMeans, col=col[3], lty=1)
                    graphics::lines(boundaries, bugsCI[, 2], col=col[3], lty=2)
                }

            } else {
                if (inla) {
                    graphics::points(boundaries, inlaCI[, 1], col=col[2], pch='-')
                    graphics::points(boundaries, inlaMeans, col=col[2], pch='-')
                    graphics::points(boundaries, inlaCI[, 2], col=col[2], pch='-')
                    for (i in 1:length(boundaries))
                        graphics::lines(rep(boundaries[i], 2), inlaCI[i, 1:2], col=col[2], lty=3)
                }
                if (posterior && !samples) {
                    graphics::points(boundaries, bugsCI[, 1], col=col[3], pch='-')
                    graphics::points(boundaries, bugsMeans, col=col[3], pch='-')
                    graphics::points(boundaries, bugsCI[, 2], col=col[3], pch='-')
                    for (i in 1:length(boundaries))
                        graphics::lines(rep(boundaries[i], 2), bugsCI[i, 1:2], col=col[3], lty=3)
                }
            }
            graphics::rug(x$modelMatrix[, covariate])
            if (legend) {
                if (pmatch("beta", var, nomatch=0))
                    est <- muEstimates
                else
                    est <- rhoEstimates
                est[1] <- FALSE  # force no prior
                legend("topright", legend=estimateNames[est], lty=1, col=col[est], bg=grDevices::hsv(s=0, a=0.5))
            }
        }
    }
}

