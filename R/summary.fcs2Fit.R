#' Detailed summary
#'
#' Provides a more detailed summary
#'
#' @param object summary object
#' @param prior from summary object
#' @param inla from summary object
#' @param posterior from summary object
#' @param allVars boolean
#'
#' @return detailed summary object
#' @seealso summary.fcs2Fit summary.fcs2Fit print.summary.fcs2Fit
#' @export
#'
#' @examples
#' \dontrun{
#' summary.fcs2Fit(object = data)
#' }
summary.fcs2Fit <-
function(object, prior = TRUE, inla = is.null(object$bugsFit), posterior = !is.null(object$bugsFit), allVars = FALSE)
{
    ret <- list()

    # call
    ret$call <- object$call

    # catch data
    ret$multiRun <- object$multiRun
    if (!object$multiRun) {
        vars <- character(0)
        nValues <- nNonZero <- numeric(0)

        if (sum(object$dataType == "run" || object$dataType == "total", na.rm=TRUE) > 0) {
            ret$allRunsTotalVar <- object$allRunsTotalVar
            vars <- object$allRunsTotalVar
            nValues <- sum(object$dataType == "run" | object$dataType == "total", na.rm=TRUE)
            nNonZero <- sum(object$modelMatrix[object$dataType == "run" | object$dataType == "total", object$allRunsTotalVar] > 0, na.rm=TRUE)
        }
        if (sum(object$dataType == "range", na.rm=TRUE) > 0) {
            ret$allRunsRangeVars <- object$allRunsRangeVars
            vars <- c(vars, object$allRunsRangeVars)
            nValues <- c(nValues, rep(sum(object$dataType == "range", na.rm=TRUE), length(object$allRunsRangeVars)))
            nNonZero <- c(nNonZero, rep(sum(object$modelMatrix[object$dataType == "range", object$allRunsRangeVars[1]] > 0, na.rm=TRUE), length(object$allRunsRangeVars)))
        }

        ret$catchData <- data.frame(values=nValues, "non-zero"=nNonZero, check.names=FALSE, row.names=vars)

    } else {
        vars <- character(0)
        nValues <- nNonZero <- numeric(0)

        if (sum(object$dataType == "run", na.rm=TRUE) > 0) {
            ret$runTotalVars <- object$runTotalVars
            vars <- object$runTotalVars
            for (i in 1:length(object$runTotalVars)) {
                nValues <- c(nValues, sum(object$modelMatrix[, object$nRunsVar] >= i & object$dataType == "run", na.rm=TRUE))
                nNonZero <- c(nNonZero, sum(object$modelMatrix[object$modelMatrix[, object$nRunsVar] >= i & object$dataType == "run", object$runTotalVars[i]] > 0, na.rm=TRUE))
            }
        }
        if (sum(object$dataType == "total", na.rm=TRUE) > 0) {
            ret$allRunsTotalVar <- object$allRunsTotalVar
            vars <- c(vars, object$allRunsTotalVar)
            nValues <- c(nValues, sum(object$dataType == "total", na.rm=TRUE))
            nNonZero <- c(nNonZero, sum(object$modelMatrix[object$dataType == "total", object$allRunsTotalVar] > 0, na.rm=TRUE))
        }
        if (sum(object$dataType == "range", na.rm=TRUE) > 0) {
            ret$allRunsRangeVars <- object$allRunsRangeVars
            vars <- c(vars, object$allRunsRangeVars)
            nValues <- c(nValues, rep(sum(object$dataType == "range", na.rm=TRUE), length(object$allRunsRangeVars)))
            nNonZero <- c(nNonZero, rep(sum(object$modelMatrix[object$dataType == "range", object$allRunsRangeVars[1]] > 0, na.rm=TRUE), length(object$allRunsRangeVars)))
        }

        ret$catchData <- data.frame(values=nValues, "non-zero"=nNonZero, check.names=FALSE, row.names=vars)
        ret$nRunsVar <- object$nRunsVar
    }

    # other data
    ret$surveyAreaVar <- object$surveyAreaVar
    ret$N <- object$N
    if (length(object$na.action) > 0 && class(object$na.action) == "omit")
        ret$nRemoved <- length(object$na.action)
    else if (sum(is.na(object$dataType)) > 0)
        ret$nMissing <- sum(is.na(object$dataType))


    # model formulae (mu and rho)
    ret$muFormula <- object$muFormula
    ret$rhoFormula <- object$rhoFormula

    # prior summaries
    if (prior) {
        params <- variable.names(object, rw=FALSE, spatial=FALSE, hyperparams="scale")
        mx <- array(dim=c(length(params), 5), dimnames=list(params, c("mean", "sd", "2.5%", "50%", "97.5%")))

        for (i in 1:nrow(mx)) {
            if (pmatch("q", params[i], nomatch=0) > 0) {
                a <- object$prior.parameters[[params[i]]]["a"]
                b <- 1 / sqrt(object$prior.parameters[[params[i]]]["b"])
                mx[i, ] <- c(a / (a + b), sqrt(a * b / (a + b + 1)) / (a + b),
                             qbeta(c(0.025, 0.5, 0.975), a, b))

            } else if (pmatch("r", params[i], nomatch=0) > 0) {
                mu <- object$prior.parameters[[params[i]]]["mu"]
                sigma <- 1 / sqrt(object$prior.parameters[[params[i]]]["tau"])
                mx[i, ] <- c(exp(mu + ((sigma ^ 2) / 2)), exp(mu + ((sigma ^ 2) / 2)) * sqrt(exp(sigma ^ 2) - 1),
                             qlnorm(c(0.025, 0.5, 0.975), mu, sigma))

            } else if(max(pmatch(c("tau", "phi"), params[i], nomatch=0)) > 0) {
                a <- object$prior.parameters[[params[i]]]["a"]
                b <- object$prior.parameters[[params[i]]]["b"]
                mx[i, ] <- c(a / b, sqrt(a) / b, qgamma(c(0.025, 0.5, 0.975), a, b))

            } else if (max(pmatch(c("beta", "gamma"), params[i], nomatch=0)) > 0) {
                mean <- object$prior.parameters[[params[i]]]["mean"]
                sd <- 1 / sqrt(object$prior.parameters[[params[i]]]["precision"])
                mx[i, ] <- c(mean, sd, qnorm(c(0.025, 0.5, 0.975), mean, sd))

            } else if (max(pmatch(c("sigma", "nu"), params[i], nomatch=0)) > 0) {
                if (pmatch("sigma", params[i], nomatch=0))
                    precVarName <- sub("sigma", "tau", params[i])
                else
                    precVarName <- sub("nu", "phi", params[i])
                a <- object$prior.parameters[[precVarName]]["a"]
                b <- object$prior.parameters[[precVarName]]["b"]
                mean <- ifelse(a > 0.5, sqrt(b) * exp(lgamma(a - 0.5) - lgamma(a)), Inf)
                sd <- ifelse(a > 1, sqrt(b / (a - 1) - mean^2), Inf)
                mx[i, ] <- c(mean, sd, qhyperprior(c(0.025, 0.5, 0.975), a, b))
            }
        }

        ret$prior <- mx
    }


    ## INLA posterior summaries
    if (inla && !is.null(object$inlaFits)) {
        ret$inla <- list()

        # abundance mu
        if (!is.null(object$inlaFits$muFit)) {
            mx <- rbind(object$inlaFits$muFit$summary.hyperpar[1, , drop=FALSE],
                        object$inlaFits$muFit$summary.fixed[, -6, drop=FALSE])

            # add significance columns if linear or individual random terms
            if (nrow(object$inlaFits$muFit$summary.fixed) > 1 || (allVars && length(object$inlaFits$muFit$summary.random) > 0))
                mx <- cbind(mx, "Signif. Pr."=NA)

            # add significance probabilities for linear terms
            if (nrow(object$inlaFits$muFit$summary.fixed) > 1) {
                for (i in 1:(nrow(object$inlaFits$muFit$summary.fixed) - 1))
                    mx[i + 2, 6] <- .significanceProbability.inla(object$inlaFits$muFit$marginals.fixed[[i + 1]])
            }

            # add individual random walk and spatial parameters
            if (allVars && length(object$inlaFits$muFit$summary.random) > 0) {
                for (i in 1:length(object$inlaFits$muFit$summary.random)) {
                    if (ncol(mx) == 6) {
                        mx <- rbind(mx, cbind(object$inlaFits$muFit$summary.random[[i]][, -c(1, 7)], "Signif. Pr."=NA))
                        for (j in 1:nrow(object$inlaFits$muFit$summary.random[[i]]))
                            mx[nrow(mx) - nrow(object$inlaFits$muFit$summary.random[[i]]) + j, 6] <-
                                    .significanceProbability.inla(object$inlaFits$muFit$marginals.random[[i]][[j]])

                    } else
                        mx <- rbind(mx, object$inlaFits$muFit$summary.random[[i]][, -c(1, 7)])
                }
            }

            # add transformed hyperparameters
            if (nrow(object$inlaFits$muFit$summary.hyperpar) > 1) {
                t <- function(object)
                    sqrt(1 / object)
                t2 <- function(object)
                    1 / object
                for (i in 1:(nrow(object$inlaFits$muFit$summary.hyperpar) - 1)) {
                    mean <- inla.expectation(t, object$inlaFits$muFit$marginals.hyperpar[[i + 1]])
                    sd <- sqrt(inla.expectation(t2, object$inlaFits$muFit$marginals.hyperpar[[i + 1]]) - mean^2)
                    if (ncol(mx) == 6)
                        mx <- rbind(mx, c(mean, sd, 1 / sqrt(object$inlaFits$muFit$summary.hyperpar[i + 1, 5:3, drop=FALSE]), "Signif. Pr."=NA))
                    else
                        mx <- rbind(mx, c(mean, sd, 1 / sqrt(object$inlaFits$muFit$summary.hyperpar[i + 1, 5:3, drop=FALSE])))
                }
            }

            # add names
            rownames(mx) <- variable.names(object, q=FALSE, r=TRUE, prevalence=FALSE, rw=ifelse(allVars, "singular", FALSE),
                                           spatial=ifelse(allVars, "singular", FALSE), hyperparams="scale")
            colnames(mx)[1:5] <- c("mean", "sd", "2.5%", "50%", "97.5%")

            # store
            ret$inla$mu <- mx

            # number fitted and DIC
            ret$inla$muN <- nrow(object$inlaFits$muFit$data)
            if (!is.null(object$inlaFits$muFit$dic))
                ret$inla$muDIC <- object$inlaFits$muFit$dic[4,]
        }


        # prevalence rho
        mx <- object$inlaFits$rhoFit$summary.fixed[, -6, drop=FALSE]

        # add significance columns if linear or individual random terms
        if (nrow(object$inlaFits$rhoFit$summary.fixed) > 1 || (allVars && length(object$inlaFits$rhoFit$summary.random) > 0))
            mx <- cbind(mx, "Signif. Pr."=NA)

        # add significance probabilities for linear terms
        if (nrow(object$inlaFits$rhoFit$summary.fixed) > 1) {
            for (i in 1:(nrow(object$inlaFits$rhoFit$summary.fixed) - 1))
                mx[i + 1, 6] <- .significanceProbability.inla(object$inlaFits$rhoFit$marginals.fixed[[i + 1]])
        }

        # add individual random walk and spatial parameters
        if (allVars && length(object$inlaFits$rhoFit$summary.random) > 0) {
            for (i in 1:length(object$inlaFits$rhoFit$summary.random)) {
                if (ncol(mx) == 6) {
                    mx <- rbind(mx, cbind(object$inlaFits$rhoFit$summary.random[[i]][, -c(1, 7)], "Signif. Pr."=NA))
                    for (j in 1:nrow(object$inlaFits$rhoFit$summary.random[[i]]))
                        mx[nrow(mx) - nrow(object$inlaFits$rhoFit$summary.random[[i]]) + j, 6] <-
                                .significanceProbability.inla(object$inlaFits$rhoFit$marginals.random[[i]][[j]])

                } else
                    mx <- rbind(mx, object$inlaFits$rhoFit$summary.random[[i]][, -c(1, 7)])
            }
        }

        # add transformed hyperparameters
        if (length(object$inlaFits$rhoFit$marginals.hyperpar) > 0) {
            t <- function(object)
                sqrt(1 / object)
            t2 <- function(object)
                1 / object
            for (i in 1:nrow(object$inlaFits$rhoFit$summary.hyperpar)) {
                mean <- inla.expectation(t, object$inlaFits$rhoFit$marginals.hyperpar[[i]])
                sd <- sqrt(inla.expectation(t2, object$inlaFits$rhoFit$marginals.hyperpar[[i]]) - mean^2)
                if (ncol(mx) == 6)
                    mx <- rbind(mx, c(mean, sd, 1 / sqrt(object$inlaFits$rhoFit$summary.hyperpar[i, 5:3, drop=FALSE]), "Signif. Pr."=NA))
                else
                    mx <- rbind(mx, c(mean, sd, 1 / sqrt(object$inlaFits$rhoFit$summary.hyperpar[i, 5:3, drop=FALSE])))
            }
        }

        # add names
        rownames(mx) <- variable.names(object, q=FALSE, r=FALSE, abundance=FALSE, rw=ifelse(allVars, "singular", FALSE),
                                       spatial=ifelse(allVars, "singular", FALSE), hyperparams="scale")
        colnames(mx)[1:5] <- c("mean", "sd", "2.5%", "50%", "97.5%")

        # store
        ret$inla$rho <- mx

        # number fitted and DIC
        ret$inla$rhoN <- nrow(object$inlaFits$rhoFit$data)
        if (!is.null(object$inlaFits$rhoFit$dic))
            ret$inla$rhoDIC <- object$inlaFits$rhoFit$dic[4,]
    }

    ## BUGS
    if (posterior && !is.null(object$bugsFit)) {
        ## BUGS posterior summaries (mean, sd, 0.025quant, median, 0.975quant)
        params <- variable.names(object, rw=ifelse(allVars, "singular", FALSE), spatial=ifelse(allVars, "singular", FALSE), hyperparams="precision")
        mx <- object$bugsFit$summary[params, -c(4, 6)]  # remove 25% and 75%, for now

        # replace hyperparameter precision means with scale means
        hyperVars <- variable.names(object, r=FALSE, q=FALSE, linear=FALSE, rw=FALSE, spatial=FALSE, hyperparams="precision")
        if (length(hyperVars) > 0) {
            for (var in hyperVars)
                mx[var, 1:5] <- c(mean(1 / sqrt(object$bugsFit$sims.list[[var]])), sd(1 / sqrt(object$bugsFit$sims.list[[var]])),
                                  1 / sqrt(mx[var, 5:3]))
            rownames(mx) <- sub("tau", "sigma", rownames(mx), fixed=TRUE)
            rownames(mx) <- sub("phi", "nu", rownames(mx), fixed=TRUE)
        }

        # add significance probabilities and codes for linear or individual random terms
        sigTerms <- c(variable.names(object, r=FALSE, q=FALSE, prevalence=FALSE, rw=ifelse(allVars, "singular", FALSE),
                                     spatial=ifelse(allVars, "singular", FALSE), hyperparams=FALSE)[-1],
                      variable.names(object, r=FALSE, q=FALSE, abundance=FALSE, rw=ifelse(allVars, "singular", FALSE),
                                     spatial=ifelse(allVars, "singular", FALSE), hyperparams=FALSE)[-1])
        if (length(sigTerms) > 0) {
            mx <- cbind(mx, "Signif. Pr."=NA)
            for (term in sigTerms)
                mx[term, ncol(mx)] <- .significanceProbability.bugs(object$bugsFit$sims.matrix[, term])
        }

        # store
        ret$posterior <- mx

        # DIC
        if (!is.null(object$bugsFit$DIC))
            ret$DIC <- object$bugsFit$DIC

        ## BUGS settings
        ret$bugsProgram <- object$bugsFit$program
        ret$n.chains <- object$bugsFit$n.chains
        ret$n.iter <- object$bugsFit$n.iter
        ret$n.burnin <- object$bugsFit$n.burnin
        ret$n.sims <- object$bugsFit$n.sims
    }

    # return summary object
    class(ret) <- "summary.fcs2Fit"
    ret
}


## print.summary.fcs2Fit
print.summary.fcs2Fit <-
function(x, signif.stars = getOption("show.signif.stars"))
{
    cat("\n")

    # data
    cat("Data:\n")
    if (!x$multiRun) {
        if (!is.null(x$allRunsTotalVar))
            cat("Catch total:    ", x$allRunsTotalVar, " (", x$catchData[x$allRunsTotalVar, 1], " values, ",
                x$catchData[x$allRunsTotalVar, 2], " non-zero)\n", sep="")
        if (!is.null(x$allRunsRangeVars))
            cat("Total min, max: ", paste(x$allRunsRangeVars, collapse=", "), " (", x$catchData[x$allRunsRangeVars[1], 1], " values, ",
                x$catchData[x$allRunsRangeVars[1], 2], " non-zero)\n", sep="")
        ### NOTE: ignoring runTotalVars[1] since with single run model this should be the same as allRunsTotalVar

    } else {
        if (!is.null(x$runTotalVars)) {
            for (i in 1:length(x$runTotalVars))
                cat("Run ", i, " total:    ", x$runTotalVars[i], " (", x$catchData[x$runTotalVars[i], 1], " values, ",
                    x$catchData[x$runTotalVars[i], 2], " non-zero)\n", sep="")
        }
        if (!is.null(x$allRunsTotalVar))
            cat("All runs total: ", x$allRunsTotalVar, " (", x$catchData[x$allRunsTotalVar, 1], " values, ",
                x$catchData[x$allRunsTotalVar, 2], " non-zero)\n", sep="")
        if (!is.null(x$allRunsRangeVars))
            cat("Total min, max: ", paste(x$allRunsRangeVars, collapse=", "), " (", x$catchData[x$allRunsRangeVars[1], 1], " values, ",
                x$catchData[x$allRunsRangeVars[1], 2], " non-zero)\n", sep="")
        cat("Number of runs:", x$nRunsVar, "\n")
    }
    cat("Survey area:   ", x$surveyAreaVar, "\n")
    cat(x$N, " values", sep='')
    if (!is.null(x$nRemoved))
        cat(" (", x$nRemoved, " removed as missing)", sep="")
    else if (!is.null(x$nMissing))
        cat(" (", x$nMissing, " missing)", sep="")
    cat("\n\n")

    # model formulae (mu and rho)
    cat("Abundance formula:\n")
    print(formula(paste("log(mu)", paste(deparse(x$muFormula), collapse=" "))), showEnv=FALSE)
    cat("\n")

    # model formulae (mu and rho)
    cat("Prevalence formula:\n")
    print(formula(paste("logit(rho)", paste(deparse(x$rhoFormula), collapse=" "))), showEnv=FALSE)
    cat("\n")

    # prior summaries
    if (!is.null(x$prior)) {
        cat("Prior:\n")
        print(x$prior)
        cat("\n")
    }


    ## INLA posterior summaries
    if (!is.null(x$inla)) {

        # abundance mu
        if (!is.null(x$inla$mu)) {
            cat("Approximate posterior from INLA summary to abundance:\n")

            # check whether summary has significance probs
            if ("Signif. Pr." %in% colnames(x$inla$mu)) {
                muSummary <- data.frame(x$inla$mu, check.names=FALSE)
                muSummary[, "Signif. Pr."] <- format(muSummary[, "Signif. Pr."])
                muSummary[is.na(x$inla$mu[, "Signif. Pr."]), "Signif. Pr."] <- ""

                # add significance stars
                if(signif.stars)
                    muSummary <- cbind(muSummary, " "=.significanceStars(x$inla$mu[, "Signif. Pr."]))

                # print summary
                print(muSummary)

            } else
                print(x$inla$mu)

            # number fitted, significance stars and DIC
            cat("---\n")
            if ("Signif. Pr." %in% colnames(x$inla$mu) && signif.stars)
                cat("Singif. codes: 0 0.001 0.01 0.05 0.1 1\n")
            cat(x$inla$muN, " values fitted (those observed or predicted to be present by INLA prevalence summary)\n", sep="")
            if (!is.null(x$inla$muDIC))
                cat("DIC: ", x$inla$muDIC, "\n", sep="")
            cat("\n")
        }


        # prevalence rho
        cat("Approximate posterior from INLA summary to prevalence:\n")

        # check whether summary has significance probs
        if ("Signif. Pr." %in% colnames(x$inla$rho)) {
            rhoSummary <- data.frame(x$inla$rho, check.names=FALSE)
            rhoSummary[, "Signif. Pr."] <- format(rhoSummary[, "Signif. Pr."])
            rhoSummary[is.na(x$inla$rho[, "Signif. Pr."]), "Signif. Pr."] <- ""

            # add significance stars
            if(signif.stars)
                rhoSummary <- cbind(rhoSummary, " "=.significanceStars(x$inla$rho[, "Signif. Pr."]))

            # print summary
            print(rhoSummary)

        } else
            print(x$inla$rho)

        # DIC
        if (!is.null(x$inla$rhoDIC) || ("Signif. Pr." %in% colnames(x$inla$rho) && signif.stars))
            cat("---\n")
        if ("Signif. Pr." %in% colnames(x$inla$rho) && signif.stars)
            cat("Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1\n")
        if (!is.null(x$inla$rhoDIC))
            cat("DIC: ", x$inla$rhoDIC, "\n", sep="")
        cat("\n")
    }


    ## BUGS
    if (!is.null(x$posterior)) {
        # BUGS posterior summaries (mean, sd, 0.025quant, median, 0.975quant)
        cat("Posterior from ", x$bugsProgram, ":\n", sep="")

        # check whether summary has significance probs
        if ("Signif. Pr." %in% colnames(x$posterior)) {
            postSummary <- data.frame(x$posterior, check.names=FALSE)
            postSummary[, "Signif. Pr."] <- format(postSummary[, "Signif. Pr."])
            postSummary[is.na(x$posterior[, "Signif. Pr."]), "Signif. Pr."] <- ""

            # add significance stars
            if(signif.stars)
                postSummary <- cbind(postSummary, " "=.significanceStars(x$posterior[, "Signif. Pr."]))

            # print summary
            print(postSummary)

        } else
            print(x$posterior)

        # DIC
        if (!is.null(x$DIC) || ("Signif. Pr." %in% colnames(x$posterior) && signif.stars))
            cat("---\n")
        if ("Signif. Pr." %in% colnames(x$posterior) && signif.stars)
            cat("Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1\n")
        if (!is.null(x$DIC))
            cat("DIC: ", x$DIC, "\n", sep="")
        cat("\n")

        # BUGS settings
        cat(x$bugsProgram, " settings:\n", sep="")
        cat(x$n.chains, " chains, each with ", x$n.iter, " iterations (first ", x$n.burnin, " discarded)\n", sep="")
        cat(x$n.sims, " Monte Carlo samples saved\n", sep="")
        cat("\n")
    }

    # invisibly return summary
    invisible(x)
}

