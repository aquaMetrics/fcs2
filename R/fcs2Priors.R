#' Prior Distributions in the FCS2 Model
#'
#' Summarises and plots the prior distribution for each parameter in the
#' \acronym{FCS2} model.  This can be used to aid the selection of parameters
#' for each prior.
#'
#' This function can also be used to translate a user-friendly specification of
#' each prior distribution to the form required by \acronym{INLA} and
#' \acronym{BUGS} to fit the model.
#'
#' Each prior distribution can be specified either by their \code{mean} and
#' standard deviation \code{sd}, \code{mean} and variance \code{var} or by
#' providing values for their parameters.
#'
#' @param parameters a named list of named vectors, each specifying the
#' parameter values for a single prior distribution.  See \sQuote{Details} for
#' how these are specified.  Alternatively, an \code{"fcs2Fit"} object can be
#' provided, as returned from \code{\link{fcs2FitModel}}.
#' @param print whether to print a summary of each prior distribution.
#' @param plot whether to plot each prior distribution.
#' @param prob a decimal between 0 and 1 giving the size of the confidence
#' interval displayed for each prior (if \code{print = TRUE}).
#' @param \dots further arguments passed to \code{\link{plot}} if plotting each
#' prior.
#' @return Invisibly returns a modified list of prior parameters with each
#' prior specified using standardised parameters, as required for
#' \acronym{INLA} and \acronym{BUGS}:
#'
#' The log-Normal prior for the shape parameter \eqn{r} is given by the mean
#' \code{mu} and precision \code{tau} of the equivalent Normal prior for
#' \eqn{log(r)}.\cr The Beta prior for the catch probability \eqn{q} is given
#' by two shape parameters \code{a} and \code{b}.\cr The Normal prior for a
#' linear variable \code{beta.}* or \code{gamma.}* is given by the \code{mean}
#' and \code{precision}.\cr The Gamma prior for a precision hyperparameter
#' \code{tau.}* or \code{phi.}* is given by the shape \code{a} and rate
#' \code{b}.
#'
#' Prior parameters specified for a scale hyperparameter \code{sigma.}* or
#' \code{nu.}* (following the \code{\link{hyperprior}} distribution) are
#' transformed to the corresponding precision parameter \code{tau.}* or
#' \code{phi.}* respectively.
#' @seealso \code{\link{fcs2FitModel}}
#' @export
#' @examples
#'
#' # Specify a range of priors by their parameters:
#' # each prior is plotted and summarised by default
#' fcs2Priors(list(r=c(mu=0, sigma=3),
#'                 q=c(a=0.5, b=0.5),
#'                 beta=c(mean=0, precision=0.5),
#'                 phi=c(a=1, b=0.3)))
#'
#' # Priors can also be specified by their mean and sd or var:
#' fcs2Priors(list(r=c(mean=3, sd=10),
#'                 q=c(mean=0.6, sd=0.2),
#'                 beta=c(mean=0, var=2),
#'                 sigma=c(mean=1.5, sd=3)))
#'
fcs2Priors <-
function(parameters = list(), print = TRUE, plot = TRUE, prob = 0.95, ...)
{
    # extract parameter list if fcs2Fit object provided
    if (class(parameters) == "fcs2Fit")
        parameters <- parameters$prior.parameters

    # add ... to parameter list
    parameters <- c(parameters, list(...))

    # set up multiple plots
    if (plot && length(parameters) > 1) {
        w <- ceiling(sqrt(length(parameters)))
        par.old <- par(mfrow=c(w - (length(parameters) <= w*(w-1)), w))
        on.exit(par(par.old))
    }

    # set colour for plots
    col <- "black"  # matches priors in plot.fcs2Fit

    for (i in 1:length(parameters)) {
        varName <- names(parameters)[i]  # variable name
        parNames <- names(parameters[[i]])  # parameter names

        if (varName == "q") {
            ## Catch probability q:

            # find parameters and mean & sd
            if (sum(ix <- pmatch(parNames, c("mean", "sd"), nomatch=0)) == 3) {
                mean <- parameters[[i]][ix == 1]
                sd <- parameters[[i]][ix == 2]
                a <- (mean ^ 2) * (1 - mean) / (sd ^ 2) - mean
                b <- a * (1 - mean) / mean

            } else if (sum(ix <- pmatch(parNames, c("mean", "var"), nomatch=0)) == 3) {
                mean <- parameters[[i]][ix == 1]
                var <- parameters[[i]][ix == 2]
                sd <- sqrt(var)
                a <- (mean ^ 2) * (1 - mean) / (sd ^ 2) - mean
                b <- a * (1 - mean) / mean

            } else if (sum(ix <- pmatch(parNames, c("a", "b"), nomatch=0)) == 3) {
                a <- parameters[[i]][ix == 1]
                b <- parameters[[i]][ix == 2]
                mean <- a / (a + b)
                sd <- sqrt(a * b / (a + b + 1)) / (a + b)

            } else {
                warning(paste("Parameters ", paste(parNames, collapse=", "), " not recognised for variable '", varName, "'", sep=''))
                next()
            }

            # set parameters
            parameters[[i]] <- c(a, b)
            names(parameters[[i]]) <- c("a", "b")
            ci <- qbeta(c((1 - prob) / 2, 1 - (1 - prob) / 2), a, b)

            # print
            if (print) {
                cat(varName, " ~ Beta(a = ", a, ", b = ", b, "):\n", sep="")
                cat("  mean   = ", mean, "\n", sep='')
                cat("  sd     = ", sd, "\n", sep='')
                cat("  ", round(prob * 100), "% CI = (", ci[1], ", ", ci[2], ")\n\n", sep='')
            }

            # plot
            if (plot) {
                x <- seq(0, 1, l=200)
                yplt <- dbeta(x, a, b)
                plot(x, yplt, t='l', ylim=c(0, max(yplt[yplt < Inf])), xlab="", ylab="Density", main=varName, col=col, lty=2)
                abline(h=0, v=c(0, 1), col="grey80")
            }

        } else if (varName == "r") {
            ## Shape parameter r:

            # find parameters and mean & sd
            if (sum(ix <- pmatch(parNames, c("mean", "sd"), nomatch=0)) == 3) {
                mean <- parameters[[i]][ix == 1]
                sd <- parameters[[i]][ix == 2]
                mu <- log(mean) - log(1 + ((sd / mean) ^ 2)) / 2
                tau <- 1 / log(1 + ((sd / mean) ^ 2))

            } else if (sum(ix <- pmatch(parNames, c("mean", "var"), nomatch=0)) == 3) {
                mean <- parameters[[i]][ix == 1]
                var <- parameters[[i]][ix == 2]
                sd <- sqrt(var)
                mu <- log(mean) - log(1 + ((sd / mean) ^ 2)) / 2
                tau <- 1 / log(1 + ((sd / mean) ^ 2))

            } else if (sum(ix <- pmatch(parNames, c("mu", "tau"), nomatch=0)) == 3) {
                mu <- parameters[[i]][ix == 1]
                tau <- parameters[[i]][ix == 2]
                mean <- exp(mu + 1 / (2 * tau))
                sd <- mean * sqrt(exp(1 / tau) - 1)

            } else if (sum(ix <- pmatch(parNames, c("mu", "sigma"), nomatch=0)) == 3 ||
                       sum(ix <- pmatch(parNames, c("location", "scale"), nomatch=0)) == 3) {
                mu <- parameters[[i]][ix == 1]
                sigma <- parameters[[i]][ix == 2]
                tau <- 1 / (sigma ^ 2)
                mean <- exp(mu + 1 / (2 * tau))
                sd <- mean * sqrt(exp(1 / tau) - 1)

            } else {
                warning(paste("Parameters ", paste(parNames, collapse=", "), " not recognised for variable '", varName, "'", sep=''))
                next()
            }

            # set parameters
            parameters[[i]] <- c(mu, tau)
            names(parameters[[i]]) <- c("mu", "tau")
            ci <- qlnorm(c((1 - prob) / 2, 1 - (1 - prob) / 2), mu, 1 / sqrt(tau))

            # print
            if (print) {
                cat(varName, " ~ log-Normal(mu = ", mu, ", sigma = ", 1 / sqrt(tau), "):\n", sep="")
                cat("  mean   = ", mean, "\n", sep='')
                cat("  sd     = ", sd, "\n", sep='')
                cat("  ", round(prob * 100), "% CI = (", ci[1], ", ", ci[2], ")\n\n", sep='')
            }

            # plot
            if (plot) {
                x <- seq(0, qlnorm(0.995, mu, 1 / sqrt(tau)), l=200)
                plot(x, dlnorm(x, mu, 1 / sqrt(tau)), t='l', xlab="", ylab="Density", main=varName, col=col, lty=2)
                abline(h=0, v=0, col="grey80")
            }

        } else if (max(pmatch(c("beta", "gamma"), varName, nomatch=0)) > 0) {
            ## Linear terms:

            # find parameters and mean & sd
            if (sum(ix <- pmatch(parNames, c("mean", "sd"), nomatch=0)) == 3) {
                mean <- parameters[[i]][ix == 1]
                sd <- parameters[[i]][ix == 2]
                prec <- 1 / (sd ^ 2)

            } else if (sum(ix <- pmatch(parNames, c("mean", "var"), nomatch=0)) == 3) {
                mean <- parameters[[i]][ix == 1]
                var <- parameters[[i]][ix == 2]
                sd <- sqrt(var)
                prec <- 1 / var

            } else if (sum(ix <- pmatch(parNames, c("mean", "precision"), nomatch=0)) == 3) {
                mean <- parameters[[i]][ix == 1]
                prec <- parameters[[i]][ix == 2]
                sd <- sqrt(1 / prec)

            } else {
                warning(paste("Parameters ", paste(parNames, collapse=", "), " not recognised for variable '", varName, "'", sep=''))
                next()
            }

            # set parameters
            parameters[[i]] <- c(mean, prec)
            names(parameters[[i]]) <- c("mean", "precision")
            ci <- qnorm(c((1 - prob) / 2, 1 - (1 - prob) / 2), mean, sd)

            # print
            if (print) {
                cat(varName, " ~ Normal(mean = ", mean, ", sd = ", sd, "):\n", sep="")
                cat("  mean   = ", mean, "\n", sep='')
                cat("  sd     = ", sd, "\n", sep='')
                cat("  ", round(prob*100), "% CI = (", ci[1], ", ", ci[2], ")\n\n", sep='')
            }

            # plot
            if (plot) {
                x <- seq(qnorm(0.005, mean, sd), qnorm(0.995, mean, sd), l=200)
                plot(x, dnorm(x, mean, sd), t='l', ylim=c(0, dnorm(mean, mean, sd)), xlab="", ylab="Density", main=varName, col=col, lty=2)
                abline(h=0, col="grey80")
            }

        } else if (max(pmatch(c("sigma", "nu", "tau", "phi"), varName, nomatch=0)) > 0) {
            ## CAR terms:

            # set names for precision and transformed scale variables
            if (max(pmatch(c("sigma", "nu"), varName, nomatch=0)) > 0) {
                scaleVarName <- varName
                if (pmatch("sigma", varName, nomatch=0))
                    precVarName <- sub("sigma", "tau", varName)
                else
                    precVarName <- sub("nu", "phi", varName)

            } else {
                if (pmatch("tau", varName, nomatch=0))
                    scaleVarName <- sub("tau", "sigma", varName)
                else
                    scaleVarName <- sub("phi", "nu", varName)
                precVarName <- varName
            }

            # find parameters and mean & sd
            if (max(pmatch(c("sigma", "nu"), varName, nomatch=0)) > 0) {
                # (prior given for scale variable)
                if (sum(ix <- pmatch(parNames, c("mean", "sd"), nomatch=0)) == 3) {
                    scaleMean <- parameters[[i]][ix == 1]
                    scaleSd <- parameters[[i]][ix == 2]

                    # solve for a by finding zero of function
                    fctn <- function(a, m, s)
                        1 - a + (m^2) * exp(2 * (lgamma(a) - lgamma(a - 0.5))) / (m^2 + s^2)
                    a <- uniroot(fctn, c(1, 1e6), m=scaleMean, s=scaleSd)$root
                    b <- (a - 1) * (scaleMean^2 + scaleSd^2)

                } else if (sum(ix <- pmatch(parNames, c("mean", "var"), nomatch=0)) == 3) {
                    scaleMean <- parameters[[i]][ix == 1]
                    scaleVar <- parameters[[i]][ix == 2]
                    scaleSd <- sqrt(scaleVar)

                    # solve for a by finding zero of function
                    fctn <- function(a, m, s)
                        1 - a + (m^2) * exp(2 * (lgamma(a) - lgamma(a - 0.5))) / (m^2 + s^2)
                    a <- uniroot(fctn, c(1, 1e6), m=scaleMean, s=scaleSd)$root
                    b <- (a - 1) * (scaleMean^2 + scaleSd^2)

                } else if (sum(ix <- pmatch(parNames, c("a", "b"), nomatch=0)) == 3 ||
                           sum(ix <- pmatch(parNames, c("shape", "scale"), nomatch=0)) == 3) {
                    a <- parameters[[i]][ix == 1]
                    b <- parameters[[i]][ix == 2]
                    scaleMean <- ifelse(a > 0.5, sqrt(b) * exp(lgamma(a - 0.5) - lgamma(a)), Inf)
                    scaleSd <- ifelse(a > 1, sqrt(b / (a - 1) - scaleMean^2), Inf)

                } else {
                    warning(paste("Parameters ", paste(parNames, collapse=", "), " not recognised for variable '", varName, "'", sep=''))
                    next()
                }
                precMean <- a / b
                precSd <- sqrt(a) / b

            } else {
                # (prior given for precision variable)
                if (sum(ix <- pmatch(parNames, c("mean", "sd"), nomatch=0)) == 3) {
                    precMean <- parameters[[i]][ix == 1]
                    precSd <- parameters[[i]][ix == 2]
                    a <- (precMean^2) / (precSd ^ 2)
                    b <- precMean / (precSd ^ 2)

                } else if (sum(ix <- pmatch(parNames, c("mean", "var"), nomatch=0)) == 3) {
                    precMean <- parameters[[i]][ix == 1]
                    precVar <- parameters[[i]][ix == 2]
                    a <- (precMean^2) / precVar
                    b <- precMean / precVar

                } else if (sum(ix <- pmatch(parNames, c("a", "b"), nomatch=0)) == 3 ||
                           sum(ix <- pmatch(parNames, c("shape", "rate"), nomatch=0)) == 3) {
                    a <- parameters[[i]][ix == 1]
                    b <- parameters[[i]][ix == 2]
                    precMean <- a / b
                    precSd <- sqrt(a) / b

                }  else if (sum(ix <- pmatch(parNames, c("shape", "scale"), nomatch=0)) == 3) {
                    a <- parameters[[i]][ix == 1]
                    b <- 1 / parameters[[i]][ix == 2]
                    precMean <- a / b
                    precSd <- sqrt(a) / b

                } else {
                    warning(paste("Parameters ", paste(parNames, collapse=", "), " not recognised for variable '", varName, "'", sep=''))
                    next()
                }
                scaleMean <- ifelse(a > 0.5, sqrt(b) * exp(lgamma(a - 0.5) - lgamma(a)), Inf)
                scaleSd <- ifelse(a > 1, sqrt(b / (a - 1) - scaleMean^2), Inf)
            }

            # set parameters
            parameters[[i]] <- c(a, b)
            names(parameters[[i]]) <- c("a", "b")
            names(parameters)[i] <- precVarName
            scaleCI <- qhyperprior(c((1 - prob) / 2, 1 - (1 - prob) / 2), a, b)

            # print
            if (print) {
                cat(precVarName, " ~ Gamma(shape = ", a, ", rate = ", b, ")\n", sep="")
                cat("  mean   = ", precMean, "\n", sep='')
                cat("  sd     = ", precSd, "\n", sep='')
                cat(scaleVarName, " = 1 / sqrt(", precVarName, "):\n", sep='')
                cat("  mean   = ", scaleMean, "\n", sep='')
                cat("  sd     = ", scaleSd, "\n", sep='')
                cat("  ", round(prob*100), "% CI = (", scaleCI[1], ", ", scaleCI[2], ")\n\n", sep='')
            }

            # plot
            if (plot) {
                x <- seq(0, qhyperprior(0.995, a, b), l=200)
                plot(x, dhyperprior(x, a, b), t='l', xlab="", ylab="Density", main=scaleVarName, col=col, lty=2)
                abline(h=0, v=0, col="grey80")
            }

        } else {
            warning(paste("Variable '", varName, "' not recognised", sep=''))
            next()
        }
    }

    invisible(parameters)
}

