#' Write BUGS Model File
#'
#' Creates a \acronym{BUGS} model file that describes the structure of the full
#' \acronym{FCS2} statistical model.  This can be used by WinBUGS or OpenBUGS
#' to fit the model.
#'
#'
#' @param filename the name or the \acronym{BUGS} model file to create.
#' @param muLinearVars a character vector of linear terms in the abundance
#' regression
#' @param muRWVars a character vector of variables with first or second-order
#' random walk terms in the abundance regression.
#' @param muSpatialVar the name of the abundance spatial variable, if present.
#' @param rhoLinearVars a character vector of linear terms in the abundance
#' regression
#' @param rhoRWVars a character vector of variables with first or second-order
#' random walk terms in the prevalence regression.
#' @param rhoSpatialVar the name of the prevalence spatial variable, if
#' present.
#' @param multiRun whether the multiple-run extension to the \acronym{FCS2}
#' model is required.
#' @param catchTotal whether catch total data are present in the model.
#' @param catchRange whether catch total range data are present in the model.
#' @param qUnif whether to use a uniform prior for the catch probability
#' \eqn{q} rather than the default beta prior.  Although the uniform is a
#' special case of the beta, this option may help \acronym{BUGS} fit the model.
#' @note The user would not usually call this function directly as it is called
#' if required by \code{\link{fcs2FitModel}}.  However, it may be useful for
#' creating a \acronym{BUGS} model file to be run by WinBUGS or OpenBUGS
#' extenally.
#' @seealso \code{\link{fcs2FitModel}}
#' @export
fcs2WriteModel <-
function(filename = "model.bug", muLinearVars = c(), muRWVars = c(), muSpatialVar = NULL,
         rhoLinearVars = c(), rhoRWVars = c(), rhoSpatialVar = NULL, multiRun = FALSE, catchTotal = TRUE, catchRange = FALSE, qUnif = FALSE)
{
    con <- file(filename, "w")

    # start file
    cat(file=con, "model {", "\n", sep="")

    ## Catch totals
    if (catchTotal) {

        # loop over number of catch total data points 'Nruns', unless no catch range data
        if (catchRange) {
            cat(file=con, "\t", "for (i in 1:Nruns) {", "\n", sep="")

        } else {
            # loop over number data points 'N'
            cat(file=con, "\t", "for (i in 1:N) {", "\n", sep="")

            # Define likelihood (using 'ones trick'):
            cat(file=con, "\t\t", "ones[i] <- 1", "\n",
                "\t\t", "ones[i] ~ dbern(likelihood[i])", "\n", sep="")
        }

        if (!multiRun) {
            # single run ZINB likelihood
            cat(file=con, "\t\t", "likelihood[i] <- (1 - rho[i]) * equals(catch[i], 0) +", "\n",
                "\t\t\t\t", "rho[i] * exp(loggam(catch[i] + r) - loggam(r) - loggam(catch[i] + 1) + r * log(r) +", "\n",
                "\t\t\t\t\t\t\t ", "catch[i] * log(area[i] * mu[i]) - (r + catch[i]) * log(r + area[i] * mu[i]))", "\n", sep='')

        } else {
            # multiple run ZINM likelihood
            cat(file=con, "\t\t", "likelihood[i] <- (1 - rho[i]) * equals(catchSum[i], 0) +", "\n",
                "\t\t\t\t", "rho[i] * exp(loggam(r + catchSum[i]) - loggam(r) - sum(logCatchFactorial[i, 1:nRuns[i]]) + r * log(r) +", "\n",
                "\t\t\t\t\t\t\t ", "sum(catchPart[i, 1:nRuns[i]]) - (r + catchSum[i]) * log(r + area[i] * mu[i] * (1 - pow(1 - q, nRuns[i]))))", "\n\n", sep='')

            # calculations over runs
            cat(file=con, "\t\t", "for (j in 1:maxNRuns) {", "\n",
                "\t\t\t", "logCatchFactorial[i, j] <- loggam(catch[i, j] + 1)", "\n",
                "\t\t\t", "catchPart[i, j] <- catch[i, j] * log(area[i] * mu[i] * q * pow(1 - q, j - 1))", "\n",
                "\t\t", "}", "\n", sep="")
        }

        # close bracket if catch range data
        if (catchRange)
            cat(file=con, "\t", "}", "\n\n", sep="")
        else
            cat(file=con, "\n", sep="")
    }


    ## Catch range data
    if (catchRange) {

        # loop over remaining number of data points, unless no catch range data
        if (catchTotal) {
            cat(file=con, "\t", "for (i in (Nruns+1):N) {", "\n", sep="")

        } else {
            # loop over number data points 'N'
            cat(file=con, "\t", "for (i in 1:N) {", "\n")

            # Define likelihood (using 'ones trick'):
            cat(file=con, "\t\t", "ones[i] <- 1", "\n",
                "\t\t", "ones[i] ~ dbern(likelihood[i])", "\n", sep="")
        }

        # ranged ZINB likelihood
        cat(file=con, "\t\t", "likelihood[i] <- (1 - rho[i]) * equals(totalCatchMin[i], 0) +", "\n",
            "\t\t\t\t", "rho[i] * exp(r * (log(r) - log(r + nbmean[i])) - loggam(r)) * sum(term[i, 1:(1 + totalCatchMax[i] - totalCatchMin[i])])", "\n\n",
            "\t\t", "logRevP[i] <- log(nbmean[i]) - log(r + nbmean[i])", "\n", sep='')

        # define mean
        if (multiRun)
            cat(file=con, "\t\t", "nbmean[i] <- area[i] * mu[i] * (1 - pow(1 - q, nRuns[i]))", "\n\n", sep="")
        else
            cat(file=con, "\t\t", "nbmean[i] <- area[i] * mu[i]", "\n\n", sep="")

        # calculations over possible total t
        cat(file=con, "\t\t", "for (t in totalCatchMin[i]:totalCatchMax[i]) {", "\n",
            "\t\t\t", "term[i, 1 + t - totalCatchMin[i]] <- exp(loggam(t + r) - loggam(t + 1) + t * logRevP[i])", "\n",
            "\t\t", "}", "\n", sep="")

        # if catch total data ...
        if (catchTotal) {
            # close bracket
            cat(file=con, "\t", "}", "\n\n", sep="")

            # loop over all points'N'
            cat(file=con, "\t", "for (i in 1:N) {", "\n")

            # Define likelihood (using 'ones trick'):
            cat(file=con, "\t\t", "ones[i] <- 1", "\n",
                "\t\t", "ones[i] ~ dbern(likelihood[i])", "\n\n", sep="")

        } else
            cat(file=con, "\n", sep="")
    }


    ## Regression for mu:
    cat(file=con, "\t\t", "log(mu[i]) <- beta.const", sep="")

    # linear mu terms
    if (length(muLinearVars) > 0) {
        for (var in muLinearVars) {
            cat(file=con, " +\n\t\t\t\t",
                "beta.", make.names(var), " * ", make.names(var), "[i]", sep='')
        }
    }

    # random walk mu terms
    if (length(muRWVars) > 0) {
        for (var in muRWVars) {
            cat(file=con, " +\n\t\t\t\t",
                "interp.lin(", make.names(var), "[i], boundaries.beta.", make.names(var), "[], beta.", make.names(var), "[])", sep='')
        }
    }

    # mu spatial term
    if (!is.null(muSpatialVar)) {
        cat(file=con, " +\n\t\t\t\t",
            "beta.", make.names(muSpatialVar), "[", make.names(muSpatialVar), "[i]]", sep='')
    }

    cat(file=con, "\n\n\t\t")


    ## Regression for rho:
    cat(file=con, "logit(rho[i]) <- gamma.const")

    # linear rho terms
    if (length(rhoLinearVars) > 0) {
        for (var in rhoLinearVars) {
            cat(file=con, " +\n\t\t\t\t",
                "gamma.", make.names(var), " * ", make.names(var), "[i]", sep='')
        }
    }

    # random walk rho terms
    if (length(rhoRWVars) > 0) {
        for (var in rhoRWVars) {
            cat(file=con, " +\n\t\t\t\t",
                "interp.lin(", make.names(var), "[i], boundaries.gamma.", make.names(var), "[], gamma.", make.names(var), "[])", sep='')
        }
    }

    # rho spatial term
    if (!is.null(rhoSpatialVar)) {
        cat(file=con, " +\n\t\t\t\t",
            "gamma.", make.names(rhoSpatialVar), "[", make.names(rhoSpatialVar), "[i]]", sep='')
    }

    cat(file=con, "\n\t")


    # end loop
    cat(file=con, "}\n\n\t")


    ## Priors:

    # prior for catch probability q
    if (multiRun) {
        if (qUnif)
            cat(file=con, "q ~ dunif(0, 1)\n\n\t")
        else
            cat(file=con, "q ~ dbeta(q.a, q.b)\n\n\t")
    }

    # prior for shape parameter r
    cat(file=con, "r ~ dlnorm(mu.r, tau.r)\n\n\t")


    ## Priors for mu:

    # constant mu term (normal)
    cat(file=con, "beta.const ~ dnorm(mean.beta.const, prec.beta.const)")

    # linear mu terms (normal)
    if (length(muLinearVars) > 0) {
        for (var in muLinearVars) {
            cat(file=con, "\n\n\t",
                "beta.", make.names(var), " ~ dnorm(mean.beta.", make.names(var), ", prec.beta.", make.names(var), ")", sep='')
        }
    }

    # random walk mu terms (gamma)
    if (length(muRWVars) > 0) {
        for (var in muRWVars) {
            cat(file=con, "\n\n\t",
                "beta.", make.names(var), "[1:nLevels.beta.", make.names(var), "] ~ car.normal(adj.beta.", make.names(var),
                        "[], weights.beta.", make.names(var), "[], num.beta.", make.names(var), "[], tau.", make.names(var), ")\n\t",
                "tau.", make.names(var), " ~ dgamma(a.tau.", make.names(var), ", b.tau.", make.names(var), ")", sep='')
        }
    }

    # mu spatial term (gamma)
    if (!is.null(muSpatialVar)) {
        cat(file=con, "\n\n\t",
            "beta.", make.names(muSpatialVar), "[1:nSpatials.", make.names(muSpatialVar), "] ~ car.normal(adj.", make.names(muSpatialVar), "[], weights.",
                     make.names(muSpatialVar), "[], num.", make.names(muSpatialVar), "[], tau.", make.names(muSpatialVar), ")\n\t",
            "tau.", make.names(muSpatialVar), " ~ dgamma(a.tau.", make.names(muSpatialVar), ", b.tau.", make.names(muSpatialVar), ")", sep='')
    }


    ## Priors for rho:

    # constant rho term (normal)
    cat(file=con, "\n\n\t",
        "gamma.const ~ dnorm(mean.gamma.const, prec.gamma.const)", sep='')

    # linear rho terms (normal)
    if (length(rhoLinearVars) > 0) {
        for (var in rhoLinearVars) {
            cat(file=con, "\n\n\t",
                "gamma.", make.names(var), " ~ dnorm(mean.gamma.", make.names(var), ", prec.gamma.", make.names(var), ")", sep='')
        }
    }

    # random walk rho terms (gamma)
    if (length(rhoRWVars) > 0) {
        for (var in rhoRWVars) {
            cat(file=con, "\n\n\t",
                "gamma.", make.names(var), "[1:nLevels.gamma.", make.names(var), "] ~ car.normal(adj.gamma.", make.names(var),
                            "[], weights.gamma.", make.names(var), "[], num.gamma.", make.names(var), "[], phi.", make.names(var), ")\n\t",
                "phi.", make.names(var), " ~ dgamma(a.phi.", make.names(var), ", b.phi.", make.names(var), ")", sep='')
        }
    }

    # rho spatial term (gamma)
    if (!is.null(rhoSpatialVar)) {
        cat(file=con, "\n\n\t",
            "gamma.", make.names(rhoSpatialVar), "[1:nSpatials.", make.names(rhoSpatialVar), "] ~ car.normal(adj.", make.names(rhoSpatialVar), "[], weights.",
                      make.names(rhoSpatialVar), "[], num.", make.names(rhoSpatialVar), "[], phi.", make.names(rhoSpatialVar), ")\n\t",
            "phi.", make.names(rhoSpatialVar), " ~ dgamma(a.phi.", make.names(rhoSpatialVar), ", b.phi.", make.names(rhoSpatialVar), ")", sep='')
    }


    # end
    cat(file=con, "\n}", sep='')

    close(con)
}

