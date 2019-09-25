#' Printing an FCS2 Model Fit
#'
#' Prints a basic summary of an \acronym{FCS2} model fit.
#'
#'
#' @param x an object of class \code{"fcs2Fit"}, usually returned by
#' \code{\link{fcs2FitModel}}.
#' @param inla whether to summarise approximate model fits obtained from
#' \acronym{INLA}, if available.
#' @param \dots Not currently used.
#' @return invisibly returns the \code{"fcs2Fit"} object \code{x}.
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{summary.fcs2Fit}} for a more detailed summary.
#' \code{\link{fcs2FitModel}} for fitting the \acronym{FCS2} model.
#' @keywords print
#' @export
print.fcs2Fit <-
function(x, inla = is.null(x$bugsFit), ...)
{
    cat("\n")

    # data
    cat("Data:\n")
    if (sum(x$dataType == "run", na.rm=TRUE) > 0)
        cat("Run totals:    ", paste(x$runTotalVars, collapse=", "), "\n")
    if (sum(x$dataType == "total", na.rm=TRUE) > 0)
        cat("All runs total:", x$allRunsTotalVar, "\n")
    if (sum(x$dataType == "range", na.rm=TRUE) > 0)
        cat("All runs range:", paste(x$allRunsRangeVars, collapse=", "), "\n")
    cat("Survey area:   ", x$surveyAreaVar, "\n\n")

    # model formulae (mu and rho)
    cat("Abundance formula:\n")
    print(formula(paste("log(mu)", paste(deparse(x$muFormula), collapse=" "))), showEnv=FALSE)
    cat("\n")

    # model formulae (mu and rho)
    cat("Prevalence formula:\n")
    print(formula(paste("logit(rho)", paste(deparse(x$rhoFormula), collapse=" "))), showEnv=FALSE)
    cat("\n")

    # posterior parameter estimates from INLA
    if (inla && !is.null(x$inlaFits)) {

        if (!is.null(x$inlaFits$muFit)) {
            cat("Approximate posterior means from INLA x to abundance:\n")
            means <- c(x$inlaFits$muFit$summary.hyperpar[1, 1], x$inlaFits$muFit$summary.fixed[, 1])

            # transform hyperparameters
            if (nrow(x$inlaFits$muFit$summary.hyperpar) > 1) {
                t <- function(x)
                    sqrt(1 / x)
                for (i in 1:(nrow(x$inlaFits$muFit$summary.hyperpar) - 1))
                    means <- c(means, inla.expectation(t, x$inlaFits$muFit$marginals.hyperpar[[i + 1]]))
            }
            names(means) <- variable.names(x, q=FALSE, r=TRUE, prevalence=FALSE, rw=FALSE, spatial=FALSE, hyperparams="scale")

            # print
            print(means)
            cat("\n")
        }

        cat("Approximate posterior means from INLA x to prevalence:\n")
        means <- x$inlaFits$rhoFit$summary.fixed[, 1]

        # transform hyperparameters
        if (length(x$inlaFits$rhoFit$marginals.hyperpar) > 0) {
            t <- function(x)
                sqrt(1 / x)
            for (i in 1:nrow(x$inlaFits$rhoFit$summary.hyperpar))
                means <- c(means, inla.expectation(t, x$inlaFits$rhoFit$marginals.hyperpar[[i]]))
        }
        names(means) <- variable.names(x, q=FALSE, r=FALSE, abundance=FALSE, rw=FALSE, spatial=FALSE, hyperparams="scale")

        # print
        print(means)
        cat("\n")
    }

    # posterior parameter estimates (singular variables only) from BUGS
    if (!is.null(x$bugsFit)) {
        cat("Posterior means from ", x$bugsFit$program, ":\n", sep="")
        params <- variable.names(x, rw=FALSE, spatial=FALSE, hyperparams="precision")
        means <- x$bugsFit$summary[params, "mean"]

        # replace hyperparameter precision means with scale means
        hyperVars <- variable.names(x, r=FALSE, q=FALSE, linear=FALSE, rw=FALSE, spatial=FALSE, hyperparams="precision")
        if (length(hyperVars) > 0) {
            for (var in hyperVars)
                means[var] <- mean(1 / sqrt(x$bugsFit$sims.list[[var]]))
            names(means) <- sub("tau", "sigma", names(means), fixed=TRUE)
            names(means) <- sub("phi", "nu", names(means), fixed=TRUE)
        }

        # print
        print(means)
        cat("\n")
    }

    invisible(x)
}

