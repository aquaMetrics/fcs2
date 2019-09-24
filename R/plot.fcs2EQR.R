#' Plot EQR Samples
#'
#' Plots a series of \acronym{EQR} variables or the corresponding probabilities
#' of each \acronym{WFD} class.
#'
#'
#' @param x an object of class \code{"fcs2EQR"} containing Monte Carlo
#' \acronym{EQR} samples, as calculated from \code{\link{fcs2SingleEQR}},
#' \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}.
#' @param survey index specifying which surveys (or water bodies etc if surveys
#' joined) to plot.
#' @param species index specifying which species to plot. This should
#' correspond to the plots defined by \code{survey}.
#' @param boundaries a vector of length 4 giving the \acronym{EQR} boundaries
#' separating the classes \emph{Bad}, \emph{Poor}, \emph{Good}, \emph{Moderate}
#' and \emph{High}.  These are used only to colour the plot.  If missing,
#' regularly spaced boundaries of \code{c(0.2, 0.4, 0.6, 0.8)} are used with a
#' warning.  If \code{NULL} and \code{type = "density"}, the plot is not
#' coloured.
#' @param type either \code{"density"} or \code{1} for a \code{\link{density}}
#' estimate of the \acronym{EQR} variable, or \code{"class"} or \code{2} to
#' plot the probability of each class (as produced by
#' \code{\link{fcs2Classify}}).
#' @param title an optional character vector giving the title for each plot.
#' Should have one entry per selection in \code{survey} and \code{species}.
#' @param adjust the bandwidth adjustment for the density estimate of the
#' \acronym{EQR}, if \code{type = "density"}.  See \code{\link{density}} for
#' further details.
#' @param \dots further arguments that are passed to either \code{\link{plot}}
#' (if \code{type = "density"}) or \code{\link{barplot}} otherwise.
#' @seealso \code{\link{fcs2SingleEQR}}, \code{\link{fcs2JointEQR}} or
#' \code{\link{fcs2JointAndSingleEQR}} for producing \code{"fcs2EQR"}
#' objects.\cr
#'
#' \code{\link{fcs2InteractivePrediction}} can be used to interactively view
#' how the \acronym{EQR} changes with covariates.
#' @keywords hplot
#' @export
plot.fcs2EQR <-
function(x, survey = 1, species = 1, boundaries, type = "density", title, adjust = 1, ...)
{
    # if boundaries missing, use regularly spaced boundaries
    if (missing(boundaries)) {
        boundaries <- seq(0.2, 0.8, by=0.2)
        warning("using regularly spaced class boundaries to colour the plot as 'boundaries' argument not supplied")
    }
    if (!is.null(boundaries))
        boundaries <- c(-0.001, boundaries, 1)

    # if single EQR, convert to multi-EQR
    if (length(dim(x)) == 2) {
        dn <- dimnames(x)
        dim(x) <- c(dim(x), 1)
        dimnames(x) <- c(dn, list(NULL))
    }

    # check survey and species not logical
    if (inherits(survey, "logical"))
        survey <- which(survey)
    if (inherits(species, "logical"))
        species <- which(species)

    # create data frame of survey and species to allow cyclic matching
    subset <- data.frame(Survey=survey, Species=species, stringsAsFactors=FALSE)

    # create titles
    if (missing(title)) {
        if (is.null(dimnames(x)[[3]])) {
            if (dim(x)[3] == 1)
                title <- paste("EQR at", colnames(x[, subset$Survey, , drop=FALSE]))
            else
                title <- paste("EQR", subset$Species, "\nat", colnames(x[, subset$Survey, , drop=FALSE]))

        } else
            title <- paste("EQR for", dimnames(x[, , subset$Species, drop=FALSE])[[3]], "\nat", colnames(x[, subset$Survey, , drop=FALSE]))
    }
    if (length(title) == 1)
        title <- rep(title, nrow(subset))

    # check eqrs available for subset, setting up subsubset
    subsubset <- (1:nrow(subset))[!is.na(x[data.frame(1, subset$Survey, subset$Species, stringsAsFactors=FALSE)])]
    if (length(subsubset) < nrow(subset))
        warning(paste(nrow(subset) - length(subsubset), "of the selected EQRs are not available"))

    # set up multiple plots
    if (length(subsubset) > 1) {
        w <- ceiling(sqrt(length(subsubset)))
        par.old <- par(mfrow=c(w - (length(subsubset) <= w*(w-1)), w))
        on.exit(par(par.old))
    }

    # check plot type
    if (type %in% c(1, "density")) {
        # density plot
        for (i in subsubset) {
            if (sd(x[, subset$Survey[i], subset$Species[i]]) == 0)
                d <- density(jitter(x[, subset$Survey[i], subset$Species[i]], factor=0.1), from=0, to=1, adjust=adjust)  # prevents bw being too large when all identical
            else
                d <- density(x[, subset$Survey[i], subset$Species[i]], from=0, to=1, adjust=adjust)
            plot(d, main=title[i], xlab="EQR", ylim=c(0, max(d$y, na.rm=TRUE)), ...)
            if (!is.null(boundaries)) {
                for (j in 1:(length(boundaries) - 1)) {
                    ix <- which(d$x >= boundaries[j] & d$x <= boundaries[j + 1])
                    if (j > 1)
                        ix <- c(ix[1] - 1, ix)  # remove gap between regions
                    polygon(d$x[c(ix, rev(ix))], c(d$y[ix], rep(0, length(ix))), border=NA, col=hsv(s=0.5, h=(j - 1 + (j == 5))/10))
                }
            }
            abline(v=c(0, 1), col="grey80")
            lines(d)
        }

    } else {
        # probability of each class
        for (i in subsubset) {
            prob <- numeric(length(boundaries) - 1)
            names(prob) <- c("Bad", "Poor", "Moderate", "Good", "High")
            for (j in 1:(length(boundaries) - 1))
                prob[j] <- mean(x[, subset$Survey[i], subset$Species[i]] > boundaries[j] & x[, subset$Survey[i], subset$Species[i]] <= boundaries[j + 1])
            barplot(prob, ylim=c(0, 1), col=hsv(s=0.5, h=c(0:3, 5)/10), ylab="Probability", main=title[i], ...)
        }
    }
}

