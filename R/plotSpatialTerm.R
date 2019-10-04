#' Plot Spatial Term
#'
#' Plots summaries of the spatial term estimates from a fitted \acronym{FCS2}
#' model on a map of the spatial regions.
#'
#'
#' @param map a \code{SpatialPolygonsDataFrame} object containing
#' polygons and a data frame corresponding to the regions used for a spatial
#' term in the \acronym{FCS2} model.  This will usually be generated from an
#' ESRI shapefile using \code{readShapePoly}.
#' @param fit an \code{"fcs2Fit"} object containing \acronym{INLA} or
#' \acronym{BUGS} model fits, as returned from \code{\link{fcs2FitModel}}.  The
#' model should include a \code{\link{spatial}} term in the abundance and/or
#' prevalence component which is based upon the region given by \code{map}.
#' @param type which summary to plot. Either \code{"mean"} for the mean or
#' \code{"sd"} for the standard deviation, or both.
#' @param abundance whether the spatial term (with these regions) is in the
#' abundance regression component and should be plotted.
#' @param prevalence whether the spatial term (with these regions) is in the
#' prevalence regression component and should be plotted.
#' @param posterior whether to plot the posterior estimates of the spatial term
#' (if available).
#' @param inla whether to plot the approximate \acronym{INLA} estimates of the
#' spatial term.
#' @param ...  further arguments passed to \code{spplot} and perhaps on
#' further to \code{levelplot}.  For example, the colour pallete can be
#' changed using \code{col.regions} which by default is set using
#' \code{\link{cm.colors}}.
#' @section Warning: This function requires the additional package \pkg{sp} to
#' provide the \code{SpatialPolygonsDataFrame} map.  This package
#' is one of the requirements of the package \pkg{maptools} which may also be
#' useful for creating \code{map} from an ESRI shapefile using the function
#' \code{readShapePoly}.
#' @seealso \code{\link{fcs2SpatialSummary}} which this function uses to attach
#' the spatial summaries to the map before plotting.
#' @keywords hplot
#' @export
plotSpatialTerm <-
function(map, fit, type = "mean", abundance = TRUE, prevalence = TRUE, posterior = !is.null(fit$bugsFit), inla = !posterior, ...)
{
    # load package 'maptools'
    if (!requireNamespace("maptools", quietly = TRUE))
        stop("`plotSpatialTerm' requires R package `maptools' - please install using `install.packages'")

    ## Add spatial terms to map
    map <- fcs2SpatialSummary(map, fit, posterior, inla, abundance, prevalence)

    ## Plot terms
    variables <- variable.names(fit, r=FALSE, q=FALSE, linear=FALSE, rw=FALSE, hyper=FALSE, abundance=abundance, prevalence=prevalence)
    which <- c()
    if (inla)
        which <- paste("inla", type, variables, sep=".")
    if (posterior)
        which <- c(which, paste(type, variables, sep="."))
    sp::spplot(map, which, ...)
}

