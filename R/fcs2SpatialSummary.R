#' Summary of Spatial Terms in FCS2 Model
#'
#' Summarises the spatial terms in a fitted \acronym{FCS2} model and attaches
#' the summary to a map of the spatial region used.
#'
#'
#' @param map a \code{sp::SpatialPolygonsDataFrame} object containing
#' polygons and a data frame corresponding to the regions used for a spatial
#' term in the \acronym{FCS2} model.  This will usually be generated from an
#' ESRI shapefile using \code{readShapePoly}.
#' @param fit an \code{"fcs2Fit"} object containing \acronym{INLA} or
#' \acronym{BUGS} model fits, as returned from \code{\link{fcs2FitModel}}.  The
#' model should include a \code{\link{spatial}} term in the abundance and/or
#' prevalence component which is based upon the region given by \code{map}.
#' @param posterior whether to summarise the posterior estimates of each region
#' (if available).
#' @param inla whether to summarise the approximate \acronym{INLA} estimates of
#' each region.
#' @param abundance whether the spatial term (with these regions) is in the
#' abundance regression component.
#' @param prevalence whether the spatial term (with these regions) is in the
#' prevalence regression component.
#' @return the \code{sp::SpatialPolygonsDataFrame} object \code{map}
#' with additional columns in the data frame component \code{map@data}
#' containing the requested summaries.  Each summary is given as two columns,
#' one for the mean of each variable named \code{mean.}* and another for the
#' standard deviation named \code{sd.}*, where * is the name of the spatial
#' variable.
#' @section Warning: This function requires the additional package
#' \pkg{maptools} which may also be useful for creating \code{map} from an ESRI
#' shapefile using the function \code{maptools::readShapePoly}.
#' @seealso \code{\link{plotSpatialTerm}} which uses this function to plot the
#' spatial term summaries on a map.
#' @export
fcs2SpatialSummary <-
function(map, fit, posterior = !is.null(fit$bugsFit), inla = !posterior, abundance = TRUE, prevalence = TRUE)
{
    # load package 'sp'
    if (!requireNamespace("sp", quietly = TRUE))
        stop("`fcs2SpatialSummary' requires R package `sp' - please install using `install.packages'")

    # INLA
    if (inla && !is.null(fit$inlaFit)) {
        if (!is.null(fit$muSpatialVar) && abundance) {
            var <- paste("beta", fit$muSpatialVar, sep='.')
            mean <- fit$inlaFits$muFit$summary.random[[make.names(fit$muSpatialVar)]]$mean
            sd <- fit$inlaFits$muFit$summary.random[[make.names(fit$muSpatialVar)]]$sd
            map@data <- data.frame(map@data, mean, sd)
            colnames(map@data)[ncol(map@data) - 1:0] <- paste("inla", c("mean", "sd"), var, sep='.')
        }

        if (!is.null(fit$rhoSpatialVar) && prevalence) {
            var <- paste("gamma", fit$rhoSpatialVar, sep='.')
            mean <- fit$inlaFits$rhoFit$summary.random[[make.names(fit$rhoSpatialVar)]]$mean
            sd <- fit$inlaFits$rhoFit$summary.random[[make.names(fit$rhoSpatialVar)]]$sd
            map@data <- data.frame(map@data, mean, sd)
            colnames(map@data)[ncol(map@data) - 1:0] <- paste("inla", c("mean", "sd"), var, sep='.')
        }
    }

    # BUGS
    if (posterior && !is.null(fit$bugsFit)) {
        if (!is.null(fit$muSpatialVar) && abundance) {
            var <- paste("beta", fit$muSpatialVar, sep='.')
            mean <- apply(fit$bugsFit$sims.list[[var]], 2, mean)
            sd <- apply(fit$bugsFit$sims.list[[var]], 2, sd)
            map@data <- data.frame(map@data, mean, sd)
            colnames(map@data)[ncol(map@data) - 1:0] <- paste(c("mean", "sd"), var, sep='.')
        }

        if (!is.null(fit$rhoSpatialVar) && prevalence) {
            var <- paste("gamma", fit$rhoSpatialVar, sep='.')
            mean <- apply(fit$bugsFit$sims.list[[var]], 2, mean)
            sd <- apply(fit$bugsFit$sims.list[[var]], 2, sd)
            map@data <- data.frame(map@data, mean, sd)
            colnames(map@data)[ncol(map@data) - 1:0] <- paste(c("mean", "sd"), var, sep='.')
        }
    }

    map
}

