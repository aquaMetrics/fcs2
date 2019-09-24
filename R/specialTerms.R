#' Special Non-Linear Terms for FCS2 Abundance and Prevalence Regressions
#'
#' Non-linear random walk and spatial terms that can be used in addition to
#' linear terms in the \acronym{FCS2} abundance and prevalence regression
#' equations. These functions can be used within the formulae when specifying the
#' abundance or prevalence regression equations via the arguments
#' \code{muFormula} and \code{rhoFormula} of \code{\link{fcs2FitModel}}.
#'
#' @aliases specialTerms rw1 rw2 spatial
#'
#' @param val a covariate to use for the non-linear term.  \code{rw1} and
#'   \code{rw2} assume a continuous variable but \code{spatial} requires an
#'   integer variable that indicates which spatial region every survey is
#'   within.
#' @param noLevels the number of discrete levels to use to represent the
#'   continuous non-linear relationship. This is ignored if \code{boundaries} is
#'   supplied.
#' @param  boundaries an optional vector specifying the location of the discrete
#'   points that are used to represent the non-linear covariate term. If
#'   missing, \code{noLevels} boundaries are created spaced regularly between
#'   the minimum and maximum values of \code{val}.
#' @param adjacency a list with three components containing adjacency
#'   information relating to the spatial region. The first component \code{num}
#'   is a vector containing the number of regions adjacent to each region. The
#'   second component \code{adj} is a vector listing the indices of each of
#'   these adjacent regions. The third component \code{sumNumNeigh} is the sum
#'   of the number of neighbours, which should equal \code{sum(num)} and
#'   \code{length(adj)}.\cr This list will usually be generated externally from
#'   \R by the adjacency tool within WinBUGS or OpenBUGS. The list should be
#'   given a name and then loaded into \R.
#' @param scale.parameters an optional vector of length 2 specifying the prior
#'   distribution of the scale hyperparameter (\eqn{\sigma} or \eqn{\nu}) that
#'   controls the variability between levels in the random walk or adjacent
#'   spatial regions.
#' @param precision.parameters an optional vector of length 2 specifying the
#'   prior distribution of the precision hyperparameter (\eqn{\tau = 1 /
#'   \sigma^2} or \eqn{\phi = 1 / \nu^2}) that controls the variability between
#'   levels in the random walk or adjacent spatial regions. This is ignored if
#'   \code{scale.parameters} is specified.
#' @param  initial.values an optional vector giving an estimate for each of the
#'   variables that make up the term, one for each random walk level or one for
#'   each spatial region. These are used as starting values for the
#'   \acronym{MCMC} chains when fitting the full model using \acronym{BUGS}.
#' @param precision.initial.value an optional starting value for the precision
#'   hyperparameter to be used when fitting the full model using \acronym{BUGS}.
#' @details These functions can be used within the \R formulae when specifying
#'   the abundance or prevalence regression equations via the arguments
#'   \code{muFormula} and \code{rhoFormula} of \code{\link{fcs2FitModel}}.
#' @seealso \code{\link{fcs2FitModel}}
#' @export

rw1 <-
function(val, noLevels = 10, boundaries, scale.parameters, precision.parameters, initial.values, precision.initial.value)
{
    ret <- list(val=val,
                name=deparse(substitute(val)),
                type="rw1")

    if (missing(boundaries))
        ret <- c(ret, noLevels=noLevels)
    else
        ret <- c(ret, list(boundaries=boundaries))

    if (!missing(scale.parameters))
        ret <- c(ret, list(scale.parameters=scale.parameters))
    else if (!missing(precision.parameters))
        ret <- c(ret, list(precision.parameters=precision.parameters))

    if (!missing(initial.values))
        ret <- c(ret, list(initial.values=initial.values))
    if (!missing(precision.initial.value))
        ret <- c(ret, precision.initial.value=precision.initial.value)

    ret
}


## rw2
rw2 <-
function(val, noLevels = 10, boundaries, scale.parameters, precision.parameters, initial.values, precision.initial.value)
{
    ret <- list(val=val,
                name=deparse(substitute(val)),
                type="rw2")

    if (missing(boundaries))
        ret <- c(ret, noLevels=noLevels)
    else
        ret <- c(ret, list(boundaries=boundaries))

    if (!missing(scale.parameters))
        ret <- c(ret, list(scale.parameters=scale.parameters))
    else if (!missing(precision.parameters))
        ret <- c(ret, list(precision.parameters=precision.parameters))

    if (!missing(initial.values))
        ret <- c(ret, list(initial.values=initial.values))
    if (!missing(precision.initial.value))
        ret <- c(ret, precision.initial.value=precision.initial.value)

    ret
}


## spatial
spatial <-
function(val, adjacency, scale.parameters, precision.parameters, initial.values, precision.initial.value)
{
    # check adjacency has required entries
    if (missing(adjacency))
        stop("adjacency information must be given for each spatial term")
    if (class(adjacency) != "list" || sum(c("num", "adj", "sumNumNeigh") %in% names(adjacency)) != 3)
        stop("'adjacency' must be a list with components 'num', 'adj' and 'sumNumNeigh'")

    ret <- list(val=as.numeric(val),
                name=deparse(substitute(val)),
                type="spatial",
                adjacency=adjacency)

    if (!missing(scale.parameters))
        ret <- c(ret, list(scale.parameters=scale.parameters))
    else if (!missing(precision.parameters))
        ret <- c(ret, list(precision.parameters=precision.parameters))

    if (!missing(initial.values))
        ret <- c(ret, list(initial.values=initial.values))
    if (!missing(precision.initial.value))
        ret <- c(ret, precision.initial.value=precision.initial.value)

    ret
}
