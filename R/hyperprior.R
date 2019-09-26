#' The Hyperprior Distribution
#'
#' @encoding UTF-8
#' @title The Hyperprior Distribution
#'
#' @name hyperprior
#' @description
#' Density, distribution function, quantile function and random generation for
#' the distribution applied as a prior for the scale hyperparameters in the
#' \acronym{FCS2} model.\cr This is the distribution of \eqn{\sigma = 1 /
#' \sqrt{\tau}} where \eqn{\tau} follows a Gamma distribution.
#'
#' @aliases dhyperprior phyperprior qhyperprior rhyperprior
#' @usage dhyperprior(x, a, b)
#' @usage phyperprior(q, a, b, lower.tail = TRUE)
#' @usage qhyperprior(p, a, b)
#' @usage rhyperprior(n, a, b)
#' @param x vector of quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of samples required.
#' @param a a shape parameter, equal to the shape parameter of the related Gamma
#'   distribution.
#' @param b a scale parameter, equal to the rate parameter of the related Gamma
#'   distribution.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P\{X \le x\}}{P{X \le x}}, otherwise, \eqn{P\{X > x\}}{P{X > x}}
#' @details The hyperprior distribution with shape parameter \eqn{a} and scale
#'   parameter \eqn{b} has density \deqn{p(x) = 2 \frac{b^a}{\Gamma(a)}
#'   \frac{1}{x^{2a + 1}} \exp\left( -\frac{b}{x^2} \right) }{% p(x) = 2 b^a /
#'   (\Gamma(a) x^(2a + 1)) exp( -b / (x^2) ) } for \eqn{x, a, b > 0}.
#'
#'   This is the distribution of \eqn{\sigma = 1 / \sqrt{\tau}}{\sigma = 1 /
#'   \sqrt(\tau)} where \eqn{\tau} follows a Gamma distribution with shape
#'   parameter \eqn{a} and rate parameter \eqn{b}.
#'
#'   The mean is given by \eqn{\Gamma(a - 0.5) \sqrt{b} / \Gamma(a)}{\Gamma(a -
#'   0.5) \sqrt(b) / \Gamma(a)} if \eqn{a > 0.5} and is infinite otherwise.\cr
#'   The variance is given by \eqn{b / (a - 1) - b \Gamma(a - 0.5)^2 /
#'   \Gamma(a)^2}{b / (a - 1) - b \Gamma(a - 0.5)^2 / (\Gamma(a)^2)} if \eqn{a >
#'   1} and is infinite otherwise.
#' @return \code{dhyperprior} gives the density, \code{phyperprior} gives the
#'   distribution function, \code{qhyperprior} gives the quantile function, and
#'   \code{rhyperprior} generates random deviates.
#' @seealso \code{\link{dgamma}} for the Gamma distribution.\cr
#'   \code{\link{specialTerms}} for the non-linear \acronym{FCS2} model terms
#'   that use this distribution as a prior for their scale parameters.
#' @keywords distribution
#' @export
dhyperprior <-
function(x, a, b)  #, log = FALSE)
{
    2 * (b ^ a) * exp(-b / (x^2)) * (x > 0) / (gamma(a) * (x ^ (2 * a + 1)) * (x > 0) + (x <= 0))
}

#' @rdname hyperprior
#' @export
phyperprior <-
function(q, a, b, lower.tail = TRUE) #, log.p = FALSE)
{
    pgamma(1 / (q^2), a, b, lower.tail=!lower.tail) * (q >= 0) + (q < 0) * (!lower.tail)
}

#' @rdname hyperprior
#' @export
qhyperprior <-
function(p, a, b)  #, lower.tail = TRUE, log.p = FALSE)
{
    1 / sqrt(qgamma(1 - p, a, b))
}

#' @rdname hyperprior
#' @export
rhyperprior <-
function(n, a, b)
{
    1 / sqrt(rgamma(n, a, b))
}


