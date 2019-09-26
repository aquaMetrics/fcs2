#' Zero-inflated Negative Binomial Distributions.
#'
#' @encoding UTF-8
#' @title Zero-inflated Negative Binomial Distributions.
#'
#' @name zinbinom
#' @aliases dzimbinom pzinbinom qzinbinom rzinbinom rzinbinom_constrained
#' @usage dzinbinom(x, size, prob, zeroprob, nbmean)
#' @usage pzinbinom(q, size, prob, zeroprob, nbmean)
#' @usage qzinbinom(p, size, prob, zeroprob, nbmean)
#' @usage rzinbinom(n, size, prob, zeroprob, nbmean)
#' @usage rzinbinom_constrained(n, size, prob, zeroprob, nbmean, min = 0, max = Inf)
#' @param x x
#' @param q q
#' @param p p
#' @param n n
#' @param size size
#' @param prob prob
#' @param zeroprob  zero prob
#' @param nbmean nbmean
#' @param min min
#' @param max max
#'
#' @return zero-inflated Negative Binomial Distributions.
#' @export
#'
#' @examples \dontrun{
#'  data$Catch <- rzinbinom(100,
#' size = 1.1, zeroprob = 0.3,
#' nbmean = 0.3 * Data$SurveyArea
#' )
#' }
dzinbinom <-
function(x, size, prob, zeroprob, nbmean)  #, log = FALSE)
{
    if (!missing(nbmean)) {
        if (!missing(prob))
            stop("'prob' and 'nbmean' both specified")
        (1 - zeroprob) * dnbinom(x, size, mu=nbmean)  +  zeroprob * (x == 0)

    } else
        (1 - zeroprob) * dnbinom(x, size, prob)  +  zeroprob * (x == 0)
}


#' @rdname zinbinom
#' @export
pzinbinom <-
function(q, size, prob, zeroprob, nbmean)  #, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(nbmean)) {
        if (!missing(prob))
            stop("'prob' and 'nbmean' both specified")
        (1 - zeroprob) * pnbinom(q, size, mu=nbmean)  +  zeroprob * (q >= 0)

    } else
        (1 - zeroprob) * pnbinom(q, size, prob)  +  zeroprob * (q >= 0)
}


#' @rdname zinbinom
#' @export
qzinbinom <-
function(p, size, prob, zeroprob, nbmean)  #, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(nbmean) && !missing(prob))
        stop("'prob' and 'nbmean' both specified")
    qnbinom(pmax((p - zeroprob) / (1 - zeroprob), 0), size, prob, nbmean)
}


#' @rdname zinbinom
#' @export
rzinbinom <-
function(n, size, prob, zeroprob, nbmean)
{
    # sample whether present
    present <- rbinom(n, 1, 1 - zeroprob)

    if (!missing(nbmean)) {
        if (!missing(prob))
            stop("'prob' and 'nbmean' both specified")
        present * rnbinom(n, size, mu=nbmean)

    } else
        present * rnbinom(n, size, prob)
}


#' @rdname zinbinom
#' @export
rzinbinom_constrained <-
function(n, size, prob, zeroprob, nbmean, min = 0, max = Inf)
{
    qzinbinom(runif(n, pzinbinom(min - 1, size, prob, zeroprob, nbmean),
                       pzinbinom(max, size, prob, zeroprob, nbmean)),
              size, prob, zeroprob, nbmean)
}

