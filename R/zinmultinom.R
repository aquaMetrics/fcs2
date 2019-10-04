#' Zero-inflated Negative Multinomial Distributions.
#'
#' @encoding UTF-8
#' @title Zero-inflated Negative Multinomial Distributions.
#'
#' @name zinmultinom
#' @aliases dzinmultinom rzinmultinom
#' @usage dzinmultinom(x, size, prob, zeroprob, nmmean)
#' @usage rzinmultinom(n, size, prob, zeroprob, nmmean)
#' @param x x
#' @param n n
#' @param size size
#' @param prob prob
#' @param zeroprob  zero prob
#' @param nmmean nmmean
#'
#' @return zero-inflated Negative Multinomial Distributions.
#' @export
dzinmultinom <-
function(x, size, prob, zeroprob, nmmean)  #, log = FALSE)
{
    # make sure x is a matrix
    if (is.null(dim(x)))
        dim(x) <- c(1, length(x))

    # find number of entries and dimensions
    n <- nrow(x)
    d <- ncol(x)

    # find prob from negative-multinomial mean
    if (!missing(nmmean)) {
        if (!missing(prob))
            stop("'prob' and 'nmmean' both specified")

        # make sure nmmean is a matrix
        if (is.null(dim(nmmean)))
            dim(nmmean) <- c(1, length(nmmean))

        prob <- nmmean / (size + rowSums(nmmean))
    }

    # make sure prob is a matrix with n rows
    if (is.null(dim(prob)))
        prob <- matrix(prob, byrow=TRUE, nrow=n, ncol=length(prob))

    # if only d probs provided, calculate p0
    if (ncol(prob) == d)
        prob <- cbind(1 - apply(prob, 1, sum), prob)

    # calculate sum of x
    sumx <- apply(x, 1, sum)

    # return density
    zeroprob * (sumx == 0) + (1 - zeroprob) * gamma(size + sumx) * (prob[, 1] ^ size) * apply(prob[, -1, drop=FALSE] ^ x, 1, prod) /
            (gamma(size) * apply(factorial(x), 1, prod))
}


## rzinmultinom
#' @rdname zinmultinom
#' @export
rzinmultinom <-
function(n, size, prob, zeroprob, nmmean)
{
    # sample whether present
    present <- rbinom(n, 1, 1 - zeroprob)

    # find prob from negative-multinomial mean
    if (!missing(nmmean)) {
        if (!missing(prob))
            stop("'prob' and 'nmmean' both specified")

        # make sure nmmean is a matrix
        if (is.null(dim(nmmean)))
            dim(nmmean) <- c(1, length(nmmean))

        prob <- nmmean / (size + rowSums(nmmean))
    }

    # make sure prob is a matrix with n rows
    if (is.null(dim(prob)))
        prob <- matrix(prob, byrow=TRUE, nrow=n, ncol=length(prob))

    # find number of dimensions
    d <- ncol(prob)

    # calculate p0 (assumes nrow(prob) is number of dimensions)
    prob <- cbind(1 - rowSums(prob), prob)

    # sample negative-multinomial
    ret <- array(dim=c(n, d))
    ret[, 1] <- rnbinom(n, size, prob[, 1] / (prob[, 1] + prob[, 2]))
    if (d > 1) {
        for (i in 2:d)
            ret[, i] <- rnbinom(n, size + rowSums(ret[, 1:(i-1), drop=FALSE]),
                                rowSums(prob[, 1:i, drop=FALSE]) / rowSums(prob[, 1:(i+1), drop=FALSE]))
    }

    # multiply by whether present
    present * ret
}
