#' Interactive FCS2 model likelihood
#'
#' Creates a plot of the \acronym{FCS2} model likelihood with interactive
#' controls that allow the model parameters to be adjusted.  The likelihood
#' represents the probability of obtaining a single fish count and the controls
#' allow you to see how these probabilities change with the parameters.  The
#' function can visualise either the zero-inflated negative binomial
#' (\acronym{ZINB}) distribution (see \code{\link{zinbinom}}) of the single catch
#' \eqn{C} in the original \acronym{EA} \acronym{FCS2} model or the
#' \acronym{ZINB} distribution of the total catch \eqn{T} over all passes from
#' the multiple-pass \acronym{FCS2} model.
#'
#'
#' @param multipass If \code{FALSE} (the default), the \acronym{ZINB}
#' likelihood of the single catch \eqn{C} from the original \acronym{EA}
#' \acronym{FCS2} model is given.  If \code{TRUE}, the \acronym{ZINB}
#' distribution of the total catch \eqn{T} over all passes from the
#' multiple-pass \acronym{FCS2} model is given.  The true likelihood for this
#' model is the multivariate \acronym{ZINM} distribution (see
#' \code{fcs2:::zinmultinom}).
#' @param eqr whether to display the \acronym{EQR} value.
#' @param boundaries a vector of length 4 giving the \acronym{EQR} boundaries
#' separating the classes \emph{Bad}, \emph{Poor}, \emph{Good}, \emph{Moderate}
#' and \emph{High}.  These are used only to colour the plot with \emph{Bad} red
#' and \emph{High} blue.  If \code{NULL} (default), the probability that
#' defines the single \acronym{EQR} is coloured blue.
#' @param r.limit,mu.limit,rho.limit,a.limit vectors of giving the lower and
#' upper limits for the interactive controls for the shape parameter \eqn{r},
#' the abundance \eqn{\mu}, the prevalence \eqn{rho} and the survey area
#' \eqn{a} respectively.
#' @param d.limit,q.limit vectors of giving the lower and upper limits for the
#' interactive controls for the number of runs \eqn{d} and the catch
#' probability \eqn{q} respectively.  These are used only if \code{multipass =
#' TRUE}.
#' @param c.limit vector of giving the lower and upper limit for the
#' interactive control for the single catch.  This is used only if \code{eqr =
#' TRUE}.
#' @section Warning: This function requires the additional package \pkg{rpanel}
#' for producing interactive controls.
#' @seealso \code{\link{fcs2InteractivePrediction}} which plot predictions from
#' a fitted \acronym{FCS2} model with interactive controls to show how the
#' predictions change when covariates are varied.
#' @keywords hplot
#' @export
#' @examples
#'
#' ## Interactive likelihood of total catch with EQR given against class boundaries
#'
#' # only run example if R is running interactively
#' if (interactive()) {
#'     # define WFD class boundaries
#'     boundaries <- c(0.001, 0.01, 0.25, 0.5)
#'
#'     # create interactive plot
#'     fcs2InteractiveLikelihood(multipass=TRUE, eqr=TRUE, boundaries=boundaries)
#' }
#'
fcs2InteractiveLikelihood <- function(multipass = FALSE, eqr = FALSE, boundaries = NULL, r.limit = c(0.1, 10), mu.limit = c(0, 10), rho.limit = c(0, 1),
                                      a.limit = c(0, 10), d.limit=c(1, 7), q.limit=c(0, 1), c.limit=c(0, 20))
{
    # load package 'rpanel'
    if (!requireNamespace("rpanel", quietly = TRUE))
        stop("`fcs2InteractiveLikelihood' requires R package `rpanel' - please install using `install.packages'")

    # calculate initial values as centre of limits
    r <- mean(r.limit)
    mu <- mean(mu.limit)
    rho <- mean(rho.limit)
    a <- mean(a.limit)
    catch <- round(mean(c.limit))

    if (!multipass) {
        # create panel with values as data
        if (eqr)
            panel <- rpanel::rp.control(sizer=r, mu=mu, rho=rho, area=a, catch=catch, title="EA FCS2 likelihood")
        else
            panel <- rpanel::rp.control(sizer=r, mu=mu, rho=rho, area=a, title="EA FCS2 likelihood")

        # create function to draw the ZINB PMF of the catch C
        draw <- function(panel) {
            nbmean <- panel$mu * panel$area
            zeroprob <- 1 - panel$rho
            xmax <- qzinbinom(0.95, panel$sizer, nbmean=nbmean, zeroprob=zeroprob)
            xplt <- 0:xmax
            prob <- dzinbinom(xplt, panel$sizer, nbmean=nbmean, zeroprob=zeroprob)

            # calculate colour
            if (is.null(boundaries)) {
                # if boundaries = NULL, colour blue and grey
                col <- grDevices::hsv(h=0.66, s=(xplt <= panel$catch) - 0.1 * (xplt < panel$catch),
                           v=0.8 + 0.2 * (xplt <= panel$catch) - 0.1 * (xplt < panel$catch))

            } else {
                # colour by WFD classification
                cumProb <- cumsum(prob)
                hue <- rep(0, length(cumProb))
                hue[cumProb > boundaries[1] & cumProb <= boundaries[2]] <- 1/10
                hue[cumProb > boundaries[2] & cumProb <= boundaries[3]] <- 2/10
                hue[cumProb > boundaries[3] & cumProb <= boundaries[4]] <- 3/10
                hue[cumProb > boundaries[4]] <- 5/10
                col <- grDevices::hsv(h=hue, s=0.3 + 0.4 * (xplt <= panel$catch) - 0.2 * (xplt < panel$catch),
                           v=0.9 + 0.1 * (xplt <= panel$catch) - 0.05 * (xplt < panel$catch))
            }

            # set title
            title <- paste("C ~ ZINB(r=", signif(panel$sizer, 3), ", m=", signif(nbmean, 3), ", z=", signif(zeroprob, 3),
                           "):\nmean=", signif(nbmean*panel$rho, 3), sep="")
            if (eqr)
                title <- paste(title, ", EQR=", signif(pzinbinom(panel$catch, panel$sizer, nbmean=nbmean, zeroprob=zeroprob), 3), sep="")

            # plot
            par.old <- graphics::par(mfrow=c(1,1), ask=FALSE)
            graphics::barplot(names=xplt, prob, xlab="Catch c", ylab="Probability",
                    main=title, col=col, space=0)
            graphics::par(par.old)

            panel
        }

        # add sliders to control parameters
        rpanel::rp.slider(panel, sizer, r.limit[1], r.limit[2], draw, "Shape r", showvalue=TRUE, resolution=signif(diff(r.limit)/100, 3))
        rpanel::rp.slider(panel, rho, rho.limit[1], rho.limit[2], draw, "Prevalence rho (probability present)", showvalue=TRUE, resolution=signif(diff(rho.limit)/100, 3))
        rpanel::rp.slider(panel, mu, mu.limit[1], mu.limit[2], draw, "Abundance mu (mean density)", showvalue=TRUE, resolution=signif(diff(mu.limit)/100, 3))
        rpanel::rp.slider(panel, area, a.limit[1], a.limit[2], draw, "Survey area a", showvalue=TRUE, resolution=signif(diff(a.limit)/100, 3))
        if (eqr)
            rpanel::rp.slider(panel, catch, c.limit[1], c.limit[2], draw, "Total catch t", showvalue=TRUE, resolution=1)

    } else {
        # calculate further initial values
        d.limit <- round(d.limit)
        d <- round(mean(d.limit))
        q <- mean(q.limit)

        # create panel with values as data
        if (eqr)
            panel <- rpanel::rp.control(sizer=r, mu=mu, rho=rho, area=a, d=d, q=q, catch=catch, title="Multiple-pass FCS2 likelihood", size=c(300, 400))
        else
            panel <- rpanel::rp.control(sizer=r, mu=mu, rho=rho, area=a, d=d, q=q, title="Multiple-pass FCS2 likelihood", size=c(300, 400))

        # create function to draw the ZINB PMF of the total catch T
        draw <- function(panel) {
            nbmean <- panel$mu * panel$area * (1 - ((1 - panel$q) ^ panel$d))
            zeroprob <- 1 - panel$rho
            xmax <- qzinbinom(0.95, panel$sizer, nbmean=nbmean, zeroprob=zeroprob)
            xplt <- 0:xmax
            prob <- dzinbinom(xplt, panel$sizer, nbmean=nbmean, zeroprob=zeroprob)

            # calculate colour
            if (is.null(boundaries)) {
                # if boundaries = NULL, colour blue and grey
                col <- grDevices::hsv(h=0.66, s=(xplt <= panel$catch) - 0.1 * (xplt < panel$catch),
                           v=0.8 + 0.2 * (xplt <= panel$catch) - 0.1 * (xplt < panel$catch))

            } else {
                # colour by WFD classification
                cumProb <- cumsum(prob)
                hue <- rep(0, length(cumProb))
                hue[cumProb > boundaries[1] & cumProb <= boundaries[2]] <- 1/10
                hue[cumProb > boundaries[2] & cumProb <= boundaries[3]] <- 2/10
                hue[cumProb > boundaries[3] & cumProb <= boundaries[4]] <- 3/10
                hue[cumProb > boundaries[4]] <- 5/10
                col <- grDevices::hsv(h=hue, s=0.3 + 0.4 * (xplt <= panel$catch) - 0.2 * (xplt < panel$catch),
                           v=0.9 + 0.1 * (xplt <= panel$catch) - 0.05 * (xplt < panel$catch))
            }

            # set title
            title <- paste("T ~ ZINB(r=", signif(panel$sizer, 3), ", m=", signif(nbmean, 3), ", z=", signif(zeroprob, 3),
                           "):\nmean=", signif(nbmean*panel$rho, 3), sep="")
            if (eqr)
                title <- paste(title, ", EQR=", signif(pzinbinom(panel$catch, panel$sizer, nbmean=nbmean, zeroprob=zeroprob), 3), sep="")

            # plot
            par.old <- graphics::par(mfrow=c(1,1), ask=FALSE)
            graphics::barplot(names=xplt, prob, xlab="Total catch t", ylab="Probability",
                    main=title, col=col, space=0)
            graphics::par(par.old)

            panel
        }

        # add sliders to control parameters
        rpanel::rp.slider(panel, sizer, r.limit[1], r.limit[2], draw, "Shape r", showvalue=TRUE, resolution=signif(diff(r.limit)/100, 3))
        rpanel::rp.slider(panel, rho, rho.limit[1], rho.limit[2], draw, "Prevalence rho (probability present)", showvalue=TRUE, resolution=signif(diff(rho.limit)/100, 3))
        rpanel::rp.slider(panel, mu, mu.limit[1], mu.limit[2], draw, "Abundance mu (mean density)", showvalue=TRUE, resolution=signif(diff(mu.limit)/100, 3))
        rpanel::rp.slider(panel, area, a.limit[1], a.limit[2], draw, "Survey area a", showvalue=TRUE, resolution=signif(diff(a.limit)/100, 3))
        rpanel::rp.slider(panel, d, d.limit[1], d.limit[2], draw, "Number of runs d", showvalue=TRUE, resolution=ceiling(diff(d.limit)/100))
        rpanel::rp.slider(panel, q, q.limit[1], q.limit[2], draw, "Catch probability q", showvalue=TRUE, resolution=signif(diff(q.limit)/100, 3))
        if (eqr)
            rpanel::rp.slider(panel, catch, c.limit[1], c.limit[2], draw, "Total catch t", showvalue=TRUE, resolution=1)
    }
}





