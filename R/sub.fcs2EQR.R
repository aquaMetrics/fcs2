#' @export
## Allow subset of eqr object to give another eqr object
`[.fcs2EQR` <- function(eqr, samples = 1:nrow(eqr), survey = 1:ncol(eqr), species, drop = TRUE)
{
    # convert to array
    eqr <- as.array(eqr)

    # return values if samples has dimensions
    if (!is.null(dim(samples))) {
        # if data frame, convert to array of indices
        if (inherits(samples, "data.frame")) {
            smpls <- array(dim=dim(samples))
            if (inherits(samples[,1], "character"))
                smpls[,1] <- match(samples[,1], rownames(eqr))
            else
                smpls[,1] <- samples[,1]
            if (inherits(samples[,2], "character"))
                smpls[,2] <- match(samples[,2], colnames(eqr))
            else
                smpls[,2] <- samples[,2]
            if (ncol(samples) > 2) {
                if (inherits(samples[,3], "character"))
                    smpls[,3] <- match(samples[,3], dimnames(eqr)[[3]])
                else
                    smpls[,3] <- samples[,3]
            }
            return(eqr[smpls])
        }

        return(eqr[samples])
    }

    # check whether multiple eqrs
    if (length(dim(eqr)) == 3) {
        # set species if necessary
        if (missing(species))
            species <- 1:(dim(eqr)[3])

        # extract subset
        eqr <- eqr[samples, survey, species, drop=FALSE]

        # if drop == TRUE and only one species, drop to matrix
        if (drop && dim(eqr)[3] == 1) {
            dnames <- dimnames(eqr)
            dim(eqr) <- dim(eqr)[1:2]
            if (!is.null(dnames))
                dimnames(eqr) <- dnames[1:2]
        }

    } else {
        # extract subset
        eqr <- eqr[samples, survey, drop=FALSE]
    }

    # set class
    class(eqr) <- "fcs2EQR"

    # return
    eqr
}

