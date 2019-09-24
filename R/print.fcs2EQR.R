#' Printing EQR Samples
#'
#' Prints a basic summary of a set of \acronym{EQR} samples.
#'
#'
#' @param x an object of class \code{"fcs2EQR"} containing Monte Carlo
#' \acronym{EQR} samples, as calculated from \code{\link{fcs2SingleEQR}},
#' \code{\link{fcs2JointEQR}} or \code{\link{fcs2JointAndSingleEQR}}.
#' @return invisibly returns the \code{"fcs2EQR"} object \code{x}.
#' @seealso \code{\link{summary.fcs2EQR}} for a more detailed summary;\cr
#' \code{\link{fcs2SingleEQR}}, \code{\link{fcs2JointEQR}} or
#' \code{\link{fcs2JointAndSingleEQR}} for generating \acronym{EQR} samples.
#' @keywords print
#' @export
print.fcs2EQR <-
function(x)
{
    cat("\n")
    cat("EQR samples:\n")

    # add number of samples
    cat(nrow(x), "Monte Carlo samples\n")

    # add survey or joinByVar name and number
    colVar <- names(dimnames(x))[2]
    if (is.null(colVar))
        colVar <- "surveys"
    else if (substr(colVar, nchar(colVar), nchar(colVar)) != "s")
        colVar <- paste(colVar, "s", sep="")
    cat(ncol(x), colVar)
    if (!is.null(na.action(x)))
        cat(" (", length(na.action(x)), " removed as missing)\n", sep="")
    else
        cat("\n")

    # if present, add species names
    if (length(dim(x)) == 3) {
        cat("Species: ", paste(dimnames(x)[[3]], collapse=", "), "\n", sep="")
    }

    cat("\n")

    invisible(x)
}

