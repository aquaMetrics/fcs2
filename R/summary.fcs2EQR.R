#' Detailed EQR summary
#'
#' Provides a more detailed EQR summary
#'
#' @param object summary object
#' @param \dots Not currently used.
#'
#' @return detailed EQR summary object
#' @export
#'
#' @examples
#' \dontrun{
#' summary.fcs2EQR(object = data)
#' }
#' @export
summary.fcs2EQR <-
function(object, ...)
{
    if (length(dim(object)) == 2)
        apply(object, 2, summary)
    else
        apply(object, 2:3, summary)
}

