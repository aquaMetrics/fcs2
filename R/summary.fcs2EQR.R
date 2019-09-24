#' Detailed EQR summary
#'
#' Provides a more detailed EQR summary
#'
#' @param object summary object
#'
#' @return detailed EQR summary object
#' @export
#'
#' @examples
#' \dontrun{
#' summary.fcs2EQR(object = data)
#' }
summary.fcs2EQR <-
function(object)
{
    if (length(dim(object)) == 2)
        apply(object, 2, summary)
    else
        apply(object, 2:3, summary)
}

