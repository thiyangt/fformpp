##' Very Simple progress bar for "for loops"
##'
##' @export
progressbar <- function(iIter, nIter)
{
    cat("\n")
    setTxtProgressBar(txtProgressBar(min = 0, max = nIter, style = 3), iIter)
    if(iIter == nIter)
    {cat("\n")}
}
