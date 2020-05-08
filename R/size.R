##' Print object.size
##'
##' @export
size <- function(x, unit = "auto")
{
    if(class(x) == "character")
    {
        x <- eval(as.name(x))
    }
    print(object.size(x), unit = unit)
}
