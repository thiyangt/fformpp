##' Add grids on plot
##'
##' @export
grid2 <- function(x.at = NA, y.at = NA, col = "black", lty="dotted", lwd = 0.5, ...)
  {
    abline(h = y.at, v = x.at, col = col, lty = lty, lwd = lwd, ...)
  }
