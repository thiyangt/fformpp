##' Trajectory MCMC.
##'
##' Details.
##' @name
##' @title
##' @param ...
##' @param interval e.g. 10%
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
MCMC.trajectory <- function(iIter, nIter, iCross, OUT.accept.probs, interval = 0.1)
{
  printIter <- c(1, seq(from = floor(nIter*interval), to = nIter, by =
                        floor(nIter*interval)))
  if(printIter[length(printIter)] != nIter)
    {
      printIter <- c(printIter, nIter)
    }

  if(iIter %in% printIter)
    {
      parsName <- names(OUT.accept.probs)

      if(iIter %in% printIter[1]) ## Print initial information
        {
          ## Print the acceptance prob.
          parsList <- NULL
          ipars <- 0
          for (i in  parsName)
            {
              subsetLen <- dim(OUT.accept.probs[[i]])[1]
              for (j in 1:subsetLen)
                {
                  ipars <- ipars +1
                  parsList[ipars] <- paste(i, "(", j, ")", sep = "")
                }
            }
          welcome <- paste("Posterior acceptance probabilities(%) in ", iCross, "/",
                           nCross, " cross-validation:\n", sep = "")
          cat(welcome)
          cat("----------------------------------------------------------------------\n")

          format.name0 <- c("Done(%)|", parsList)
          format.name <- format(c("Done(%)|", parsList), width = 13, justify = "right")
          cat(format.name, "\n")

        }
      else
        {
          curIter <- which(printIter == iIter)
          preIter <- printIter[curIter-1]
          probs.tmp <- NULL

          for (i in  parsName)
            {
              ## probs.cur.tmp <- as.vector(rowMeans(OUT.accept.probs[[i]][, preIter:iIter, iCross, drop = FALSE]))
              probs.cur.tmp <- as.vector(rowMeans(OUT.accept.probs[[i]][, 1:iIter, iCross, drop = FALSE]))
              probs.tmp <- c(probs.tmp, probs.cur.tmp)
            }
          doneinPer <- iIter/nIter*100
          format.output <- format(c(paste(doneinPer, "|"), round(probs.tmp*100,
                                                                 digits=2)), nsmall=2,
                                  width  =  15, justify = "right") # FIXME: change width

          cat(format.output,"\n")
          if(iIter == nIter)
            {
              cat("----------------------------------------------------------------------\n\n")
            }
        }

    }
}
