##' Organize the subsets of the parameters by taking away the fixed parameters.
##'
##' <details>
##' @title
##' @param p
##' @param splineArgs
##' @param Params_Fixed
##' @param Params_subsetsArgs
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
Params.subsets <- function(p, splineArgs, Params_Fixed, Params_subsetsArgs)
  {
    out <- list()
##----------------------------------------------------------------------------------------
    ## knots
##----------------------------------------------------------------------------------------
    s.dim <- splineArgs$thinplate.s.dim
    a.locate <- splineArgs$thinplate.a.locate

    knots.maxlen <- prod(s.dim) + sum(a.locate) # The maximum number of knots parameters
                                        # maybe used

    ## The index matrix
    knots.idx <- knots_mat2list(1:knots.maxlen, splineArgs)
    knots.comp <- names(knots.idx)

    knots.subsets <- list()

    ## Split the surface part(if exists)
    if("thinplate.s" %in% knots.comp)
      {
        knots.s.partiArgs <- Params_subsetsArgs[["knots"]][["thinplate.s"]]

        knots.s.idx <- knots.idx[["thinplate.s"]]
        knots.s <- 1:s.dim[1]
        fixed.s <- Params_Fixed[["knots"]]$thinplate.s
        knots.s.remain <- knots.s[! knots.s  %in%  fixed.s]
        length.s.remain <- length(knots.s.remain)

        if(length.s.remain != 0) # not all surface knots are fixed
          {
            ## The full index for the unfixed knots
            knots.s.remain.idx <- knots.s.idx[as.vector(knots.s.remain), , drop = FALSE]
            length.knots.s.remain.idx <- length(knots.s.remain.idx)

            ## Partition the index
            s.subsets.idx <- data.partition(length.knots.s.remain.idx, knots.s.partiArgs)

            for(i in 1:length(s.subsets.idx))
              {
                ## surface knots are labeled row-wise
                knots.subsets[["thinplate.s"]][[i]] <-
                  knots.s.remain.idx[as.vector(s.subsets.idx[[i]])]
              }
          }
        else
          {
            knots.subsets[["thinplate.s"]] <- NULL
          }

      }
    ## Split the additive part(if exists)
    if(("thinplate.a" %in% knots.comp))
      {
        knots.a.partiArgs <- Params_subsetsArgs[["knots"]][["thinplate.a"]]

        knots.a.idx <- knots.idx[["thinplate.a"]]
        knots.a <- 1:sum(a.locate)
        fixed.a <- Params_Fixed[["knots"]]$thinplate.a
        knots.a.remain <- knots.a[! knots.a  %in%  fixed.a]
        length.a.remain <- length(knots.a.remain)
        if(length.a.remain  != 0)
          {
            ## The full index for the unfixed knots
            knots.a.remain.idx <- knots.a.idx[as.vector(knots.a.remain)]
            length.knots.a.remain.idx <- length(knots.a.remain.idx)

            a.subsets.idx <- data.partition(length.a.remain, knots.a.partiArgs)

            for(i in 1:length(a.subsets.idx))
              {
                knots.subsets[["thinplate.a"]][[i]] <-
                  knots.a.remain.idx[as.vector(a.subsets.idx[[i]])]
              }
          }
        else
          {
            knots.subsets[["thinplate.a"]] <- NULL
          }
      }

    ## Merge or split
    ## Special case when the subsets in the surface and additive components are
    ## the of the same length
    knots.split <- Params_subsetsArgs[["knots"]]$split
    no.subsets <- unlist(lapply(knots.subsets, length))

    if(length(no.subsets) == 2 &&
       no.subsets[1] == no.subsets[2] &&
       knots.split == FALSE) # merge
      {

        knots.N.subsets <- no.subsets[1]
        out.knots <- list()
        for(i in 1:knots.N.subsets)
          {
            out.knots[[i]] <- c(knots.subsets[["thinplate.s"]][[i]],
                                knots.subsets[["thinplate.a"]][[i]])
          }
      }
    else # split
      {
        out.knots <- unlist(knots.subsets, recursive = FALSE)
        names(out.knots) <- NULL
      }

    ## output
    if(length(out.knots) == 0)
      {
        out[["knots"]] <- list(NULL)
      }
    else
      {
        out[["knots"]] <- out.knots
      }

##----------------------------------------------------------------------------------------
    ## shrinkages
##----------------------------------------------------------------------------------------

    shrinkages.fixed <- Params_Fixed[["shrinkages"]]
    model.comp <- splineArgs[["comp"]][ "intercept" != splineArgs[["comp"]] ]
    ncomp <- length(model.comp)
    shrinkages.idx <- 1:(p*ncomp)

    shrinkages.remain <- shrinkages.idx[! shrinkages.idx %in% shrinkages.fixed]
    shrinkages.remain.len <- length(shrinkages.remain)

    if(shrinkages.remain.len != 0)
      {
        shrinkages.parti.idx <- data.partition(shrinkages.remain.len,
                                               Params_subsetsArgs[["shrinkages"]])
        out.shrinkages <- list()

        for(i in 1:length(shrinkages.parti.idx))
          {
            out.shrinkages[[i]] <- shrinkages.idx[as.vector(shrinkages.remain)][as.vector(shrinkages.parti.idx[[i]])]
          }
        out[["shrinkages"]] <- out.shrinkages
      }
    else
      {
        out[["shrinkages"]] <- list(NULL)
      }
##----------------------------------------------------------------------------------------
    ## coefficients
##----------------------------------------------------------------------------------------

    out[["coefficients"]] <- list(NULL) # been integrated out
##----------------------------------------------------------------------------------------
    ## covariance
##----------------------------------------------------------------------------------------
    ## FIXME: PARTITION this
    out[["covariance"]] <- list(1:(p*(p+1)/2))
    return(out)
  }

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
## Params.subsets(p, splineArgs, Params_Fixed, Params_subsetsArgs)
