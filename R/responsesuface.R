##' @export
responsesurface <- function(Y, x, logpost.fun.name, crossvaid.struc, splineArgs,
                 priorArgs, OUT.Params, Params_Transform, burn.in, LPDS.sampleProp) {

    B = fitted.model[["B"]]


    for(i in 1:nObs)
    {

        for(j in 1:100)
        {

            Params.j <- lapply(OUT.Params, function(x) apply(x[, , j, iCross, drop =
                                                                                  FALSE], c(1, 2), "["))
            caller.log.like <- call(logpost.fun.name,Y = Y.iTesting, x = x.iTesting,
                                    Params = Params.j, callParam = list(id =
                                                                            "likelihood"), priorArgs =
                                                                                               priorArgs, splineArgs = splineArgs, Params_Transform
                                    = Params_Transform)
            log.like <- eval(caller.log.like)
        }
    }

    y = exp(log.like)


    ## Calculate Y = XB



    ## algorithm


    ## Y.....
}
