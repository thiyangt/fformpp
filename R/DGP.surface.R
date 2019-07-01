##' DGP Surface nested
##'
##' A DGP process where the true model nests in the fitted model
##' @title
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Mon Apr 11 09:18:38 CEST 2011;
##'       Current:       Tue Jul 31 08:54:18 CEST 2012.
##' @param n
##' @param p
##' @param q.o
##' @param q.s
##' @param Sigma
##' @param splineArgs
##' @param splineArgs.crl
##' @param sPlotDGP
##' @export
DGP.surface <- function(n, p, q.o, q.s, Sigma, splineArgs, splineArgs.crl,
                        PlotDGP = FALSE)
  {
    ## Initial DGP status
    DGP.OK <- FALSE
    nRun <- 0

    ## What criterion used
    ## TODO: Hard coded,  use as an input argument
    check.crit <- TRUE

    ## If check.cirt == TRUE,  which options to check
    check.pval <- FALSE # If check p-value
    check.nlf <- TRUE # If check nonlinear factor
    check.man <- FALSE # If check manually with plots

    while(DGP.OK ==  FALSE)
      {
        ## nCenters <- ceiling(sqrt(q.s))
        nCenters <- q.s
        weights <- runif(nCenters, 0, 1) # random weights

        means <- matrix(runif(nCenters*q.o, -1, 1), q.o, nCenters, byrow = TRUE)
        sigma0 <- diag(q.o)/10

        sigmas <- array(sigma0, c(q.o, q.o, nCenters))

        ## x.gen <- matrix(rnorm(n*q.o), n, q.o)

        x.gen <- rmixnorm(n, means, sigmas, weights)

        ## knots.s.idx <- sample(1:n, q.s)
        ## knots.s <- x.gen[knots.s.idx, , drop = FALSE] #
        ## TODO: MUST be a bug in rdist() function.
        knots.s <- matrix(runif(q.s*q.o, min(x.gen), max(x.gen)), nrow = q.s)

        knots.gen <- list(thinplate.s = knots.s)

        X.desi <- d.matrix(x = x.gen, knots = knots.gen, splineArgs = splineArgs)

        q <- dim(X.desi)[2]

        ## Generate coefficients matrix
        seq0 <- c(-1, 1)
        b0 <- rep(seq0, floor(q/2))
        if(q%%2 == 1) b0 <- c(b0, seq0[1])
        B <- matrix(b0, q, p)
        if(p >= 2)
          {
            c0 <- rep(c(1, -1), floor(p/2))
            if(p%%2 == 1) c0 <- c(c0, 1)
            C0 <- matrix(c0, q, p, byrow = TRUE)
            B <- B*C0
          }
        ## B <- matrix(runif(q*p), q, p)

        SurfaceMean <- X.desi %*% B

        Errors <-  rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)

        y.gen <- SurfaceMean + Errors

        ## Generate new x for predictions TODO: Ellpses Mardia 39
        x.testing <- matrix(runif(nPred*q.o, min(x.gen), max(x.gen)), nPred)

        ## Signal to noise,  make sure the noise are not too big
        Errors.sd <- apply(Errors, 2, sd)
        Sig2Noise <- mean(abs(SurfaceMean/Errors.sd))

###----------------------------------------------------------------------------
### ## DGP Check
###----------------------------------------------------------------------------
        ## Regression with additive knots To control the DGP has a relevantly
        ## good nonlinear surface, we tried to fit the results with only one
        ## additive model. Just make sure the surface is not too flat.

        if(check.crit == FALSE) # don't check anything. quit
          {
            DGP.OK <- TRUE
            NonlinFactor <- NA
            Sig2Noise <- NA
          }
        else
          {
            ## Option to check p-value
            if(check.pval == TRUE)
              {
                knots.ctrl <- make.knots(x = x.gen, method = "k-means", splineArgs.ctrl)
                X.desi.ctrl <- d.matrix(x = x.gen, knots = knots.ctrl, splineArgs =
                                        splineArgs.ctrl)
                lmgen <- summary(lm(y.gen~x.gen)) # regression with only linear part
                lmgen.ctrl <- summary(lm(y.gen~0+X.desi.ctrl)) # regression with one knot

                pval <- matrix(NA, (dim(x.gen)[2]+1), p)
                pval.ctrl <- matrix(NA, dim(X.desi.ctrl)[2], p)
                R2 <- matrix(NA, p, 1)

                if(p>1)
                  {
                    for(i in 1:p)
                      {
                        pval[, i] <- lmgen[[i]][[4]][, 4]
                        pval.ctrl[, i] <- lmgen.ctrl[[i]][[4]][, 4]
                        R2[i] <- lmgen[[i]]$adj.r.squared
                      }
                  }
                else
                  {
                    pval[, 1] <- lmgen[[4]][, 4]
                    pval.ctrl[, 1] <- lmgen.ctrl[[4]][, 4]
                    R2 <- lmgen$adj.r.squared
                  }
                if(all(R2 <= 0.1) && all(pval.ctrl <= 0.1))
                  {
                    DGP.OK <- TRUE # PASS
                    NonlinFactor <- NA
                    Sig2Noise <- NA
                  }
              }

            ## Option to check nonlinear factor
            if(check.nlf == TRUE)
              {
                ## Compute the nonlinear factor
                ## TODO: Think about the intercept seriously

                B.lin <- matrix(lm(y.gen~x.gen)$coef, , p)
                y4ctrl <- cbind(1, x.gen)%*%B.lin
                NLMean <- matrix(y4ctrl-SurfaceMean, , p)
                NonlinFactor <- apply(NLMean, 2, sd) # TODO: relevant measure for
                                        # different dataset,  correlations

                if(all(NonlinFactor > 3)) # ad-hoc let the NL in 10-100. when q =
                                        # 10,  TODO: remember to change this later.
                  {
                    DGP.OK <- TRUE
                  }
              }

            ## Option to check manually via plot
            if(check.man == TRUE)
              {
                plot(x.gen, y.gen)
                points(sort(x.gen), SurfaceMean[order(x.gen)], type = "l", col = "red")
                answer <- substr(readline("Is this OK (yes/no)?  "), 1L, 1L)
                if(tolower(answer) %in% c("yes", "y"))
                  {
                    DGP.OK <- TRUE
                    NonlinFactor <- NA
                    Sig2Noise <- NA
                  }
                else
                  {
                    DGP.OK <- FALSE
                  }
              }
          }

        nRun <- nRun + 1
      }

    out <- list(x = x.gen, Y = y.gen, knots = knots.gen, Errors = Errors, B = B, SurfaceMean =
                SurfaceMean, nRun = nRun, NonlinFactor = NonlinFactor, Sig2Noise =
                Sig2Noise, x.testing = x.testing)

    ## Plotting option
    if(PlotDGP && q.o == 2 && p == 1)
      {
        require(rgl)
        plot3d(x = x.gen[, 1], y = x.gen[, 2], z = y.gen, col = "red",
               xlab = "X1", ylab = "X2", zlab = "y")
      }

    return(out)
  }

##----------------------------------------------------------------------------------------
## TESTS:
##----------------------------------------------------------------------------------------
## n <- 200
## p <- 1
## q.o <- 4
## q.s <- 5

## splineArgs <- list(comp = c("intercept", "covariates", "thinplate.s"), # the
##                                         # components of the design matrix.
##                    thinplate.s.dim = c(q.s, q.o), # the dimension of the knots for surface.
##                    thinplate.a.locate = c(0, 0, 0, 0)) # no. of knots used in each
##                                         # covariates for the additive part. zero means no
##                                         # knots for that covariates
