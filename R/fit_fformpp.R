#' Function to fit Efficient Bayesian Surface Regression Models.
#'
#' This function is written based on the example codes in movingknots
#' package by Dr. Feng Li. Please see https://github.com/feng-li/movingknots for more details.
#'
#' @param feamat matrix of features, rows corresponds to each time series,
#' columns coresponds to features
#' @param accmat matrix of forecast errors from each method,
#' rows represent time series, columns represent forecast algorithms
#' @param sknots the dimension of knots for surface, default 2
#' @param aknots no. of knots used in each covariates for the additive part, default 2
#' @param fix.s number of knots to be fixed in the surface components, 0 means all are
#' updated
#' @param fix.a number of knots to be fixed in the additive components, default is 0 which means all are updated
#' @param fix.shrinkage, number of shrinkage covariates not to be updated, defalut is, 1:p,
#' @param fix.covariance, number of knots to be fixed in the covariance, default is 0, all are updated
#' @param fix.coefficients, number of knots to be fixed in the coefficients, default is 0, all are updated
#' @param n.iter number of ierations
#' @param knot.moving.algorithm, select either "KStepNewton" or "Random-Walk", to fasten the
#' code use Random-Walk
#' @param ptype For fixing gprior, This could be c("X'X", "identity", "identity") or c("identity", "identity", "identity")
#' @param prior.knots this could be n, log(n) or 1 # to set priors for knots
#' @return returns a list contanining the fitted model and arguments for splines
#' @export
fit_fformpp <- function(feamat, accmat, sknots=2, aknots=2,
                        fix.s=0, fix.a=0, fix.shrinkage, fix.covariance=0,
                        fix.coefficients=0, n.iter=100,
                        knot.moving.algorithm="Random-Walk",
                        ptype=c("identity", "identity", "identity"),
                        prior.knots){
  ## processing data to arrange according to the format of movingknots functionalities
  colnames(feamat) <- NULL

  ## preparation of "Y" matrix
  Y <- accmat
  colnames(Y) <- NULL

  ## MCMC TRAJECTORY
  track.MCMC = TRUE

  ## standardizing the data
  data <- StdData(feamat, method = "norm-0-1")
  x <<- data[["data"]]

  ## no. of observations
  n <- dim(Y)[1]

  ## no. of dimensions
  p <- dim(Y)[2]

  ## no. of original covariates
  m <- dim(x)[2]

  ##----------------------------------------------------------------------------------------
  ## Model configurations
  ##----------------------------------------------------------------------------------------
  ## MODEL NAME
  Model_Name <- "linear"

  ## ARGUMENTS FOR SPLINES
  splineArgs <- list(
    ## the components of the design matrix.
    comp = c("intercept", "covariates", "thinplate.s", "thinplate.a"),
    ## the dimension of the knots for surface.
    thinplate.s.dim = c(sknots, m), # m is the number of features
    ## no. of knots used in each covariates for the additive part. zero means no knots for
    ## that covariates
    thinplate.a.locate = rep(aknots, m))

  ## PARAMETERS UPDATED USING GIBBS
  ## You have to change this when "splineArgs$comp" has
  ## changed. Coefficients are updated by directly sampling
  Params4Gibbs <- c("knots", "shrinkages", "covariance")

  ## FIXED PARAMETERS
  Params_Fixed <- list(
    ## which knots from which part of model are not updated.
    "knots" = list(thinplate.s = fix.s, thinplate.a = fix.a),
    "shrinkages" = fix.shrinkage, # the shrinkages for covariates not updated
    "covariance"  = fix.covariance,   # zero means all are updated
    "coefficients" = fix.coefficients)

  ## ARGUMENTS FOR PARTITION PARAMETERS (BATCHES UPDATE)
  ## The split argument is only used when surface and additive subsets are of the
  ## same length

  if(knot.moving.algorithm=="KStepNewton"|knot.moving.algorithm=="SGLD"){
    Params_subsetsArgs <- list(
      "knots" = list(thinplate.s = list(N.subsets = 1, partiMethod = "systematic"),
                     thinplate.a = list(N.subsets = 1, partiMethod = "systematic"), split = FALSE),

      "shrinkages" = list(N.subsets = 1, partiMethod = "systematic"),
      "covariance"  = list(N.subsets = 1, partiMethod = "systematic"),
      "coefficients" = list(N.subsets = 1, partiMethod = "systematic"))
  } else {
    Params_subsetsArgs <- list("knots" = list(
      thinplate.s = list(
        N.subsets = 20,
        partiMethod = "systematic"),

      thinplate.a = list(
        N.subsets = 4,
        partiMethod = "systematic"),
      split = FALSE),

      "shrinkages" = list(N.subsets = 1, partiMethod = "systematic"),
      "covariance"  = list(N.subsets = 1, partiMethod = "systematic"),
      "coefficients" = list(N.subsets = 1, partiMethod = "systematic"))
  }

  ##----------------------------------------------------------------------------------------
  ## Parameters settings
  ##----------------------------------------------------------------------------------------

  ## TRANSFORMATION FUNCTION
  Params_Transform <- list("knots" = "identity",
                           "shrinkages" = "log",
                           "covariance" = "identity",
                           "coefficients" = "identity")

  ## HESSIAN METHODS
  hessMethods <- list("knots" = "outer",
                      "shrinkages" = "outer",
                      "covariance" = NA,
                      "coefficients" = NA)

  ## Propose method in Metropolis-Hasting
  if(knot.moving.algorithm=="Random-Walk"){
  propMethods <- list("knots" = "Random-Walk",
                      "shrinkages" = "KStepNewton",
                      "covariance" = "Inverse-Wishart", # random MH without K-step Newton
                      "coefficients" = NA) }
  if(knot.moving.algorithm=="SGLD"){
    propMethods <- list("knots" = "SGLD",
                        "shrinkages" = "SGLD",
                        "covariance" = "Inverse-Wishart", # random MH without K-step Newton
                        "coefficients" = NA) }


  ##----------------------------------------------------------------------------------------
  ## MCMC configurations
  ##----------------------------------------------------------------------------------------

  ## NO. OF ITERATIONS
  nIter <- n.iter

  ## BURN-IN
  burn.in <- 0.2  # [0, 1) If 0: use all MCMC results.

  if(knot.moving.algorithm=="Random-Walk"){
  ## LPDS SAMPLE SIZE
  LPDS.sampleProp <- 0.05 # Sample proportion to the total posterior after burn-in.

  ## CROSS-VALIDATION
  crossValidArgs <- list(N.subsets = 0, # No. of folds. If 0:, no cross-validation.
                           partiMethod = "systematic", # How to partition the data
                           full.run = FALSE)     # Also include a full run.

  ## NO. OF FINTE NEWTON MOVE FOR EACH PARAMETERS
  nNewtonSteps <- list("knots" = 1,
                       "shrinkages" = 1,
                       "covariance" = NA, # random MH
                       "coefficients" = NA) # integrated out

  ## THE DF. FOR A MULTIVARIATE T-PROPOSAL IN MH ALGORITHM.
  MH.prop.df <- list("knots" = 5,
                     "shrinkages" = 5,
                     "covariance" = NA,
                     "coefficients" = NA)}

  if(knot.moving.algorithm=="SGLD"){
    ## LPDS SAMPLE SIZE
    LPDS.sampleProp <- 1 # Sample proportion to the total posterior after burn-in.

    ## CROSS-VALIDATION
    crossValidArgs <- list(N.subsets = 5, # No. of folds. If 0:, no cross-validation.
                           partiMethod = "systematic" # How to partition the data
    )

    algArgs = list(knots = list(minibatchProp = 0.1, nEpoch= 2, calMHAccRate = FALSE, # Welling & Teh (2011), p 3.
                                stepsizeSeq = make_stepsize(
                                  steprange = c(0.01, 0.0001), n = 20 * nIter, # Welling & Teh (2011), p 5.
                                  args = list(method = "exp-decay", lambda = 0.55))),
                   shrinkages = list(minibatchProp = 0.1, nEpoch= 2, calMHAccRate = FALSE,
                                     stepsizeSeq = make_stepsize(
                                       steprange = c(0.01, 0.0001), n = 20 * nIter,
                                       args = list(method = "exp-decay", lambda = 0.45))),
                   covariance = NA,
                   coefficients = NA)
    nInner = 20 # 1/minibatchProp * nEpoch
  }


  ##----------------------------------------------------------------------------------------
  ## Set up Priors
  ##----------------------------------------------------------------------------------------

  ## TODO: The prior should be set in the transformed scale when the linkages is not
  ## "identity". Write a general function to handle this.

  ## Regression
  knots.location.gen <- fformpp::make.knots(x = x, method = "k-means", splineArgs)

  X.init <- fformpp::d.matrix(x, knots = knots.location.gen, splineArgs)
  lm.init <- stats::lm(Y~0+X.init)
  S0.init <- matrix(stats::var(lm.init$residual), p, p)
  q <- dim(X.init)[2]

  ## P MATRIX TYPE
  ## P.type <- c("identity", "identity", "identity") # can be "identity" or "X'X"
  P.type <- ptype # can be "identity" or "X'X"

  ## PRIOR FOR COVARIANCE
  covariance.priType <- "Inverse-Wishart"
  covariance.df0 <- 10
  covariance.S0 <- S0.init # p-by-p, see Mardia p.158

  ## PRIOR FOR COEFFICIENTS
  coefficients.priType <- "mvnorm"
  coefficients.mu0 <- matrix(0, q*p, 1)  # mean of B|Sigma, assume no covariates in.

  ## PRIOR FOR KNOTS
  knots.priType <- "mvnorm"
  knots.mu0 <- knots_list2mat(knots.location.gen) # mean from k-means
  knots.Sigma0 <- make.knotsPriVar(x, splineArgs) # the covariance for each knots came from x'x
  knots.c <- prior.knots # The shrinkage

  ## PRIOR FOR SHRINKAGES

  ## how many components does the model have
  model.comp.len <- length(splineArgs[["comp"]][ "intercept" != splineArgs[["comp"]] ])
  # how many components does the model have
  shrinkages.pri.trans <- convert.densParams(mean = n/2, var = (n/2)^2, linkage =
                                               Params_Transform[["shrinkages"]]) # assume
  # normal prior with "mean" and "var"
  shrinkages.priType <- "mvnorm"
  shrinkages.mu0 <- matrix(rep(shrinkages.pri.trans[1], p*model.comp.len)) # The mean of
  # shrinkage,  "n" is unit information
  # prior. (n*(X'X)^(-1))
  shrinkages.Sigma0 <- diag(rep(shrinkages.pri.trans[2], p), p*model.comp.len) # The variance
  # for the shrinkage parameter.
  shrinkages.c <- n # The shrinkage

  ## Organize the arguments
  priorArgs <- list(P.type = P.type,

                    knots.priType = knots.priType,
                    knots.mu0 = knots.mu0, # prior for knots
                    knots.Sigma0 = knots.Sigma0,
                    knots.c = knots.c,

                    shrinkages.priType = shrinkages.priType,
                    shrinkages.mu0 = shrinkages.mu0, # prior for shrinkages
                    shrinkages.Sigma0 = shrinkages.Sigma0,
                    shrinkages.c  = shrinkages.c,

                    coefficients.priType = coefficients.priType,
                    coefficients.mu0 = coefficients.mu0, # prior for coefficients

                    covariance.priType = covariance.priType,
                    covariance.df0 = covariance.df0, # prior for covariance
                    covariance.S0 = covariance.S0)

  ##----------------------------------------------------------------------------------------
  ## Initial values
  ##----------------------------------------------------------------------------------------
  ## TODO: The initial values should be transformed into the new scale according to the
  ## linkages if it is not "identity"

  ## INITIAL KNOTS LOCATIONS, "list"
  INIT.knots <- knots.location.gen

  ## INITIAL SHRINKAGE FOR MODEL COVARIANCE "matrix"
  INIT.shrinkages <- shrinkages.mu0

  ## INITIAL COVARIANCE "matrix"
  INIT.covariance <- covariance.S0

  ##########################################################################################
  ##                                   System settings
  ##########################################################################################

  ##----------------------------------------------------------------------------------------
  ## Initialize the data
  ##----------------------------------------------------------------------------------------
  ## Gradient function name
  gradhess.fun.name <- tolower(paste(Model_Name, "gradhess", sep = "_"))

  ## Log posterior function name
  logpost.fun.name <-  tolower(paste(Model_Name, "logpost", sep = "_"))

  ##----------------------------------------------------------------------------------------
  ## Set up cross validation etc
  ##----------------------------------------------------------------------------------------

  ## The training($training) and testing($testing) structure.
  ## If no cross-validation, $training is also $testing.
  ## If full run is required, the last list in $training and $testing is for a full run.
  crossvalid.struc <<- fformpp::set.crossvalid(nObs = n, crossValidArgs = crossValidArgs)

  ## No. of total runs
  nCross <<- length(crossvalid.struc$training)

  ## No. of training obs. in each data subset.
  nTraining <- unlist(lapply(crossvalid.struc$training, length))

  ## Params
  Params <- list("knots" = fformpp::knots_list2mat(INIT.knots),
                 "shrinkages" = INIT.shrinkages,
                 "covariance" = fformpp::vech(INIT.covariance),
                 "coefficients" = matrix(NA, q, p))

  ## The parameters subset structures.
  Params.sub.struc <- Params.subsets(p, splineArgs, Params_Fixed, Params_subsetsArgs)

  ##----------------------------------------------------------------------------------------
  ## Construct the output formats
  ##----------------------------------------------------------------------------------------

  ## NOTATIONS TO USE
  ## The output is alway with "OUT.XXX"
  ## The last dimension is always for the i:th cross-validation subsets.

  ## Accept probabilities for MH.
  OUT.accept.probs <- mapply(function(x) array(NA, c(length(x), nIter, nCross)),
                             Params.sub.struc, SIMPLIFY = FALSE)

  ## Parameters updates in each MH step
  INIT.knots.mat <- knots_list2mat(INIT.knots)

  OUT.Params <- list("knots" = array(INIT.knots.mat, c(length(INIT.knots.mat), 1, nIter, nCross)),
                     "shrinkages" = array(INIT.shrinkages, c(p*model.comp.len, 1, nIter, nCross)),
                     "coefficients" = array(NA, c(q, p, nIter, nCross)),
                     "covariance" = array(fformpp::vech(INIT.covariance), c((p+1)*p/2, 1, nIter, nCross)))

  ##########################################################################################
  ##                                 Testings
  ##########################################################################################
  ## See the "tests" folder and tests at end of each function.

  ##########################################################################################
  ##                                   Main algorithm
  ##########################################################################################

  ##----------------------------------------------------------------------------------------
  ## Stabilize the initial values
  ##----------------------------------------------------------------------------------------
  ## see "tests/test.init.BFGS.R" file

  ##----------------------------------------------------------------------------------------
  ## MovingKnots MCMC
  ##----------------------------------------------------------------------------------------
  if(knot.moving.algorithm=="Random-Walk"){
  OUT.FITTED <- MovingKnots_MCMC_rw(gradhess.fun.name = gradhess.fun.name,
                                 logpost.fun.name =  logpost.fun.name,
                                 nNewtonSteps =  nNewtonSteps,
                                 nIter = nIter,
                                 Params = Params,
                                 Params4Gibbs = Params4Gibbs,
                                 Params.sub.struc =  Params.sub.struc,
                                 hessMethods = hessMethods,
                                 Y = Y,
                                 x0 = x,
                                 splineArgs = splineArgs,
                                 priorArgs = priorArgs,
                                 MH.prop.df = MH.prop.df,
                                 Params_Transform = Params_Transform,
                                 propMethods = propMethods,
                                 crossvalid.struc = crossvalid.struc,
                                 OUT.Params = OUT.Params,
                                 OUT.accept.probs = OUT.accept.probs,
                                 burn.in = burn.in,
                                 LPDS.sampleProp = LPDS.sampleProp,
                                 track.MCMC = track.MCMC)}

  if(knot.moving.algorithm=="SGLD"){
    OUT.FITTED <- MovingKnots_MCMC_sgld(gradhess.fun.name = gradhess.fun.name,
                                                logpost.fun.name =  logpost.fun.name,
                                                nIter = nIter,
                                                Params = Params,
                                                Params4Gibbs = Params4Gibbs,
                                                Params.sub.struc =  Params.sub.struc,
                                                Y = Y,
                                                x0 = x,
                                                splineArgs = splineArgs,
                                                priorArgs = priorArgs,
                                                algArgs = algArgs,
                                                Params_Transform = Params_Transform,
                                                propMethods = propMethods,
                                                crossvalid.struc = crossvalid.struc,
                                                OUT.Params = OUT.Params,
                                                OUT.accept.probs = OUT.accept.probs,
                                                burn.in = burn.in,
                                                LPDS.sampleProp = LPDS.sampleProp,
                                                track.MCMC = track.MCMC)}

  ##----------------------------------------------------------------------------------------
  ## Save outputs to files
  ##----------------------------------------------------------------------------------------

  return(list(out.fitted=OUT.FITTED, spline.args=splineArgs))

}

#'@examples
#'require(Mcomp)
#'require(seer)
#'data(M1)
#'yearly_m1 <- subset(M1, "yearly")
#'feamat <- as.matrix(cal_features(yearly_m1, database="M1", h=6, highfreq=FALSE))
#'acccal <- fcast_accuracy(yearly_m1, models= c("arima","ets","rw","rwd", "theta", "nn"), database ="M1", cal_MASE, h=6, length_out = 1, fcast_save = FALSE)
#'accmat <- as.matrix(acccal$accuracy)
#'fformpp <- fit_fformpp(feamat, accmat, sknots=2, aknots=2,
#' fix.s=0, fix.a=0, fix.shrinkage=1:6, fix.covariance=0,
#' fix.coefficients=0, n.iter=100,
#' knot.moving.algorithm="Random-Walk",
#' ptype=c("identity", "identity", "identity"),
#' prior.knots=181)
#'
