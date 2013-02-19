##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param setupfile "character".
##'        The path where the setup script is located.
##'
##' @return "MCMC-details"
##'
##'        This will return the MCMC details into the global environment.
##'
##' @references Li, F., 2012 Covariate-dependent copula
##'
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##'
##' @note Created: Thu Feb 02 13:33:06 CET 2012;
##'       Current: Wed Mar 07 16:56:30 CET 2012.
CplMain <- function(configfile)
{

###----------------------------------------------------------------------------
### Load user setup file
###----------------------------------------------------------------------------
  source(configfile, local = TRUE)

###----------------------------------------------------------------------------
### INITIALIZE THE STORAGE AND DATA STRUCTURE
###----------------------------------------------------------------------------
  ## Generating the numerical tabular for the inverse Kendall's tau

  tauTabular <- kendalltauTabular(CplNM = CplNM, tol = 1e-3)

  ## Indices for training and testing sample according to cross-validation
  crossValidIdx <- set.crossvalid(nObs,crossValidArgs)
  nCrossFold <- length(crossValidIdx[["training"]])

  ## Split the data into folds for cross validation,  if no cross-validation,
  ## only one fold used
  MdlTraining.X <- vector("list", nCrossFold)
  MdlTraining.Y <- vector("list", nCrossFold)

  MdlTesting.X <- vector("list", nCrossFold)
  MdlTesting.Y <- vector("list", nCrossFold)

  MdlMCMC.beta <- vector("list", nCrossFold)
  MdlMCMC.betaIdx <- vector("list", nCrossFold)
  MdlMCMC.par <- vector("list", nCrossFold)
  MdlMCMC.AccProb <- vector("list", nCrossFold)

###----------------------------------------------------------------------------
### Initialize the MCMC
### This section should be in a function so that it can be called easily for
### parallel computing and other purpose. The input should be the
### cross-validation training and testing indices.
### TODO:crossValidIdx should be reorganized
###----------------------------------------------------------------------------

  iCross <- 1
  Training.Idx <- crossValidIdx[["training"]][[iCross]]
  Testing.Idx <- crossValidIdx[["testing"]][[iCross]]

  nTraining <- length(Training.Idx)

  ## envCurr <- environment()
  ## MCMCFun <- function(Traning.Idx, Testing.Idx, parent.env = envCurr)
  ## {

  ## Extract the training and testing data according to cross-validation
  subsetFun <- function(x, idx)x[idx, , drop = FALSE]
  MdlTraining.X <- rapply(object=Mdl.X,
                          f = subsetFun,
                          idx = Training.Idx,
                          how = "replace")
  MdlTraining.Y <- rapply(object=Mdl.Y,
                          f = subsetFun,
                          idx = Training.Idx,
                          how = "replace")
  MdlTesting.X <- rapply(object=Mdl.X,
                         f = subsetFun,
                         idx = Testing.Idx,
                         how = "replace")
  MdlTesting.Y <- rapply(object=Mdl.Y,
                         f = subsetFun,
                         idx = Testing.Idx,
                         how = "replace")
###----------------------------------------------------------------------------
### STABILIZE THE INITIAL VALUES VIA NEWTON ITERATIONS
###----------------------------------------------------------------------------

  ## Switch all the updating indicators ON
  parUpdate <- MCMCUpdate

  ## Initialize "staticArgs" structure
  u <- matrix(NA, dim(MdlTraining.Y[[1]])[1], length(MdlTraining.Y),
              dimnames = list(NULL, names(MdlTraining.Y)))

  staticArgs <- list(Mdl.logPri =  MdlDataStruc,
                     Mdl.par = MdlDataStruc,
                     Mdl.u = u,
                     Mdl.d = u,
                     tauTabular = tauTabular)

  ## Assign the initial values
  initParOut <- initPar(
      varSelArgs = varSelArgs,
      betaInit = betaInit,
      Mdl.X = MdlTraining.X,
      Mdl.Y = MdlTraining.Y)
  Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
  Mdl.beta <- initParOut[["Mdl.beta"]]

  ## Generate initial values that does not let log posterior be -Inf.
  ## Loop and count how many times tried for generating initial values
  optimInit <- FALSE
  betaTest <- NULL

  ## source("/home/fli/workspace/copula/code/inst/scripts/Plot-tau.R")

  if(optimInit == TRUE &&
     any(tolower(unlist(betaInit)) == "random"))
    {
      InitGood <- FALSE
      nLoopInit <- 0
      while(InitGood == FALSE)
        {
          ## Reassign the initial values
          initParOut <- initPar(
              varSelArgs = varSelArgs,
              betaInit = betaInit,
              Mdl.X = MdlTraining.X,
              Mdl.Y = MdlTraining.Y)
          Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
          Mdl.beta <- initParOut[["Mdl.beta"]]

          ## Dry run to obtain the initial "staticArgs"
          ## NOTE: As this is the very first run, the "parUpdate" should be all on for
          ## this time.

          staticArgs <- logPost(
              CplNM = CplNM,
              Mdl.Y = MdlTraining.Y,
              Mdl.X = MdlTraining.X,
              Mdl.beta = Mdl.beta,
              Mdl.betaIdx = Mdl.betaIdx,
              Mdl.parLink = Mdl.parLink,
              varSelArgs = varSelArgs,
              MargisTypes = MargisTypes,
              priArgs = priArgs,
              parUpdate = rapply(parUpdate, function(x) TRUE, how = "replace"),
              staticArgs = staticArgs,
              staticArgsOnly = TRUE,
              parUpdate4Pri = parUpdate)[["staticArgs"]]

          ## Optimize the initial values via BFGS.
          ## NOTE: The variable selection indicators are fixed (not optimized)
          betaVecInit <- parCplSwap(
              betaInput = Mdl.beta,
              Mdl.beta = Mdl.beta,
              Mdl.betaIdx = Mdl.betaIdx,
              parUpdate = parUpdate)

          betaVecOptim <- try(optim(
              par = betaVecInit,
              fn = logPostOptim,
              control = list(fnscale = -1, maxit = 1000),
              method = "L-BFGS-B",
              CplNM = CplNM,
              Mdl.Y = MdlTraining.Y,
              Mdl.X = MdlTraining.X,
              Mdl.beta = Mdl.beta,
              Mdl.betaIdx = Mdl.betaIdx,
              Mdl.parLink = Mdl.parLink,
              varSelArgs = varSelArgs,
              MargisTypes = MargisTypes,
              priArgs = priArgs,
              parUpdate = parUpdate,
              staticArgs = staticArgs,
              parUpdate4Pri = parUpdate), silent = TRUE)

          if(is(betaVecOptim, "try-error") == TRUE ||
             betaVecOptim$convergence != 0L)
            {
              InitGood <- FALSE
            }
          else
            {
              InitGood <- TRUE
              Mdl.beta <- parCplSwap(
                  betaInput = betaVecOptim[["par"]],
                  Mdl.beta = Mdl.beta,
                  Mdl.betaIdx = Mdl.betaIdx,
                  parUpdate = parUpdate)

              ## betaTest <- rbind(betaTest, betaVecOptim[["par"]])
              print(betaVecOptim)
            }

          nLoopInit <- nLoopInit +1

          ## Too many failures,  abort!
          if(nLoopInit >= 10)
            {
              ## InitGood <- TRUE
              warning(paste(
                  " The initializing algorithm failed more that",
                  nLoopInit, "times.\n",
                  "Continue without initial value optimization."), immediate. = TRUE)
              break
            }

        }
    }

  ## Dry run to obtain staticArgs for the initial values
  ## Again this time all the parameters should be updated.
  staticArgs <- logPost(
      CplNM = CplNM,
      Mdl.Y = MdlTraining.Y,
      Mdl.X = MdlTraining.X,
      Mdl.beta = Mdl.beta,
      Mdl.betaIdx = Mdl.betaIdx,
      Mdl.parLink = Mdl.parLink,
      varSelArgs = varSelArgs,
      MargisTypes = MargisTypes,
      priArgs = priArgs,
      parUpdate = rapply(parUpdate, function(x) TRUE, how = "replace"),
      staticArgs = staticArgs,
      staticArgsOnly = TRUE)[["staticArgs"]]

###----------------------------------------------------------------------------
### THE METROPOLIS-HASTINGS WITHIN GIBBS
###----------------------------------------------------------------------------

  ## Allocate the storage for the final parameters in current fold
  MCMC.beta <- MdlDataStruc
  MCMC.betaIdx <- MdlDataStruc
  MCMC.par <- MdlDataStruc
  MCMC.AccProb <- MdlDataStruc
  for(i in names(MdlDataStruc))
    {
      for(j in names(MdlDataStruc[[i]]))
        {
          ncolX.ij <- ncol(MdlTraining.X[[i]][[j]])
          namesX.ij <- colnames(MdlTraining.X[[i]][[j]])

          ## The MCMC storage
          MCMC.betaIdx[[i]][[j]] <- matrix(
              Mdl.betaIdx[[i]][[j]], nIter, ncolX.ij, byrow = TRUE,
              dimnames = list(NULL, namesX.ij))
          MCMC.beta[[i]][[j]] <- matrix(
              Mdl.beta[[i]][[j]], nIter, ncolX.ij, byrow = TRUE,
              dimnames = list(NULL, namesX.ij))
          MCMC.par[[i]][[j]] <- matrix(
              staticArgs[["Mdl.par"]][[i]][[j]], nIter, nTraining, byrow = TRUE)
          ## The Metropolis-Hasting acceptance rate
          MCMC.AccProb[[i]][[j]] <- matrix(NA, c(nIter, 1))
        }
    }

  ## Switch all the updating indicators OFF
  parUpdate <- rapply(parUpdate, function(x) FALSE, how = "replace")

  ## The updating matrix,  ordered
  UpdateMat <- parCplCaller(parUpdate = MCMCUpdate, parUpdateOrder = MCMCUpdateOrder)

  ## The MCMC loops
  for(iIter in 1:nIter)
    {
      ## The Gibbs loop
      for(iUpdate in 1:dim(UpdateMat)[1])
        {
          CompCaller <- UpdateMat[iUpdate, 1]
          parCaller <- UpdateMat[iUpdate, 2]

          ## Switch current updating parameter indicator on
          parUpdate[[CompCaller]][[parCaller]] <- TRUE

          ## Call the mail Metropolis--Hastings
          algmArgs <- propArgs[[CompCaller]][[parCaller]][["algorithm"]]

          if(tolower(algmArgs[["type"]]) == "gnewtonmove")
            {
              ## staticArgs <- list(u = u, Mdl.par = Mdl.par)
              MHOut <- MHWithGNewton(
                  CplNM = CplNM,
                  propArgs = propArgs,
                  varSelArgs = varSelArgs,
                  priArgs = priArgs,
                  parUpdate = parUpdate,
                  Mdl.Y = MdlTraining.Y,
                  Mdl.X = MdlTraining.X,
                  Mdl.beta = Mdl.beta,
                  Mdl.betaIdx = Mdl.betaIdx,
                  MargisTypes = MargisTypes,
                  Mdl.parLink = Mdl.parLink,
                  staticArgs = staticArgs)

              if(MHOut$errorFlag == FALSE)
                {
                  ## Update the MH results to the current parameter structure
                  staticArgs <- MHOut[["staticArgs"]]
                  Mdl.beta[[CompCaller]][[parCaller]] <- MHOut[["beta"]]
                  Mdl.betaIdx[[CompCaller]][[parCaller]] <- MHOut[["betaIdx"]]

                  ## Export the parameters in each cross-validation fold
                  MCMC.beta[[CompCaller]][[parCaller]][iIter, ] <- MHOut[["beta"]]
                  MCMC.betaIdx[[CompCaller]][[parCaller]][iIter, ] <- MHOut[["betaIdx"]]
                  MCMC.AccProb[[CompCaller]][[parCaller]][iIter,] <- MHOut[["accept.prob"]]
                  MCMC.par[[CompCaller]][[parCaller]][iIter, ] <-
                    staticArgs[["Mdl.par"]][[CompCaller]][[parCaller]]
                }
            }
          else
            {
              stop("Unknown proposal algorithm!")
            }
          ## Switch current updating parameter indicator off
          parUpdate[[CompCaller]][[parCaller]] <- FALSE
        }

      ## MCMC trajectory
      MCMC.trajectory(iIter = iIter, nIter = nIter, interval = 0.1,
                      MCMC.beta = MCMC.beta,
                      MCMC.betaIdx = MCMC.betaIdx,
                      MCMC.par = MCMC.par,
                      MCMC.AccProb = MCMC.AccProb)
    }
  out <- list(MCMC.beta = MCMC.beta,
              MCMC.betaIdx = MCMC.betaIdx,
              MCMC.par = MCMC.par,
              MCMC.AccProb = MCMC.AccProb)

###----------------------------------------------------------------------------
### RUN THE MCMC
### Cross-validation is independent for all folds. In order to parallel it via
### mcmlapply() function in the parallel package in R, we write the sequential
### code with malapply().
### ---------------------------------------------------------------------------

  ## Temporally Disabled for debugging mode
  parallel <- FALSE
  if(parallel == TRUE)
    {
      require(parallel)
      ## Use mcmlapply function
    }
  else
    {
      ## Sequential loops over folds.
      ## mapply(FUN = MCMCFun)
    }

  ## Fetch everything at current environment to a list
  out <- as.list(environment())
  list2env(out, envir = .GlobalEnv)
}
