##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param CplConfigFile "character".
##'        The path where the setup script is located.
##'
##' @param Mdl.Idx.training "vector"
##'
##' @return "MCMC-details"
##'
##'        This will return the MCMC details into the global environment.
##'
##' @references Li, F., 2014 Covariate-dependent copula
##'
##' @author Feng Li, Central University of Finance and Economics.
##'
##' @note Created: Thu Feb 02 13:33:06 CET 2012; Current: Fri Mar 27 12:08:32 CST 2015.
CplMain <- function(Mdl.Idx.training, CplConfigFile)
{
###----------------------------------------------------------------------------
### LOAD USER SETUP FILE
###----------------------------------------------------------------------------
  ## Load dependences
  require("mvtnorm", quietly = TRUE)

  ## Load the sourceDir tool
  R_CPL_LIB_ROOT_DIR <- Sys.getenv("R_CPL_LIB_ROOT_DIR")
  if(length(R_CPL_LIB_ROOT_DIR) == 0L)
    {
      stop("R_CPL_LIB_ROOT_DIR is not set properly!")
    }

  sys.source(file.path(R_CPL_LIB_ROOT_DIR, "R/flutils/sourceDir.R"),
             envir = .GlobalEnv)

  ## Load the whole library
  Cpl.source <- sourceDir(file.path(R_CPL_LIB_ROOT_DIR, "R"),
                          byte.compile = 0,
                          recursive = TRUE,
                          ignore.error = TRUE)

  ## Source the configuration file for the model
  source(CplConfigFile, local = TRUE)


###----------------------------------------------------------------------------
### INITIALIZE THE DATA STRUCTURE AND INITIAL VALUES
###----------------------------------------------------------------------------

  ## Extract the training and testing data according to cross-validation
  subsetFun <- function(x, idx)x[idx, , drop = FALSE]

  ## Mdl.Idx.training <- crossValidIdx[["training"]][[iCross]] # obtained from inputs
  nTraining <- length(Mdl.Idx.training)
  Mdl.X.training <- rapply(object=Mdl.X, f = subsetFun,
                           idx = Mdl.Idx.training, how = "replace")
  Mdl.Y.training <- rapply(object=Mdl.Y, f = subsetFun,
                           idx = Mdl.Idx.training, how = "replace")

  ## Assign the initial values
  initParOut <- initPar(varSelArgs = varSelArgs, betaInit = betaInit,
                        Mdl.X = Mdl.X.training, Mdl.Y = Mdl.Y.training,
                        Mdl.parLink = Mdl.parLink)

  Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
  Mdl.beta <- initParOut[["Mdl.beta"]]

###----------------------------------------------------------------------------
### STABILIZE THE INITIAL VALUES VIA NEWTON ITERATIONS
###----------------------------------------------------------------------------

  ## Generate initial values that does not let log posterior be -Inf.
  ## Loop and count how many times tried for generating initial values
  optimInit <- FALSE

  if(optimInit == TRUE &&
     any(tolower(unlist(betaInit)) == "random"))
    {
      require(optimx)

      Mdl.Idx.training.sample <- Mdl.Idx.training[seq(1, nTraining, length.out = 30)]
      Mdl.X.training.sample <- rapply(object=Mdl.X, f = subsetFun,
                                      idx = Mdl.Idx.training.sample, how = "replace")
      Mdl.Y.training.sample <- rapply(object=Mdl.Y, f = subsetFun,
                                      idx = Mdl.Idx.training.sample, how = "replace")

      cat("Optimizing initial values, may take a few minutes...\n\n")

      for(CompCaller in names(Mdl.beta))
        {
          ## If nothing to update, optimization inside this components skipped.
          if(all(unlist(MCMCUpdate[[CompCaller]]) == FALSE)) next

          cat("\nInitializing model component:", CompCaller, "...\n")
          InitGoodCompCurr <- FALSE
          nLoopInit <- 0
          maxLoopInit <- 3

          ## Only current component is updated.
          parUpdate <- rapply(MCMCUpdate, function(x) FALSE,
                              how = "replace")
          parUpdate[[CompCaller]] <- MCMCUpdate[[CompCaller]]

          while(InitGoodCompCurr == FALSE)
            {
              ## Reassign the initial values
              initParOut.CompCurr <- initPar(
                      varSelArgs = varSelArgs,
                      betaInit = betaInit,
                      Mdl.X = Mdl.X.training.sample,
                      Mdl.Y = Mdl.Y.training.sample,
                      Mdl.parLink = Mdl.parLink)
              Mdl.betaIdx[[CompCaller]] <- initParOut.CompCurr[["Mdl.betaIdx"]][[CompCaller]]
              Mdl.beta[[CompCaller]] <- initParOut.CompCurr[["Mdl.beta"]][[CompCaller]]

              ## Dry run to obtain the initial "staticCache" NOTE: As this is the
              ## very first run, the "parUpdate" should be all on for this time.

              staticCache.sample <- logPost(
                      CplNM = CplNM,
                      Mdl.Y = Mdl.Y.training.sample,
                      Mdl.X = Mdl.X.training.sample,
                      Mdl.beta = Mdl.beta,
                      Mdl.betaIdx = Mdl.betaIdx,
                      Mdl.parLink = Mdl.parLink,
                      varSelArgs = varSelArgs,
                      MargisTypes = MargisTypes,
                      priArgs = priArgs,
                      parUpdate = MCMCUpdate,
                      MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

              ## Optimize the initial values via BFGS. NOTE: The variable selection
              ## indicators are fixed (not optimized) loop over all the marginal
              ## models and copula via TWO STAGE OPTIMIZATION

              betaVecInitComp <- parCplSwap(
                      betaInput = Mdl.beta,
                      Mdl.beta = Mdl.beta,
                      Mdl.betaIdx = Mdl.betaIdx,
                      parUpdate = parUpdate)

              ## Optimize the initial values
              betaVecOptimComp <- try(optimx(
                      par = betaVecInitComp,
                      fn = logPostOptim,
                      control = list(maximize = TRUE,
                        ## all.methods = TRUE,
                        maxit = 100),
                      method = "BFGS",
                      hessian = FALSE,
                      CplNM = CplNM,
                      Mdl.Y = Mdl.Y.training.sample,
                      Mdl.X = Mdl.X.training.sample,
                      Mdl.beta = Mdl.beta,
                      Mdl.betaIdx = Mdl.betaIdx,
                      Mdl.parLink = Mdl.parLink,
                      varSelArgs = varSelArgs,
                      MargisTypes = MargisTypes,
                      priArgs = priArgs,
                      staticCache = staticCache.sample,
                      parUpdate = parUpdate,
                      MCMCUpdateStrategy = "twostage"
                      ), silent = FALSE)

              if(is(betaVecOptimComp, "try-error") == TRUE)
                {# It does not have to be converged.
                  cat("Initializing algorithm failed,  retry...\n")
                  InitGoodCompCurr <- FALSE
                  ## break
                }
              else
                {
                  InitGoodCompCurr <- TRUE
                  Mdl.beta <- parCplSwap(
                          betaInput = as.numeric(betaVecOptimComp[1,
                            1:length(betaVecOptimComp)]),
                          Mdl.beta = Mdl.beta,
                          Mdl.betaIdx = Mdl.betaIdx,
                          parUpdate = parUpdate)
                }
              nLoopInit <- nLoopInit +1
              if((nLoopInit >= maxLoopInit) & InitGoodCompCurr  == FALSE)
                {
                  ## InitGood <- TRUE
                  ## Too many failures,  abort!
                  cat("The initializing algorithm failed more that", nLoopInit, "times.\n")
                  cat("Trying to continue without initial value optimization in this component.\n\n")
                  break
                }
            }
          cat("The initial values for beta coefficients are:\n(conditional on variable selection indicators)\n")
          print(Mdl.beta[[CompCaller]])
        }
    }

###----------------------------------------------------------------------------
###  ALLOCATE THE STORAGE
###----------------------------------------------------------------------------

  ## The final parameters in current fold
  MCMC.beta <- MdlDataStruc
  MCMC.betaIdx <- MdlDataStruc
  MCMC.par <- MdlDataStruc
  MCMC.AccProb <- MdlDataStruc

  if(!exists("MCMC.density"))
    {
      MCMC.density <- list()
      MCMC.density[["d"]] <- array(NA, c(nTraining, length(MargisNM),  nIter))
      MCMC.density[["u"]] <- array(NA, c(nTraining, length(MargisNM),  nIter))
      ## FIXME: This is really big ~ 1G
    }

  for(i in names(MdlDataStruc))
    {
      for(j in names(MdlDataStruc[[i]]))
        {
          ncolX.ij <- ncol(Mdl.X.training[[i]][[j]])
          nPar.ij <- Mdl.parLink[[i]][[j]][["nPar"]]
          namesX.ij <- rep(colnames(Mdl.X.training[[i]][[j]]), nPar.ij)

          ## The MCMC storage
          MCMC.betaIdx[[i]][[j]] <- matrix(NA, nIter, ncolX.ij*nPar.ij,
                                           dimnames = list(NULL, namesX.ij))
          MCMC.beta[[i]][[j]] <- matrix(NA, nIter, ncolX.ij*nPar.ij,
                                        dimnames = list(NULL, namesX.ij))
          MCMC.par[[i]][[j]] <- array(NA, c(nIter, nTraining, nPar.ij))

          ## The Metropolis-Hasting acceptance rate
          MCMC.AccProb[[i]][[j]] <- matrix(NA, nIter, 1)
        }
    }

###----------------------------------------------------------------------------
### THE GIBBS SAMPLER (WITH METROPOLIS-HASTINGS)
###----------------------------------------------------------------------------

  cat("Posterior sampling using Metropolis-Hastings within Gibbs\n")


## Dry run to obtain staticcache for the initial values. Again this time all the
## parameters should be updated.

  staticCache <- logPost(
          CplNM = CplNM,
          Mdl.Y = Mdl.Y.training,
          Mdl.X = Mdl.X.training,
          Mdl.beta = Mdl.beta,
          Mdl.betaIdx = Mdl.betaIdx,
          Mdl.parLink = Mdl.parLink,
          varSelArgs = varSelArgs,
          MargisTypes = MargisTypes,
          priArgs = priArgs,
          parUpdate = MCMCUpdate,
          MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

  ## Clear all warnings during initial value optimization. NOTE: not working
  ## warningsClear(envir = environment())

  ## The updating matrix
  UpdateMat <- parCplRepCaller(CplNM = CplNM,
                               parUpdate = MCMCUpdate,
                               parUpdateOrder = MCMCUpdateOrder)
  nInner <- nrow(UpdateMat)

  Starting.time <- Sys.time()

  for(iUpdate in 1:(nInner*nIter))
    {
      iInner <- ifelse((iUpdate%%nInner) == 0, nInner, iUpdate%%nInner)
      iIter <- floor((iUpdate-1)/nInner)+1


      ##print(c(iUpdate, iInner, iIter))

      chainCaller <- UpdateMat[iInner, ]
      CompCaller <- UpdateMat[iInner, 1] ## SP100, SP600,  BB7
      parCaller <- UpdateMat[iInner, 2] ## mean,  df...

      ## Switch all the updating indicators OFF and switch current updating parameter
      ## indicator ON.
      parUpdate <- rapply(MCMCUpdate, function(x) FALSE, how = "replace")
      parUpdate[[CompCaller]][[parCaller]] <- TRUE

      ## Call the Metropolis-Hastings algorithm
      MHOut <- MetropolisHastings(
              CplNM = CplNM,
              propArgs = propArgs,
              varSelArgs = varSelArgs,
              priArgs = priArgs,
              parUpdate = parUpdate,
              Mdl.Y = Mdl.Y.training,
              Mdl.X = Mdl.X.training,
              Mdl.beta = Mdl.beta,
              Mdl.betaIdx = Mdl.betaIdx,
              MargisTypes = MargisTypes,
              Mdl.parLink = Mdl.parLink,
              staticCache = staticCache,
              MCMCUpdateStrategy = MCMCUpdateStrategy)

      if(MHOut$errorFlag == FALSE)
        {
          ## Update the MH results to the current parameter structure
          staticCache <- MHOut[["staticCache"]]
          Mdl.beta[[CompCaller]][[parCaller]] <- MHOut[["beta"]]
          Mdl.betaIdx[[CompCaller]][[parCaller]] <- MHOut[["betaIdx"]]

          ## Export the parameters in each cross-validation fold
          MCMC.beta[[CompCaller]][[parCaller]][iIter, ] <- MHOut[["beta"]]
          MCMC.betaIdx[[CompCaller]][[parCaller]][iIter, ] <- MHOut[["betaIdx"]]
          MCMC.AccProb[[CompCaller]][[parCaller]][iIter,] <- MHOut[["accept.prob"]]
          MCMC.par[[CompCaller]][[parCaller]][iIter, ,] <-
            staticCache[["Mdl.par"]][[CompCaller]][[parCaller]]

          ## Save the marginal densities.
          ## FIXME: This is not updated for every iteration.
          MCMC.density[["u"]][, , iIter] <-
            MHOut[["staticCache"]][["Mdl.u"]]
          MCMC.density[["d"]][, , iIter] <-
            MHOut[["staticCache"]][["Mdl.d"]]
        }
      else
        {
          ## Set acceptance probability to zero if this
          ## iteration fails.
          MCMC.AccProb[[CompCaller]][[parCaller]][iIter,] <- 0
        }

      ## MCMC trajectory
      if(track.MCMC == TRUE && iInner == nInner)
        {
          CplMCMC.summary(iIter = iIter, nIter = nIter,
                          interval = 0.1, burnin = burnin,
                          MCMC.beta = MCMC.beta,
                          MCMC.betaIdx = MCMC.betaIdx,
                          MCMC.par = MCMC.par,
                          MCMC.AccProb = MCMC.AccProb,
                          Starting.time = Starting.time,
                          MCMCUpdate = MCMCUpdate)
        }

    }


  ## Fetch everything at current environment to a list
  ## list2env(out, envir = .GlobalEnv)
  out <- as.list(environment())
  return(out)
}
