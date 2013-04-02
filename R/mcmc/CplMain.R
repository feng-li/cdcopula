##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param CplConfigFile "character".
##'        The path where the setup script is located.
##'
##' @param Training.Idx "vector"
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
CplMain <- function(Training.Idx, CplConfigFile)
{
###----------------------------------------------------------------------------
### Load user setup file
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
### INITIALIZE THE STORAGE AND DATA STRUCTURE
###----------------------------------------------------------------------------

###----------------------------------------------------------------------------
### Initialize the MCMC
###----------------------------------------------------------------------------

  ## Extract the training and testing data according to cross-validation
  subsetFun <- function(x, idx)x[idx, , drop = FALSE]

  ## Training.Idx <- crossValidIdx[["training"]][[iCross]] # obtained from inputs
  nTraining <- length(Training.Idx)
  MdlTraining.X <- rapply(object=Mdl.X,
                          f = subsetFun,
                          idx = Training.Idx,
                          how = "replace")
  MdlTraining.Y <- rapply(object=Mdl.Y,
                          f = subsetFun,
                          idx = Training.Idx,
                          how = "replace")

  ## Switch all the updating indicators ON
  parUpdate <- MCMCUpdate

  ## Assign the initial values
  initParOut <- initPar(
      varSelArgs = varSelArgs,
      betaInit = betaInit,
      Mdl.X = MdlTraining.X,
      Mdl.Y = MdlTraining.Y)

  Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
  Mdl.beta <- initParOut[["Mdl.beta"]]

###----------------------------------------------------------------------------
### STABILIZE THE INITIAL VALUES VIA NEWTON ITERATIONS
###----------------------------------------------------------------------------

  ## Generate initial values that does not let log posterior be -Inf.
  ## Loop and count how many times tried for generating initial values
  optimInit <- TRUE

  ## source("/home/fli/workspace/copula/code/inst/scripts/Plot-tau.R")

  if(optimInit == TRUE &&
     any(tolower(unlist(betaInit)) == "random"))
    {
      InitGood <- FALSE
      nLoopInit <- 0

      cat("Optimizing initial values, may take a few minutes...\n\n")
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

          ## Dry run to obtain the initial "staticCache"
          ## NOTE: As this is the very first run, the "parUpdate" should be all on for
          ## this time.

          staticCache <- logPost(
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
              call.out = "staticCache")[["staticCache"]]

          ## Optimize the initial values via BFGS.
          ## NOTE: The variable selection indicators are fixed (not optimized)

          ## loop over all the marginal models and copula via two stage
          ## optimization

          for(iComp in names(Mdl.beta))
            {

              ## Only current component is updated.
              parUpdateComp <- rapply(parUpdate, function(x) FALSE, how = "replace")
              parUpdateComp[[iComp]] <- parUpdate[[iComp]]


              betaVecInitComp <- parCplSwap(
                  betaInput = Mdl.beta,
                  Mdl.beta = Mdl.beta,
                  Mdl.betaIdx = Mdl.betaIdx,
                  parUpdate = parUpdateComp)

              ## Optimize the initial values
              betaVecOptimComp <- try(optim(
                  par = betaVecInitComp,
                  fn = logPostOptim,
                  control = list(fnscale = -1, maxit = 100),
                  method = "BFGS",
                  CplNM = CplNM,
                  Mdl.Y = MdlTraining.Y,
                  Mdl.X = MdlTraining.X,
                  Mdl.beta = Mdl.beta,
                  Mdl.betaIdx = Mdl.betaIdx,
                  Mdl.parLink = Mdl.parLink,
                  varSelArgs = varSelArgs,
                  MargisTypes = MargisTypes,
                  priArgs = priArgs,
                  staticCache = staticCache,
                  parUpdate = parUpdateComp,
                  split = TRUE,
                  ), silent = FALSE)


              if(is(betaVecOptim, "try-error") == TRUE) # It does not have to be converged.
                {
                  cat("Initializing algorithm failed,  retry again...\n")
                  InitGood <- FALSE
                  break
                }
              else
                {
                  InitGood <- TRUE
                  Mdl.beta <- parCplSwap(
                      betaInput = betaVecOptimComp[["par"]],
                      Mdl.beta = Mdl.beta,
                      Mdl.betaIdx = Mdl.betaIdx,
                      parUpdate = parUpdateComp)
                }

            }


          nLoopInit <- nLoopInit +1

          ## Too many failures,  abort!
          if(nLoopInit >= 1)
            {
              ## InitGood <- TRUE
              cat("The initializing algorithm failed more that", nLoopInit, "times.\n")
              cat("Trying to continue without initial value optimization.\n\n")
              break
            }


          ## Update with BFGS method,  fast but might collapse
          ## betaVecOptim <- try(optim(
          ##     par = betaVecInit,
          ##     fn = logPostOptim,
          ##     control = list(fnscale = -1, maxit = 100),
          ##     method = "BFGS",
          ##     CplNM = CplNM,
          ##     Mdl.Y = MdlTraining.Y,
          ##     Mdl.X = MdlTraining.X,
          ##     Mdl.beta = Mdl.beta,
          ##     Mdl.betaIdx = Mdl.betaIdx,
          ##     Mdl.parLink = Mdl.parLink,
          ##     varSelArgs = varSelArgs,
          ##     MargisTypes = MargisTypes,
          ##     priArgs = priArgs,
          ##     staticCache = staticCache,
          ##     parUpdate = parUpdate), silent = FALSE)

        }
    }

  ## Dry run to obtain staticCache for the initial values
  ## Again this time all the parameters should be updated.
  staticCache <- logPost(
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
      call.out = "staticCache")[["staticCache"]]

###----------------------------------------------------------------------------
### THE METROPOLIS-HASTINGS WITHIN GIBBS
###----------------------------------------------------------------------------
  cat("Posterior sampling using Metropolis-Hastings within Gibbs\n")

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
              staticCache[["Mdl.par"]][[i]][[j]], nIter, nTraining, byrow = TRUE)
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
              ## staticCache <- list(u = u, Mdl.par = Mdl.par)
              MHOut <- MHWithGNewtonMove(
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
                  staticCache = staticCache)

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
                  MCMC.par[[CompCaller]][[parCaller]][iIter, ] <-
                    staticCache[["Mdl.par"]][[CompCaller]][[parCaller]]


                  ## if(CompCaller  == "SP100"& parCaller  == "df"&
                  ##    iIter>1  & MHOut[["accept.prob"]]>0.3 &
                  ##    all(MCMC.par[[CompCaller]][[parCaller]][iIter, ] ==
                  ##    MCMC.par[[CompCaller]][[parCaller]][iIter-1, ])) browser()
                  ## print(MCMC.par[[CompCaller]][[parCaller]][iIter, 1])

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
      CplMCMC.summary(iIter = iIter, nIter = nIter, interval = 0.1, burnin = burnin,
                      MCMC.beta = MCMC.beta,
                      MCMC.betaIdx = MCMC.betaIdx,
                      MCMC.par = MCMC.par,
                      MCMC.AccProb = MCMC.AccProb,
                      MCMCUpdate = MCMCUpdate)
    }

  ## Fetch everything at current environment to a list
  ## list2env(out, envir = .GlobalEnv)
  out <- as.list(environment())
  return(out)
}
