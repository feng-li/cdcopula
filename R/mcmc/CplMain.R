##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param MdlConfigFile "character".
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
CplMain <- function(Mdl.Idx.training, MdlConfigFile)
{

###----------------------------------------------------------------------------
### DEBUGGING
###----------------------------------------------------------------------------

  DEBUGGING <- FALSE
  if(DEBUGGING)
  {## Turn warnings into error
    options(warn = 100)
  }
###----------------------------------------------------------------------------
### LOAD USER SETUP FILE
###----------------------------------------------------------------------------
  ## Load dependences
  require("mvtnorm", quietly = TRUE)
  require("numDeriv", quietly = TRUE)

  ## Load the sourceDir tool
  R_CPL_LIB_ROOT_DIR <- Sys.getenv("R_CPL_LIB_ROOT_DIR")
  if(nchar(R_CPL_LIB_ROOT_DIR) == 0L)
  {
    stop("R_CPL_LIB_ROOT_DIR is not set properly!")
  }

  R_CPL_NPARALLEL <- Sys.getenv("R_CPL_NPARALLEL")
  if(nchar(R_CPL_NPARALLEL) == 0L)
  {
    stop("R_CPL_NPARALLEL is not set properly!")
  }


  sys.source(file.path(R_CPL_LIB_ROOT_DIR, "R/flutils/sourceDir.R"),
             envir = .GlobalEnv)

  ## Load the whole library
  Cpl.source <- sourceDir(file.path(R_CPL_LIB_ROOT_DIR, "R"),
                          byte.compile = 0,
                          recursive = TRUE,
                          ignore.error = TRUE)

  ## Source the configuration file for the model
  source(MdlConfigFile, local = TRUE)

  ## Parallel Setting up
  if(as.numeric(R_CPL_NPARALLEL)>1)
  {
    require("parallel", quietly = TRUE)
    require("snow", quietly = TRUE)
    cl4MCMC <- makeCluster(as.numeric(R_CPL_NPARALLEL), type = "MPI")
    setDefaultCluster(cl4MCMC)
    ce4MCMC <- clusterEvalQ(cl4MCMC,
    {sourceDir(file.path(Sys.getenv("R_CPL_LIB_ROOT_DIR"), "R"),
               byte.compile = 0,
               recursive = TRUE,
               ignore.error = TRUE)
    })
  }

  ## Generating simple model information
  Starting.time <- Sys.time()
  ModelDescription <- paste(c(MargisNM[-length(MargisNM)],"+",  MargisType, "+" , "nObs",
                              nObs, "nCross", nCross,  "+",
                              format(Starting.time, "%Y%m%d@%H.%M")),
                            collapse = "")
###----------------------------------------------------------------------------
### INITIALIZE THE DATA STRUCTURE AND INITIAL VALUES
###----------------------------------------------------------------------------
  ## Extract the training and testing data according to cross-validation
  subsetFun <- function(x, idx)x[idx, , drop = FALSE]

  nTraining <- length(Mdl.Idx.training)
  Mdl.Y.training <- rapply(object=Mdl.Y, f = subsetFun,
                           idx = Mdl.Idx.training, how = "replace")

  if(tolower(MargisType[length(MargisType)]) %in% c("gogarch", "dccgarch"))
  {## Special case when a foreign model is introduced
    browser()
    out <-ForeignModelEval(model  = MargisType[length(MargisType)],
                           spec = ForeignModelSpec,
                           data = Mdl.Y.training)
    print(out) # Model summary

    return(out)
  }

  if(any(rapply(Mdl.X, class) != "matrix"))
  { ## Evaluating Foreign marginal models.
    cat("Evaluating foreign marginal models...\n")

    Mdl.X.Fit <- MargiModelForeignEval(MargisNM = MargisNM,
                                       MargisType = MargisType,
                                       MargisForeignConfig = Mdl.X,
                                       Mdl.Y = Mdl.Y.training)

    Mdl.X.training <- c(Mdl.X.Fit[["Mdl.X"]],
                        rapply(object=Mdl.X[MargisNM[length(MargisNM)]],
                               f = subsetFun,
                               idx = Mdl.Idx.training,
                               how = "replace"))
    Mdl.ForeignFit <- Mdl.X.Fit[["Mdl.ForeignFit"]]
  }
  else
  {## Native marginal model structure
    Mdl.X.training <- rapply(object=Mdl.X, f = subsetFun,
                             idx = Mdl.Idx.training, how = "replace")
  }
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
  optimInit <- TRUE
  if(optimInit == TRUE &&
     any(tolower(unlist(betaInit)) == "random"))
  {
    require("optimx")

    if(nTraining<50)
    {
      nOptim.Sample <- nTraining
    }
    else
    {
      nOptim.Sample <- max(min(floor(nTraining*.10), 500), 50)
    }

    Mdl.Idx.training.sample <- seq(1, nTraining, length.out = nOptim.Sample)
    Mdl.X.training.sample <- rapply(object=Mdl.X.training, f = subsetFun,
                                    idx = Mdl.Idx.training.sample, how = "replace")
    Mdl.Y.training.sample <- rapply(object=Mdl.Y.training, f = subsetFun,
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
        initParOut.CompCurr <- initPar(varSelArgs = varSelArgs,
                                       betaInit = betaInit,
                                       Mdl.X = Mdl.X.training.sample,
                                       Mdl.Y = Mdl.Y.training.sample,
                                       Mdl.parLink = Mdl.parLink)
        Mdl.betaIdx[[CompCaller]] <- initParOut.CompCurr[["Mdl.betaIdx"]][[CompCaller]]
        Mdl.beta[[CompCaller]] <- initParOut.CompCurr[["Mdl.beta"]][[CompCaller]]

        ## Dry run to obtain the initial "staticCache" NOTE: As this is the
        ## very first run, the "parUpdate" should be all on for this time.

        staticCache.sample <- logPost(MargisType = MargisType,
                                      Mdl.Y = Mdl.Y.training.sample,
                                      Mdl.X = Mdl.X.training.sample,
                                      Mdl.beta = Mdl.beta,
                                      Mdl.betaIdx = Mdl.betaIdx,
                                      Mdl.parLink = Mdl.parLink,
                                      varSelArgs = varSelArgs,
                                      priArgs = priArgs,
                                      parUpdate = MCMCUpdate,
                                      MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

        ## Optimize the initial values via BFGS. NOTE: The variable selection
        ## indicators are fixed (not optimized) loop over all the marginal
        ## models and copula via TWO STAGE OPTIMIZATION

        betaVecInitComp <- parCplSwap(betaInput = Mdl.beta,
                                      Mdl.beta = Mdl.beta,
                                      Mdl.betaIdx = Mdl.betaIdx,
                                      parUpdate = parUpdate)

        ## Optimize the initial values
        betaVecOptimComp <- try(optimx(par = betaVecInitComp,
                                       fn = logPostOptim,
                                       control = list(maximize = TRUE,
                                                      ## all.methods = TRUE,
                                                      maxit = 100),
                                       method = "BFGS",
                                       hessian = FALSE,
                                       MargisType = MargisType,
                                       Mdl.Y = Mdl.Y.training.sample,
                                       Mdl.X = Mdl.X.training.sample,
                                       Mdl.beta = Mdl.beta,
                                       Mdl.betaIdx = Mdl.betaIdx,
                                       Mdl.parLink = Mdl.parLink,
                                       varSelArgs = varSelArgs,
                                       priArgs = priArgs,
                                       staticCache = staticCache.sample,
                                       parUpdate = parUpdate,
                                       MCMCUpdateStrategy = "twostage"
                                       ), silent = FALSE)

        if(is(betaVecOptimComp, "try-error") == TRUE
           || any(is.na(as.numeric(betaVecOptimComp[1, 1:length(betaVecOptimComp)]))))
        {# It does not have to be converged.
          cat("Initializing algorithm failed,  retry...\n")
          InitGoodCompCurr <- FALSE
          ## break
        }
        else
        {
          InitGoodCompCurr <- TRUE
          Mdl.beta <- parCplSwap(betaInput = as.numeric(betaVecOptimComp[1, 1:length(betaVecOptimComp)]),
                                 Mdl.beta = Mdl.beta,
                                 Mdl.betaIdx = Mdl.betaIdx,
                                 parUpdate = parUpdate)
        }
        nLoopInit <- nLoopInit +1
        if((nLoopInit >= maxLoopInit) & InitGoodCompCurr  == FALSE)
        {## Too many failures,  abort!
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
  MCMC.beta <- MCMCUpdate
  MCMC.betaIdx <- MCMCUpdate
  MCMC.par <- MCMCUpdate
  MCMC.AccProb <- MCMCUpdate

  if(!exists("MCMC.density"))
  {
    ## MCMC.density <- list()
    ## MCMC.density[["d"]] <- array(NA, c(nTraining, length(MargisNM) +1,  MCMC.nIter))
    ## MCMC.density[["u"]] <- array(NA, c(nTraining, length(MargisNM),  MCMC.nIter))
    ## FIXME: This is really big ~ 1G
  }

  for(CompCaller in names(MCMCUpdate))
  {
    for(parCaller in names(MCMCUpdate[[CompCaller]]))
    {
      ncolX.ij <- ncol(Mdl.X.training[[CompCaller]][[parCaller]])
      nPar.ij <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]
      namesX.ij <- rep(colnames(Mdl.X.training[[CompCaller]][[parCaller]]), nPar.ij)

      if(CompCaller %in% MargisNM[length(MargisNM)])
      {
        ## browser()
        nDim <- length(MargisNM)-1
        namesParFull.ij <- matrix(paste(matrix(1:nDim, nDim, nDim),
                                        matrix(1:nDim, nDim, nDim, byrow = TRUE), sep = "."), nDim)
        namesPar.ij <- namesParFull.ij[lower.tri(namesParFull.ij, diag = FALSE)]
      }
      else
      {
        namesPar.ij <- "1.1"
      }


      nMCMC <- ifelse(MCMCUpdate[[CompCaller]][[parCaller]], MCMC.nIter, 1)

      ## The MCMC storage
      MCMC.betaIdx[[CompCaller]][[parCaller]] <- matrix(Mdl.betaIdx[[CompCaller]][[parCaller]],
                                                        nMCMC, ncolX.ij*nPar.ij,
                                                        dimnames = list(NULL, namesX.ij),
                                                        byrow = TRUE)
      MCMC.beta[[CompCaller]][[parCaller]] <- matrix(Mdl.beta[[CompCaller]][[parCaller]],
                                                     nMCMC, ncolX.ij*nPar.ij,
                                                     dimnames = list(NULL, namesX.ij),
                                                     byrow = TRUE)

      MCMC.par[[CompCaller]][[parCaller]] <- array(NA, c(nMCMC, nTraining, nPar.ij),
                                                   dimnames = list(NULL, NULL, namesPar.ij))

      ## The Metropolis-Hasting acceptance rate
      MCMC.AccProb[[CompCaller]][[parCaller]] <- matrix(NA, nMCMC, 1)
    }
  }


###----------------------------------------------------------------------------
### THE GIBBS SAMPLER (WITH METROPOLIS-HASTINGS)
###----------------------------------------------------------------------------

  cat("Posterior sampling using Metropolis-Hastings within Gibbs\n")


  ## Dry run to obtain staticcache for the initial values. Again this time all the
  ## parameters should be updated.

  staticCache <- logPost(MargisType = MargisType,
                         Mdl.Y = Mdl.Y.training,
                         Mdl.X = Mdl.X.training,
                         Mdl.beta = Mdl.beta,
                         Mdl.betaIdx = Mdl.betaIdx,
                         Mdl.parLink = Mdl.parLink,
                         varSelArgs = varSelArgs,
                         priArgs = priArgs,
                         parUpdate = MCMCUpdate,
                         MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

  ## Clear all warnings during initial value optimization. NOTE: not working
  ## warningsClear(envir = environment())

  ## The updating matrix
  UpdateMat <- parCplRepCaller(parUpdate = MCMCUpdate, parUpdateOrder = MCMCUpdateOrder)
  nInner <- nrow(UpdateMat)

  for(iUpdate in 1:(nInner*MCMC.nIter))
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
    MHOut <- MetropolisHastings(MargisType = MargisType,
                                propArgs = propArgs,
                                varSelArgs = varSelArgs,
                                priArgs = priArgs,
                                parUpdate = parUpdate,
                                Mdl.Y = Mdl.Y.training,
                                Mdl.X = Mdl.X.training,
                                Mdl.beta = Mdl.beta,
                                Mdl.betaIdx = Mdl.betaIdx,
                                Mdl.parLink = Mdl.parLink,
                                staticCache = staticCache,
                                MCMCUpdateStrategy = MCMCUpdateStrategy)

    if(MHOut$errorFlag == FALSE)
    { ## Update the MH results to the current parameter structure
      staticCache <- MHOut[["staticCache"]]
      Mdl.beta[[CompCaller]][[parCaller]] <- MHOut[["beta"]]
      Mdl.betaIdx[[CompCaller]][[parCaller]] <- MHOut[["betaIdx"]]

      AccProbCurr <- MHOut[["accept.prob"]]
    }
    else
    { ## Keep everything unchanged. Only set acceptance probability to zero if this
      ## iteration fails.
      AccProbCurr <- 0
    }

    ## Export the parameters in each cross-validation fold
    MCMC.beta[[CompCaller]][[parCaller]][iIter, ] <- Mdl.beta[[CompCaller]][[parCaller]]
    MCMC.betaIdx[[CompCaller]][[parCaller]][iIter, ] <- Mdl.betaIdx[[CompCaller]][[parCaller]]
    MCMC.par[[CompCaller]][[parCaller]][iIter, ,] <- staticCache[["Mdl.par"]][[CompCaller]][[parCaller]]

    MCMC.AccProb[[CompCaller]][[parCaller]][iIter,] <- AccProbCurr

    ## Save the marginal densities.  FIXME: This is not updated for every iteration.
    ## MCMC.density[["u"]][, , iIter] <- staticCache[["Mdl.u"]]
    ## MCMC.density[["d"]][, , iIter] <- staticCache[["Mdl.d"]]

    ## MCMC trajectory
    if(MCMC.track == TRUE && iInner == nInner)
    {
      CplMCMC.summary(iIter = iIter, MCMC.nIter = MCMC.nIter,
                      interval = 0.01, MCMC.burninProp = MCMC.burninProp,
                      OUT.MCMC = list(MCMC.beta = MCMC.beta,
                                      MCMC.betaIdx = MCMC.betaIdx,
                                      MCMC.par = MCMC.par,
                                      MCMC.AccProb = MCMC.AccProb,
                                      MCMCUpdate = MCMCUpdate,
                                      Starting.time = Starting.time))
    }

  }


  if(R_CPL_NPARALLEL>1)
  {
    stopCluster(cl4MCMC)
  }
  ## Fetch everything at current environment to a list
  ## list2env(out, envir = .GlobalEnv)
  ## GENERATING SHORT MODEL DESCRIPTION
  ModelDescription <- paste(c(MargisNM[-length(MargisNM)],"+",  MargisType, "+" ,
                              "nObs", nObs, "nCross", nCross, "MCMC.nIter", MCMC.nIter, "+",
                              format(Starting.time, "%Y%m%d@%H.%M")),
                            collapse = "")

  gc()

  OUT.MCMC = list(MCMC.beta = MCMC.beta,
                  MCMC.betaIdx = MCMC.betaIdx,
                  MCMC.par = MCMC.par,
                  MCMC.AccProb = MCMC.AccProb,
                  MCMCUpdate = MCMCUpdate,
                  Starting.time = Starting.time)

  out <- as.list(environment())
  return(out)
}
