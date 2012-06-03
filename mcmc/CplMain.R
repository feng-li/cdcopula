##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param setupfile "character".
##'        The path where the setup script is located.
##'
##' @return
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
  tauTabular <- kendalltauTabular(CplNM = CplNM, tol = 0.01)

  ## Indices for training and testing sample according to cross-validation
  nObs <- length(Mdl.Y[[1]])
  crossValidIdx <- set.crossvalid(nObs,crossValidArgs)
  nCrossFold <- length(crossValidIdx[["training"]])

  ## Split the data into folds for cross validation,  if no cross validation,
  ## only one fold used
  MdlTraining.X <- vector("list", nCrossFold)
  MdlTraining.Y <- vector("list", nCrossFold)

  MdlTesting.X <- vector("list", nCrossFold)
  MdlTesting.Y <- vector("list", nCrossFold)

  MdlMCMC.beta <- vector("list", nCrossFold)
  MdlMCMC.betaIdx <- vector("list", nCrossFold)
  MdlMCMC.par <- vector("list", nCrossFold)
  MdlMCMC.AccPbeta <- vector("list", nCrossFold)

###----------------------------------------------------------------------------
### Initialize the MCMC
### This section should be in a function so that it can be called easily for
### parallel computing and other purpose. The input should be the
### cross-validation training and testing indices.
### TODO:crossValidIdx should be reorganized
###----------------------------------------------------------------------------
  ## browser()

  iCross <- 1
  Training.Idx <- crossValidIdx[["training"]][[iCross]]
  Testing.Idx <- crossValidIdx[["testing"]][[iCross]]
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

  ## Allocate the storage for the final parameters
  MCMC.beta <- MdlDataStruc
  MCMC.betaIdx <- MdlDataStruc
  MCMC.par <- MdlDataStruc
  MCMC.AccPbeta <- MdlDataStruc
  for(i in names(MdlDataStruc))
    {
      for(j in names(MdlDataStruc[[i]]))
        {
          ncolX.ij <- ncol(MdlTraining.X[[i]][[j]])
          ## The MCMC storage
          MCMC.betaIdx[[i]][[j]] <- array(NA, c(ncolX.ij, nIter))
          MCMC.beta[[i]][[j]] <- array(NA, c(ncolX.ij, nIter))
          MCMC.par[[i]][[j]] <- array(NA, c(nObs, nIter))
          ## The Metropolis-Hasting acceptance rate
          MCMC.AccPbeta[[i]][[j]] <- array(NA, c(nIter, 1))
        }
    }

  ## Assign the initial values
  initParOut <- initPar(varSelArgs = varSelArgs,
                        betaInit = betaInit,
                        Mdl.X = MdlTraining.X)
  Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
  Mdl.beta <- initParOut[["Mdl.beta"]]

  ## Switch all the updating indicators ON
  parUpdate <- MCMCUpdate

  ## Initialize "staticArgs"
  ## Dry run the logPost once to obtain it.
  ## FIXME: Using optimization routine to get good initial values
  u <- matrix(NA, dim(MdlTraining.Y[[1]])[1], length(MdlTraining.Y),
              dimnames = list(NULL, names(MdlTraining.Y)))
  staticArgs <- list(Mdl.logPri =  MdlDataStruc,
                     Mdl.par = MdlDataStruc,
                     Mdl.u = u,
                     Mdl.d = u,
                     tauTabular = tauTabular)
  staticArgs <- logPost(CplNM = CplNM,
                        Mdl.Y = MdlTraining.Y,
                        Mdl.X = MdlTraining.X,
                        Mdl.beta = Mdl.beta,
                        Mdl.betaIdx = Mdl.betaIdx,
                        Mdl.parLink = Mdl.parLink,
                        varSelArgs = varSelArgs,
                        MargisTypes = MargisTypes,
                        priArgs = priArgs,
                        parUpdate = parUpdate,
                        staticArgs = staticArgs)[["staticArgs"]]

###----------------------------------------------------------------------------
### STABILIZE THE INITIAL VALUES VIA NEWTON ITERATIONS
###----------------------------------------------------------------------------
  ## We use a simple two stage optimizations for the marginal and copula model.
  ##  optim()
  browser()

###----------------------------------------------------------------------------
### THE METROPOLIS-HASTINGS WITHIN GIBBS
###----------------------------------------------------------------------------
  ## Switch all the updating indicators OFF
  parUpdate <- rapply(parUpdate, function(x) FALSE, how = "replace")
  CompNM <- names(MdlDataStruc)

  ## The MCMC loops
  for(iIter in 1:nIter)
    {
      ## The Gibbs loop
      for(CompCurr in CompNM)
        {
          parNM <- names(MdlDataStruc[[CompCurr]])
          for(parCurr in parNM)
            {
              if (MCMCUpdate[[CompCurr]][[parCurr]] == TRUE)
                {
                  ## Switch current updating parameter indicator on
                  parUpdate[[CompCurr]][[parCurr]] <- TRUE

                  ## Call the mail Metropolis--Hastings
                  algmArgs <- propArgs[[CompCurr]][[parCurr]][["algorithm"]]

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
                    }
                  else
                    {
                      stop("Unknown proposal algorithm!")
                    }
                }
              ## Switch current updating parameter indicator off
              parUpdate[[CompCurr]][[parCurr]] <- FALSE
            }
        }
    }

  ## }
###----------------------------------------------------------------------------
### RUN THE MCMC
### Cross-validation is independent for all folds. In order to parallel it via
### mcmlapply() function in the parallel package in R, we write the sequential
### code with malapply().
### ----------------------------------------------------------------------------


  parallel <- FALSE
  if(parallel == TRUE)
    {
      ## Use mcmlapply function
    }
  else
    {
      ## Sequential loops over folds.
      ## mapply(FUN = MCMCFun)
    }

###----------------------------------------------------------------------------
### POSTERIOR INFERENCE, PREDICTION ETC
###----------------------------------------------------------------------------
      ## The final update the parameters in each fold
      MdlMCMC.beta[[iCross]] <- Mdl.betaIdx.iCross
      MdlMCMC.betaIdx[[iCross]] <- Mdl.beta.iCross
      MdlMCMC.par[[iCross]] <- Mdl.par.iCross
      MdlMCMC.AccPbeta[[iCross]] <- Mdl.AccPbeta.iCross


  ## Fetch everything at current environment to a list
  ## Use list2env() to extract the entries
  out <- as.list(environment())
  ## return(out)
}
