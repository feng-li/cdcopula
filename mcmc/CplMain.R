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
  tauTabular <- kendalltauTabular(CplNM = CplNM, tol = 0.005)

  ## Indices for training and testing sample according to cross-validation
  nObs <- length(Mdl.Y[[1]])
  crossValidIdx <- set.crossvalid(nObs,crossValidArgs)
  nCrossFold <- length(crossValidIdx[["training"]])

  ## Split the data into folds for cross validation,  if no cross validation,
  ## only one fold used
  subsetFun <- function(x, idx)x[idx, , drop = FALSE]
  MdlTraining.X <- vector("list", nCrossFold)
  MdlTraining.Y <- vector("list", nCrossFold)
  MdlTesting.X <- vector("list", nCrossFold)
  MdlTesting.Y <- vector("list", nCrossFold)

  MdlMCMC.beta <- vector("list", nCrossFold)
  MdlMCMC.betaIdx <- vector("list", nCrossFold)
  MdlMCMC.par <- vector("list", nCrossFold)
  MdlMCMC.AccPbeta <- vector("list", nCrossFold)


  for(iCross in 1:nCrossFold)
    {
      Training.Idx.iCross <- crossValidIdx[["training"]][[iCross]]
      Testing.Idx.iCross <- crossValidIdx[["testing"]][[iCross]]
      ##Extract the training and testing data according to cross-validation
      MdlTraining.X[[iCross]] <- rapply(object=Mdl.X,
                                        f = subsetFun,
                                        idx = Training.Idx.iCross,
                                        how = "replace")
      MdlTraining.Y[[iCross]] <- rapply(object=Mdl.X,
                                        f = subsetFun,
                                        idx = Training.Idx.iCross,
                                        how = "replace")
      MdlTesting.X[[iCross]] <- rapply(object=Mdl.X,
                                       f = subsetFun,
                                       idx = Testing.Idx.iCross,
                                       how = "replace")
      MdlTesting.Y[[iCross]] <- rapply(object=Mdl.X,
                                       f = subsetFun,
                                       idx = Testing.Idx.iCross,
                                       how = "replace")
    }


###----------------------------------------------------------------------------
### Initialize the MCMC
### This section should be in a function so that it can be called easily for
### parallel computing and other purpose.
###----------------------------------------------------------------------------

  ## Allocate the storage for the final parameters and split the covariates
  ## according to the variable selection settings
  Mdl.beta.iCross <- MdlDataStruc
  Mdl.betaIdx.iCross <- MdlDataStruc
  Mdl.par.iCross <- MdlDataStruc
  Mdl.AccPbeta .iCross <- MdlDataStruc
  for(i in names(MdlDataStruc))
    {
      for(j in names(MdlDataStruc[[i]]))
        {
          ncolX.ij <- ncol(Mdl.X[[i]][[j]])
          ## The MCMC storage
          Mdl.betaIdx.iCross[[i]][[j]] <- array(NA, c(ncolX.ij, nIter))
          Mdl.beta.iCross[[i]][[j]] <- array(NA, c(ncolX.ij, nIter))
          Mdl.par.iCross[[i]][[j]] <- array(NA, c(nObs, nIter))
          ## The Metropolis-Hasting acceptance rate
          Mdl.AccPbeta.iCross[[i]][[j]] <- array(NA, c(nIter, 1))
        }
    }


  ## Assign the initial values
  initParOut <- initPar(varSelArgs, betaInit, Mdl.X.iCross)
  Mdl.betaIdx <- initParOut[["Mdl.betaIdx.iCross"]]
  Mdl.beta <- initParOut[["Mdl.beta.iCross"]]

  ## Switch all the updating indicators ON
  parUpdate <- MCMCUpdate

  ## Allocate and initialize the static argument
  u <- matrix(NA, nObs, length(Mdl.Y), dimnames = list(NULL, names(Mdl.Y)))
  staticArgs <- list(Mdl.logPri =  MdlDataStruc,
                     Mdl.par = MdlDataStruc,
                     Mdl.u = u,
                     Mdl.d = u,
                     tauTabular = tauTabular)

  ## Dry run to initialize "staticArgs"
  ## FIXME: Using optimization routine to get good initial values
  logPostOut <- logPost(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx, Mdl.parLink,
                        varSelArgs, MargisTypes, priArgs, parUpdate, staticArgs)
  staticArgs <- logPostOut[["staticArgs"]]  .

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
                  algArgs <- propArgs[[CompCurr]][[parCurr]][["algorithm"]]

                  if(tolower(algArgs[["type"]]) == "gnewtonmove")
                    {
                      ## Cross validation subsets
                      subset <- crossValidIdx[["training"]][[iCross]]
                      ## staticArgs <- list(u = u, Mdl.par = Mdl.par)

                      MHOut <- MHWithGNewton(
                                 CplNM = CplNM,
                                 propArgs = propArgs,
                                 varSelArgs = varSelArgs,
                                 priArgs = priArgs,
                                 parUpdate = parUpdate,
                                 Mdl.Y = Mdl.Y,
                                 Mdl.X = Mdl.X,
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

  ## The final update the parameters in each fold
  MdlMCMC.beta[[iCross]] <- Mdl.betaIdx.iCross
  MdlMCMC.betaIdx[[iCross]] <- Mdl.beta.iCross
  MdlMCMC.par[[iCross]] <- Mdl.par.iCross
  MdlMCMC.AccPbeta[[iCross]] <- Mdl.AccPbeta.iCross

###----------------------------------------------------------------------------
### RUN THE MCMC
### Cross-validation is independent for all folds. In order to parallel it via
### mcmlapply() function in the parallel package in R, we write the sequential
### code with malapply().
### ----------------------------------------------------------------------------
  browser()

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


  ## Fetch everything at current environment to a list
  ## Use list2env() to extract the entries
  out <- as.list(environment())
  ## return(out)
}
