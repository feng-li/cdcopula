##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param setupfile "character".
##'        The path where the setup script is located.
##' 
##' @return
##' 
##' @references
##' 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' 
##' @note Created: Thu Feb 02 13:33:06 CET 2012;
##'       Current: Wed Feb 08 15:35:52 CET 2012.
CplMHMain <- function(setupfile)
{
###----------------------------------------------------------------------------
### Load user setup file
###----------------------------------------------------------------------------  
  source(setupfile, local = TRUE)
  
###----------------------------------------------------------------------------
### INITIALIZE THE STORAGE AND DATA STRUCTURE 
###----------------------------------------------------------------------------
  
  ## Storage for the beta coefficients
  MdlData.beta <- MdlDataStruc
  MdlDataAccP.beta <- MdlDataStruc

  ## Storage for the variable selection index vector
  MdlData.betaIdx <- MdlDataStruc
  MdlDataAccP.betaIdx <- MdlDataStruc

  ## Storage for the parameters
  MdlData.par <- MdlDataStruc  

  ## Allocate the storage
  for(i in 1:length(MdlDataStruc))
    {
      for(j in 1:length(MdlDataStruc[[i]]))
        {
          ncolX.ij <- ncol(Mdl.X[[i]][[j]])

          ## The MCMC storage
          MdlData.betaIdx[[i]][[j]] <- array(NA, c(ncolX.ij, nIter, nCrossFold))
          MdlData.beta[[i]][[j]] <- array(NA, c(ncolX.ij, nIter, nCrossFold))
          MdlData.par[[i]][[j]] <- array(NA, c(nObs, nIter, nCrossFold))

          ## The Metropolis-Hasting acceptance rate
          MdlDataAccP.beta[[i]][[j]] <- array(NA, c(nIter, nCrossFold))
        } 
    }

###----------------------------------------------------------------------------
### Initialize the MCMC  
###----------------------------------------------------------------------------

  ## Assign the initial values FIXME:
  initParOut <- initPar(varSelArgs, betaInit, Mdl.X)
  Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
  Mdl.beta <- initParOut[["Mdl.beta"]]

  ## Switch all the updating indicators ON
  parUpdate <- MCMCUpdate
  
  ## Initialize the prior
  priCurr <- MdlDataStruc
  priCurr <- logPriors(Mdl.X = Mdl.X,
                       Mdl.parLink = Mdl.parLink,
                       Mdl.beta = Mdl.beta,
                       Mdl.betaIdx = Mdl.betaIdx,
                       varSelArgs = varSelArgs,
                       priArgs = priArgs,
                       priCurr = priCurr,
                       parUpdate = MCMCUpdate)

  ## Switch all the updating indicators OFF
  parUpdate <- rapply(parUpdate, function(x) FALSE, how = "replace")


  ## Update static argument
  staticArgs <- list()
  staticArgs[["priCurr"]] <- priCurr
  staticArgs[["u"]] <- NA
  staticArgs[["Mdl.par"]] <- Mdl.par
  staticArgs[["tauTabular"]] <- tauTabular
  
###----------------------------------------------------------------------------
### RUN THE MCMC 
###----------------------------------------------------------------------------
  for(iCross in nCrossFold) ## TODO: Parallel computing 
    {
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
                      if(tolower(algArgs[["type"]]) == "kstepsnewtonmove")
                        {
                          ## Cross validation subsets 
                          subset <- crossValidIdx[["training"]][[iCross]]
                          ## staticArgs <- list(u = u, Mdl.par = Mdl.par) 

                          out <- MHPropWithKStepNewton(
                                                       CplNM = CplNM,
                                                       propArgs = propArgs,
                                                       varSelArgs = varSelArgs,
                                                       priArgs = priArgs,
                                                       betaIdxProp = betaIdxProp,
                                                       parUpdate = parUpdate,
                                                       Mdl.Y = Mdl.Y,
                                                       Mdl.X = Mdl.X,
                                                       Mdl.beta = Mdl.beta,
                                                       Mdl.betaIdx = Mdl.betaIdx,
                                                       MargisTypes = MargisTypes, 
                                                       staticArgs = staticArgs
                                                       )   

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
    }

  ## Fetch everything at current environment to a list
  ## Use list2env() to extract the entries
  out <- as.list(environment())
  ## return(out)
}
