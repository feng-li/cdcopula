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
CplMain <- function(setupfile)
{
###----------------------------------------------------------------------------
### Load user setup file
###----------------------------------------------------------------------------  
  source(setupfile, local = TRUE)
  
###----------------------------------------------------------------------------
### INITIALIZE THE STORAGE AND DATA STRUCTURE 
###----------------------------------------------------------------------------

  ## Storage for the final parameters
  MdlData.beta <- MdlDataStruc
  MdlData.betaIdx <- MdlDataStruc
  MdlData.par <- MdlDataStruc  
  MdlDataAccP.beta <- MdlDataStruc

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

  ## Assign the initial values FIXME: nObs is not right 
  initParOut <- initPar(varSelArgs, betaInit, Mdl.X)
  Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
  Mdl.beta <- initParOut[["Mdl.beta"]]
  
  ## Switch all the updating indicators ON
  parUpdate <- MCMCUpdate
    
  ## Allocate and initialize the static argument
  u <- matrix(NA, nObs, length(Mdl.Y), dimnames = list(NULL, names(Mdl.Y)))
  
  staticArgs <- list(Mdl.logPri =  MdlDataStruc,
                     Mdl.par = MdlDataStruc,
                     Mdl.u = u,
                     Mdl.d = u, 
                     tauTabular = tauTabular)

  ## FIXME: Using optimization routine to get good initial values
  
  logPostOut <- logPost(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx, Mdl.parLink,
                  varSelArgs, MargisTypes, priArgs, parUpdate, staticArgs)  

  staticArgs <- logPostOut[["staticArgs"]]  # Dry run to initialize staticArgs.
  
  ## Switch all the updating indicators OFF
  parUpdate <- rapply(parUpdate, function(x) FALSE, how = "replace")
  
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
                                                       parUpdate = parUpdate,
                                                       Mdl.Y = Mdl.Y,
                                                       Mdl.X = Mdl.X,
                                                       Mdl.beta = Mdl.beta,
                                                       Mdl.betaIdx = Mdl.betaIdx,
                                                       MargisTypes = MargisTypes, 
                                                       Mdl.parLink = Mdl.parLink, 
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
