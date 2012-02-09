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
  browser()

###----------------------------------------------------------------------------
### INITIALIZE THE STORAGE AND DATA STRUCTURE 
###----------------------------------------------------------------------------
  
  ## Storage for the beta coefficients
  MdlData.beta <- MdlDataStruc
  MdlCurr.beta <- MdlDataStruc
  MdlDataAccP.beta <- MdlDataStruc

  ## Storage for the variable selection index vector
  MdlData.betaIdx <- MdlDataStruc
  MdlCurr.betaIdx <- MdlDataStruc
  MdlDataAccP.betaIdx <- MdlDataStruc

  ## Storage for the parameters
  MdlData.par <- MdlDataStruc  
  MdlCurr.par <- MdlDataStruc  
  MdlDataAccP.par <- MdlDataStruc

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
          MdlDataAccP.betaIdx[[i]][[j]] <- array(NA, c(nIter, nCrossFold))
          MdlDataAccP.beta[[i]][[j]] <- array(NA, c(nIter, nCrossFold))
          MdlDataAccP.par[[i]][[j]] <- array(NA, c(nIter, nCrossFold))
        } 
    }

###----------------------------------------------------------------------------
### Initialize the MCMC  
###----------------------------------------------------------------------------

  ## Assign the initial values FIXME:
  initParOut <- initPar(varSelArgs, betaInit, Mdl.X)
  Mdl.betaIdx <- initParOut[["MdlCurr.betaIdx"]]
  Mdl.beta <- initParOut[["MdlCurr.beta"]]

  ## Switch all the updating indicators ON
  parCurrUpdate <- MdlDataStruc
  parCurrUpdate <- rapply(parCurrUpdate, function(x) TRUE, how = "replace")

  ## Initialize the prior
  priCurr <- MdlDataStruc
  priCurr <- logPriors(Mdl.X, Mdl.parLink, Mdl.beta, Mdl.betaIdx,
                       varSelArgs, priArgs, priCurr, parCurrUpdate)

  ## Switch all the updating indicators OFF
  parCurrUpdate <- rapply(parCurrUpdate, function(x) FALSE, how = "replace")

###----------------------------------------------------------------------------
### RUN THE MCMC 
###----------------------------------------------------------------------------

  for(iCross in nCrossFold) ## TODO: Parallel computing 
    {
      CompNM <- names(MdlDataStruc)
      ## The MCMC loops
      for(iterCurr in 1:nIter)
        { 
          ## The Gibbs loop
          for(CompCurr in CompNM)
            {
              parNM <- names(MdlDataStruc[[CompCurr]])
              for(parCurr in parNM)
                {
                  if (parCurrUpdate[[CompCurr]][[parCurr]] == TRUE)
                    {
                      ## Switch current updating parameter indicator on
                      parCurrUpdate[[CompCurr]][[parCurr]] <- TRUE
                      
                      ## This is the mail MH file,                        
                      if(tolower(propMethod) == "kstepnewton")
                        {
                          out <- MHPropWithKStepNewton()
                        }
                      else
                        {
                          stop("Wrong argument for propDensity!")
                        }                      
                    }
                  ## Switch current updating parameter indicator off
                  parCurrUpdate[[CompCurr]][[parCurr]] <- FALSE
                }
            }
        }
    }
  out <- list
  return(NA)
}
