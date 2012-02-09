##' This is the prior settings for the Covariates dependent copula model
##'
##' The detailed documentation is in the main setting file for each parameter.
##' @title Copula model prior settings.
##' @param Mdl.X 
##' @param Mdl.parLink 
##' @param MdlCurr.beta 
##' @param MdlCurr.betaIdx 
##' @param varSelArgs 
##' @param priArgs "list"
##' @param priCurr "list"
##' @param parCurrUpdate 
##' @return "list" synced 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Dec 15 10:45:56 CET 2011;
##'       Current: Thu Dec 15 10:46:05 CET 2011.
logPriors <- function(Mdl.X, Mdl.parLink, Mdl.beta, Mdl.betaIdx,
                      varSelArgs, priArgs, priCurr, parUpdate)
{
  ## Loop over all updated parameter candidates
  
  CompNM <- names(priArgs)
  for(i in CompNM)
    {        
###----------------------------------------------------------------------------
### Only update priors for parameters that need to update.
###----------------------------------------------------------------------------
      parUpdateIdx <- which(parUpdate[[i]] == TRUE) 
      for(j in parUpdateIdx)
        {

          ## Initial the storage structure for current log prior
          outCurr <-  priArgs[[i]][[j]]
          
###----------------------------------------------------------------------------
### Prior for variable selection indicators
###----------------------------------------------------------------------------
          priArgsCurr <- priArgs[[i]][[j]][["indicators"]]
          betaIdxCurr <- Mdl.betaIdx[[i]][[j]]

          if(tolower(priArgsCurr[["type"]]) == "bern") # Bernoulli prior 
            {
              ## The probability is set for each variable involved in the
              ## variable selection procedure via "valSelArgs" 
              ## The probability is recycled if necessary
              prob <- priArgsCurr[["prob"]]
              
              candIdx <- varSelArgs[[i]][[j]][["cand"]] # Variable section candidates
              candLen <- length(candIdx)
              probVec <- matrix(prob, candLen, 1)

              idx1 <- (betaIdxCurr[candIdx] == 1)
              idx0 <- (betaIdxCurr[candIdx] == 0)
              logDens <- sum(log(probVec[idx1])) + sum(log(1-probVec[idx0]))
              outCurr[["indicators"]] <- logDens
            }

###----------------------------------------------------------------------------
### Prior for coefficients 
###----------------------------------------------------------------------------

### intercept as special case
          priArgsCurr <- priArgs[[i]][[j]][["beta"]][["intercept"]]
          xCurr <- Mdl.beta[[i]][[j]][1] # the intercept
          linkCurr <- Mdl.parLink[[i]][[j]]

          if(tolower(priArgsCurr[["type"]]) == "custom")
            {
              ## Call the any2any() function                
              densOutput <- any2any(densArgs = priArgsCurr, linkType = linkCurr)

              mean <- densOutput$mean
              variance <- densOutput$variance
              shrinkage <- priArgsCurr[["output"]][["shrinkage"]]

              logDens <- dnorm(xCurr, mean = mean,
                               sd = sqrt(variance*shrinkage), log = TRUE)
              outCurr[["beta"]][["intercept"]] <- logDens
            }
          
          
### coefficients conditional on variable selection indicators
          priArgsCurr <- priArgs[[i]][[j]][["beta"]][["slopes"]]
          xCurr <- Mdl.beta[[i]][[j]][-1] # Slopes(taking away intercept)
          betaIdxNoIntCurr <- Mdl.betaIdx[[i]][[j]][-1] # Variable section
                                        # indicator without intercept

          if(tolower(priArgsCurr[["type"]]) == "cond-mvnorm")
            {
              ## Normal distribution condition normal The full beta
              ## vector is assumed as normal. Since variable selection is
              ## included in the MCMC, The final proposed beta are those
              ## non zeros. We need to using the conditional normal
              ## density See Mardia p. 63

              ## Subtract the prior information for the full beta
              mean <- priArgsCurr[["mean"]] # mean of density
              covariance <- priArgsCurr[["covariance"]] # Covariates
              shrinkage <- priArgsCurr[["shrinkage"]] # Shrinkage

              ## Split the beta vector by zero and nonzero.
              Idx1 <- which(betaIdxNoIntCurr == 1)
              Idx0 <- which(betaIdxNoIntCurr == 0)
              
              Idx0Len <- length(Idx0)
              betaLen <- length(betaIdxNoIntCurr)                

              ## The mean vector (recycled if necessary)
              meanVec <- matrix(mean, betaLen, 1) 
              
              ## The covariance matrix for the whole beta vector
              if(tolower(covariance) == "g-prior")
                {
                  X <- Mdl.X[[i]][[j]][, -1, drop = FALSE]
                  coVar <- qr.solve(crossprod(X)) 
                }
              else if(tolower(covariance) == "identity")
                {
                  coVar <- diag(length(Idx1))
                }


              ## Consider three situations:
              if(Idx0Len == 0 || Idx0Len == betaLen)
                {
                  ## 1. all are selected 2. non are selected:
                  ## Switch to unconditional prior.
                  outCurr[["beta"]][["slopes"]] <-
                    dmvnorm(t(xCurr), meanVec, coVar*shrinkage, log = TRUE)
                }
              else if(Idx0Len > 0 && Idx0Len < betaLen)
                {
                  ## 2. some are selected
                  ## The conditional prior
                  
                  A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
                  condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
                  condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]
                  
                  outCurr[["beta"]][["slopes"]] <-
                    dmvnorm(matrix(xCurr[Idx1], 1),
                            condMean, condCovar*shrinkage, log = TRUE)
                }
              else
                {
                  stop("Debug me: Unknown situation for conditional priors.")
                }
            }
          
###----------------------------------------------------------------------------
### Update the output for prior
###----------------------------------------------------------------------------            
          priCurr[[i]][[j]] <- outCurr
          
        }
    }
  return(priCurr)
}
