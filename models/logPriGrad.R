##' A collection of gradient for the copula model
##'
##' @name logPriGrad
##' @title Gradient for priors
##' 
##' @param Mdl.X 
##' @param Mdl.parLink 
##' @param Mdl.beta 
##' @param Mdl.betaIdx 
##' @param varSelArgs 
##' @param priArgs 
##' @param priCurr 
##' @param chainCaller 
##' @return "list". The gradient and hessian matrix, see below.
##' \item   {gradObsPri}
##'         {"matrix". One-column.}
##' 
##' \item   {hessObsPri}
##'         {"matrix". A squre matrix. Dimension same as prior_type$Sigma0.}
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Mar 30 09:33:23 CEST 2010;
##'       Current: Mon Feb 27 15:22:27 CET 2012.
logPriGrad <- function(Mdl.X, Mdl.beta, Mdl.betaIdx, Mdl.parLink,
                       varSelArgs, priArgs, chainCaller)
  {
    ## Only update priors for parameters that need to update.
    ## Initial the storage structure for current log prior
    CompCurr <- chainCaller[[1]]
    parCurr <- chainCaller[[2]]



    
    ## Gradient for the intercept as a special case
    priArgsCurr <- priArgs[[CompCurr]][[parCurr]][["beta"]][["intercept"]]
    xCurr <- Mdl.beta[[CompCurr]][[parCurr]][1] # the intercept
    linkCurr <- Mdl.parLink[[CompCurr]][[parCurr]]
    
    if(tolower(priArgsCurr[["type"]]) == "custom")
      {
        ## Call the any2any() function                
        densOutput <- any2any(densArgs = priArgsCurr, linkType = linkCurr)
        mean <- densOutput$mean
        variance <- densOutput$variance
        shrinkage <- priArgsCurr[["output"]][["shrinkage"]]
        
        GradInitOut <- - 1/(variance*shrinkage)*(xCurr-mean)
      }
                        
    ## Gradient for the coefficients conditional on variable selection indicators
    priArgsCurr <- priArgs[[CompCurr]][[parCurr]][["beta"]][["slopes"]]
    xCurr <- Mdl.beta[[CompCurr]][[parCurr]][-1] # Slopes(taking away intercept)
    betaIdxNoIntCurr <- Mdl.betaIdx[[CompCurr]][[parCurr]][-1] # Variable section
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
            X <- Mdl.X[[CompCurr]][[parCurr]][, -1, drop = FALSE]
            coVar <- qr.solve(crossprod(X)) 
          }
        else if(tolower(covariance) == "identity")
          {
            coVar <- diag(length(Idx1))
          }

        ## Consider three situations for gradient:
        if(Idx0Len == 0 || Idx0Len == betaLen)
          {
            ## 1. all are selected,  2. non are selected:
            ## FIXME: The situation 2. need to reconsider.
            ## Switch to unconditional prior.
            outCurr <- - solve(condCovar*shrinkage)*(xCurr-condMean)
          }
        else if(Idx0Len > 0 && Idx0Len < betaLen)
          {
            ## 2. some are selected
            ## The conditional prior
            A <- coVar[Idx1, Idx0, drop = FALSE]%*%
              solve(coVar[Idx0, Idx0, drop = FALSE])
            condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
            condCovar <- coVar[Idx1, Idx1, drop = FALSE] -
              A%*%coVar[Idx0, Idx1, drop = FALSE]
            
            GradSlopOut <- - solve(condCovar*shrinkage)*(xCurr-condMean) 
          }
        else
          {
            stop("Debug me: Unknown situation for conditional priors.")
          }
        ## The final gradient output. 
        out <- matrix(c(GradInitOut, GradSlopOut))
      }
    return(out)
  }
