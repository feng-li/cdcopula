##' This is the Newton move with dimension changes for the copula model. 
##'
##' @title Generalized Newton method for the copula model
##' @param propArgs 
##' @param varSelArgs 
##' @param priArgs 
##' @param betaIdxProp 
##' @param parUpdate 
##' @param CplNM 
##' @param Mdl.Y 
##' @param Mdl.X 
##' @param Mdl.beta 
##' @param Mdl.betaIdx 
##' @param Mdl.parLink 
##' @param MargisTypes 
##' @param staticArgs 
##' @param param.cur "matrix".
##'         The initial values for the Newton update.
##' @return "list". See bellow.
##' \item   {gradObs}
##'         {"matrix". The gradient}
##' \item   {HessObsInv}
##'         {"matrix". The inverse Hessian matrix.}
##' \item   {param}
##'         {"matrix". The updated parameters after K step Newton integrations.}
##' @references Li Villani Kohn 2010
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Wed Sep 29 17:18:22 CEST 2010;
##'       Current: Mon Mar 05 10:33:29 CET 2012.
GNewtonMove <- function(propArgs, varSelArgs, priArgs, betaIdxProp,
                             parUpdate, CplNM, Mdl.Y, Mdl.X, Mdl.beta,
                             Mdl.betaIdx, Mdl.parLink, MargisTypes, staticArgs)    
{
  ## The updating component parameter chain
  cp <- parCaller(parUpdate)
  CompCurr <- cp[1]
  parCurr <- cp[2]
  
  ## The current parameters 
  X <- Mdl.X[[CompCurr]][[parCurr]]
  betaCurr <- Mdl.beta[[CompCurr]][[parCurr]] # p-by-1
  betaIdxCurr <- Mdl.betaIdx[[CompCurr]][[parCurr]] # p-by-1
  
  ## Finite Newton move. K steps to approach the mode plus one more step to
  ## update the gradient and Hessian at i:th step.
  kSteps <- propArgs[[CompCurr]][[parCurr]][["algorithm"]][["ksteps"]]
  hessMethod <- propArgs[[CompCurr]][[parCurr]][["algorithm"]][["hess"]]

  ##-----------------------Initialize the Newton move---------------------------

  ## Initialize the Newton move 
  Mdl.betaIdx[[CompCurr]][[parCurr]] <- betaIdxProp

  staticArgsOld <- staticArgs
  
  ##--------------------The k-step Generalized Newton Move---------------------- 
  for(iStep in 1:(kSteps+1))
    {
      ## Obtain the gradient and Hessian information
      gradHess.prop <- logPostGradHess(CplNM = CplNM,
                                       Mdl.Y = Mdl.Y,
                                       Mdl.X = Mdl.X,
                                       Mdl.parLink = Mdl.parLink,
                                       Mdl.beta = Mdl.beta,
                                       MargisTypes = MargisTypes, 
                                       Mdl.betaIdx = Mdl.betaIdx,
                                       parUpdate = parUpdate,
                                       priArgs = priArgs, 
                                       staticArgs = staticArgs)

      ## Gradient and Hessian for the likelihood
      logLikGrad.prop <- gradHess.prop[["logLikGradObs"]] # n-by-pp
      logLikHess.prop <- hessApprox(logLikGrad.prop, hessMethod)

      ## The alternative variable selection indices for the covariates
      betaIdxProp2 <- which(betaIdxProp == 1)
      betaIdxCurr2 <- which(betaIdxCurr == 1)

      ## The selected covariates
      X.prop <- X[ , betaIdxProp2, drop = FALSE] # n-by-pp
      X.curr <- X[ , betaIdxCurr2, drop = FALSE] # n-by-pc
    
      ## Gradient Hessian for the prior *including non selected covariates*
      ## NOTE: The Hessian matrix of the prior is also approximated, we should
      ## use the explicit Hessian whenever possible.
      
      logPriGrad.prop <- gradHess.prop[["logPriGradHessObs"]][["gradObs"]] # pp-by-1
      logPriHess.prop <- gradHess.prop[["logPriGradHessObs"]][["HessObs"]] # p-by-p

      ## The gradient and Hessian subsets due to variable selection
      logPriGrad.pp <- logPriGrad.prop[betaIdxProp2] # scaler
      logPriHess.pp <- logPriHess.prop[betaIdxProp2, betaIdxProp2,
                                       drop = FALSE] # pp-by-pp
      logPriHess.pc <- logPriHess.prop[betaIdxProp2, betaIdxCurr2,
                                       drop = FALSE] # pp-by-pc

      ## The gradient and Hessian in the general Newton's update
      gradObs.prop <- matrix(rowSums(Md(t(X.prop), logLikGrad.prop))+
                             logPriGrad.pp) # pp-by-1
      HessObs.pp <- tMdN(X.prop, logLikHess.prop, X.prop)+
        logPriHess.pp # pp-by-pp
      HessObs.pc <- tMdN(X.prop, logLikHess.prop, X.curr)+
        logPriHess.pc # pp-by-pc
      HessObsInv.pp <- solve(HessObs.pp) # pp-by-pp

      ## The general Newton's Update
      if((iStep <= kSteps)) 
        {
          ## update the proposed parameters via the general Newton formula
          param <- HessObsInv.pp%*%(HessObs.pc%*%param - gradObs.prop)
          
          ## Update the parameter with current updated results.
          betaIdxCurr <- betaIdxProp 
          paramTmp <- matrix(0, length(betaIdxCurr), 1) # Temporary with all zeros. 
          paramTmp[betaIdxCurr] <- param
          Mdl.beta[[CompCurr]][[parCurr]] <- paramTmp
          staticArgs <- gradHess.prop[["staticArgs"]]
        }
      else if(iStep == (kSteps+1)) # (k+1):th step.  Make a output 
        {
          out <- list(gradObs = gradObs.prop,
                      HessObs = HessObs.pp,
                      HessObsInv = HessObsInv.pp,
                      param = param, 
                      staticArgs = staticArgs) 
        }
       
    }
  return(out) 
}
