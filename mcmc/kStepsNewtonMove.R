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
##' @param staticArgs 
##' @param param.cur "matrix".
##'         The initial values for the Newton update.
##' @return "list". See bellow.
##' \item   {gradObs.cur}
##'         {"matrix". The gradient}
##' \item   {invHessObs.cur}
##'         {"matrix". The inverse Hessian matrix.}
##' \item   {param.cur}
##'         {"matrix". The updated paramters after K step Newton integrations.}
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Wed Sep 29 17:18:22 CEST 2010;
##'       Current: Thu Feb 09 14:59:27 CET 2012.
kStepsNewtonMove <- function(propArgs, varSelArgs, priArgs, betaIdxProp,
                             parUpdate, CplNM, Mdl.Y, Mdl.X, Mdl.beta,
                             Mdl.betaIdx, Mdl.parLink, MargisTypes, staticArgs)    
{

  ## The updating component parameter chain
  cp <- parCaller(parUpdate)
  CompCurr <- cp[1]
  parCurr <- cp[2]
  
  ## The current parameters 
  X <- Mdl.X[[CompCurr]][[parCurr]]
 
  betaCurr <- Mdl.beta[[CompCurr]][[parCurr]]
  betaIdxCurr <- Mdl.betaIdx[[CompCurr]][[parCurr]]

  ## Initialize the Newton move 
  Mdl.betaIdx[[CompCurr]][[parCurr]] <- betaIdxProp
  
  ## Finite Newton move. K steps to approach the mode plus one more step to
  ## update the gradient and Hessian at i:th step.
  kSteps <- propArgs[[CompCurr]][[parCurr]][["algorithm"]][["ksteps"]]
  hessMethod <- propArgs[[CompCurr]][[parCurr]][["algorithm"]][["hess"]]

  staticArgsOld <- staticArgs

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
      logLikGrad.prop <- gradHess.prop[["logLikGrad"]] # n-by-pp
      logLikHess.prop <- hessApprox(logLikGrad.prop, hessMethod)

      ## The alternative variable selection indices for the covariates
      betaIdxProp2 <- which(betaIdxProp == 1)
      betaIdxCurr2 <- which(betaIdxCurr == 1)

      ## The selected covariates
      X.prop <- X[ , betaIdxProp2, drop = FALSE] # n-by-pp
      X.curr <- X[ , betaIdxCurr2, drop = FALSE] # n-by-pc

      browser()
      
      ## Gradient Hessian for the prior *including non selected covariates*
      ## NOTE: The Hessian matrix of the prior is also approximated, we should
      ## use the explicit Hessian whenever possible.
      
      logPriGrad.prop <- gradHess.prop[["logPriGrad"]] # pp-by-1
      logPriHess.prop <- diag(hessApprox(logPriGrad.prop, hessMethod),
                            length(betaIdxProp))

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
      HessObsInv.pp <- solve(HessObs.pp)
      
      if((iStep <= kSteps)) ## The general Newton's Update
        {
          ## update the proposed parameters
          param <- HessObsInv.pp%*%(HessObs.pc%*%param -
                                        gradObs.prop)
          
          ## Update the parameter with current updated results.
          Mdl.beta[[CompCurr]][[parCurr]] <- param
          betaIdxCurr <- betaIdxProp
        }
      else if(iStep == (kSteps+1)) # (k+1):th step.  Make a output 
        {
          out <- list(gradObs = gradObs.prop,
                      HessObs = HessObs.prop,
                      HessObsInv = HessObsInv,
                      param = param) 
        }
       
    }
  return(out) 
}
