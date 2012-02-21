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
                             Mdl.betaIdx, staticArgs)    
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

      browser()
      
      ## Obtain the gradient and Hessian information
      gradHess.prop <- logPostGrad(CplNM = CplNM,
                                   Mdl.Y = Mdl.Y,
                                   Mdl.X = Mdl.X,
                                   Mdl.parLink = Mdl.parLink,
                                   Mdl.beta = Mdl.beta,
                                   Mdl.betaIdx = Mdl.betaIdx,
                                   parUpdate = parUpdate,
                                   staticArgs = staticArgs)
      
      ## Gradient and Hessian for the likelihood
      logLikGrad.prop <- gradHess.prop[["logLikGrad"]] # n-by-pp
      logLikHess.prop <- hessApprox(logLikGrad.pror, hessMethod)

      ## Gradient Hessian for the prior
      logPriGrad.prop <- gradHess.prop[["logPriGrad"]] # pp-by-1
      logLikHess.prop <- hessApprox(logPriGrad.pror, hessMethod)

      X.prop <- X[, betaIdxProp, drop = FALSE] # n-by-pc
      X.curr <- X[, betaIdxCurr, drop = FALSE] # n-by-pp
      
      ## The gradient in the general Newton's update
      gradObs.prop <- Md(t(X.prop), logLikGrad.prop) + logPriGrad.prop
      
      ## The Hessian in the general Newton's update
      HessObs.pp <- tMdN(X.prop, logLikHess.prop, X.prop) + logPriHess.prop
      HessObs.pc <- tMdN(X.prop, logLikHess.prop, X.curr) + logPriHess.prop
      
      if((iStep <= kSteps)) ## The general Newton's Update
        {
          ## update the parameters
          param <- solve(HessObs.pp)%*%(HessObs.pc%*%param -
                                        gradObs.prop)
          
          ## Update the parameter with current updated results.
          Mdl.beta[[CompCurr]][[parCurr]] <- param
          betaIdxCurr <- betaIdxProp
          
        }
      else if(iStep == (kSteps+1)) # (k+1):th step.  Make a output 
        {
          out <- list(gradObs.cur = gradObs.cur,
                      hessObs.cur = hessObs.cur,
                      invHessObs.cur = invHessObs.cur,
                      param = param.cur) 
        }
       
    }
  return(out) 
}
