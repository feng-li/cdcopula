##' Metropolis–Hastings algorithm with K-step Newton method with variable
##' selection.
##'
##' @title 
##' @param param.cur 
##' @param gradhess.fun.name 
##' @param logpost.fun.name 
##' @param nNewtonStep 
##' @param Params 
##' @param hessMethod 
##' @param Y 
##' @param x 
##' @param callParam 
##' @param splineArgs 
##' @param priorArgs 
##' @param prop.df 
##' @param Params_Transform 
##' @return 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Thu Feb 17 14:03:14 CET 2011;
##'       Current: Wed Feb 01 16:09:04 CET 2012.
##' DEPENDS: mvtnorm
##' TODO: replace the old multivariate t functions to mvtnorm functions
MHPropWithKStepNewton <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                                  Mdl.parLink 
                                  parUpdate,
                                  priorArgs, varSelArgs, propArgs, staticArgs)   
{

###----------------------------------------------------------------------------
### Variable Selection proposal 
###----------------------------------------------------------------------------
  
  ## The updating component parameter chain
  cp <- parCaller(parUpdate)
  CompCurr <- cp[1]
  parCurr <- cp[2]

  ## The current variable selection indicators
  betaIdxCurr <- Mdl.betaIdx[[CompCurr]][[parCurr]]

  ## Randomly propose a subset TODO: A general proposal function
  varSelCand <- varSelArgs[[CompCurr]][[parCurr]]$cand
  IdxArgs <- propArgs[[CompCurr]][[parCurr]][["indicators"]]
  
  betaIdxProp <- betaIdxCurr # The proposal base line

  ## Binomial proposal a small subset
  betaIdxPropCand <- which(rbinom(n = length(varSelCand) , size = 1,
                             prob = IdxArgs$prob) == 1)

  ## Special case to make sure at least one variable is proposed a change
  if(length(betaIdxPropCand) == 0)
    {
      betaIdxPropCand <- sample(varSelCand, 1)
    }

  ## The proposed variable selection indicators 
  betaIdxProp[betaIdxPropCand] <- 1 - betaIdxCurr[betaIdxPropCand]
  
###----------------------------------------------------------------------------
### Make good proposal via K-steps Newton's method
###----------------------------------------------------------------------------  
  
  ## Newton method to approach the posterior for the current draw 
  KStepNewton1 <- kStepsNewtonMove(propArgs = propArgs,
                                   varSelArgs = varSelArgs,
                                   priArgs = priArgs,
                                   betaIdxProp = betaIdxProp,
                                   parUpdate = parUpdate,
                                   CplNM = CplNM,
                                   Mdl.Y = Mdl.Y,
                                   Mdl.X = Mdl.X,
                                   Mdl.parLink = Mdl.parLink,
                                   Mdl.beta = Mdl.beta,
                                   Mdl.betaId = Mdl.betaIdx,
                                   staticArgs)   
 
  ## The information for proposed density via K-step Newton's method
  param.cur.prop <- KStepNewton1$param # mean information 
  HessObs.cur.prop <- KStepNewton1$hessObs # variance information
  invHessObs.cur.prop <- KStepNewton1$invHessObs

  ## Check if it is a good proposal
  if(any(is.na(HessObs.cur.prop)) ||
     is(try(chol(-invHessObs.cur.prop), silent=T), "try-error")) # Something is wrong.
    {
       logjump.cur2prop <- NaN
       param.prop <- NaN
    }
  else # Continues with the Metropolis-Hastings
    {
      ## Propose a draw from multivariate t-distribution based on the proposed information
      param.prop <- RndMultiT(param.cur.prop, -invHessObs.cur.prop, prop.df)
      
      ## The jump density from the proposed draw 
      logjump.cur2prop <- DensMultiT(param.prop, param.cur.prop, -HessObs.cur.prop, prop.df)
                                        
    }
  
###----------------------------------------------------------------------------
### 
###----------------------------------------------------------------------------  

  if(any(is.na(param.prop))) # Previous proposal unsuccessful
    {
      HessObs.prop.prop <- NaN
    }
  else # all are right
    {
      ## Newton method to approach the posterior for the proposed draw 
      KStepNewton2 <- KStepNewtonMove()  

      ## The information for proposed density via K-step Newton's method
      param.prop.prop <- KStepNewton2$param.cur 
      HessObs.prop.prop <- KStepNewton2$hessObs.cur
      invHessObs.prop.prop <- KStepNewton2$invHessObs.cur
    }

  if(any(is.na(HessObs.prop.prop))) # Something is wrong at KStepNewton2,  reject it.
    {
      logpost.cur <- NaN
      logpost.prop <- NaN
      logjump.prop2cur <- NaN
    }
  else # all are right
    {
      # The jump density from propose draw to current draw.
      logjump.prop2cur <- DensMultiT(param.cur, param.prop.prop, -invHessObs.prop.prop,
                                     prop.df)  
      Params.prop <- Params
      Params.prop[[callParam$id]][callParam$subset] <- param.prop
      Params.cur <- Params
      Params.cur[[callParam$id]][callParam$subset] <- param.cur

      ## The log posterior for the proposed draw
      caller.prop <- call(logpost.fun.name,Y = Y, x = x, Params = Params.prop,callParam =
                          callParam ,  priorArgs = priorArgs, splineArgs = splineArgs,
                          Params_Transform = Params_Transform)
      logpost.prop <-  eval(caller.prop) 
      
      ## The log posterior for the current draw
      caller.cur <- call(logpost.fun.name, Y = Y,  x = x, Params = Params.cur,callParam =
                         callParam, priorArgs = priorArgs, splineArgs = splineArgs,
                         Params_Transform = Params_Transform) 
      logpost.cur <- eval(caller.cur)
    }
  
###----------------------------------------------------------------------------
### Compute the MH ratio and make decision to update or keep current draw. 
###----------------------------------------------------------------------------  

  ## compute the MH ratio.
  log.r <- logpost.prop - logpost.cur + logjump.prop2cur - logjump.cur2prop
  r <- exp(log.r) 

  ## the acceptance probability
  accept.prob <- min(1, r) 

  if(!is.na(accept.prob) && runif(1) < accept.prob) 
    {
      param.out <- param.prop  # keep update
    }
  else 
    {
      param.out <- param.cur # keep current
      accept.prob <- 0 # set acceptance probs to zero.
    }
  
###----------------------------------------------------------------------------
### The output
###----------------------------------------------------------------------------  
  out <- list(param.out = param.out, accept.prob = accept.prob)

  ##cat("prop:", param.prop, "cur:", param.cur, "\n")
  return(out)
}
