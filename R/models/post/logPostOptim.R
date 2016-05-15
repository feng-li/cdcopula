##' logPostOptim
##'
##' A wrapper of the log posterior with first argument as an vector
##' @param betaVec "vector".
##' @param CplNM "character".
##' @param Mdl.Y "list".
##' @param Mdl.X "list".
##' @param Mdl.beta
##' @param Mdl.betaIdx "list".
##' @param Mdl.parLink "list".
##' @param MCMC.varSelArgs "list".
##' @param Mdl.MargisType "list".
##' @param Mdl.priArgs "list".
##' @param parUpdate "list".
##' @param staticCache "list".
##' @return "scaler" The log posterior.
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Sun Jun 03 19:13:54 CEST 2012;
##'       Current: Sun Jun 03 19:14:00 CEST 2012.
logPostOptim <- function(betaVec, Mdl.MargisType, Mdl.Y, Mdl.X, Mdl.beta,
                         Mdl.betaIdx,Mdl.parLink,MCMC.varSelArgs,
                         Mdl.priArgs,parUpdate, staticCache, MCMC.UpdateStrategy)
{
  ## a wrapper of the log posterior function that can be used for directly
  ## optimization via Newton's method
  Mdl.beta <- parCplSwap(betaInput = betaVec,
                         Mdl.beta = Mdl.beta,
                         Mdl.betaIdx = Mdl.betaIdx,
                         parUpdate = parUpdate)
  logPostOut <- logPost(Mdl.MargisType = Mdl.MargisType,
                        Mdl.Y = Mdl.Y,
                        Mdl.X = Mdl.X,
                        Mdl.beta = Mdl.beta,
                        Mdl.betaIdx = Mdl.betaIdx,
                        Mdl.parLink = Mdl.parLink,
                        MCMC.varSelArgs = MCMC.varSelArgs,
                        Mdl.priArgs = Mdl.priArgs,
                        staticCache = staticCache,
                        parUpdate = parUpdate,
                        MCMC.UpdateStrategy = MCMC.UpdateStrategy)

  out <- logPostOut[["Mdl.logPost"]]
  return(out)
}
