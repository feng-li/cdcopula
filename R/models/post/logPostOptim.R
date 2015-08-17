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
##' @param varSelArgs "list".
##' @param MargisType "list".
##' @param priArgs "list".
##' @param parUpdate "list".
##' @param staticCache "list".
##' @return "scaler" The log posterior.
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Sun Jun 03 19:13:54 CEST 2012;
##'       Current: Sun Jun 03 19:14:00 CEST 2012.
logPostOptim <- function(betaVec, MargisType, Mdl.Y, Mdl.X, Mdl.beta,
                         Mdl.betaIdx,Mdl.parLink,varSelArgs,
                         priArgs,parUpdate, staticCache, MCMCUpdateStrategy)
{
  ## a wrapper of the log posterior function that can be used for directly
  ## optimization via Newton's method
  Mdl.beta <- parCplSwap(betaInput = betaVec,
                         Mdl.beta = Mdl.beta,
                         Mdl.betaIdx = Mdl.betaIdx,
                         parUpdate = parUpdate)

  logPostOut <- logPost(MargisType = MargisType,
                        Mdl.Y = Mdl.Y,
                        Mdl.X = Mdl.X,
                        Mdl.beta = Mdl.beta,
                        Mdl.betaIdx = Mdl.betaIdx,
                        Mdl.parLink = Mdl.parLink,
                        varSelArgs = varSelArgs,
                        priArgs = priArgs,
                        staticCache = staticCache,
                        parUpdate = parUpdate,
                        MCMCUpdateStrategy = MCMCUpdateStrategy)

  out <- logPostOut[["Mdl.logPost"]]

  return(out)
}
