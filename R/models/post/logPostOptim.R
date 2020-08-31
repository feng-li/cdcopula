#' logPostOptim
#'
#' A wrapper of the log posterior with first argument as an vector
#' @param betaVec "vector".
#' @param CplNM "character".
#' @param Mdl.Y "list".
#' @param Mdl.X "list".
#' @param Mdl.beta NA
#' @param Mdl.betaIdx "list".
#' @param Mdl.parLink "list".
#' @param Mdl.varSelArgs "list".
#' @param Mdl.MargisType "list".
#' @param Mdl.priArgs "list".
#' @param parUpdate "list".
#' @param staticCache "list".
#' @return "scaler" The log posterior.
#' @references Li 2012
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @export
logPostOptim <- function(betaVec, Mdl.MargisType, Mdl.Y, Mdl.X, Mdl.beta,
                         Mdl.betaIdx,Mdl.parLink,Mdl.varSelArgs, Mdl.algorithm,
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
                          Mdl.varSelArgs = Mdl.varSelArgs,
                          Mdl.priArgs = Mdl.priArgs,
                          staticCache = staticCache,
                          parUpdate = parUpdate,
                          Mdl.algorithm = Mdl.algorithm, # full run?
                          MCMC.UpdateStrategy = MCMC.UpdateStrategy)

    out <- logPostOut[["Mdl.logPost"]]
    return(out)
}
