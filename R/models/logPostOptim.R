##' A wrapper of the log likelihood with first argument as an vector
##'
##'
##' @title log posterior function
##' @param betaVec "vector"
##' @param CplNM
##' @param Mdl.Y
##' @param Mdl.X
##' @param Mdl.betaIdx
##' @param Mdl.parLink
##' @param varSelArgs
##' @param MargisTypes
##' @param priArgs
##' @param parUpdate
##' @param staticArgs
##' @return "scaler" The log posterior.
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Sun Jun 03 19:13:54 CEST 2012;
##'       Current: Sun Jun 03 19:14:00 CEST 2012.
logPostOptim <- function(betaVec, CplNM, Mdl.Y, Mdl.X, Mdl.beta,
                         Mdl.betaIdx,Mdl.parLink,varSelArgs, MargisTypes,
                         priArgs,parUpdate,staticArgs)
  {
    ## a wrapper of the log posterior function that can be used for directly
    ## optimization via Newton's method
    Mdl.beta <- parSwap(
        betaInput = betaVec,
        Mdl.beta = Mdl.beta,
        Mdl.betaIdx = Mdl.betaIdx,
        parUpdate = parUpdate)

    ## Update logPost
    out <- logPost(
        CplNM = CplNM,
        Mdl.Y = Mdl.Y,
        Mdl.X = Mdl.X,
        Mdl.beta = Mdl.beta,
        Mdl.betaIdx = Mdl.betaIdx,
        Mdl.parLink = Mdl.parLink,
        varSelArgs = varSelArgs,
        MargisTypes = MargisTypes,
        priArgs = priArgs,
        parUpdate = parUpdate,
        staticArgs = staticArgs)[["Mdl.logPost"]]

    return(out)
  }
