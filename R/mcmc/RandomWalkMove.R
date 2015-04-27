randomwalk <- function(propArgs, varSelArgs, priArgs, betaIdxProp, parUpdate,
                       CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                       Mdl.parLink, MargisTypes, staticCache, MCMCUpdateStrategy)
  {

    if(errorFlag)
      {
        out <- list(errorFlag = errorFlag)
      }
    else
      {
        out <- list(param = param,
                    gradObs = gradObs.pp,
                    HessObs = HessObs.pp,
                    HessObsInv = HessObs.pp.Inv,
                    staticCache = staticCache.curr,
                    errorFlag = errorFlag)
        ## print(gradObs.pp)
      }
  }
