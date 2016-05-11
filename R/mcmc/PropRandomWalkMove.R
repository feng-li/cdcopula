randomwalk <- function(MCMC.propArgs, Mdl.varSelArgs, Mdl.priArgs, betaIdxProp, parUpdate,
                       CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                       Mdl.parLink, Mdl.MargisType, staticCache, MCMC.UpdateStrategy)
{
  errorFlag <- NA
    if(errorFlag)
      {
        out <- list(errorFlag = errorFlag)
      }
    else
      {
        out <- list(param = NA,
                    gradObs = NA,
                    HessObs = NA,
                    HessObsInv = NA,
                    staticCache = NA,
                    errorFlag = NA)
        ## print(gradObs.pp)
      }
  }
