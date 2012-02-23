##' The gradient for the log posterior in the copula model
##'
##' This function usually does not need to modify. Note that this function does
##' not give the complete posteriors but its components,  i.e. gradient for the
##' likelihood function, gradient for the prior,  gradient for the linkage. 
##' @title The gradient for the log posterior
##' @param CplNM 
##' @param Mdl.Y 
##' @param Mdl.par 
##' @param parUpdate 
##' @param u 
##' @return 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Feb 02 22:45:42 CET 2012;
##'       Current: Thu Feb 02 22:45:48 CET 2012.
logPostGrad <- function(CplNM, MargisTypes, Mdl.Y, Mdl.X, Mdl.parLink,
                        Mdl.beta, Mdl.betaIdx, parUpdate, staticArgs)
{

  ## The updating chain
  chainCaller <- parCaller(parUpdate)
  Mdl.par <- staticArgs$Mdl.par
  
###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE LIKELIHOOD
###----------------------------------------------------------------------------  
  if(tolower(chainCaller[1]) != tolower(CplNM))
    {
      cplCaller <- chainCaller[1] # The Copula parameter caller is the
                                        # marginal CDF

      ## Update the current CDF margin for the marginal model.
      parMargis <- Mdl.par[names(MargisTypes)]
      
      staticArgs[["Mdl.u"]] <-
        MargisModels(Mdl.Y = Mdl.Y,
                     MargisTypes = MargisTypes,
                     parMargis = parMargis,
                     whichMargis = cplCaller,
                     staticArgs = staticArgs)[["Mdl.u"]]

      ## Gradient Fraction in the marginal component. n-by-1
      FracGradOut <-
        MargisGrad(parMargis, Mdl.Y, MargisTypes, chainCaller)
    }
  else
    {
      cplCaller <- chainCaller[2]
      
      ## Gradient Fraction in the copula component.
      FracGradOut <- 1
    }

  ## The gradient for the copula function. n-by-1
  logCplGradOut <- logCplGrad(CplNM = CplNM,
                              u = staticArgs$Mdl.u,
                              parCpl = Mdl.par[[CplNM]],
                              cplCaller = cplCaller,
                              staticArgs = staticArgs)

  ## The gradient for the link function n-by-1
  LinkGradOut <-
    parMeanFunGrad(par = Mdl.par[[chainCaller[1]]][[chainCaller[2]]],
                   link = Mdl.parLink[[chainCaller[1]]][[chainCaller[2]]])

  ## The gradient and Hessian for the likelihood
  logLikGradOut <- (logCplGradOut*FracGradOut)*LinkGradOut

###----------------------------------------------------------------------------
### GRADIENT IN THE PRIOR COMPONENT
###----------------------------------------------------------------------------

  ## p-by-1
  logPriGradOut <- logPriGrad(Mdl.X = Mdl.X,
                              Mdl.parLink = Mdl.parLink,
                              MdlCurr.beta = MdlCurr.beta,
                              MdlCurr.betaIdx = MdlCurr.betaIdx,
                              varSelArgs = varSelArgs,
                              priArgs = priArgs,
                              parUpdate = parUpdate, 
                              chainCaller = chainCaller)
  
###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------    
  
  out <- list(logLikGrad = logLikGradOut,
              logPriGrad = logPriGradOut, 
              staticArgs = staticArgs) 
  return(out)
}
