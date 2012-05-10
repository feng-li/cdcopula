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
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Feb 02 22:45:42 CET 2012;
##'       Current: Thu Feb 02 22:45:48 CET 2012.
logPostGradHess <- function(CplNM, MargisTypes, Mdl.Y, Mdl.X, Mdl.parLink,
                            Mdl.beta, Mdl.betaIdx, parUpdate, priArgs, staticArgs)
{
  ## The updating chain
  chainCaller <- parCaller(parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]
  Mdl.par <- staticArgs$Mdl.par

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE LIKELIHOOD
###----------------------------------------------------------------------------
  if(tolower(CompCaller) != tolower(CplNM))
    {
      ## Gradient Fraction in the marginal component. n-by-1
      parMargis <- Mdl.par[names(MargisTypes)]
      FracGradObs <- MargisGrad(
                       parMargis = parMargis,
                       Mdl.Y = Mdl.Y,
                       MargisTypes = MargisTypes,
                       chainCaller = chainCaller)

      cplCaller <- "u" # The Copula parameter caller is the marginal CDF

      staticArgs[["Mdl.u"]] <- MargisModels(
                                 Mdl.Y = Mdl.Y,
                                 MargisTypes = MargisTypes,
                                 parMargis = parMargis,
                                 whichMargis = CompCaller,
                                 staticArgs = staticArgs)[["Mdl.u"]]
    }
  else
    {
      ## Only update the gradient for copula parameters
      ## Gradient Fraction in the copula component.
      FracGradObs <- 1
      cplCaller <- chainCaller[2]
    }

  ## The gradient for the copula function. n-by-1
  logCplGradObs <- logCplGrad(
                     CplNM = CplNM,
                     u = staticArgs$Mdl.u,
                     parCpl = Mdl.par[[CplNM]],
                     cplCaller = cplCaller,
                     staticArgs = staticArgs)

  ## The gradient for the link function n-by-1
  LinkGradObs <-  parMeanFunGrad(
                    par = Mdl.par[[chainCaller[1]]][[chainCaller[2]]],
                    link = Mdl.parLink[[chainCaller[1]]][[chainCaller[2]]])

  ## The gradient and Hessian for the likelihood
  logLikGradObs <- (logCplGradObs*FracGradObs)*LinkGradObs

###----------------------------------------------------------------------------
### GRADIENT AND HESSIAN IN THE PRIOR COMPONENT
###----------------------------------------------------------------------------

  logPriGradHessObs <- logPriGradHess(Mdl.X = Mdl.X,
                                      Mdl.parLink = Mdl.parLink,
                                      Mdl.beta = Mdl.beta,
                                      Mdl.betaIdx = Mdl.betaIdx,
                                      varSelArgs = varSelArgs,
                                      priArgs = priArgs,
                                      chainCaller = chainCaller)

###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------

  out <- list(logLikGradObs = logLikGradObs,
              logPriGradHessObs = logPriGradHessObs,
              staticArgs = staticArgs)
  return(out)
}
