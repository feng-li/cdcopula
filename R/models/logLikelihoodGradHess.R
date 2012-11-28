##' The gradient for the log posterior in the copula model
##'
##' This function usually does not need to modify. Note that this function does
##' not give the complete posteriors but its components,  i.e. gradient for the
##' likelihood function, gradient for the prior,  gradient for the linkage.
##' @title The gradient for the log posterior
##' @param CplNM
##' @param MargisTypes
##' @param Mdl.Y
##' @param Mdl.X
##' @param Mdl.parLink
##' @param Mdl.beta
##' @param Mdl.betaIdx
##' @param parUpdate
##' @param priArgs
##' @param staticArgs
##' @return
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Feb 02 22:45:42 CET 2012;
##'       Current: Thu Feb 02 22:45:48 CET 2012.
logLikelihoodGradHess <- function(CplNM,
                            MargisTypes,
                            Mdl.Y,
                            Mdl.X,
                            Mdl.parLink,
                            Mdl.beta,
                            Mdl.betaIdx,
                            parUpdate,
                            priArgs,
                            varSelArgs,
                            staticArgs,
                            gradMethods = c("analytic", "numeric"))
{
  ## The updating chain
  chainCaller <- parCaller(parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]
  Mdl.par <- staticArgs[["Mdl.par"]]
  ## if(parCaller == "tau") browser()

  if(tolower(CompCaller) != tolower(CplNM))
    {
      yCurr <- Mdl.Y[[CompCaller]]
      parCurr <- Mdl.par[[CompCaller]]
      typeCurr <- MargisTypes[[CompCaller]]

      ## Gradient Fraction in the marginal component. n-by-1
      if("analytic" %in% tolower(gradMethods))
        {
          FracGradObs <- MargiModelGrad(
              par = parCurr,
              y = yCurr,
              type = typeCurr,
              parCaller = parCaller)
        }
      if("numeric" %in% tolower(gradMethods))
        {
          ## A simple wrapper that calculates the numerical gradient for given
          ## parameters

          ## NOTE: that this part depends on Chris Sims' numgrad() function. It
          ## is under the flutils module.

          MargiModelGradNumFun <- function(
              x, parCaller, parCurr, yCurr, typeCurr)
            {
              parCurr[[parCaller]] <- x
              MargiLogLikObs <- MargiModel(
                  par = parCurr,
                  y = yCurr,
                  type = typeCurr)$u
              out <- MargiLogLikObs
              return(out)
            }

          FracGradObs <- matrix(NA, length(yCurr), 1)
          for(i in 1:length(yCurr))
            {
              FracGradObs[i] <- numgrad(
                  fcn = MargiModelGradNumFun,
                  x = parCurr[[parCaller]][i],
                  parCaller = parCaller,
                  parCurr = lapply(parCurr, function(x, i)x[i], i = i),
                  yCurr = yCurr[i],
                  typeCurr = typeCurr)$g
            }

        }

      ## The Copula parameter caller is the marginal CDF, i.e. u1,  u2, ...
      cplCaller <- paste("u", which(CompCaller%in%names(MargisTypes)), sep = "")

      staticArgs[["Mdl.u"]][, CompCaller] <- MargiModel(
          y = yCurr,
          type = typeCurr,
          par = parCurr)[["u"]]
    }
  else
    {
      ## Only update the gradient for copula parameters
      ## Gradient Fraction in the copula component.
      FracGradObs <- 1
      cplCaller <- parCaller
    }

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE LIKELIHOOD
###----------------------------------------------------------------------------

  if("analytic" %in% tolower(gradMethods))
    {
      ## The gradient for the copula function. n-by-1
      logCplGradObs <- logCplGrad(
          CplNM = CplNM,
          u = staticArgs$Mdl.u,
          parCpl = Mdl.par[[CplNM]],
          cplCaller = cplCaller,
          staticArgs = staticArgs)

      ## The gradient for the link function n-by-1
      LinkGradObs <- parCplMeanFunGrad(
          CplNM = CplNM,
          Mdl.par = Mdl.par,
          Mdl.parLink = Mdl.parLink,
          chainCaller = chainCaller)
    }

  if("numeric" %in% tolower(gradMethods))
    {
      ## The gradient for the copula function. scaler input and output
      logCplGradNumFun <- function(x, u,  CompCaller, parCaller, cplCaller,
                                   CplNM, parCpl, staticArgs)
        {
          if(tolower(cplCaller) %in% c("u1", "u2"))
            {
              ## Calling the marginal CDF u1, u2
              ## u <- staticArgs$Mdl.u
              u[, CompCaller] <- x
            }
          else
            {
              ## Calling copula parameters
              parCpl[[parCaller]] <- x
            }

          out <- logCplLik(u = u, CplNM = CplNM,
                           parCpl = parCpl,
                           staticArgs = staticArgs,
                           logLik = FALSE)
          return(out)
        }

      nObs <- length(Mdl.Y[[1]])
      logCplGradObs2 <- matrix(NA, nObs, 1)
      for(i in 1:nObs)
        {
          if(tolower(cplCaller) %in% c("u1", "u2"))
            {
              ## Calling the marginal CDF u1, u2
              xCurr <- staticArgs$Mdl.u[i, CompCaller]
            }
          else
            {
              ## Calling copula parameters
              xCurr <- Mdl.par[[CompCaller]][[parCaller]][i]
            }

          logCplGradObs2[i] <- numgrad(
              fcn = logCplGradNumFun,
              x = xCurr,
              u = staticArgs$Mdl.u[i, , drop = FALSE],
              CompCaller = CompCaller,
              parCaller = parCaller,
              cplCaller = cplCaller,
              CplNM =  CplNM,
              parCpl = lapply(Mdl.par[[CplNM]], function(x, i)x[i], i = i),
              staticArgs = staticArgs)$g
        }
      ## The gradient for the link function n-by-1

    }

  ## The gradient for the likelihood,  n-by-1
  logLikGradObs <- (logCplGradObs*FracGradObs)*LinkGradObs


###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------

  out <- list(logLikGradObs = logLikGradObs,
              staticArgs = staticArgs)

  return(out)
}
