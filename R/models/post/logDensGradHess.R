##' The gradient for the log posterior in the copula model
##'
##' This function usually does not need to modify. Note that this function does not give
##' the complete posteriors but its components, i.e. gradient for the likelihood function,
##' gradient for the prior, gradient for the linkage.
##' @param CplNM "character".
##'
##' @param MargisType "list".
##'
##' @param Mdl.Y "list".
##'
##' @param Mdl.parLink "list".
##'
##' @param parUpdate "list".
##'
##' @param varSelArgs "list".
##'
##' @param staticCache "list".
##'
##' @param gradMethods "character" If is "numeric", the numeric gradient is returned; else
##' return the analytical gradient.
##'
##' @return "list".
##'
##'        errorFlag: If anything went wrong, quit with TRUE.
##'
##'        logLilGradObs: Gradient for the log likelihood
##'
##' @references Li 2012
##' @author Feng Li, Central University of Finance and Economics.
##' @note Created: Thu Feb 02 22:45:42 CET 2012; Current: Mon Dec 22 20:25:44 CST 2014
logDensGradHess <- function(MargisType, Mdl.Y, Mdl.parLink, parUpdate,
                            gradMethods = c("analytic","numeric")[1],
                            staticCache, MCMCUpdateStrategy)
{
  ## The updating chain
  CplNM <- MargisType[length(MargisType)]
  chainCaller <- parCplRepCaller(parUpdate)

  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  Mdl.par <- staticCache[["Mdl.par"]]

  CompNM <- names(parUpdate)
  MargisNM <- CompNM[(CompNM  != CplNM)]
  names(MargisType) <- MargisNM

  R_CPL_NPARALLEL <- as.numeric(Sys.getenv("R_CPL_NPARALLEL"))
###----------------------------------------------------------------------------
### SPLIT THE GRADIENT INTO COPULA AND MARGINAL ACCORDING TO MCMC STRATEGY
###----------------------------------------------------------------------------
  if(CompCaller != CplNM)
  {
    if(tolower(MCMCUpdateStrategy) == "joint")
    {
      evalCpl <- TRUE
      cplCaller <- paste("u", which(MargisNM%in%CompCaller), sep = "")

      evalMargi <- TRUE
      densCaller <- c("u", "d")
    }
    else if(tolower(MCMCUpdateStrategy) == "twostage")
    {
      ## Stage one of the two stage approach
      evalCpl <- FALSE
      cplCaller <- NA

      evalMargi <- TRUE
      densCaller <- c("d")
    }
    else if(tolower(MCMCUpdateStrategy) == "margin")
    {
      evalCpl <- FALSE
      cplCaller <- NA

      evalMargi <- TRUE
      densCaller <- c("d")
    }
    else
    {
      stop(paste("MCMC update strategy:", MCMCUpdateStrategy,
                 "not implemented!"))
    }
  }
  else
  {
    evalCpl <- TRUE
    cplCaller <- parCaller

    evalMargi <- FALSE
    densCaller <- NA
  }
###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE MARGINAL LIKELIHOOD
###----------------------------------------------------------------------------
  if(evalMargi == TRUE)
  {
    yCurr <- Mdl.Y[[CompCaller]]
    parCurr <- Mdl.par[[CompCaller]]
    typeCurr <- MargisType[CompCaller]

    ## Gradient Fraction in the marginal component. n-by-1
    if("analytic" %in% tolower(gradMethods))
    {
      if(R_CPL_NPARALLEL>1)
      {
        MargiGradFUN.NM <- "MargiModelGradParallel"
      }
      else
      {
        MargiGradFUN.NM <- "MargiModelGrad"
      }
      MargiGradObs.ana.caller <- call(MargiGradFUN.NM, par = parCurr, y = yCurr,
                                      type = typeCurr, parCaller = parCaller,
                                      densCaller = densCaller)
      MargiGradObs.ana <- eval(MargiGradObs.ana.caller)
      MargiGradObs.u <- MargiGradObs.ana[, "u"]
      MargiGradObs.d <- MargiGradObs.ana[, "d"]
    }

    ## Numerical Gradients, check with analytical gradients.
    if("numeric" %in% tolower(gradMethods))
    {
      MargiGradObs.num <- MargiModelGradNum(y = yCurr, par = parCurr, type = typeCurr,
                                            parCaller = parCaller,
                                            densCaller = densCaller)
      MargiGradObs.u <- MargiGradObs.num[, "u"]
      MargiGradObs.d <- MargiGradObs.num[, "d"]

      ## DEBUG: Check if any gradient component is not correctly computed.  To check the
      ## overall gradient chain, look at the "PropGNewtonMove()" function. Below evaluates
      ## if the numeric and analytic gradients are consistent
      ## try(plot(sort(MargiGradObs.ana), MargiGradObs.num[order(MargiGradObs.ana)], type
      ## = "l", pch = 20, main = chainCaller), silent = TRUE)
    }
  }
  else
  { ## Only update the gradient for copula parameters
    ## Gradient Fraction in the copula component.
    MargiGradObs.u <- 1
    MargiGradObs.d <- 0
  }

  ## Error checking
  if(any(is.na(MargiGradObs.u)) || any(is.infinite(MargiGradObs.u)))
  {
    return(list(errorFlag = TRUE))
  }

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE COPULA LIKELIHOOD
###----------------------------------------------------------------------------

  if(evalCpl == TRUE)
  {
    if("analytic" %in% tolower(gradMethods))
    {
      if(R_CPL_NPARALLEL>1)
      {
        CplGradFUN.NM <- "logCplRepGradParallel"
      }
      else
      {
        CplGradFUN.NM <- "logCplRepGrad"
      }

      logCplGrad.caller <- call(CplGradFUN.NM, CplNM = CplNM,
                                u = staticCache$Mdl.u,
                                parCplRep = Mdl.par[[CplNM]],
                                parCaller = cplCaller)

      ## The gradient for the copula function. n-by-1
      logCplGradObs.ana <- eval(logCplGrad.caller)
      logCplGradObs <- logCplGradObs.ana
    }
    if("numeric" %in% tolower(gradMethods))
    {
      logCplGradObs.num <- logCplRepGradNum(CplNM = CplNM, u = staticCache$Mdl.u,
                                            parCplRep = Mdl.par[[CplNM]],
                                            parCaller = cplCaller)
      logCplGradObs <- logCplGradObs.num
    }
  }
  else
  {
    logCplGradObs <- 1
  }

  ## Error checking
  if(length(logCplGradObs) == 0) browser()
  if(any(is.na(logCplGradObs)) || any(is.infinite(logCplGradObs)))
  {
    return(list(errorFlag = TRUE))
  }


###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------
  ## The gradient for the full likelihood,  n-by-1
  logDensGradObs <- (logCplGradObs*MargiGradObs.u) + MargiGradObs.d #*LinkGradObs

  ## The output
  out <- list(logGradObs = logDensGradObs, # n-by-1
              logHessObs = NA,
              errorFlag = FALSE)
  return(out)
}
