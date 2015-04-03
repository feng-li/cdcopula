##' The gradient for the log posterior in the copula model
##'
##' This function usually does not need to modify. Note that this function does
##' not give the complete posteriors but its components,  i.e. gradient for the
##' likelihood function, gradient for the prior,  gradient for the linkage.
##' @param CplNM "character".
##'
##' @param MargisTypes "list".
##'
##' @param Mdl.Y "list".
##'
##' @param Mdl.X "list".
##'
##' @param Mdl.parLink "list".
##'
##' @param Mdl.beta "list".
##'
##' @param Mdl.betaIdx "list".
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
logLikelihoodGradHess <- function(CplNM, MargisTypes, Mdl.Y, Mdl.X, Mdl.parLink,
                                  Mdl.beta, Mdl.betaIdx, parUpdate, varSelArgs,
                                  staticCache,
                                  gradMethods = c("analytic", "numeric")[1:2],
                                  MCMCUpdateStrategy)
{
  ## The updating chain
  chainCaller <- parCplRepCaller(CplNM = CplNM, parUpdate)

  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  Mdl.par <- staticCache[["Mdl.par"]]

  CompNM <- names(Mdl.beta)
  MargisNM <- CompNM[(CompNM  != CplNM)]
  names(MargisTypes) <- MargisNM

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
          densCaller <- c("u")
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
      typeCurr <- MargisTypes[CompCaller]

      ## Gradient Fraction in the marginal component. n-by-1
      if("analytic" %in% tolower(gradMethods))
        {
          MargiGradObs.ana <- MargiModelGrad(
                  par = parCurr,
                  y = yCurr,
                  type = typeCurr,
                  parCaller = parCaller,
                  densCaller = densCaller)
          MargiGradObs <- MargiGradObs.ana
        }
      if("numeric" %in% tolower(gradMethods))
        {
          ## A simple wrapper that calculates the numerical gradient for given
          ## parameters NOTE: The numeric gradient depends on numDeriv package

          require(numDeriv)
          MargiModelGradNumFun <- function(x, parCaller, parCurr, yCurr, typeCurr)
            {
              parCurr[[parCaller]] <- x
              MargiLogLikObs <- MargiModel(
                      par = parCurr,
                      y = yCurr,
                      type = typeCurr)$u
              out <- MargiLogLikObs
              return(out)
            }

          nObs <- length(Mdl.Y[[1]])
          MargiGradObs.num <- matrix(NA, nObs, 1)
          for(i in 1:nObs)
            {
              gradTry <- try(grad(
                      func = MargiModelGradNumFun,
                      x = parCurr[[parCaller]][i],
                      parCaller = parCaller,
                      parCurr = lapply(parCurr, function(x, i)x[i], i = i),
                      yCurr = yCurr[i],
                      typeCurr = typeCurr), silent = TRUE)

              if(is(gradTry, "try-error"))
                {
                  MargiGradObs.num[i] <- NA
                }
              else
                {
                  MargiGradObs.num[i] <- gradTry
                }
            }
          MargiGradObs <- MargiGradObs.num

          ## Evaluate if the numeric and analytic gradients are consistent
          try(plot(sort(MargiGradObs.ana),
               MargiGradObs.num[order(MargiGradObs.ana)],
               type = "l", pch = 20, main = chainCaller), silent = TRUE)

        }


      staticCache[["Mdl.u"]][, CompCaller] <- MargiModel(
              y = yCurr,
              type = typeCurr,
              par = parCurr)[["u"]]
    }
  else
    { ## Only update the gradient for copula parameters
      ## Gradient Fraction in the copula component.
      MargiGradObs <- 1
    }

  ## Error checking
  if(any(is.na(MargiGradObs)) || any(is.infinite(MargiGradObs)))
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

          ## The gradient for the copula function. n-by-1
          logCplGradObs.ana <- logCplGrad(
                  CplNM = CplNM,
                  u = staticCache$Mdl.u,
                  parCplRep = Mdl.par[[CplNM]],
                  cplCaller = cplCaller,
                  Mdl.X = Mdl.X,
                  Mdl.beta = Mdl.beta)
          logCplGradObs <- logCplGradObs.ana
        }
      if("numeric" %in% tolower(gradMethods))
        {
          ## The gradient for the copula function. scaler input and output NOTE: The
          ## numerical gradient may not work well if the tabular version of Kendall's tau
          ## is used (due to the precision).

          require(numDeriv)
          logCplGradNumFun <- function(x, u, iRun, CompCaller, parCaller, cplCaller,
                                       CplNM, parCplRep, staticCache)
            {
              if(tolower(cplCaller) %in% paste("u", 1:ncol(u), sep = ""))
                {
                  u[iRun] <- x
                }
              else
                { ## Calling copula parameters
                  parCplRep[[parCaller]][iRun] <- x
                }

              nObs <- nrow(u)
              iRunInRow = iRun%%nObs
              iRunInRow[iRunInRow == 0] <- nObs

              out <- logCplLik(
                      CplNM = CplNM,
                      u = u[iRunInRow, , drop = FALSE],
                      parCplRep = lapply(parCplRep, function(x) x[iRunInRow, , drop = FALSE]),
                      sum = FALSE)
              return(out)
            }

          nDim <- length(Mdl.Y)
          if(tolower(cplCaller) %in% paste("u", 1:nDim, sep = ""))
            { ## Calling the marginal CDF u1, u2
              nDimGrad <- 1
            }
          else
            { ## Calling copula parameters
              nDimGrad <- nDim
            }

          nObs <- length(Mdl.Y[[1]])
          logCplGradObs.num <- matrix(NA, nObs, nDimGrad)

          for(iRun in 1:(nObs*nDimGrad))
            {
              if(tolower(cplCaller) %in% paste("u", 1:nDim, sep = ""))
                { ## Calling the marginal CDF u_i
                  xCurr <- staticCache$Mdl.u[iRun, CompCaller]
                }
              else
                { ## Calling copula parameters
                  xCurr <- Mdl.par[[CompCaller]][[parCaller]][iRun]
                }

              gradTry <- try(grad(
                      func = logCplGradNumFun,
                      x = xCurr,
                      u = staticCache$Mdl.u,
                      iRun = iRun,
                      CompCaller = CompCaller,
                      parCaller = parCaller,
                      cplCaller = cplCaller,
                      CplNM =  CplNM,
                      parCplRep = Mdl.par[[CplNM]],
                      staticCache = staticCache), silent = TRUE)

              if(is(gradTry, "try-error"))
                {
                  logCplGradObs.num[iRun] <- NA
                }
              else
                {
                  logCplGradObs.num[iRun] <- gradTry
                }
            }

          ## The Gradient For The Link Function n-by-1
          logCplGradObs <- logCplGradObs.num

          ## Debugging plot for the numerical and analytical gradients
          if(cplCaller == "u1") try(plot(sort(logCplGradObs.ana),
                 logCplGradObs.num[order(logCplGradObs.ana)],
                 type = "l", pch = 20, main = chainCaller), silent = TRUE)

        }
    }
  else
    {
      logCplGradObs <- 1
    }

  ## Error checking
  if(any(is.na(logCplGradObs)) || any(is.infinite(logCplGradObs)))
    {
      return(list(errorFlag = TRUE))
    }

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE LINK FUNCTION
###----------------------------------------------------------------------------
  ## The gradient for the link function n-by-1
  LinkGradObs <- parCplMeanFunGrad(
          CplNM = CplNM,
          Mdl.par = Mdl.par,
          Mdl.parLink = Mdl.parLink,
          chainCaller = chainCaller)

  ## Error checking
  if(any(is.na(LinkGradObs)) || any(is.infinite(LinkGradObs)))
    {
      return(list(errorFlag = TRUE))
    }

###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------

  ## The gradient for the likelihood,  n-by-1
  logLikGradObs <- (logCplGradObs*MargiGradObs)*LinkGradObs

  ## The output
  out <- list(logLikGradObs = logLikGradObs,
              logLikHessObs = NA,
              errorFlag = FALSE)


  return(out)
}
