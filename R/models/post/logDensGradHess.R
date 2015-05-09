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
                            gradMethods = c("analytic","numeric")[1:2],
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
          MargiGradObs.ana <- MargiModelGrad(
                  par = parCurr,
                  y = yCurr,
                  type = typeCurr,
                  parCaller = parCaller,
                  densCaller = densCaller)
          MargiGradObs.u <- MargiGradObs.ana[["u"]]
          MargiGradObs.d <- MargiGradObs.ana[["d"]]

        }
      if("numeric" %in% tolower(gradMethods))
        {
          ## A simple wrapper that calculates the numerical gradient for given parameters
          ## NOTE: The numeric gradient depends on numDeriv package
          MargiModelGradNumFun.subtask <- function(subtask, data.parent.env, data.global.env)
            {
              ## This order is very important. The data in local environment should always
              ## be nested inside global environment.
              list2env(data.parent.env, envir = environment())

              nSubTask <- length(subtask)
              out <- matrix(NA, nSubTask, 1)
              require("numDeriv")
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

              for(i in 1:nSubTask)
                {
                  gradTry <-  try(
                          grad(
                                  func = MargiModelGradNumFun,
                                  x = parCurr[[parCaller]][i],
                                  parCaller = parCaller,
                                  parCurr = lapply(parCurr, function(x, i)x[i], i = i),
                                  yCurr = yCurr[subtask][i],
                                  typeCurr = typeCurr), silent = TRUE)

                  if(is(gradTry, "try-error"))
                    {
                      out[i] <- NA
                    }
                  else
                    {
                      out[i] <- gradTry
                    }
                }
              return(out)
            }

          nObs <- length(Mdl.Y[[1]])
          nCores <- detectCores()
          tasks <- data.partition(nObs, list(N.subsets = nCores, partiMethod = "ordered"))
          data.current.env <- as.list(environment())

          MargiGradObs.num.uLst <- lapply(X = tasks,
                                      FUN = MargiModelGradNumFun.subtask,
                                      data.parent.env = data.current.env)
          MargiGradObs.u <- unlist(MargiGradObs.num.uLst)
          MargiGradObs.d <- logDensGradHessNum(MargisType, Mdl.Y, Mdl.parLink, parUpdate,
                                               staticCache, MCMCUpdateStrategy = "twostage")$logGradObs

          ## DEBUG: Check if any gradient component is not correctly computed.  To check
          ## the overall gradient chain, look at the "PropGNewtonMove()" function. Below
          ## evaluates if the numeric and analytic gradients are consistent
          ## try(plot(sort(MargiGradObs.ana),
          ##          MargiGradObs.num[order(MargiGradObs.ana)],
          ##          type = "l", pch = 20, main = chainCaller), silent = TRUE)

        }

      staticCache[["Mdl.u"]][, CompCaller] <- MargiModel(
              y = yCurr,
              type = typeCurr,
              par = parCurr)[["u"]]
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

          ## The gradient for the copula function. n-by-1
          logCplGradObs.ana <- logCplGrad(
                  CplNM = CplNM,
                  u = staticCache$Mdl.u,
                  parCplRep = Mdl.par[[CplNM]],
                  cplCaller = cplCaller)
          logCplGradObs <- logCplGradObs.ana
        }
      if("numeric" %in% tolower(gradMethods))
        {
          ## The gradient for the copula function. scaler input and output NOTE: The
          ## numerical gradient may not work well if the tabular version of Kendall's tau
          ## is used (due to the precision).
          require("numDeriv")
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
          if(cplCaller == "u1")
            {
              ## try(plot(sort(logCplGradObs.ana),
              ##          logCplGradObs.num[order(logCplGradObs.ana)],
              ##          type = "p", pch = 20, main = chainCaller), silent = TRUE)
            }
          ## browser()
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
