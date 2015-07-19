## The gradient for the copula function. scaler input and output NOTE: The numerical
## gradient may not work well if the tabular version of Kendall's tau is used (due to the
## precision).
logCplRepGradNum <- function(CplNM, u, parCplRep, parCaller)
{
  require("numDeriv")
  logCplGradNumFun <- function(x, u, iRun, parCaller, CplNM, parCplRep)
  {
    if(tolower(parCaller) %in% paste("u", 1:ncol(u), sep = ""))
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

    out <- logCplLik(CplNM = CplNM,
                     u = u[iRunInRow, , drop = FALSE],
                     parCplRep = lapply(parCplRep, function(x) x[iRunInRow, , drop = FALSE]),
                     sum = FALSE)
    return(out)
  }

  nDim <- ncol(u)
  nObs <- nrow(u)

  if(tolower(parCaller) %in% paste("u", 1:nDim, sep = ""))
  { ## Calling the marginal CDF u1, u2
    nDimGrad <- 1
  }
  else
  { ## Calling copula parameters
    nDimGrad <- ncol(parCplRep[[parCaller]])
  }

  logCplGradObs.num <- matrix(NA, nObs, nDimGrad)

  for(iRun in 1:(nObs*nDimGrad))
  {
    if(tolower(parCaller) %in% paste("u", 1:nDim, sep = ""))
    { ## Calling the marginal CDF u_i
      xCurr <- u[iRun, parCaller]
    }
    else
    { ## Calling copula parameters
      xCurr <- parCplRep[[parCaller]][iRun]
    }

    gradTry <- try(grad(func = logCplGradNumFun,
                        x = xCurr,
                        u = u,
                        iRun = iRun,
                        parCaller = parCaller,
                        CplNM =  CplNM,
                        parCplRep = parCplRep), silent = TRUE)

    if(is(gradTry, "try-error"))
    {
      logCplGradObs.num[iRun] <- NA
    }
    else
    {
      logCplGradObs.num[iRun] <- gradTry
    }
  }
  return(logCplGradObs.num)
}
