## The gradient for the copula function. scaler input and output NOTE: The numerical
## gradient may not work well if the tabular version of Kendall's tau is used (due to the
## precision) in BB7 copula.
#' @export
logCplRepGradNum <- function(CplNM, Mdl.u, parCplRep, parCaller)
{

    if(all(sapply(Mdl.u, ncol)  == 1))
    {
        u <- do.call(cbind, Mdl.u)
    }
    else
    {
        stop("Mixed margins are not ready yet.")
    }



  require("numDeriv", quietly = TRUE)
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

    nObs <- nrow(u[[1]])
    iRunInRow = iRun%%nObs
    iRunInRow[iRunInRow == 0] <- nObs

    out <- logCplLik(CplNM = CplNM,
                     Mdl.u = as.list(u[iRunInRow]),
                     parCplRep = lapply(parCplRep, function(x) x[iRunInRow, , drop = FALSE]),
                     sum = FALSE)
    return(out)
  }

  nDim <- ncol(u[[1]])
  nObs <- nrow(u[[1]])

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
      which.u <- as.numeric(substr(parCaller, 2, nchar(parCaller)))
      xCurr <- u[iRun, which.u]
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
