## APPROACH TWO: This version calculates the numerical gradient (for log
## posterior) via the full log likelihood function
logDensGradHessNum <- function(MargisType, Mdl.Y, Mdl.parLink, parUpdate,
                               staticCache, MCMCUpdateStrategy)
{
  ## The updating chain
  chainCaller <- parCplRepCaller(parUpdate)


  Mdl.par <- staticCache[["Mdl.par"]]
  Mdl.u <- staticCache[["Mdl.u"]]
  Mdl.d <- staticCache[["Mdl.d"]]

  require("parallel")
  nSubTasks <- detectCores()
  dataSubIdxLst <- data.partition(
          nObs = nrow(Mdl.u), args = list(N.subsets = nSubTasks, partiMethod = "ordered"))

  ## Uncomment this can use in-node parallelism
  logDensGradObs.Lst <- parLapply(
          cl, dataSubIdxLst,
          logDensGradNum,
          Mdl.Y = Mdl.Y,
          Mdl.u = Mdl.u,
          Mdl.d = Mdl.d,
          Mdl.par = Mdl.par,
          parUpdate = parUpdate,
          MargisType = MargisType,
          MCMCUpdateStrategy = MCMCUpdateStrategy)

  ## logDensGradObs.Lst <- lapply(
  ##         dataSubIdxLst,
  ##         logDensGradNum,
  ##         Mdl.Y = Mdl.Y,
  ##         Mdl.u = Mdl.u,
  ##         Mdl.d = Mdl.d,
  ##         Mdl.par = Mdl.par,
  ##         parUpdate = parUpdate,
  ##         MargisType = MargisType,
  ##         MCMCUpdateStrategy = MCMCUpdateStrategy)


  ## stopCluster(cl)

  logGradObs <- do.call(rbind, logDensGradObs.Lst)

  return(list(logGradObs = logGradObs,
              logHessObs = NA))
}

logDensGradNum <- function(dataSubIdx, MargisType, Mdl.Y,Mdl.u, Mdl.d, Mdl.par,parUpdate,
                           MCMCUpdateStrategy)
  {
    require("numDeriv")

    ## The updating chain
    chainCaller <- parCplRepCaller(parUpdate)
    CompCaller <- chainCaller[1]
    parCaller <- chainCaller[2]

    Par <- Mdl.par[[CompCaller]][[parCaller]]

    nSubObs <- length(dataSubIdx)
    nCol <- ncol(Par)
    out <- matrix(NA, nSubObs, nCol)

    subfun <- function(x, iSubObs, dataSubIdx) x[dataSubIdx[iSubObs], , drop = FALSE]
    subfundeep <- function(x, iSubObs, dataSubIdx)
      {
        lapply(x, subfun, iSubObs = iSubObs, dataSubIdx = dataSubIdx)
      }


    for(iSubObs in 1:nSubObs)
      {
        for(jPar in 1:nCol)
          {
            ## cat(iSubObs, dataSubIdx[iSubObs], "\n")
            gradTry <- try(
                    grad(func = logDensOptim,
                         x = Par[dataSubIdx[iSubObs], jPar],
                         jPar = jPar,
                         MargisType = MargisType,
                         Mdl.Y = lapply(Mdl.Y, subfun, iSubObs = iSubObs,
                             dataSubIdx = dataSubIdx),
                         Mdl.par = rapply(Mdl.par, subfun, iSubObs = iSubObs,
                             dataSubIdx = dataSubIdx, how  = "replace"),
                         Mdl.u = Mdl.u[dataSubIdx[iSubObs], , drop = FALSE],
                         Mdl.d = Mdl.d[dataSubIdx[iSubObs], , drop = FALSE],
                         parUpdate = parUpdate,
                         MCMCUpdateStrategy = MCMCUpdateStrategy)
                   ,silent = TRUE)

            if(is(gradTry, "try-error"))
              {
                out[iSubObs, jPar] <- NA
              }
            else
              {
                out[iSubObs, jPar] <- gradTry
              }

          }
      }
    return(out)
  }
