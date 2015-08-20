##' Gradient for reparameterized log copula function
##'
##'
##' @title Log copula gradient
##' @param CplNM
##' @param u
##' @param parCpl
##' @param parCaller
##' @return
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Fri May 11 12:42:20 CEST 2012; Current: Fri Mar 27 17:47:58 CST 2015.
logCplRepGrad <- function(CplNM, u, parCplRep, parCaller)
{
  if(tolower(CplNM) == "bb7")
  {
    ## The name of marginal model
    MargisNM <- dimnames(u)[[2]]
    nObs <- dim(u)[1]

    parCpl <- parCplRep2Std(CplNM, parCplRep)

    if(tolower(parCaller) == "lambdal")
    {
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = c("delta")) # n-by-1

      lambdaGrad.par <- lambdaGrad(CplNM = CplNM, parCpl = parCpl,
                                       parCaller = "delta")

      out <- logCplGrad.par[["delta"]]*(1/lambdaGrad.par[["delta"]])
    }
    else if(tolower(parCaller) == "lambdau")
    {
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = c("theta")) # n-by-1

      lambdaGrad.par <- lambdaGrad(CplNM = CplNM, parCpl = parCpl,
                                   parCaller = "theta")

      out <- logCplGrad.par[["theta"]]*(1/lambdaGrad.par[["theta"]])
    }

    else if(tolower(parCaller) == "tau")
    {
      ## logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
      ##                              parCpl = parCpl, parCaller = c("theta", "delta")) # n-by-1

      ## ## Gradient w.r.t. tau
      ## kendalltauGrad.par <- kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
      ##                                      parCaller = c("theta", "delta"))

      ## ## The chain gradient for complex functions
      ## out <- (logCplGrad.par[["theta"]]*(1/kendalltauGrad.par[["theta"]]) +
      ##         logCplGrad.par[["delta"]]*(1/kendalltauGrad.par[["delta"]]))



      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = c("theta")) # n-by-1

      ## Gradient w.r.t. tau
      kendalltauGrad.par <- kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
                                           parCaller = c("theta"))

      ## The chain gradient for complex functions
      out <- (logCplGrad.par[["theta"]]*(1/kendalltauGrad.par[["theta"]]))

    }
    else
    {
      ## Gradient w.r.t u_i
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = parCaller) # n-by-1
      out <- logCplGrad.par[["u"]]
    }
  }
  else if(tolower(CplNM) == "mvt")
  {
    ## The name of marginal model
    MargisNM <- dimnames(u)[[2]]
    nObs <- dim(u)[1]

    parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)
    df <- parCpl[["df"]] # n-by-1
    ## rho <- parCpl[["rho"]] # n-by-lq

    u.quantile <- qt(u, df)
    if(tolower(parCaller) == "lambdal")
    { ## CopulaDensity-MVT.nb

      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = c("df")) # n-by-1

      lambdaGrad.par <- lambdaGrad(CplNM = CplNM, parCpl = parCpl,
                                   parCaller = c("df"))

      ## The chain gradient
      out <- (logCplGrad[["df"]]*(1/lambdaGrad.par[["df"]])) # n-by-lq

    }
    else if(tolower(parCaller) == "tau")
    {
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = c("rho")) # n-by-lq

      kendalltauGrad.par <- kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
                                           parCaller = "rho")

      out <- 2*logCplGrad.par[["rho"]]*(1/kendalltauGrad.par[["rho"]]) # n-by-lq

      ## FIXME: for some reason, the analytical result is always 1/2 of the numerical
      ## result. need further verification.

    }
    else
    {
      ## The gradient with respect to u_i
      ## Reorder the parameters.
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = parCaller) # n-by-lq

      out <- logCplGrad.par[["u"]]
    }
  }
  else if(tolower(CplNM) == "gumbel")
  {
    parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)

    if(tolower(parCaller) == "tau")
    {## browser()
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = c("delta")) # n-by-lq

      kendalltauGrad.par <- kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
                                           parCaller = "delta")
      out <- logCplGrad.par[["delta"]]*(1/kendalltauGrad.par[["delta"]])
    }
    else
    {
      ## The gradient with respect to u_i
      ## Reorder the parameters.
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = parCaller) # n-by-lq

      out <- logCplGrad.par[["u"]]
    }
  }
  else if(tolower(CplNM) == "clayton")
  {
    parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)
    if(tolower(parCaller) == "tau")
    {## browser()
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = c("delta")) # n-by-lq

      kendalltauGrad.par <- kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
                                           parCaller = "delta")
      out <- logCplGrad.par[["delta"]]*(1/kendalltauGrad.par[["delta"]])
    }
    else
    {
      ## The gradient with respect to u_i
      ## Reorder the parameters.
      logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                   parCpl = parCpl, parCaller = parCaller) # n-by-lq

      out <- logCplGrad.par[["u"]]
    }
  }
  else
  {
    stop("No such copula implemented!")
  }

  return(out)
}



## Parallel version to gradients
logCplRepGradParallel <- function(CplNM, u, parCplRep, parCaller)
{
  cl <- parallel:::defaultCluster()
  nObs <- nrow(u)

  dataSubIdxLst <- data.partition(nObs = nObs,
                                  list(N.subsets = length(cl), partiMethod = "ordered"))

  subfun <- function(index, data)data[index, , drop = FALSE]

  u.Lst <- lapply(dataSubIdxLst, subfun, data = u)

  splitlist <- function(data, index) lapply(data, subfun, index = index)
  parCplRep.Lst <- rapply(dataSubIdxLst, splitlist, data = parCplRep, how = "replace")

  out.Lst <- clusterMap(cl, logCplRepGrad,
                        u = u.Lst,
                        parCplRep = parCplRep.Lst,
                        MoreArgs = list(CplNM = CplNM,
                                        parCaller = parCaller))
  out <- do.call(rbind, out.Lst)

  return(out)
}
