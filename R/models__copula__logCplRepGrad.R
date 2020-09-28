#' Gradient for reparameterized log copula function
#'
#'
#' Log copula gradient
#' @param CplNM NA
#' @param Mdl.u  NA
#' @param parCplRep NA
#' @param parCaller  NA
#' @return NA
#' @references Li 2012
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @export
logCplRepGrad <- function(CplNM, Mdl.u, parCplRep, parCaller)
{

    CplNM0 <- strsplit(CplNM, split = "_")[[1]][1]

    if(all(sapply(Mdl.u, ncol)  == 1))
    {
        u <- do.call(cbind, Mdl.u)
    }
    else
    {
        stop("Mixed margins are not ready yet.")
    }

    if(tolower(CplNM0) %in% c("bb7", "sjc"))
    {
        ## The name of marginal model
        Mdl.MargisNM <- dimnames(u)[[2]]
        nObs <- nrow(Mdl.u[[1]])

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
    else if(tolower(CplNM0) == "mvt")
    {
        ## The name of marginal model
        Mdl.MargisNM <- dimnames(u)[[2]]
        nObs <- dim(u)[1]

        parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)
        df <- parCpl[["df"]] # n-by-1

        if(tolower(parCaller) %in% c("lambdal", "df"))
        { ## CopulaDensity-MVT.nb

            logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                         parCpl = parCpl, parCaller = c("df")) # n-by-1

            if(tolower(parCaller) == "df")
            {
                ## No reparameterization
                parRepGrad <- 1
            }
            else
            {
                lambdaGrad.par <- lambdaGrad(CplNM = CplNM, parCpl = parCpl,
                                             parCaller = c("df"))
                parRepGrad <- (1/lambdaGrad.par[["df"]])
            }

            ## The chain gradient
            out <- (logCplGrad.par[["df"]]*parRepGrad) # n-by-lq
        }
        else if(tolower(parCaller) %in% c("tau", "rho"))
        {
            logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                         parCpl = parCpl, parCaller = c("rho")) # n-by-lq

            if(tolower(parCaller) == "rho")
            {
                parRepGrad <- 1
            }
            else
            {
                kendalltauGrad.par <- kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
                                                     parCaller = "rho")
                parRepGrad <- (1/kendalltauGrad.par[["rho"]])
            }

            out <- logCplGrad.par[["rho"]]*parRepGrad # n-by-lq

            ## FIXME: for some reason, the analytical result is always 1/2 of the numerical
            ## result. need further verification.
        }
        else
        {
            ## The gradient with respect to u_i
            ## Reorder the parameters.
            logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u, parCpl = parCpl,
                                         parCaller = parCaller) # n-by-lq

            out <- logCplGrad.par[["u"]]
        }
    }
    else if(tolower(CplNM0) == "gumbel")
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
        else if(tolower(parCaller) == "lambdau")
        {## browser()

            logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                         parCpl = parCpl, parCaller = c("delta")) # n-by-1

            lambdaGrad.par <- lambdaGrad(CplNM = CplNM, parCpl = parCpl,
                                         parCaller = "delta")

            out <- logCplGrad.par[["delta"]]*(1/lambdaGrad.par[["delta"]])

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
    else if(tolower(CplNM0) == "clayton")
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
        else if(tolower(parCaller) == "lambdal")
        {
            ## browser()
            logCplGrad.par <- logCplGrad(CplNM = CplNM, u = u,
                                         parCpl = parCpl, parCaller = c("delta")) # n-by-1

            lambdaGrad.par <- lambdaGrad(CplNM = CplNM, parCpl = parCpl,
                                         parCaller = "delta")

            out <- logCplGrad.par[["delta"]]*(1/lambdaGrad.par[["delta"]])

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
#' @export
logCplRepGradParallel <- function(CplNM, Mdl.u, parCplRep, parCaller)
{
    cl <- parallel:::defaultCluster()

    nObs <- nrow(Mdl.u[[1]])
    # u = do.call(cbind, Mdl.u)

    dataSubIdxLst <- data.partition(nObs = nObs,
                                    list(N.subsets = length(cl), partiMethod = "ordered"))
    # browser()
    # subfun <- function(index, data)data[index, , drop = FALSE]
    # u.Lst <- lapply(dataSubIdxLst, subfun, data = u)

    splitlist <- function(data, index) lapply(data, subfun, index = index)
    Mdl.u.Lst <- rapply(dataSubIdxLst, splitlist, data = Mdl.u, how = "replace")
    parCplRep.Lst <- rapply(dataSubIdxLst, splitlist, data = parCplRep, how = "replace")

    ## Parallel code is far slow than serial code due to communications.
    ## system.time(out0 <- logCplRepGrad(u = u, parCplRep = parCplRep,
    ##                                   CplNM = CplNM, parCaller = parCaller))
    out.Lst <- clusterMap(cl, logCplRepGrad,
                          Mdl.u = Mdl.u.Lst,
                          parCplRep = parCplRep.Lst,
                          MoreArgs = list(CplNM = CplNM, parCaller = parCaller))

    out <- do.call(rbind, out.Lst)

    return(out)
}
