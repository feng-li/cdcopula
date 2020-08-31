#' This function make the log-likelihood function of the copula component
#'
#' Compute the copula likelihood of copula function. This function only update the copula
#' density part. The whole log-likelihood should consist of two parts: copula and
#' marginal parts.
#' @param Mdl.u NA
#' @param CplNM  NA
#' @param parCplRep NA
#' @param sum  NA
#' @param logLik "logical"
#' @return "matrix";
#' @references Joe 1997, p. 153
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @export
logCplLik <- function(Mdl.u, CplNM, parCplRep, sum = TRUE)
{
    parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)

    logCplDensObs <- dCplMixed(CplNM = CplNM, Mdl.u = Mdl.u, parCpl = parCpl, log = TRUE)
    ## browser()
    ## The output
    if(sum)
    {
        ## The sum of log copula density,  scaler
        out <- sum(logCplDensObs)
    }
    else
    {
        ## The log copula density,  vector
        out <- logCplDensObs
    }
    return(out)
}

#' @export
logCplLikParallel <- function(Mdl.u, CplNM, parCplRep, sum = TRUE)
{
    parCplRep <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)

    cl <- parallel:::defaultCluster()
    nObs <- nrow(Mdl.u[[1]])
    u = do.call(cbind, Mdl.u)
    dataSubIdxLst <- data.partition(nObs = nObs,
                                    list(N.subsets = length(cl),
                                         partiMethod = "ordered"))

    subfun <- function(index, data)data[index, , drop = FALSE]
    u.Lst <- lapply(dataSubIdxLst, subfun, data = u)

    splitlist <- function(data, index) lapply(data, subfun, index = index)
    parCpl.Lst <- rapply(dataSubIdxLst, splitlist,
                         data = parCplRep, how = "replace")

    ## Parallel is huge slow due to communication
    ## system.time(out0 <- MargiModelGrad(y = y, par = par, type = type,
    ##                                    parCaller = parCaller, densCaller = densCaller))
    logCplDensObs.Lst <- clusterMap(cl = cl, fun = dCpl,
                                    u = u.Lst,
                                    parCpl = parCpl.Lst,
                                    MoreArgs = list(CplNM = CplNM, log = TRUE))
    logCplDensObs <- do.call(rbind, logCplDensObs.Lst)

    if(sum)
    {
        ## The sum of log copula density,  scaler
        out <- sum(logCplDensObs)
    }
    else
    {
        ## The log copula density,  vector
        out <- logCplDensObs
    }
    return(out)
}
