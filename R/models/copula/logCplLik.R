##' This function make the log-likelihood function of the copula component
##'
##' This function only update the copula density part. The whole log-likelihood
##' should consist of two parts: copula and marginal parts.
##' @title Compute the copula likelihood of copula function
##' @param u
##' @param CplNM
##' @param parCpl
##' @param logLik "logical"
##' @return "matrix";
##' @references Joe 1997, p. 153
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Oct 20 18:15:13 CEST 2011;
##'       Current: Sun Nov 08 20:10:52 CST 2015.
logCplLik <- function(u, CplNM, parCplRep, sum = TRUE)
{
    parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)

    logCplDensObs <- dCpl(CplNM = CplNM, u = u, parCpl = parCpl, log = TRUE)
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

logCplLikParallel <- function(u, CplNM, parCplRep, sum = TRUE)
{
    parCplRep <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)

    cl <- parallel:::defaultCluster()
    nObs <- nrow(u)
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
