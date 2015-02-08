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
##'       Current: Mon May 21 14:37:01 CEST 2012.
logCplLik <- function(u, CplNM, parCplRep, sum = TRUE)
{
    browser()
    parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)

    logCplDensObs <- dCpl(CplNM = CplNM, u = u,
                          parCpl = parCpl, log = TRUE)


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
