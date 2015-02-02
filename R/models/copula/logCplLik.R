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
logCplLik <- function(u, CplNM, parCpl, sum = TRUE)
{

###----------------------------------------------------------------------------
### Copula likelihood numerical correction if u -> 0 or 1
###----------------------------------------------------------------------------

    ## Fix u on the cliff if any u -> 0 or u -> 1.
    ## Thanks to the advice from M. Smith

    ## Debugging symbol: if the warning should be printed out immediately.
    immediate. <- FALSE


###----------------------------------------------------------------------------
### Compute the copula likelihood
###----------------------------------------------------------------------------

    ## The sum of log copula density
    if(tolower(CplNM) == "bb7")
        {
            ## Subtract the parameters list.
            tau <- as.vector(parCpl[["tau"]])
            lambdaL <- as.vector(parCpl[["lambdaL"]])
            lambdaU <- as.vector(kendalltauInv(
                CplNM = CplNM, parRepCpl = parCpl))

            delta <- -log(2)/log(lambdaL)
            theta <- log(2)/log(2-lambdaU)

            logCplDensObs <- cCpl(CplNM = CplNM, u = u,
                                  theta = theta, delta = delta, log = TRUE)


        }
    else if(tolower(CplNM) == "gaussian")
        {
            ## The Gaussian copula
        }
    else if((tolower(CplNM) == "mvt"))
        {
            tau <- as.vector(parCpl[["tau"]])
            lambda <- as.vector(parCpl[["lambda"]])

            rho <- sin(tau*pi/2)
            df <- as.vector(lambdaInv(
                CplNM = CplNM, parRepCpl = parCpl))

            logCplDensObs <- cCpl(CplNM = CplNM, u = u,
                                  df = df, rho = rho, log = TRUE)
        }

    else
        {
            stop("The copula is not implemented yet!")
        }

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
