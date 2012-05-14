##' This function is used to compute the log gradient with respect to the
##' parameters in the marginal part. of copula model.
##'
##' This should be exactly the same as the usual way we do.
##' @title Log likelihood for marginal densities
##' @param parMargis
##' @param Mdl.Y
##' @param MargisTypes
##' @param chainCaller
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Fri Oct 21 18:17:01 CEST 2011;
##'       Current: Mon May 14 19:07:18 CEST 2012.
MargiModelGrad <- function(y, par, type, parCaller)
  {

    if(tolower(type) == "gaussian")
      {
        ## Subtract parameters and data
        mu <- par[["mu"]]
        sigma <- par[["sigma"]]

        ## Calculate the log density TODO: call from log likelihood?
        logMargiDens <- dnorm(y, mean = mu, sd = sigma, log = TRUE)

        if(tolower(parCaller) == "mu")
          {
            out <- -exp(logMargiDens)
          }
        else if(tolower(parCaller) == "sigma")
          {
            ## Calculate the fractions in the log form and transform back to
            ## avoid overflow/underflow.
            logGradFrac <- log(sigma^2-(x-mu)^2) - log(sigma) - log(x-mu)
            out <- exp(logMargiDens + logGradFrac)
          }
      }
    else if(tolower(margiType) == "student-t")
      {
        ## The marginal likelihood

      }
    else
      {
        stop("This marginal density is not implemented!")
      }

    ## The output
    return(out)
  }
