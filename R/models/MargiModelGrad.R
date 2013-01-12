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
        phi <- par[["phi"]]

        ## Calculate the log density TODO: call from log likelihood?
        logMargiDens <- dnorm(y, mean = mu, sd = phi, log = TRUE)

        if(tolower(parCaller) == "mu")
          {
            out <- -exp(logMargiDens)
          }
        else if(tolower(parCaller) == "phi")
          {
            ## GradFrac <- (phi^2-(y-mu)^2)/((y-mu)*phi)
            ## out <- exp(logMargiDens)*GradFrac
            out <- -(y-mu)/phi*exp(logMargiDens)
          }
      }
    else if(tolower(margiType) == "splitt")
      {
        ## The marginal likelihood
        ## Literal translation from GSMMatlab code AsymStudT

        mu = par[["mu"]]  # location parameter
        df = par[["df"]] # Degrees of freedom
        phi = par[["phi"]] # Scaling Parameter
        lmd = par[["lmd"]] # Skewness Parameter

        ## PDF
        ## logMargiDens <- dsplitt(x = y, mu = mu, df = df, phi = phi, lmd = lmd,
        ##                         log = FALSE)

        if(tolower(parCaller) == "mu")
          {
            I0 <- (y<=mu)
            I <- (!I0)
            sign <- 1*I0 + lmd*I

            out <- -2*sign*sqrt(1/((y-mu)^2+sign^2*df*phi^2))*
              (sign^2*df*phi^2/
               ((y-mu)^2+sign^2*df*phi^2))^(df/2)/
                ((1+lmd)*beta(df/2, 1/2))

          }
        else if(tolower(parCaller) == "df")
          {
            stop("Not implemented yet")
          }
        else if(tolower(parCaller) == "phi")
          {
            I0 <- (y<=mu)
            I <- (!I0)
            sign <- 1*I0 + lmd*I

            out <- -2*sign*(y-mu)*sqrt(1/((y-mu)^2+sign^2*df*phi^2))*
              (sign^2*df*phi^2/((y-mu)^2+sign^2*df*phi^2))^(df/2)/
                ((1+lmd)*phi*beta(df/2, 1/2))

          }
        else if(tolower(parCaller) == "lmd")
          {
            stop("Not implemented yet")
          }
        else
          {
            stop("No such parameter!")
          }

      }
    else
      {
        stop("This marginal density is not implemented!")
      }

    ## The output
    return(out)

  }
