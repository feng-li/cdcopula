##' This function is used to compute the log gradient with respect to the
##' parameters in the marginal part. of copula model.
##'
##' This should be exactly the same as the usual way we do.
##' @title Log likelihood for marginal densities
##' @param y "vector"
##' @param par "list"
##' @param type "list" Type of the marginal model
##' @param parCaller "character string" Indicates which parameter is called.
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
    else if(tolower(type) == "splitt")
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

            browser()
          }
        else if(tolower(parCaller) == "df")
          {
            I0 <- (y<=mu)
            I <- (!I0)
            sign <- 1*I0 + lmd*I
            sign2 <- I0*(-1) + I*1

            A <- cbind(1/2, df/2, df/2)
            B <- cbind(1+df/2, 1+df/2)
            Z <- (df*phi^2*sign^2)/((y-mu)^2+sign^2*df*phi^2)

            out <- (sign/(2*(1+lmd)*df^2*beta(df/2, 1/2)))*
              (
                  sign2*4*Z^(df/2)*
                  ghypergeo(A, B, Z)+
                  (
                      df*(
                          -2*(y-mu)*sqrt(1/((y-mu)^2+sign^2*df*phi^2))*Z^(df/2)-
                          sign2*df*ibeta(x = Z, a = df/2, b = 1/2)*
                          (log(Z)-digamma(df/2)+digamma((1+df)/2))
                          )
                      )
                  )
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
            I0 <- (y<=mu)
            I <- (!I0)

            ## Reserve the storage
            out <- mu
            out[0:length(mu)] <- NA

            if(any(I0))
              {
                y0 <- y[I0]
                mu0 <- mu[I0]
                phi0 <- phi[I0]
                df0 <- df[I0]
                lmd0 <- lmd[I0]

                A0 <- df0*phi0^2/((y0-mu0)^2+df0*phi0^2)
                out0 <- -ibeta(A0, df0/2, 1/2, reg = TRUE)/(1+lmd0)^2
                out[I0] <- out0
              }
            if(any(I))
              {
                y1 <- y[1]
                mu1 <- mu[I]
                phi1 <- phi[I]
                df1 <- df[I]
                lmd1 <- lmd[I]

                B1 <- ((y1-mu1)^2+df1*phi1^2*lmd1^2)
                A1 <- df1*phi1^2*lmd1^2/B1

                out1 <- -(2*(1+lmd1)*(y1-mu1)*sqrt(1/B1)*
                    A1^(df1/2)+ibeta(A1, df1/2, 1/2))/((1+lmd)^2*beta(df1/2, 1/2))
                out[I] <- out1
              }
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
