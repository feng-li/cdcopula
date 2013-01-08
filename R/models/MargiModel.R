##' Calculate the CDF and PDF of the marginal distribution.
##'
##' The detailed description can be found in the main setting file for each
##' input variable.
##' @param y "list".
##'        The response variables for the marginal model.
##' @param type "vector" with character input.
##'        The type of the marginal model.
##' @param par "list".
##'        The parameters input for the marginal model.
##' @return "matrix" with marginal names attributed.
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Jan 17 19:27:25 CET 2012;
##'       Current: Tue Jan 17 19:27:30 CET 2012.
MargiModel <- function(y, type, par)
  {
    if(tolower(type) == "gaussian")
      {
        ## The mean and standard deviation for Gaussian density
        mu <- par[["mu"]] # scaler
        sigma <- par[["sigma"]]   # scaler

        ## The percentile representation
        u <- pnorm(y, mean = mu, sd = sigma, log = FALSE)

        ## The quantile representation
        d <- dnorm(y, mean = mu, sd = sigma, log = TRUE)

        ## The output
        out <- list(u = u, d = d)
      }
    else if (tolower(margiType) == "splitt")
      {
        ## The marginal likelihood
        ## Literal translation from GSMMatlab code AsymStudT

        mu = par[["mu"]]  # location parameter
        df = par[["df"]] # Degrees of freedom
        phi = par[["phi"]]; # Scaling Parameter
        lmd = par[["lmd"]]; # Skewness Parameter

        ## CDF
        u <- psplitt(x = y, mu = mu, df = df, phi = phi, lmd = lmd,
                     log = FALSE)

        ## PDF
        d <- dsplitt(x = y, mu = mu, df = df, phi = phi, lmd = lmd,
                     log = FALSE)

        out <- list(u = u, d = d)
      }
    else
      {
        stop("This type of margin is not implemented.")
      }


    ## The output
    return(out)
  }
