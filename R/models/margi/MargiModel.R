##' CDF and PDF of the marginal distribution.
##'
##' The detailed description can be found in the main setting file for each
##' input variable.
##' @param y "vector".
##'        The response variables for the marginal model.
##' @param type "vector" with character input.
##'        The type of the marginal model.
##' @param par "list".
##'        The parameters input for the marginal model.
##' @return "list".
##'
##'        Return the percentile and log density of the marginal models.
##'
##'
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Jan 17 19:27:25 CET 2012;
##'       Current: Tue Jan 17 19:27:30 CET 2012.
MargiModel <- function(y, type, par)
  {
    if(tolower(type) == "gaussian")
      {
        ## The mean and standard deviation for Gaussian density
        mu <- par[["mu"]] # mean parameter
        phi <- par[["phi"]]   # standard deviation

        ## The percentile representation
        u <- pnorm(y, mean = mu, sd = phi, log = FALSE)

        ## The quantile representation
        d <- dnorm(y, mean = mu, sd = phi, log = TRUE)

        ## The output
        out <- list(u = u, d = d)
      }
    else if (tolower(type) == "splitt")
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
                     log = TRUE)

        out <- list(u = u, d = d)
      }
    else
      {
        stop("This type of margin is not implemented.")
      }


    ## The output
    return(out)
  }
