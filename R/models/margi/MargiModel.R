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
##'       Current: Mon Jan 05 16:24:51 CST 2015.
MargiModel <- function(y, type, par, densCaller = c("u", "d"))
{
  ## The out storage
  out <- list()

  if(tolower(type) %in% c("gaussian", "garch-normal"))
  {
    ## The mean and standard deviation for Gaussian density
    mu <- par[["mu"]] # mean parameter
    phi <- par[["phi"]]   # standard deviation

    ## The percentile representation
    if("u" %in% tolower(densCaller))
    {
      u <- pnorm(y, mean = mu, sd = phi, log = FALSE)
      out[["u"]] <- u
    }
    ## The quantile representation
    if("d" %in% tolower(densCaller))
    {
      d <- dnorm(y, mean = mu, sd = phi, log = TRUE)
      out[["d"]] <- d
    }
  }
  else if (tolower(type)  == "splitt")
  {
    ## The marginal likelihood
    ## Literal translation from GSMMatlab code AsymStudT

    mu <- par[["mu"]]  # location parameter
    df <- par[["df"]] # Degrees of freedom
    phi <- par[["phi"]]; # Scaling Parameter
    lmd <- par[["lmd"]]; # Skewness Parameter

    ## CDF
    if("u" %in% tolower(densCaller))
    {
      u <- psplitt(x = y, mu = mu, df = df, phi = phi, lmd = lmd,
                   log = FALSE)
      out[["u"]] <- u
    }
    ## PDF

    if("d" %in% tolower(densCaller))
    {
      d <- dsplitt(x = y, mu = mu, df = df, phi = phi, lmd = lmd,
                   log = TRUE)

      out[["d"]] <- d
    }
  }
  else
  {
    stop("This type of margin is not implemented.")
  }


  ## The output
  return(out)
}
