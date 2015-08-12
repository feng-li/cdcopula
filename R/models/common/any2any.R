##' Approximate translate a distribution to another one through link functions.
##'
##' User need to input distribution and output distribution and linkage
##' type.
##' @title Translate distributions.
##' @param densArgs "list" arguments for densities
##' @param linkType "character" Type of linkage.
##' @return "list" features of the destiny distribution,  e.g. mean, and variance.
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Jan 12 00:05:00 CET 2012;
##'       Current: Mon May 14 15:15:39 CEST 2012.
any2any <- function(densArgs, linkArgs)
{
  ## Input distribution type
  inType <- densArgs[["input"]][["type"]]

  ## Output distribution type
  outType <- densArgs[["output"]][["type"]]

  linkType <- linkArgs[["type"]]
###----------------------------------------------------------------------------
### The input distribution
###----------------------------------------------------------------------------

  if(tolower(inType) == "beta")
  {
    ## Beta distribution
    ## Subtract the parameter input in a convenient way
    mean <- densArgs[["input"]][["mean"]] # mean of beta density
    variance <- densArgs[["input"]][["variance"]] # Covariates
    ## shrinkage <- densArgs[["shrinkage"]] # Shrinkage

    if(variance >= mean*(1-mean))
    {
      stop("Mean and variance should be in intervals (0, 1) and (0, mean(1-mean)),  respectively.")
    }

    ## Transform to standard parametrization
    alpha <- - mean*(mean^2 -mean + variance)/(variance)
    beta <- (mean-1)^2*mean/(variance) + mean -1

    if(tolower(linkType) == "logit")
    {
      ## Assume this normal, See Beta-rep.nb
      meanLinked<- digamma(alpha)-digamma(beta)
      varLinked <- trigamma(alpha)+trigamma(beta)
    }
    else
    {
      stop("This link function has not been implemented yet!")
    }
  }
  else if(tolower(inType) == "gbeta")
  {
    ## The generalized Beta distribution in interval [a, b]
    ## Subtract the parameter input in a convenient way
    mean <- densArgs[["input"]][["mean"]] # mean of beta density
    variance <- densArgs[["input"]][["variance"]] # Covariates
    a <- densArgs[["input"]][["a"]]
    b <- densArgs[["input"]][["b"]]

    mean0 <- (mean-a)/(b-a)
    variance0 <- variance/(b-a)^2

    if(variance0 >= mean0*(1-mean0))
    {
      stop("Mean and variance should be in intervals \n (a, b) and (0, (mean-a)(b-mean)),  respectively.")
    }

    ## shrinkage <- densArgs[["shrinkage"]] # Shrinkage

    ## Transform to standard parametrization
    alpha0 <- -mean0*(mean0^2 -mean0 + variance0)/(variance0)
    beta0 <- (mean0-1)^2*mean0/(variance0) + mean0 -1

    if(tolower(linkType) == "glogit")
    {
      ## Assume this normal, See Beta-rep.nb
      meanLinked<- digamma(alpha0)-digamma(beta0)
      varLinked <- trigamma(alpha0)+trigamma(beta0)
    }
    else
    {
      stop("This link function has not been implemented yet!")
    }
  }
  else if(tolower(inType) == "lognorm")
  {
    ## Lognormal distribution.
    ## Subtract the parameter input in a convenient way
    mean <- densArgs[["input"]][["mean"]] # mean of beta density
    variance <- densArgs[["input"]][["variance"]] # Covariates

    ## Transform to standard parametrization
    sigma2 <- log(1+variance/mean^2)
    mu <- log(mean) - sigma2/2

    if(tolower(linkType) == "log")
    {
      meanLinked<- mu
      varLinked <- sigma2
    }
    else
    {
      stop("Not implemented yet!")
    }

  }
  else if(tolower(inType) == "glognorm")
  {
    ## Lognormal distribution.
    ## Subtract the parameter input in a convenient way
    mean <- densArgs[["input"]][["mean"]] # mean of beta density
    variance <- densArgs[["input"]][["variance"]] # Covariates
    a <- densArgs[["input"]][["a"]]

    mean0 <- mean-a
    variance0 <- variance

    ## Transform to standard parametrization
    sigma2 <- log(1+variance0/mean0^2)
    mu <- log(mean0) - sigma2/2

    if(tolower(linkType) == "glog")
    {
      meanLinked<- mu
      varLinked <- sigma2
    }
    else
    {
      stop("Not implemented yet!")
    }
  }
  else if(tolower(inType) == "norm")
  {
    ## Lognormal distribution.
    ## Subtract the parameter input in a convenient way
    mean <- densArgs[["input"]][["mean"]] # mean of beta density
    variance <- densArgs[["input"]][["variance"]] # Covariates
    ## shrinkage <- densArgs[["shrinkage"]] # Shrinkage

    ## Transform to standard parametrization
    mu <- mean
    sigma2 <- variance

    if(tolower(linkType) == "identity")
    {
      meanLinked<- mu
      varLinked <- sigma2
    }
    else
    {
      stop("Not implemented yet!")
    }

  }

###----------------------------------------------------------------------------
### The output distribution
###----------------------------------------------------------------------------

  if(tolower(outType) == "norm")
  {
    out <- list(mean = meanLinked,
                variance = varLinked)
  }
  else
  {
    stop("Not implemented yet!")
  }

###----------------------------------------------------------------------------
### The final output
###----------------------------------------------------------------------------
  return(out)

}
