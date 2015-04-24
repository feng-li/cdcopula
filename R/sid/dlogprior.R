##' A collection of log densities for common priors.
##'
##' The parameters after "..." should be matched exactly.
##'
##' @name log_prior
##' @title logarithm density for priors
##'
##' @param beta "matrix" or "scaler".
##'         The parameter that need to be added with a prior. The prior density
##'         is calculated conditional on beta. beta can be a scaler of
##'         one-column matrix.
##' @param priorArgs "list".
##'         when priorArgs$type: is set to "mvnorm", you have to provide
##'         priorArgs$mean: "matrix", the mean of parameter, mu0 should be
##'         always an one-column matrix; priorArgs$covariance: "matrix", the
##'         covariance matrix. A g-prior can be constructed by setting it to
##'         X'X, where X is the covariates matrix.; priorArgs$shrinkage:
##'         "numeric", the shrinkage for the covariance.
##' @param link
##' @return "scalor". The log density of priors.
##'
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Mar 30 09:33:23 CEST 2010;
##'       Current:       Thu Dec 15 16:08:56 CET  2011.
##'
##' DEPENDS: mvtnorm::dmvnorm
##' TODO:
dlogprior <- function(beta, priArgs, link = NA)
{
  if (tolower(priArgs[["type"]]) == "mvnorm")
    {
      ## beta ~ N(mean, shrinkage*covariance)
      mean <- priArgs$mean

      covariance <- priArgs$covariance

      shrinkage <- priArgs$shrinkage
      out <- dmvnorm(t(beta), mean, covariance*shrinkage, log = TRUE)

    }
  else if (tolower(priArgs[["type"]]) == "norm-via-beta")
    {
      mean <- priArgs$mean # mean of beta density
      variance <- priArgs$variance # Covariates of beta density
      shrinkage <- priArgs$shrinkage # Shrinkage of beta density

      ## Transform to standard parametrization
      alpha <- - mean*(mean^2 -mean + variance*shrinkage)/
                 (covariance*shrinkage)
      beta <- (mean-1)^2*mean/(variance*shrinkage) + mean -1

      ## ## Assume this normal,
      if(tolower(link) == "logit")
        {
          ## See Beta-rep.nb
          meanLinked<- digamma(alpha)-digamma(beta)
          varLinked <- trigamma(alpha)+trigamma(beta)
          ## skewnessLinked <- (psigamma(alpha, 2)-psigamma(beta, 2))/
          ##   (trigamma(alpha)+trigamma(beta))^(3/2)
        }
      out <- sum(dnorm(beta, meanLinked, varLinked, log = TRUE))
    }
  else if(tolower(priArgs[["type"]]) == "bern")
    {
      ## The probability is recycled if necessary
      prob <- priArgs[["prob"]]
      betaLen <- length(beta)
      probVec <- matrix(prob, betaLen, 1)

      idx1 <- (beta == 1)
      idx0 <- (beta == 0)
      logDens <- sum(log(probVec[idx1])) + sum(log(1-probVec[idx0]))
      out <- logDens
    }

  else
    {
      stop("The prior is not implemented yet!")
    }

  return(out)
}

###----------------------------------------------------------------------------
### TESTING
###----------------------------------------------------------------------------
## a <- 1.1
## b <- 2.2
## n <- 10000
## mv <- rbeta(n, a, b)
## mv2 <- log(mv/(1-mv))
## mean(mv2)
## var(mv2)

## digamma(a)-digamma(b)
## trigamma(a)+trigamma(b)
