##' This function make the log-likelihood function of the copula component
##'
##' This function only update the copula density part. The whole log-likelihood
##' should consist of two parts: copula and marginal parts.
##' @title Compute the copula likelihood of copula function
##' @param u
##' @param CplNM
##' @param parCpl
##' @param staticArgs
##' @param logLik "logical"
##' @return "matrix";
##' @references Joe 1997, p. 153
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Oct 20 18:15:13 CEST 2011;
##'       Current: Mon May 21 14:37:01 CEST 2012.
logCplLik <- function(u, CplNM, parCpl, staticArgs, logLik = TRUE)
{

###----------------------------------------------------------------------------
### Copula likelihood numerical correction if u -> 0 or 1
###----------------------------------------------------------------------------

  ## Fix u on the cliff if any u -> 0 or u -> 1.
  ## Thanks to the advice from M. Smith

  tol <- .Machine$double.eps*1e3
  u.bad1 <- (u > 1-tol)
  u.bad0 <- (u < 0+tol)

  if(any(u.bad1))
    {
      u[u.bad1] <- 1-tol
    }
  if(any(u.bad0))
    {
      u[u.bad0] <- 0 +tol
    }

###----------------------------------------------------------------------------
### Compute the copula likelihood
###----------------------------------------------------------------------------

  ## The sum of log copula density
  if(tolower(CplNM) == "bb7")
    {
      ## The name of marginal model
      MargisNM <- names(u)

      ## Subtract the parameters list.
      tau <- as.vector(parCpl[["tau"]])
      lambdaL <- as.vector(parCpl[["lambdaL"]])
      lambdaU <- as.vector(kendalltauInv(CplNM = CplNM, parRepCpl = parCpl,
                               tauTabular = staticArgs[["tauTabular"]]))
      ## The standard copula parameters (recycled if necessary, should not have
      ## dimension attributed).

      delta <- -log(2)/log(lambdaL)
      theta <- log(2)/log(2-lambdaU)

      ## temporal data
      ## L1 <- 1-(1-u1)^theta
      ## L2 <- 1-(1-u2)^theta
      TC1 <- 1-(1-u)^theta

      ## L3 <- (1-u1)^(-1+theta)
      ## L4 <- (1-u2)^(-1+theta)
      TC2 <- (1-u)^(-1+theta)
      ## L3 <- TC2[, 1]
      ## L4 <- TC2[, 2]

      ## L5 <- -1+L1^(-delta)+L2^(-delta)
      L5 <- rowSums(TC1^(-delta)) - 1

      L6 <- 1-L5^(-1/delta) # FIXME: log(L6)->Inf when u->1,  v->1.

      logCplDensObs <- (-1-delta)*rowSums(log(TC1))+
        rowSums(log(TC2))-
          2*(1+delta)/delta*log(L5)+
            (-2+1/theta)*log(L6)+
              log(-1+theta+L5^(1/delta)*L6*(1+delta)*theta)
    }
  else if(tolower(CplNM) == "gaussian")
    {
      ## The Gaussian copula
    }
  else
    {
      stop("The copula is not implemented yet!")
    }


  ## The output
  if(logLik)
    {
      ## The sum of log marginal density,  scaler
      out <- sum(logCplDensObs)
    }
  else
    {
      ## The log marginal density,  vector
      out <- logCplDensObs
    }
  return(out)
}
