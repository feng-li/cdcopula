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

  tol <- .Machine$double.eps*1e8
  u.bad1 <- (u > 1-tol)
  u.bad0 <- (u < 0+tol)

  if(any(u.bad1))
    {
      u[u.bad1] <- 1-tol
      warning("u is to close to 1. Adjusted...",
              immediate. = TRUE)
    }
  if(any(u.bad0))
    {
      u[u.bad0] <- 0 +tol
      warning("u is to close to 0. Adjusted...",
              immediate. = TRUE)
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
      lambdaU <- as.vector(kendalltauInv(
          CplNM = CplNM, parRepCpl = parCpl,
          tauTabular = staticArgs[["tauTabular"]]))

      print(lambdaU)

      ## The standard copula parameters (recycled if necessary, should not have
      ## dimension attributed).

      delta <- -log(2)/log(lambdaL)
      theta <- log(2)/log(2-lambdaU)

      ## temporal data
      TC1 <- 1-(1-u)^theta

      ## Numeric check if L6 is too close to zero.
      ## tol <- .Machine$double.eps*1e8
      ## TC1.bad0 <- (TC1<tol)
      ## if(any(TC1.bad0))
      ##   {
      ##     TC1[TC1.bad0] <- tol
      ##     warning("Numerical unstable on TC1,  adjusted on the cliff...",
      ##             immediate. = TRUE)

      ##   }
      ## TC1.bad1 <- (TC1>(1-tol))
      ## if(any(TC1.bad1))
      ##   {
      ##     TC1[TC1.bad1] <- 1-tol
      ##     warning("Numerical unstable on TC1,  adjusted on the cliff...",
      ##             immediate. = TRUE)
      ##   }

      TC2 <- (1-u)^(-1+theta)


      ## Numeric check if L6 is too close to zero.
      ## tol <- .Machine$double.eps*1e8
      ## TC2.bad0 <- (TC2<tol)
      ## if(any(TC2.bad0))
      ##   {
      ##     TC2[TC2.bad0] <- tol
      ##     warning("Numerical unstable on TC2,  adjusted on the cliff...",
      ##             immediate. = TRUE)

      ##   }
      ## TC2.bad1 <- (TC2>(1-tol))
      ## if(any(TC2.bad1))
      ##   {
      ##     TC2[TC2.bad1] <- 1-tol
      ##     warning("Numerical unstable on TC2,  adjusted on the cliff...",
      ##             immediate. = TRUE)
      ##   }

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
