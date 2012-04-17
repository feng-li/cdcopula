##' Generate random variables using copulas
##'
##' The copula expression can be find at Trivedi and Zimmer 2005
##' @title Random *uniform* variables generator by bivariate copula
##' @param n
##' @param copula
##' @param parCpl "list" Any additional parameters input for the copula. For
##' example in the t copula, par$df: the degrees of freedom is needed.
##' @return "list" The random variables with some detailed information
##' @references Trivedi and Zimmer 2005
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     DEPENDS: mvtnorm
##'     Created: Mon Sep 26 09:43:12 CEST 2011;
##'     Current: Mon Apr 16 13:40:47 CEST 2012.
ruCpl <- function(n, parCpl, copula, exArgs)
{
  if(tolower(copula) == "bb7") # Joe 1997
    {
      ## Subtract the parameters list.
      tau <- parCpl[["tau"]]
      lambdaL <- parCpl[["lambdaL"]]

      ## Transform the parameters into the standard form
      ## FIXME: Consider to speed it up if it is really slow
      parOut <- kendalltauInv(CplNM = copula, parRepCpl = parCpl,
                              tauTabular = exArgs[["tauTabular"]])

      ## The standard copula parameters (recycled if necessary).
      delta <- as.vector(parOut[["delta"]])
      theta <- as.vector(parOut[["theta"]])

      ## Parameters healthy conditions TODO: more serious check
      ## hCond <- (theta[1] >= 1) && (theta[2]>0)
      ## if(!hCond) stop("BB7 copula should have: theta >= 1, delta >0.")

      p <- 2 #Hard code for bivariate copula
      v <- matrix(runif(n*p, 0, 1), n, p)
      u <- matrix(NA, n, p)
      u[, 1] <- v[, 1]

      ## The conditional method for sampling copula
      TC1 <- 1-(1-v)^theta
      L1 <- TC[, 1, drop = FALSE]
      TC2 <- (1-v)^(-1+theta)
      L5 <- -1 + rowSums(TC1^(-delta))
      L6 <- 1-L5^(-1/delta) # FIXME: log(L6)->Inf when u->1,  v->1.

      logv <- -(1+delta)*log(L1) - (1+1/delta)*log(L5) +
        (-1+1/theta)*log(L6) + (-1+theta)*log(1-v[, 1, drop = FALSE])

      u[, 2] <- exp(logv)

      ## The theoretical upper and lower tail dependence parameters
      lambdaU <- 2-2^(1/theta[1])
      lambdaL <- 2^(-1/theta[2])

      ## Kendall's tau,  empirical
      emptau <- cor(u, method = "kendall")[2]

      ## Kendall's tau,  theoretical

      tol <- 0.0001
      if(theta<(2-tol))
        {
          theotau <- 1-2/(delta*(2-theta))+4/(theta^2*delta)*beta(delta+2, 2/theta-1)
        }
      else if (theta>(2+tol))
        {
          theotau <- 1-2/(delta*(2-theta))-
            4*pi/(theta^2*delta*(2+delta)*sin(2*pi/theta)*beta(2-2/theta, 1+delta+2/theta))
        }
      else
        {
          theotau <- 1 + 1/delta*(1+digamma(1)-digamma(2+delta))
        }

      out <- list(u = u, lambdaU = lambdaU, lambdaL = lambdaL,
                  emptau = emptau, theotau = theotau)
    }
  else if(tolower(copula) == "gaussian")
    {
      ## Elliptical sampling,  Trivedi and Zimmer 2005 (p. 109)
      p <- 2 #Hard code for bivariate copula
      v <- matrix(rnorm(n*p, 0, 1), n, p)
      v1 <- v[, 1]
      v2 <- v[, 2]

      y <- matrix(NA, n, p)
      y[, 1] <- v1
      y[, 2] <- v1*theta + v2*sqrt(1-theta^2)
      u <- pnorm(y)
      out <- u
    }
  else if(tolower(copula) == "mvt") ## the multivariate t Demarta Mcneil (2005)
    {
      p <- (1+sqrt(8*length(theta)+1))/2 # p is obtained from the lower
                                        # triangular matrix.
      df <- par[["df"]]
      ## theta is the vector for the lower triangular correlation matrix. P is
      ## the full correlation matrix.  TODO: Better construction of the
      ## correlation matrix.
      P.mat <- matrix(1, p, p) # The correlation matrix
      P.mat[lower.tri(P.mat)] <- theta
      P.mat[upper.tri(P.mat)] <- t(P.mat)[upper.tri(P.mat)]
      corr <- P.mat # The covariance matrix with scale 1.

      x <- matrix(rmvt(n*p, df = df, sigma = corr), n, p)
      u <- pt(x, df = df) # The percentile

      ## Kendall's tau,  empirical
      emptau <- cor(u, method = "kendall")[2]

      ## Kendall's tau,  theoretical
      theotau <- 2/pi*asin(corr)

      out <- list(u = u, emptau = emptau, theotau = theotau)

    }
  else if(tolower(copula) == "fgm")
    {
      if(!(theta >= -1 && theta <= 1))
        {stop("FGM copula should have -1 <= theta <= 1.")}

      ## Conditional sampling Trivedi and Zimmer 2005(p. 108)
      p <- 2 #Hard code for bivariate copula
      v <- matrix(runif(n*p, 0, 1), n, p)
      v1 <- v[, 1]
      v2 <- v[, 2]

      u <- matrix(NA, n, p)
      u1 <- v1
      u[, 1] <- u1

      ## A <- theta*(2*u1-1)
      ## B <- (1-A)^2 + 4*v2*A
      ## u2 <- 2*v2/(sqrt(B)-A) TODO: Wrong formula in book.
      u2 <- v2 + (2*u1-1)*(v2-1)*v2*theta

      u[, 2] <- u2
      out <- u
    }
  else if(tolower(copula) == "gumbel")
    {
      ## Mixture of power simulation,  Trivedi and Zimmer 2005(p. 110)
     alpha <- 1/theta
     gamma <- rps(n, alpha)
     v <-matrix(runif(2*n, 0, 1), n, 2)
     t <- -1/gamma*log(v)
     u <- exp(-t^(1/theta))
     out <- u
    }
  else stop("Not implemented for given copula name.")
  return(out)
}
