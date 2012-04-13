##' Generate random variables using copulas
##'
##' The copula expression can be find at Trivedi and Zimmer 2005
##' @title Random *uniform* variables generator by bivariate copula
##' @param n
##' @param theta
##' @param copula
##' @param par "list"
##'         Any additional parameters input for the copula. For example in the
##'         t copula, par$df: the degrees of freedom is needed.
##' @return
##' @references
##'     Trivedi and Zimmer 2005
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     DEPENDS: mvtnorm
##'     Created: Mon Sep 26 09:43:12 CEST 2011;
##'     Current: Mon Sep 26 09:43:21 CEST 2011.
ruCpl <- function(n, parCpl, copula)
{
  if(tolower(copula) == "bb7") # Joe 1997
    {
      ## The copula parameters are recycled if necessary
      theta <- matrix(parCpl[[1]], n, 1)
      delta <- matrix(parCpl[[2]], n, 1)

      ## Parameters healthy conditions
      hCond <- (theta[1] >= 1) && (theta[2]>0)
      if(!hCond) stop("BB7 copula should have: theta >= 1, delta >0.")

      p <- 2 #Hard code for bivariate copula
      v <- matrix(runif(n*p, 0, 1), n, p)
      u <- matrix(NA, n, p)
      u[, 1] <- v[, 1]

      ## No closed form for inverse of conditional function
      ## TODO: Potential overflow when theta is big. Consider using Chris Sims solver.
      for(i in 1:n)
        {
          out.cur <- uniroot(function(x, theta, v0)
                             {

                               theta1 = theta[1]
                               theta2 = theta[2]
                               u1 <- v0[1]
                               v2 <- v0[2]
                               u2 <- x

                               L1 <- 1-(1-u1)^theta1
                               L2 <- 1-(1-u2)^theta1
                               L5 <- -1 + L1^(-theta2) + L2^(-theta2)
                               L6 <- 1- L5^(-1/theta2)
                               out <- L1^(-1-theta2)*L5^(-1-1/theta2)*
                                 L6^(-1+1/theta1)*(1-u1)^(-1+theta1)-v2
                               return(out)
                             },
                             interval = c(0, 1),
                             theta = theta, v0 = v[i, ])
          u[i, 2] <- out.cur[["root"]]
        }

      ## The theoretical upper and lower tail dependence parameters
      lambdaU <- 2-2^(1/theta[1])
      lambdaL <- 2^(-1/theta[2])

      ## Kendall's tau,  empirical
      emptau <- cor(u, method = "kendall")[2]

      ## Kendall's tau,  theoretical
      theta1 = theta[1] # theta
      theta2 = theta[2] # delta

      tol <- 0.0001
      if(theta1<(2-tol))
        {
          theotau <- 1-2/(theta2*(2-theta1))+4/(theta1^2*theta2)*beta(theta2+2, 2/theta1-1)
        }
      else if (theta1>(2+tol))
        {
          theotau <- 1-2/(theta2*(2-theta1))-
            4*pi/(theta1^2*theta2*(2+theta2)*sin(2*pi/theta1)*beta(2-2/theta1, 1+theta2+2/theta1))
        }
      else
        {
          theotau <- 1 + 1/theta2*(1+digamma(1)-digamma(2+theta2))
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
  else if(tolower(copula) == "mvt") ## the multivariate t demarta mcneil (2005)
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
      out <- u
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
  else stop("Given copula name is not implemented.")
  return(out)
}
