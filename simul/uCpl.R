##' Compute the Copula function for given margins
##'
##' The function is used for bivariate copula C(U1<u2, U2<u1)
##' @title Copula distribution
##' @param u 
##' @param theta 
##' @param copula 
##' @param par "list"
##'     Any additional parameters needed in the copula. In the
##' t-copula, par$df: indicates the degrees of freedom.
##' @return 
##' @references 
##' @author
##'     Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     DEPENDS: nvtnorm
##'     Created: Mon Sep 26 13:54:13 CEST 2011;
##'     Current: Mon Sep 26 13:55:12 CEST 2011.
uCpl <- function(u, theta, copula, par)
{
  if(tolower(copula) == "bb7")
    {
      p <- 2
      theta1 <- theta[1]
      theta2 <- theta[2]
      u1 <- u[, 1]
      u2 <- u[, 2]
      percentile <- 1 - (1 - ((1 - (1 - u)^theta1)^(-theta2) + (1 - (1 - 
         v)^theta1)^(-theta2) - 1)^(-1/theta2))^(1/theta1)
     out <- matrix(percentile)
    }
  else if(tolower(copula) == "gaussian")
    {
      u.quantile <- qnorm(u)
      p <- 2 # Hard code 
      corr <- matrix(c(1, theta, theta, 1), 2, 2) # The theta is the
                                        # correlation
      percentile <- apply(u.quantile, 1, function(x) pmvnorm(lower = c(-Inf, -Inf),
                          upper = x, corr = corr)[1])
      ## TODO: apply pmvnorm is very slow      
      out <- matrix(percentile)
    }
  else if(tolower(copula) == "mvt")
    {
      ##
      p <- ncol(u)
      df <- par[["df"]]
      ## theta is the vector for the lower triangular correlation matrix. P is
      ## the full correlation matrix.  TODO: Better construction of the
      ## correlation matrix.
      P.mat <- matrix(1, p, p) # The correlation matrix
      P.mat[lower.tri(P.mat)] <- theta
      P.mat[upper.tri(P.mat)] <- t(P.mat)[upper.tri(P.mat)]
      corr <- P.mat # The covariance matrix with scale 1.
      u.quantile <- qt(u, df = df)
      percentile <- apply(u.quantile, 1, function(x) pmvt(lower = rep(-Inf, p),
                          upper = x, corr = corr, df = df)[1])
      ## TODO: apply pmvt is very slow      
      out <- matrix(percentile)
    }
  
  else if(tolower(copula) == "fgm")
    {
      u1 <- u[, 1]
      u2 <- u[, 2]
      percentile <- u1*u2*(1+theta*(1-u1)*(1-u2))
      out <- matrix(percentile)
    }
  else if(tolower(copula) == "gumbel")
    {
      u.tilde <- -log(u)
      percentile <- exp(-rowSums(u.tilde^theta)^(1/theta))
      out <- matrix(percentile)
    }
  else stop("Copula name is not implemented.")
  return(out)  
}
