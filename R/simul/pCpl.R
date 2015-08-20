##' Compute the Copula function for given margins
##'
##' The function is used for bivariate copula C(U1<u2, U2<u1)
##' @title Copula distribution
##' @param u "matrix".
##' @param CplNM  "character".
##'
##'        The copula name.
##'
##' @param parCpl "list"
##'     Any additional parameters needed in the copula. In the
##' t-copula, parCpl$df: indicates the degrees of freedom.
##' @param copula "character"
##' @return "vector"
##' @references Nelsen 2006
##' @author
##'     Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     DEPENDS: mvtnorm
##'     Created: Mon Sep 26 13:54:13 CEST 2011;
##'     Current: Tue Jul 14 13:06:55 CST 2015.
pCpl <- function(u, CplNM, parCpl, log = FALSE)
{
  if(tolower(CplNM) == "bb7")
  {
    p <- 2
    theta <- parCpl[["theta"]]
    delta <- parCpl["delta"]
    u1 <- u[, 1]
    u2 <- u[, 2]
    percentile <- (1 - (1 - ((1 - (1 - u1)^theta)^(-delta) +
                                                  (1 - (1 -  u2)^theta)^(-delta) - 1)^(-1/delta))^(1/theta))
    out <- matrix(percentile)
  }
  else if(tolower(CplNM) == "gaussian")
  {
    theta <- parCpl[["theta"]]
    u.quantile <- qnorm(u)
    p <- 2 # Hard code
    corr <- matrix(c(1, theta, theta, 1), 2, 2) # The theta is the
                                        # correlation
    percentile <- apply(u.quantile, 1, function(x) pmvnorm(lower = c(-Inf, -Inf),
                                                           upper = x, corr = corr)[1])
    ## TODO: apply pmvnorm is very slow
    out <- matrix(percentile)
  }
  else if(tolower(CplNM) == "mvt")
  {
    p <- ncol(u)
    df <- parCpl[["df"]]
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
  else if(tolower(CplNM) == "fgm")
  {
    u1 <- u[, 1]
    u2 <- u[, 2]
    percentile.log <- log(u1)+log(u2)+log(1+theta*(1-u1)*(1-u2))
    out.log <- matrix(percentile.log)
    if(log)
    {
      out <- out.log
    }
    else
    {
      out <- exp(out.log)
    }

  }
  else if(tolower(CplNM) == "gumbel")
  {# Joe 1997. p.142 Family B6
    delta <- as.vector(parCpl[["delta"]])
    u.tilde <- -log(u)
    percentile.log <- -rowSums(u.tilde^delta)^(1/delta)
    out.log <- matrix(percentile.log)

    if(log)
    {
      out <- out.log
    }
    else
    {
      out <- exp(out.log)
    }
  }
  else if(tolower(CplNM) == "frechet-upper")
  {
    out <- apply(u, 1, min)
  }
  else if(tolower(CplNM) == "frechet-lower")
  {
    d <- ncol(u)
    z <- cbind((rowSums(u)+1-d), 0)
    out <- apply(z, 1, max)
  }
  else
  {
    stop("Given copula is not implemented.")
  }

  return(out)
}
