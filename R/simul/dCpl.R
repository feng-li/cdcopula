##' Compute the copula density
##'
##' @title Calculate the density function for give bivariate copula.
##'
##' @param u "matrix"
##'     Each column of the matrix is the value in the unit interval of the real
##' line for the corresponding margin of the copula function.
##' @param theta "list";
##'     The parameter list supplied to the copula function.
##' @param copula "character";
##'     This is the option for different types of copulas used in the function.
##' Currently "gaussian" for Gaussian copula.
##' @param par "list";
##'     Any further parameters need to pass to the copula. Which is equivalent
##' to the ... argument but less possible to make error.
##' @return "list" see below.
##'     \item {copula} {"matrix"; The copula function with the same length of
##' that in each column entries in "u".}
##'     \item {tao} {The theoretical Kendall's tao.}
##'     \item {rho} {The theoretical Spearman's rho.}
##' @references
##'     Trivedi and Zimmer (2007).
##' @author
##'     Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     Created: Mon Sep 12 13:36:10 CEST 2011;
##'     Current: Tue Sep 13 11:38:16 CEST 2011.
##' TODO: the first argument maybe should be x instead of u,
##'       Return log density instead of current form.
dCpl <- function(CplNM, u, parCpl, log = TRUE)
{
  if(tolower(CplNM) == "bb7")
  {
    ## theta,  delta
    theta <- parCpl[["theta"]]
    delta <- parCpl[["delta"]]

    TC1 <- 1-(1-u)^theta
    TC2 <- (1-u)^(-1+theta)

    L5 <- rowSums(TC1^(-delta)) - 1
    L6 <- 1-L5^(-1/delta) # FIXME: log(L6)->Inf when u->1,  v->1.

    logCplDensObs <- (-1-delta)*rowSums(log(TC1))+
                     rowSums(log(TC2))-
                     2*(1+delta)/delta*log(L5)+
                     (-2+1/theta)*log(L6)+
                     log(-1+theta+L5^(1/delta)*L6*(1+delta)*theta)

    out.log <- matrix(logCplDensObs)
  }
  else if(tolower(CplNM) == "gaussian")
  {
    rho <- parCpl[["theta"]] # n-by-lq

    ## The quantile for normal CDF
    u.quantile <- qnorm(u) # n-by-q

    nObs <- nrow(u.quantile)
    ## The CplNM density function C_12(u1, u2)

    dmvNormVecFun <- function(i, x, rho)
    {
      Sigma = vech2m(rho[i, ], diag = FALSE)
      out <- dmvnorm(x = x[i, , drop = FALSE],
                     sigma = Sigma,
                     log = TRUE)
      return(out)
    }
    logDensUpper <- apply(matrix(1:nObs),1,dmvNormVecFun,
                          x = u.quantile, rho = rho)

    logDensLower <- apply(dnorm(u.quantile, log = TRUE), 1, sum)

    logDens <- logDensUpper-logDensLower

    ## The output
    out.log <- matrix(logDens)
  }
  else if(tolower(CplNM) == "mvt") # The multivariate t-copula
  {## Demarta & McNeil (2005),  The t copula and related copulas

    require("mvtnorm")

    ## df, corr
    df <- parCpl[["df"]] # n-by-1
    rho <- parCpl[["rho"]] # n-by-lq

    ## The quantile for *univariate* t
    u.quantile <- qt(u, df = df)
    ## u.quantile <- X


    ## The copula density function C_12(u1, u2)
    nObs <- length(df)
    dmvtVecFun <- function(i, x, rho, df)
    {
      Sigma = vech2m(rho[i, ], diag = FALSE)
      out <- dmvt(x = x[i, , drop = FALSE],
                  sigma = Sigma,
                  df = df[i], log = TRUE)

      return(out)
    }
    logDensUpper <- apply(matrix(1:nObs),1,dmvtVecFun,
                          x = u.quantile, rho = rho, df = df)

    logDensLower <- apply(dt(u.quantile, df = df, log = TRUE), 1, sum)

    logDens <- logDensUpper-logDensLower

    ## The output
    out.log <- matrix(logDens)
  }
  else if(tolower(CplNM) == "fgm")
  {
    u1 <- u[, 1]
    u2 <- u[, 2]
    density <- 1+theta*(1-2*u1)*(1-2*u2)
    out <- density
  }
  else if(tolower(CplNM) == "gumbel")
  {# Joe 1997. p.142 Family B6
    delta <- as.vector(parCpl[["delta"]])
    pctl.log <- pCpl(u = u, parCpl = parCpl, CplNM = "gumbel", log = TRUE)
    u.tilde <- -log(u)

    u.tildeProd <- u.tilde[, 1]*u.tilde[, 2]
    u.tildeSumdelta <- rowSums(u.tilde^delta)

    out.log <- (pctl.log-log(u[, 1])-log(u[, 2])+
                (delta-1)*(log(u.tilde[, 1])+log(u.tilde[, 2]))-
                (2-1/delta)*log(u.tildeSumdelta)+
                log(u.tildeSumdelta^(1/delta)+delta-1))
  }
  else if(tolower(CplNM)  == "frank")
  {
    ## The quantile for normal CDF
    u1 <- u[, 1]
    u2 <- u[, 2]
    ## The copula function
    density <- -1/theta*log(1+(exp(-theta*u1)-1)*(exp(-theta*u2)-1)/
                            (exp(-theta)-1))
    ## The output
    out <- density
  }
  else
  {
    stop("Given copula name is not implemented.")
  }



  if(log)
  {
    out <- out.log
  }
  else
  {
    out <- exp(out.log)
  }

  return(out)
}
