##' Empirical lower tail dependence coefficients.
##'
##' Based on nonparametric methods. Be aware it can be biased and
##' inconsistent. See the referenced paper for details.
##' @param x "numeric"
##'
##'        Random sample of length n.
##'
##' @param y "numeric".
##'
##'        Same as x.
##'
##' @param k "scaler".
##'
##'        Number of observations in the lower tail.
##'
##' @param method "character".
##'
##'        If "simple", use e.q (6) in Dobric and Schimd (2005). If "ols" use
##' e.q (7). If "mixture", the result are base on the idea that the copula can
##' be approximated by a mixture of the comonotonicity and independence copulas
##' (third method in Dobric and Schimd).
##'
##' @return "scaler".
##'
##'         The empirical tail dependence coefficient.
##'
##' @references Dobric and Schmid (2005)
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Jan 08 13:02:11 CET 2013;
##'
##'       Current: Tue Jan 08 13:02:16 CET 2013.
empltdc <- function(x, y, k, method)
  {
    n <- length(x)
    if(tolower(method) == "simple")
      {
        u <- k/n
        v <- k/n
        out <- (k/n)^(-1)*hatCpl(u = u, v = v, x = x, y = y)
      }
    else if(tolower(method) == "ols")
      {
        u <- (1:k)/n
        v <- (1:k)/n
        out <- (sum(((1:k)/n)^2))^(-1)*
          sum((1:k)/n*hatCpl(u = u, v = v, x = x, y = y))
      }
    else if(tolower(method) == "mixture")
      {
        u <- (1:k)/n
        v <- (1:k)/n
        out <- sum((hatCpl(u = u, v = v, x = x, y = y)-((1:k)/n)^2)*
                   ((1:k)/n-((1:k)/n)^2))/
                     (sum(((1:k)/n-((1:k)/n)^2)^2))
      }
    return(out)
  }
