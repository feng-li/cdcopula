##' Empirical copula estimation
##'
##' Non-parametric method of estimating copula
##' @param u "numeric vector"
##' @param v "numeric vector"
##'
##'        Note that u and v should of the same length and both u and v should
##'        be in the interval of [0, 1].
##'
##' @param x "numeric"
##' @param y "numeric"
##' @return "numeric"
##'
##'         The numeric estimation of copula for data x and y.
##'
##' @references Dobric and Schmid 2005
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Jan 08 17:37:50 CET 2013;
##'       Current: Tue Jan 08 17:37:55 CET 2013.
hatCpl <- function(u, v, x, y)
  {

    n <- length(x)
    i <- round(u*n)
    j <- round(v*n)

    x.ord <- order(x)
    y.ord <- order(y)

    if(length(i) != length(j))
      {
        stop("i and j should have the same length")
      }

    n.ij <- length(i)

    x.mat <- matrix(x, n, n.ij)
    y.mat <- matrix(y, n, n.ij)

    x.ord.mat <- matrix(x[x.ord][i], n, n.ij, byrow = TRUE)
    y.ord.mat <- matrix(y[y.ord][j], n, n.ij, byrow = TRUE)

    x.sign.mat <- (x.mat <=  x.ord.mat)
    y.sign.mat <- (y.mat <=  y.ord.mat)

    sign.mat <- matrix(0, n, n.ij)
    sign.mat[as.vector(x.sign.mat & y.sign.mat)] <- 1
    out <- colMeans(sign.mat)
    return(out)
  }
