##' Simple data generating process for the deterministic linear model Y = XB
##' given Y and coefficients.
##'
##' Remember there is no error term in the model. This function is primarily
##' used in the covariate-dependent type of models to generate desired
##' features.
##' @title Deterministic linear model DGP.
##' @param Y "scaler"
##' @param beta "scaler"
##' @param Xlim "vector"
##' @param intercept "logical"
##'        The intercept be considered in the DGP if the value for the argument
##' is TRUE.
##' @return "matrix"
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Wed Mar 14 10:46:20 CET 2012;
##'       Current: Wed Mar 14 10:46:25 CET 2012.
DGPlm <- function(Y, beta, Xlim, intercept = TRUE)
  {
    n <- length(Y)
    q <- length(beta)
    B <- matrix(beta) # Currently only allow univariate response

    if(intercept)
      {
        q0 <- q-1
        X0 <- matrix(runif(n*q0, min = Xlim[1], max = Xlim[2]), n, q0)
        X <- cbind(1, X0)
        q1 <- 2
      }
    else
      {
        q0 <- q
        X0 <- matrix(runif(n*q0, min = Xlim[1], max = Xlim[2]), n, q0)
        X <- X0
        q1 <- 1
      }

    ## Randomly select some positions to be determined.
    IdxLast <- sample(x = q1:q, size = n, replace = TRUE)
    IdxLast1 <- whichInd(arr.ind = cbind(1:n, IdxLast), dims = c(n, NA))

    ## Make it zero temporally.
    X[IdxLast1] <- 0

    ## Determine the last hole
    ## TODO: What if the last hole is terrible(far away from the X domain)?
    Y1 <- X%*%B
    XLast <- (Y-Y1)/beta[IdxLast]

    ## Fix the hole
    X[IdxLast1] <- XLast

    return(X)
  }
###----------------------------------------------------------------------------
### TESTING
###----------------------------------------------------------------------------

## Y <- matrix(rnorm(20, 1))
## beta <- runif(6)
## Xlim <- c(-3, 3)

## DGPlm(Y, beta, Xlim, intercept = TRUE)
