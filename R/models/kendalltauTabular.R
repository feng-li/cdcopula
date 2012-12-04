##' Generating tabular data for the Kendall's tau with respect to its input
##' arguments.
##'
##' This file is used for the looking up dictionary if the inverse Kendall's is
##' difficult to calculate.
##' @title Kendall's tau tabular
##' @param CplNM "character" Copula name
##' @param tol "numeric" Scaler show how precise of the dictionary is needed.
##' @return "list" All needed information for the copula Kendall's tau inverse.
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Dec 26 17:17:38 CET 2011;
##'       Current: Thu Apr 19 18:04:39 CEST 2012.
##' TODO: This is a very standard way of creating tables. It should be possible
##' to make a general function.
kendalltauTabular <- function(CplNM, tol = 1e-4)
{
  if(tolower(CplNM) == "bb7")
    {
      ## The dictionary lookup method The input argument. We choose to use the
      ## lower and upper tail dependence because they are fixed in [0, 1] for
      ## BB7 The code is only used once during the initialization.  If need
      ## more precisions is needed , we consider using iterative way to handle
      ## the memory problem.

      lambdaLGrid <- seq(0+tol, 1-tol, tol)
      lambdaUGrid <- lambdaLGrid

      nGrid <- length(lambdaLGrid)
      tauMat <- matrix(NA, nGrid, nGrid)

      ## Big table takes huge amount of memory. We split the calculation if we
      ## require a very precise table.
      ## Split the calculations
      MaxLenCurr <- min(nGrid, 2000)
      LoopIdx <- c(seq(1, nGrid, MaxLenCurr), nGrid)
      LoopIdx[1] <- 0
      tauIdxCurr0 <- 0

      nLoops <- length(LoopIdx)-1
      for(j in 1:nLoops)
      {
        IdxCurr0 <- LoopIdx[j]+1
        IdxCurr1 <- LoopIdx[j+1]

        lambdaL <- rep(lambdaLGrid, times = IdxCurr1-IdxCurr0+1)
        lambdaU <- rep(lambdaUGrid[IdxCurr0:IdxCurr1], each = nGrid)

        delta <- -log(2)/log(lambdaL)
        theta <- log(2)/log(2-lambdaU)
        parCpl <- list(theta = theta, delta = delta)
        tauCurr <- kendalltau(CplNM, parCpl)

        tauMat[, IdxCurr0:IdxCurr1] <- tauCurr
      }

      out <- list(tauMat = tauMat,
                  nGrid = nGrid,
                  tol = tol,
                  lambdaUGrid = lambdaUGrid)
    }
  return(out)
}

###----------------------------------------------------------------------------
### TESTING
###----------------------------------------------------------------------------
## out <- list(tauMat = tauMat, nGrid = nGrid, tol = tol,
##             lambda = lambda, tau = tau, lambdaLGrid = lambdaLGrid,
##             lambdaUGrid = lambdaUGrid)

## Approximate the inverse function with additive spline
## theta ~ AdditiveSplineMolde(tau, lambdaL)
## x <- cbind(tau, lambda[, 1])
## nKnots <- round(length(lambdaL)*0.5) # sqrt(n) knots.
## splineArgs <- list(comp = c("intercept", "covariates", "thinplate.a"),
##                    thinplate.a.locate = rep(nKnots, ncol(x)))
## knots <- make.knots(x = x, method = "k-means",
##                     splineArgs = splineArgs)
## Xdesign <- d.matrix(x = x, knots = knots, splineArgs = splineArgs)
## the residual is big
