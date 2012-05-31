##' Calculate the inverse function of the Kendall's tau
##'
##' The dictionary looking up method is used. This is a lot faster within 4
##' decimal numbers. If more precision is needed, consider using iteration
##' method or spline approximation.
##' @title Inverse Kendall's tau
##' @param CplNM "character" The copula name
##' @param parRepCpl "list" The parameter list in the reparmeterized copula.
##' @param tauTabular "list" The dictionary provide. See function
##' kendalltauTabular().
##' @return "list" The traditional parameter list for the copula.
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note DEPENDS: flutils
##'       Created: Sat Dec 31 02:05:14 CET 2011;
##'       Current: Mon Apr 16 14:40:16 CEST 2012.
##' TODO: Split the output by using a caller.


kendalltauInv <- function(CplNM, parCpl, method = "tabular", ...)
  {
    if(tolower(method) == "tabular")
      {
        out <- kendalltauInv.tab(CplNM = CplNM,
                                 parRepCpl = parRepCpl,
                                 tauTabular = tauTabular)

      }
    else if(tolower(method) == "iterative")
      {
        out <- kendalltauInv.iter(CplNM = CplNM,
                                  parRepCpl = parRepCpl,
                                  tauTabular = tauTabular)
      }
      return(out)
  }

###----------------------------------------------------------------------------
### The tabular looking up method (faster)
###----------------------------------------------------------------------------
kendalltauInv.tab <- function(CplNM, parRepCpl, tauTabular)
  {
    if(tolower(CplNM) == "bb7")
      {
        out <- vector("list", length(parRepCpl))
        lambdaL <- parRepCpl[["lambdaL"]]
        tau <- parRepCpl[["tau"]]

        ## The dictionary look up method for the upper tail dependence given
        ## lower tail dependence and Kendall's tau.
        tol <- tauTabular$tol
        nGrid <- tauTabular$nGrid
        tauMat <- tauTabular$tauMat
        lambdaUGrid <- tauTabular$lambdaUGrid

        ## The lower tail dependence indices.
        lambdaLIdxRaw <- lambdaL/tol
        lambdaLIdxFloor <- round(lambdaLIdxRaw)

        ## Extra work to avoid under and over flow
        lambdaLIdxFloor[lambdaLIdxFloor < 1] <- 1
        lambdaLIdxFloor[lambdaLIdxFloor > nGrid] <- nGrid

        tauMatTabFloor <- tauMat[lambdaLIdxFloor, ,drop = FALSE]

        ## Find the indices of the closed values close to tau's left and right side
        nObs <- length(tau)
        tauTest <- matrix(tau, nObs, nGrid)

        tauFloorDev0 <- abs(tauTest-tauMatTabFloor)

        ## The indices of lambdaU
        ## TODO: This is the bottom neck of speed.
        lambdaUFloorIdx0 <- max.col(tauFloorDev0)

        lambdaUFloor0 <- lambdaUGrid[lambdaUFloorIdx0]
        out <- lambdaUFloor0
      }
    return(out)
  }


###----------------------------------------------------------------------------
### The iterative method
###----------------------------------------------------------------------------
kendalltauInv0 <- function(CplNM, parRepCpl, parCaller = "theta")
  {
    ## TODO: The max interval could not handle Inf in the uniroot function.
    ## TODO: Consider the error handle. i.e., In theory (see the Appendix in
    ## the paper) you can't have small tau with big delta (lower tail
    ## dependent) which yields non accurate root.
    if(tolower(CplNM) == "bb7")
      {

        if(tolower(parCaller) == "delta")
          {
            lambdaU <- parRepCpl[["lambdaU"]]
            tau <- parRepCpl[["tau"]]
            theta <- log(2)/(log(2)-lambdaU)

            out <- theta
            thetaLen <- length(theta)
            out[1:thetaLen] <- NA

            for(i in 1:thetaLen)
              {
                tauCurr <- tau[i]
                deltaCurr <- delta[i]
                outRootCurr <-
                  uniroot(function(x, thetaCurr, deltaCurr, tauCurr)
                          {
                            kendalltau(CplNM = CplNM,
                                       parCpl = list(theta = thetaCurr,
                                         delta = x))-tauCurr
                          },
                          interval = c(1e-13, 1000), CplNM = CplNM,
                          thetaCurr = thetaCurr, tauCurr =
                          tauCurr)
                ##print(outRootCurr)
                out[i] <- outRootCurr$root
              }

          }
        else if(tolower(parCaller) == "theta")
          {
            lambdaL <- parRepCpl[["lambdaL"]]
            tau <- parRepCpl[["tau"]]

            delta <- -log(2)/log(lambdaL)

            out <- delta
            deltaLen <- length(delta)
            out[1:deltaLen] <- NA

            for(i in 1:deltaLen)
              {
                tauCurr <- tau[i]
                deltaCurr <- delta[i]
                outRootCurr <-
                  uniroot(function(x, CplNM, deltaCurr, tauCurr)
                          {
                            kendalltau(CplNM = CplNM,
                                       parCpl = list(theta = x,
                                         delta = deltaCurr))-tauCurr
                          },
                          interval = c(1, 1000), CplNM = CplNM,
                          deltaCurr = deltaCurr,
                          tauCurr = tauCurr)
                ##print(outRootCurr)
                out[i] <- outRootCurr$root
              }
          }
      }
    return(out)
  }
###----------------------------------------------------------------------------
### TESTING
###----------------------------------------------------------------------------
## nObs <- 5000 ## about 1.7sec
## tauTabular <- kendalltauTabular("BB7", tol = 0.001)

## tau <- runif(nObs, 0.001, 0.999)
## lambdaLMax <- 2^(1/2-1/(2*tau))

## lambdaL <- NA
## for(i in 1:nObs) lambdaL[i] <- runif(1, 0.002, lambdaLMax[i])

## parRepCpl <- list(lambdaL = lambdaL, tau = tau)

## a <- proc.time()
## lambdaU <- kendalltauInv(CplNM = "BB7", parRepCpl = parRepCpl,
##                         tauTabular = tauTabular)
## print(proc.time()-a)

## delta <- -log(2)/log(lambdaL)
## theta <- log(2)/log(2-lambdaU)
## tauEst <- kendalltau("BB7", parCpl = list(theta = theta, delta = delta))
