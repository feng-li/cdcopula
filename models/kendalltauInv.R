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
kendalltauInv <- function(CplNM, parRepCpl, tauTabular)
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
        lambdaLIdxFloor <- floor(lambdaLIdxRaw)
        lambdaLIdxCeiling <- ceiling(lambdaLIdxRaw)

        ## Extra work to avoid under and over flow
        lambdaLIdxFloor[lambdaLIdxFloor == 0] <- 1
        lambdaLIdxCeiling[lambdaLIdxCeiling > nGrid] <- nGrid

        tauMatTabFloor <- tauMat[lambdaLIdxFloor, ,drop = FALSE]
        tauMatTabCeiling <- tauMat[lambdaLIdxCeiling, ,drop = FALSE]

        ## Find the indices of the closed values close to tau's left and right side
        nObs <- length(tau)
        tauTest <- matrix(tau, nObs, nGrid)

        tauFloorDev0 <- tauTest-tauMatTabFloor
        tauFloorDev0[tauFloorDev0>0] <- -Inf

        tauFloorDev1 <- tauMatTabFloor-tauTest
        tauFloorDev1[tauFloorDev1>0] <- -Inf

        tauCeilingDev0 <- tauTest-tauMatTabCeiling
        tauCeilingDev0[tauCeilingDev0>0] <- -Inf

        tauCeilingDev1 <- tauMatTabCeiling - tauTest
        tauCeilingDev1[tauCeilingDev1>0] <- -Inf

        ## The indices of lambdaU
        ## TODO: This is the bottom neck of speed.
        lambdaUFloorIdx0 <- max.col(tauFloorDev0)
        lambdaUFloorIdx1 <- max.col(tauFloorDev1)

        lambdaUCeilingIdx0 <- max.col(tauCeilingDev0)
        lambdaUCeilingIdx1 <- max.col(tauCeilingDev1)

        ## lambdaUFloorIdx0 <- apply(tauFloorDev0, 1, which.max)
        ## lambdaUFloorIdx1 <- apply(tauFloorDev1, 1, which.max)

        ## lambdaUCeilingIdx0 <- apply(tauCeilingDev0, 1, which.max)
        ## lambdaUCeilingIdx1 <- apply(tauCeilingDev0, 1, which.max)

        lambdaUFloor0 <- lambdaUGrid[lambdaUFloorIdx0]
        lambdaUFloor1 <- lambdaUGrid[lambdaUFloorIdx1]

        lambdaUCeiling0 <- lambdaUGrid[lambdaUCeilingIdx0]
        lambdaUCeiling1 <- lambdaUGrid[lambdaUCeilingIdx1]

        ## And the corresponding tau
        tauFloorIdx0 <- whichInd(cbind(1:nObs, lambdaUFloorIdx0), c(nObs, nGrid))
        tauFloorIdx1 <- whichInd(cbind(1:nObs, lambdaUFloorIdx1), c(nObs, nGrid))

        tauCeilingIdx0 <- whichInd(cbind(1:nObs, lambdaUCeilingIdx0), c(nObs, nGrid))
        tauCeilingIdx1 <- whichInd(cbind(1:nObs, lambdaUCeilingIdx1), c(nObs, nGrid))

        tauFloor0 <- tauMatTabCeiling[tauFloorIdx0]
        tauFloor1 <- tauMatTabCeiling[tauFloorIdx1]

        tauCeiling0 <- tauMatTabFloor[tauCeilingIdx0]
        tauCeiling1 <- tauMatTabFloor[tauCeilingIdx1]

        ## Smooth the final results.
        ## (tauFloor -tau)/(tau-tauCeiling) = (lambdaUFloor-lambdaU)/(lambdaU-lambdaUCeiling)
        tauRatio0 <- exp(log(tauFloor0 -tau)-log(tau-tauFloor1))
        lambdaU0 <- (lambdaUFloor0 + lambdaUFloor1*tauRatio0)/(1+tauRatio0)

        tauRatio1 <- exp(log(tauCeiling0 -tau)-log(tau-tauCeiling1))
        lambdaU1 <- (lambdaUCeiling0 + lambdaUCeiling1*tauRatio1)/(1+tauRatio1)

        lambdaU <- (lambdaU0+lambdaU1)/2

        ## The output in the traditional form
        ## delta <- -log(2)/log(lambdaL)
        ## theta <- log(2)/log(2-lambdaU)

        out <- lambdaU
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


## The iterative method
## kendalltauInv0 <- function(CplNM, parRepCpl, parCaller = "theta")
##   {
##     ## TODO: The max interval could not handle Inf in the uniroot function.
##     ## TODO: Consider the error handle. i.e., In theory (see the Appendix in
##     ## the paper) you can't have small tau with big delta (lower tail
##     ## dependent) which yields non accurate root.
##     if(tolower(CplNM) == "bb7")
##       {

##         if(tolower(parCaller) == "delta")
##           {

##             theta <- parRepCpl[["theta"]]
##             tau <- parRepCpl[["tau"]]

##             out <- theta
##             thetaLen <- length(theta)
##             out[1:thetaLen] <- NA

##             for(i in 1:thetaLen)
##               {
##                 tauCurr <- tau[i]
##                 deltaCurr <- delta[i]
##                 outRootCurr <-
##                   uniroot(function(x, thetaCurr, deltaCurr, tauCurr)
##                           {
##                             kendalltau(CplNM = CplNM,
##                                        parCpl = list(theta = thetaCurr,
##                                          delta = x))-tauCurr
##                           },
##                           interval = c(1e-13, 1000), CplNM = CplNM,
##                           thetaCurr = thetaCurr, tauCurr =
##                           tauCurr)
##                 ##print(outRootCurr)
##                 out[i] <- outRootCurr$root
##               }
##           }
##         else if(tolower(parCaller) == "theta")
##           {
##             lambdaL <- parRepCpl[["lambdaL"]]
##             tau <- parRepCpl[["tau"]]

##             delta <- -log(2)/log(lambdaL)

##             out <- delta
##             deltaLen <- length(delta)
##             out[1:deltaLen] <- NA

##             for(i in 1:deltaLen)
##               {
##                 tauCurr <- tau[i]
##                 deltaCurr <- delta[i]
##                 outRootCurr <-
##                   uniroot(function(x, CplNM, deltaCurr, tauCurr)
##                           {
##                             kendalltau(CplNM = CplNM,
##                                        parCpl = list(theta = x,
##                                          delta = deltaCurr))-tauCurr
##                           },
##                           interval = c(1, 1000), CplNM = CplNM,
##                           deltaCurr = deltaCurr,
##                           tauCurr = tauCurr)
##                 ##print(outRootCurr)
##                 out[i] <- outRootCurr$root
##               }
##           }
##       }
##     return(out)
##   }
