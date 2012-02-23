##' Calculate the inverse function of the Kendall's tau
##'
##' The dictionary looking up method is used. This is a lot faster within 4
##' decimal numbers. If more precision is needed, consider using iteration
##' method or spline approximation.
##' @title Inverse Kendall's tau
##' @param CplNM 
##' @param parRepCpl 
##' @param tauTabular "list" The dictionary provide. See function
##' kendalltauTabular(). 
##' @return 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Sat Dec 31 02:05:14 CET 2011;
##'       Current: Sat Dec 31 02:05:21 CET 2011.
##' TODO: Merge iterative and dictionary as two options
##' TODO: Split the output by using a caller.
kendalltauInv <- function(CplNM, parRepCpl, tauTabular)
  {
    ## The dictionary look up method
    if(tolower(CplNM) == "bb7") 
      {
        out <- vector("list", length(parRepCpl))
        lambdaL <- parRepCpl[["lambdaL"]]
        tau <- parRepCpl[["tau"]]

        ## The first parameter: use direct method
        delta <- -log(2)/log(lambdaL)

        ## The second parameter
        tol <- tauTabular$tol
        nGrid <- tauTabular$nGrid
        tauMat <- tauTabular$tauMat
        lambdaUGrid <- tauTabular$lambdaUGrid

        ## The idices.
        lambdaLIdx <- round(lambdaL/tol)
        
        ## Extra work to avoid under and over flow
        lambdaLIdx[lambdaLIdx == 0] <- 1
        lambdaLIdx[lambdaLIdx == (nGrid+1)] <- nGrid
            
        tauMatTab <- tauMat[lambdaLIdx, ,drop = FALSE]
        nObs <- length(tau)
        tauTest <- matrix(tau, nObs, nGrid)
        tauAbsDev <- -abs(tauMatTab-tauTest)
       
        ## lambdaUIdx <- apply(tauAbsDev, 1, which.min) # max.col is slightly faster.
        lambdaUIdx <- max.col(tauAbsDev)
        
        ## The upper tail dependence via dictionary search.
        lambdaU <- matrix(lambdaUGrid[lambdaUIdx])
        
        ## The cross ponding theta
        theta <- log(2)/log(2-lambdaU)

        ## The output
        out <- list(theta = theta, delta = delta)

      }
    return(out)
  }

###----------------------------------------------------------------------------
### TESTING CODE
###----------------------------------------------------------------------------

## The iterative method 
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
            
            theta <- parRepCpl[["theta"]]
            tau <- parRepCpl[["tau"]]
              
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
## nObs <- 1
## a <- proc.time()
## tauTabular <- kendalltauTabular("BB7", tol = 0.0005)
## print(proc.time()-a)

## tau <- rep(c(0.5, 0.7, 0.9, 0.3), nObs)
## lambdaL <- rep(c(0.5, 0.3, 0.7, 0.1), nObs)

## parRepCpl <- list(lambdaL = lambdaL, tau = tau)

## a <- proc.time()
## ## Rprof()
## parCpl <- kendalltauInv(CplNM = "BB7", parRepCpl = parRepCpl,
##                         tauTabular = tauTabular)
## ## Rprof(NULL)
## ## summaryRprof()

## print(proc.time()-a)

## tau0 <- kendalltau("BB7", parCpl = parCpl)

## a <- proc.time()
## theta0 <- kendalltauInv0("bb7", parRepCpl = parRepCpl, parCaller = "theta")
## print(proc.time()-a)

