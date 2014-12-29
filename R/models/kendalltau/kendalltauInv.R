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
##' TODO: Split the output by using a caller. The inverse look will fail when lambdaL>8.5
##' due to the reason of non-monotonic.
kendalltauInv <- function(
    CplNM,
    parRepCpl,
    method = c("tabular", "iterative")[1])
  {
    if(tolower(method) == "tabular")
      {
        ## If the tauTabular not exist, create it
        if(!exists("tauTabular", envir = .GlobalEnv))
          {
            tauTabular <<- kendalltauTabular(CplNM = CplNM, tol = 1e-3)
          }
        out <- kendalltauInv.tab(CplNM = CplNM,
                                 parRepCpl = parRepCpl,
                                 tauTabular = tauTabular)
      }
    else if(tolower(method) == "iterative")
      {
        out <- kendalltauInv.iter(CplNM = CplNM,
                                  parRepCpl = parRepCpl,
                                  parCaller = "theta")
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
        ## out <- vector("list", length(parRepCpl))
        lambdaL <- parRepCpl[["lambdaL"]]
        tau <- parRepCpl[["tau"]]


        if(length(lambdaL) !=length(tau))
          {
            stop("The input parameters should be of the same length.")
          }

        ## The dictionary look up method for the upper tail dependence given
        ## lower tail dependence and Kendall's tau.
        tol <- tauTabular$tol
        nGridL <- tauTabular$nGridL
        nGridU <- tauTabular$nGridU

        tauMat <- tauTabular$tauMat
        lambdaUGrid <- tauTabular$lambdaUGrid

        ## The lower tail dependence indices.
        lambdaLIdxRaw <- lambdaL/tol
        lambdaLIdxFloor <- round(lambdaLIdxRaw)

        ## Extra work to avoid under and over flow
        lambdaLIdxFloor1 <- (lambdaLIdxFloor < 1)
        lambdaLIdxFloor2 <- (lambdaLIdxFloor > nGridL)
        if(any(lambdaLIdxFloor1))
          {
            lambdaLIdxFloor[lambdaLIdxFloor1] <- 1
          }
        if(any(lambdaLIdxFloor2))
          {
            lambdaLIdxFloor[lambdaLIdxFloor2] <- nGridL
          }
        tauMatTabFloor <- tauMat[lambdaLIdxFloor, ,drop = FALSE]

        ## Find the indices of the closed values close to tau's left and right side
        nObs <- length(tau)
        tauTest <- matrix(tau, nObs, nGridU)

        tauFloorDev0 <- - abs(tauTest-tauMatTabFloor)


        ## The indices of lambdaU
        ## This is the bottom neck of speed.
        lambdaUFloorIdx0 <- max.col(tauFloorDev0)
        lambdaUFloor0 <- lambdaUGrid[lambdaUFloorIdx0]

        ## Make sure the output format is same as the input
        out <- tau
        out[1:length(out)] <- lambdaUFloor0
      }
    return(out)
  }

###----------------------------------------------------------------------------
### The iterative method
###----------------------------------------------------------------------------
kendalltauInv.iter <- function(CplNM, parRepCpl, parCaller = "theta")
  {
    ## TODO: The max interval could not handle Inf in the uniroot function.
    ## TODO: Consider the error handle. i.e., In theory (see the Appendix in
    ## the paper) you can't have small tau with big delta (lower tail
    ## dependent) which yields non accurate root.
    if(tolower(CplNM) == "bb7")
      {
        if(tolower(parCaller) == "theta")
          {
            lambdaL <- parRepCpl[["lambdaL"]]
            tau <- parRepCpl[["tau"]]
            delta <- -log(2)/log(lambdaL)

            theta.interval <- c(1, 1000)
            ## theta.interval <- c(1, 1.9)
            warning("Iteration method set restricted theta interval!")

            out.theta <- delta
            parLen <- length(delta)
            out.theta <- delta
            out.theta[1:parLen] <- NA

            for(i in 1:parLen)
              {
                tauCurr <- tau[i]
                deltaCurr <- delta[i]
                outRootCurr <-
                  try(uniroot(function(x, CplNM, deltaCurr, tauCurr)
                          {
                            kendalltau(CplNM = CplNM,
                                       parCpl = list(theta = x,
                                         delta = deltaCurr))-tauCurr
                          },
                          interval = theta.interval, CplNM = CplNM,
                          deltaCurr = deltaCurr,
                          tauCurr = tauCurr), silent = TRUE)

                if(is(outRootCurr, "try-error"))
                  {out.theta[i] <- NA}
                else
                  {
                    out.theta[i] <- outRootCurr$root
                  }
              }
            out <- 2 - 2^(1/out.theta)
          }
        else if(tolower(parCaller) == "delta")
          {
            lambdaU <- parRepCpl[["lambdaU"]]
            tau <- parRepCpl[["tau"]]
            theta <- log(2)/(log(2)-lambdaU)

            out.delta <- theta
            parLen <- length(theta)
            out.delta[1:parLen] <- NA

            for(i in 1:thetaLen)
              {
                tauCurr <- tau[i]
                deltaCurr <- delta[i]
                outRootCurr <-
                  try(uniroot(
                      function(x, thetaCurr, deltaCurr, tauCurr)
                      {
                        kendalltau(CplNM = CplNM,
                                   parCpl = list(theta = thetaCurr,
                                     delta = x))-tauCurr
                      },
                      interval = c(1e-10, 1000), CplNM = CplNM,
                      thetaCurr = thetaCurr, tauCurr =
                      tauCurr), silent = TRUE)
                ##print(outRootCurr)
                if(is(outRootCurr, "try-error"))
                  {out.delta[i] <- NA}
                else
                  {
                    out.delta[i] <- outRootCurr$root
                  }

              }

            out <- 2^(1/out.delta)
          }

      }
    return(out)
  }

###----------------------------------------------------------------------------
### TESTING
###----------------------------------------------------------------------------
## grd <- mesh.grid(seq(0.001, 0.999, 0.01))

## lambdaL <- grd[, 1]# runif(nObs)
## lambdaU <- grd[, 2]# runif(nObs)
## delta <- -log(2)/log(lambdaL)
## theta <- log(2)/log(2-lambdaU)
## tau <- kendalltau(CplNM = "BB7",
##                   parCpl = list(theta = theta, delta = delta))
