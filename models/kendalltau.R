##' The theoretical Kendall's tau for different copulas
##'
##' Most Kendall's tau has explicit solutions. The Kendall's tau for elliptical
##' class are all the same.
##' @title Theoretical Kendall's tau
##' @param CplNM "character" Copula name.
##' @param parCpl "list" Parameters in the copula
##' @return "vector" for the theoretical Kendall's tau.
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Apr 17 18:45:16 CEST 2012;
##'       Current: Tue Apr 17 18:45:23 CEST 2012.
kendalltau <- function(CplNM, parCpl)
  {
    if(tolower(CplNM) == "bb7")
      {
        theta <- parCpl[["theta"]]
        delta <- parCpl[["delta"]]
        nObs <- length(theta)

        ## The storage
        out <- theta
        out[0:nObs] <- NA

        ## Healthy condition and condition for the stepwise Kendall's tau
        deltaHcond <- (delta > 0)
        Idx12 <- which(theta >= 1 & theta < 2 & deltaHcond)
        IdxLarge <- which(theta > 2 & deltaHcond)
        Idx2 <- which(theta  == 2 & deltaHcond)

        ## The Kendall's tau
        if(length(Idx12)>0)
          {
            thetaCurr <- theta[Idx12]
            deltaCurr <- delta[Idx12]
            out[Idx12] <- 1-2/(deltaCurr*(2-thetaCurr)) +
              4*beta(deltaCurr+2, 2/thetaCurr-1)/(thetaCurr^2*deltaCurr)
          }
        if(length(IdxLarge)>0)
          {
            thetaCurr <- theta[IdxLarge]
            deltaCurr <- delta[IdxLarge]
            out[IdxLarge] <- 1-2/(deltaCurr*(2-thetaCurr)) -
              4*pi/(thetaCurr^2*deltaCurr*(2+deltaCurr)*sin(2*pi/thetaCurr)*
                    beta(1+deltaCurr+2/thetaCurr, 2-2/thetaCurr))
          }
        if(length(Idx2)>0)
          {
            thetaCurr <- theta[Idx2]
            deltaCurr <- delta[Idx2]
            out[Idx2] <- 1- (digamma(2+deltaCurr)-digamma(1)-1)/deltaCurr
          }
      }
    else if(tolower(CplNM) == "fgm")
      {
        theta <- parCpl[["theta"]]
        out <- 2/9*theta
      }
    else if (tolower(CplNM) %in% c("gaussian", "student-t"))
      {
        rho <- parCpl[["rho"]]
        out <- 2/pi*asin(rho)
      }


    return(out)
  }


## kendalltau0 <- function(CplNM, parCpl)
##   {
##     if(tolower(CplNM) == "bb7")
##       {
##         theta <- parCpl[["theta"]]
##         delta <- parCpl[["delta"]]

##         ## Healthy condition of parameters.
##         deltaHcond <- (delta > 0)

##         ## The gradient w.r.t. theta conditional on delta
##         if(theta >= 1 && theta < 2 && deltaHcond)
##           {
##             out <- 1-2/(delta*(2-theta)) +
##               4*beta(delta+2, 2/theta-1)/(theta^2*delta)
##           }
##         else if (theta > 2 && deltaHcond)
##           {
##             out <- 1-2/(delta*(2-theta)) -
##               4*pi/(theta^2*delta*(2+delta)*sin(2*pi/theta)*
##                     beta(1+delta+2/theta, 2-2/theta))
##           }
##         else if(theta  == 2 && deltaHcond)
##           {
##             ## You computer will never hit this branch. But we keep it anyway.
##             out <- 1- (digamma(2+delta)-digamma(1)-1)/delta
##           }
##         else # bad condition met.
##           {
##             out <- NA
##           }
##       }
##     return(out)
##   }



###----------------------------------------------------------------------------
### TESTING
###----------------------------------------------------------------------------
## n <- 100
## out <- numeric(n)
## i <- 0
## thetaCurr <- seq(1, 20, length.out = n)
## for(i in 1:n)
## {
##   out[i] <- kendalltau("BB7", delta = 1, thetaCurr = thetaCurr[i])
## }

## plot(thetaCurr, out)


## CplNM <- "bb7"
## parCpl <- list(delta = 5, theta = 3)
## Rprof()
## a <- proc.time()
## for(i in 1:1e5)
##   {
##     if(tolower(CplNM) == "bb7")
##       {
##         theta <- parCpl[["theta"]]
##         delta <- parCpl[["delta"]]

##         ## Healthy condition of parameters.
##         deltaHcond <- (delta > 0)

##         ## The gradient w.r.t. theta conditional on delta
##         if(theta >= 1 && theta < 2 && deltaHcond)
##           {
##             out <- 1-2/(delta*(2-theta)) +
##               4*beta(delta+2, 2/theta-1)/(theta^2*delta)
##           }
##         else if (theta > 2 && deltaHcond)
##           {
##             out <- 1-2/(delta*(2-theta)) -
##               4*pi/(theta^2*delta*(2+delta)*sin(2*pi/theta)*
##                     beta(1+delta+2/theta, 2-2/theta))
##           }
##         else if(theta  == 2 && deltaHcond)
##           {
##             ## You computer will never hit this branch. But we keep it anyway.
##             out <- 1- (digamma(2+delta)-digamma(1)-1)/delta
##           }
##         else # bad condition met.
##           {
##             out <- NA
##           }
##       }
##   }
## print(proc.time()-a)
## Rprof(NULL)
## summaryRprof()


## beta <- function(a, b)
##   {
##     gamma(a)*gamma(b)/gamma(a+b)
##   }


## The easy way to obtain the inverse function. But it seems very difficult
## with a long scaler input.


## tol <- 0.01
## lambdaL <- seq(0.01, 0.99, tol)
## lambdaU <- seq(0.01, 0.99, tol)
## lambda <- mesh.grid(lambdaL, lambdaU)
## delta <- -log(2)/log(lambda[, 1])
## theta <- log(2)/log(2-lambda[, 2])
## parCpl <- list(theta = theta, delta = delta)
## tau0 <- kendalltau0("bb7", parCpl)
## out <- cbind(lambda, tau0)



## CplNM <- "BB7"
## nObs <- 50
## tau <- runif(nObs)
## lambdaL <- runif(nObs)
## parRepCpl <- list(tau = tau, lambdaL = lambdaL)

## tauTol <- 0.01
## tauTabular <- kendalltauTabular(CplNM, tol = tauTol)

## parCpl <- list(delta = delta, theta = NA)



## parCplOut <- kendalltauInv(CplNM, parCpl, tau, tauTabular, parCaller = "theta")





## delta0 <- 1.6

## a <- proc.time()
## index <- (abs(tau0 - testtau)<0.001 & abs(parout[, 2]-delta0)<0.001)
## res <- parout[index, ]
## print(proc.time() -a)
