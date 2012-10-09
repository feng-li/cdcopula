##' TODO: Vectorized this
kendalltauGrad <- function(CplNM, theta, delta, caller)
  {
    if(tolower(CplNM) == "bb7")
      {
        ## The storage
        nObs <- length(theta)
        out <- theta
        out[0:nObs] <- NA

        ## Healthy condition for the stepwise Kendall's tau
        tol = 0.001
        deltaHcond <- (delta > 0+tol)
        deltaHcondSmall <- (delta <=  0+tol & delta > 0) # A special case when delta is very close to zero (almost no tail dependence)

        Idx12 <- which(theta >= 1 & theta < 2-tol & deltaHcond)
        IdxLarge <- which(theta > 2+tol & deltaHcond)
        Idx2 <- which(abs(theta-2)<tol & deltaHcond) # You will never reach this
        IdxDeltaSmall <- which(deltaHcondSmall) # sometimes you hit this,  use
                                        # the asymptotic results.

        if(tolower(caller) == "theta")
          {
            ## The gradient w.r.t. theta conditional on delta
            if(length(Idx12)>0)
              {
                thetaCurr <- theta[Idx12]
                deltaCurr <- delta[Idx12]
                out[Idx12] <- -2/((thetaCurr-2)^2*deltaCurr) -
                  8*beta(2+deltaCurr, 2/thetaCurr-1)*
                  (thetaCurr + digamma(2/thetaCurr-1) -
                   digamma(2/thetaCurr+deltaCurr+1))/
                    (thetaCurr^4*deltaCurr)
              }
            if(length(IdxLarge)>0)
              {
                thetaCurr <- theta[IdxLarge]
                deltaCurr <- delta[IdxLarge]
                out[IdxLarge] <- (-2*(2+deltaCurr)*thetaCurr^4*
                                  beta(1+deltaCurr+2/thetaCurr, 2-2/thetaCurr)-
                                  8*pi^2*(thetaCurr-2)^2*cos(2*pi/thetaCurr)/sin(2*pi/thetaCurr)^2-
                                  8*pi*(thetaCurr-2)^2*
                                  (digamma(1+deltaCurr+2/thetaCurr)-
                                   digamma(2-2/thetaCurr)-thetaCurr)/sin(2*pi/thetaCurr))/
                                    (deltaCurr*(2+deltaCurr)*(thetaCurr-2)^2*thetaCurr^4*
                                     beta(1+deltaCurr+2/thetaCurr, 2-2/thetaCurr))
              }
            if(length(Idx2)>0)
              {
                thetaCurr <- theta[Idx2]
                deltaCurr <- delta[Idx2]

                ## You computer will never hit this branch. But we keep it anyway.
                out[Idx2] <- -(12+24*digamma(1)+6*digamma(1)^2+pi^2-
                               12*(2+digamma(1))*digamma(2+deltaCurr)+
                               6*digamma(2+deltaCurr)^2-6*trigamma(2+deltaCurr))/
                                 (24*deltaCurr)
              }
            if(length(IdxDeltaSmall)>0)
              {
                out[IdxDeltaSmall] <- 2*(1-harmonic(2/theta))/(theta-2)^2-
                  4*trigamma(2/theta+1)/((theta-2)*theta^2)
              }

          }
        else if(tolower(caller) == "delta")
          {
            ## The gradient w.r.t. theta conditional on delta
            if(length(Idx12)>0)
              {
                thetaCurr <- theta[Idx12]
                deltaCurr <- delta[Idx12]
                out[Idx12] <- -2/((thetaCurr-2)*deltaCurr^2) +
                  4*beta(2+deltaCurr, 2/thetaCurr-1)*
                    (digamma(2+deltaCurr) -
                     digamma(2/thetaCurr+deltaCurr+1)-1/deltaCurr)/
                       (thetaCurr^2*deltaCurr)
              }
            if(length(IdxLarge)>0)
              {
                thetaCurr <- theta[IdxLarge]
                deltaCurr <- delta[IdxLarge]
                out[IdxLarge] <- -2/((thetaCurr-2)*deltaCurr^2)-
                  4*pi*(digamma(3+deltaCurr)-
                        digamma(2/thetaCurr+deltaCurr+1)-
                        2*(1+deltaCurr)/(2*deltaCurr+deltaCurr^2))/
                          ((2+deltaCurr)*deltaCurr*thetaCurr^2*
                           sin(2*pi/thetaCurr)*
                           beta(1+deltaCurr+2/thetaCurr, 2-2/thetaCurr))
              }
            if(length(Idx2)>0)
              {
                thetaCurr <- theta[Idx2]
                deltaCurr <- delta[Idx2]

                ## You computer will never hit this branch. But we keep it anyway.
                out[Idx2] <- (digamma(2+deltaCurr)-
                              deltaCurr*trigamma(2+deltaCurr)
                              -digamma(1)-1)/deltaCurr^2
              }
          }
        else
          {
            stop("No such copula implementation for Kendall's tau.")
          }
      }

    return(out)
  }
###----------------------------------------------------------------------------
### TESTING: PASSED
###----------------------------------------------------------------------------
## n = 10000
## out <- numeric(n)
## theta <- seq(1, 20, length.out = n)
## for(i in 1:n)
##   out[i] <- kendalltauGrad("BB7", delta = 20, theta  = theta[i])
## plot(theta,out,  pch = ".")
