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
        deltaHcond <- (delta > 0)
        Idx12 <- which(theta >= 1 & theta < 2-tol & deltaHcond)
        IdxLarge <- which(theta > 2+tol & deltaHcond)
        Idx2 <- which(abs(theta-2)<tol & deltaHcond) # You will never reach this

        if(tolower(caller) == "theta")
          {

            ## The gradient w.r.t. theta conditional on delta
            if(length(Idx12)>0)
              {
                thetaCurr <- theta[Idx12]
                deltaCurr <- delta[Idx12]
                out[Idx12] <- -2/((theta-2)^2*delta) - 8*beta(2+delta, 2/theta-1)*
                  (theta + digamma(2/theta-1) - digamma(2/theta+delta+1))/
                    (theta^4*delta)
              }
            if(length(IdxLarge)>0)
              {
                thetaCurr <- theta[IdxLarge]
                deltaCurr <- delta[IdxLarge]
                out[IdxLarge] <- (-2*(2+delta)*theta^4*beta(1+delta+2/theta, 2-2/theta)-
                                  8*pi^2*(theta-2)^2*cos(2*pi/theta)/sin(2*pi/theta)^2-
                                  8*pi*(theta-2)^2*
                                  (digamma(1+delta+2/theta)-digamma(2-2/theta)-theta)/
                                  sin(2*pi/theta))/
                                    (delta*(2+delta)*(theta-2)^2*theta^4*
                                     beta(1+delta+2/theta, 2-2/theta))
              }
            if(length(Idx2)>0)
              {
                thetaCurr <- theta[Idx2]
                deltaCurr <- delta[Idx2]

                ## You computer will never hit this branch. But we keep it anyway.
                out[Idx2] <- -(12+24*digamma(1)+6*digamma(1)^2+pi^2-
                               12*(2+digamma(1))*digamma(2+delta)+
                               6*digamma(2+delta)^2-6*trigamma(2+delta))/
                                 (24*delta)
              }
          }
        else if(tolower(caller) == "delta")
          {
            ## Healthy condition of parameters.
            deltaHcond <- delta > 0

            ## The gradient w.r.t. theta conditional on delta
            if(length(Idx12)>0)
              {
                thetaCurr <- theta[Idx12]
                deltaCurr <- delta[Idx12]
                out[Idx12] <- -2/((theta-2)*delta^2) + 4*beta(2+delta, 2/theta-1)*
                  (digamma(2+delta) - digamma(2/theta+delta+1)-1/delta)/
                    (theta^2*delta)
              }
            if(length(IdxLarge)>0)
              {
                thetaCurr <- theta[IdxLarge]
                deltaCurr <- delta[IdxLarge]
                out[IdxLarge] <- -2/((theta-2)*delta^2)-
                  4*pi*(digamma(3+delta)-digamma(2/theta+delta+1)-2*(1+delta)/(2*delta+delta^2))/
                    ((2+delta)*delta*theta^2*sin(2*pi/theta)*beta(1+delta+2/theta, 2-2/theta))
              }
            if(length(Idx2)>0)
              {
                thetaCurr <- theta[Idx2]
                deltaCurr <- delta[Idx2]

                ## You computer will never hit this branch. But we keep it anyway.
                out[Idx2] <- (digamma(2+delta)-delta*trigamma(2+delta)-digamma(1)-1)/delta^2
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
