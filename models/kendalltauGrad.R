##' TODO: Vectorized this
kendalltauGrad <- function(CplNM, theta, delta, caller)
  {
    if(tolower(CplNM) == "bb7")
      {
        if(tolower(caller) == "theta")
          {
            ## Healthy condition of parameters.
            deltaHcond <- (delta > 0)

            ## The gradient w.r.t. theta conditional on delta
            if(theta >= 1 && theta < 2 && deltaHcond)
              {
                out <- -2/((theta-2)^2*delta) - 8*beta(2+delta, 2/theta-1)*
                  (theta + digamma(2/theta-1) - digamma(2/theta+delta+1))/
                    (theta^4*delta)
              }
            else if (theta > 2 && deltaHcond)
              {
                out <- (-2*(2+delta)*theta^4*beta(1+delta+2/theta, 2-2/theta)-
                         8*pi^2*(theta-2)^2*cos(2*pi/theta)/sin(2*pi/theta)^2-
                         8*pi*(theta-2)^2*
                         (digamma(1+delta+2/theta)-digamma(2-2/theta)-theta)/
                         sin(2*pi/theta))/
                        (delta*(2+delta)*(theta-2)^2*theta^4*
                         beta(1+delta+2/theta, 2-2/theta))
              }
            else if(theta  == 2 && deltaHcond)
              {
                ## You computer will never hit this branch. But we keep it anyway.
                out <- -(12+24*digamma(1)+6*digamma(1)^2+pi^2-
                          12*(2+digamma(1))*digamma(2+delta)+
                          6*digamma(2+delta)^2-6*trigamma(2+delta))/
                        (24*delta)
              }
            else # bad condition met.
              {
                out <- NA
              }
          }
        else if(tolower(caller) == "delta")
          {
            ## Healthy condition of parameters.
            deltaHcond <- delta > 0

            ## The gradient w.r.t. theta conditional on delta
            if(theta >= 1 && theta < 2 && deltaHcond)
              {
                out <- -2/((theta-2)*delta^2) + 4*beta(2+delta, 2/theta-1)*
                  (digamma(2+delta) - digamma(2/theta+delta+1)-1/delta)/
                    (theta^2*delta)
              }
            else if (theta > 2 && deltaHcond)
              {
                out <- -2/((theta-2)*delta^2)-
                  4*pi*(digamma(3+delta)-digamma(2/theta+delta+1)-2*(1+delta)/(2*delta+delta^2))/
                    ((2+delta)*delta*theta^2*sin(2*pi/theta)*beta(1+delta+2/theta, 2-2/theta))
              }
            else if(theta  == 2 && deltaHcond)
              {
                ## You computer will never hit this branch. But we keep it anyway.
                out <- (digamma(2+delta)-delta*trigamma(2+delta)-digamma(1)-1)/delta^2
              }
            else # bad condition met.
              {
                out <- NA
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
