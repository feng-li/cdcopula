kendalltauGrad <- function(CplNM, parCpl, parCaller)
{
  out <- list()
  if(tolower(CplNM) == "bb7")
  {
    ## The storage
    theta <- parCpl[["theta"]]
    delta <- parCpl[["delta"]]

    nObs <- length(theta)

    ## Healthy condition for the stepwise Kendall's tau
    tol = 0.001
    deltaHcond <- (delta > 0+tol)
    ## A special case when delta is very close to zero (almost no tail
    ## dependence)
    deltaHcondSmall <- (delta <=  0+tol & delta > 0)

    Idx12 <- which(theta >= 1 & theta < 2-tol & deltaHcond)
    IdxLarge <- which(theta > 2+tol & deltaHcond)
    Idx2 <- which(abs(theta-2)<tol & deltaHcond) # You will never reach this
    IdxDeltaSmall <- which(deltaHcondSmall) # sometimes you hit this,  use
                                        # the asymptotic results.

    if("theta" %in% tolower(parCaller))
    {
      ## The gradient w.r.t. theta conditional on delta
      ## Conditional on the lower tail dependency lambdaL(delta),  the
      ## gradient for the Kendall's tau with respect to theta should be
      ## the same as usual.
      out.theta <- theta
      out.theta[0:nObs] <- NA


      if(length(Idx12)>0)
      {
        thetaCurr <- theta[Idx12]
        deltaCurr <- delta[Idx12]

############################################
        ## Debugging
        ## thetaCurr <- 1.5
        ## deltaCurr <- 2.5
        ## PASSED with analytical expression in Mathematica
############################################

        out.theta[Idx12] <- -2/((thetaCurr-2)^2*deltaCurr) -
                       8*beta(2+deltaCurr, 2/thetaCurr-1)*
                       (thetaCurr + digamma(2/thetaCurr-1) -
                        digamma(2/thetaCurr+deltaCurr+1))/
                       (thetaCurr^4*deltaCurr)
      }
      if(length(IdxLarge)>0)
      {
############################################
        ## Debugging
        ## thetaCurr <- 3.5
        ## deltaCurr <- 4.5
        ## PASSED with analytical expression in Mathematica
############################################

        thetaCurr <- theta[IdxLarge]
        deltaCurr <- delta[IdxLarge]
        out.theta[IdxLarge] <- (-2*(2+deltaCurr)*thetaCurr^4*
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

############################################
        ## Debugging
        ## thetaCurr <- 2
        ## deltaCurr <- 4.5
        ## PASSED with analytical expression in Mathematica
############################################

        ## You computer will never hit this branch. But we keep it anyway.
        out.theta[Idx2] <- -(12+24*digamma(1)+6*digamma(1)^2+pi^2-
                                                                12*(2+digamma(1))*digamma(2+deltaCurr)+
                                                                6*digamma(2+deltaCurr)^2-6*trigamma(2+deltaCurr))/
                            (24*deltaCurr)
      }
      if(length(IdxDeltaSmall)>0)
      {
############################################
        ## Debugging
        ## thetaCurr <- 3
        ## deltaCurr <- 0.0005
        ## PASSED with analytical expression in Mathematica
        ## NOTE: The two if conditions are essentially the same, interesting:)
############################################

        thetaCurr.A <- theta[IdxDeltaSmall & Idx12]
        thetaCurr.B <- theta[IdxDeltaSmall & IdxLarge]

        if(length(thetaCurr.A)>0)
        {
          out.theta[IdxDeltaSmall & Idx12] <- -
                    (thetaCurr.A^2*(4-2*digamma(1)*(thetaCurr.A-2)
                      +(thetaCurr.A-6)*thetaCurr.A) +
                                 2*(thetaCurr.A-2)*thetaCurr.A^2*digamma(2/thetaCurr.A-1)+
                                                               4*(thetaCurr.A-2)^2*trigamma((2+thetaCurr.A)/thetaCurr.A)
                    )/
                    (
                      (thetaCurr.A-2)^3*thetaCurr.A^2
                    )
        }
        if(length(thetaCurr.B)>0)
        {
          out.theta[IdxDeltaSmall & IdxLarge] <- 2*(1-harmonic(2/thetaCurr.B))/(thetaCurr.B-2)^2-
                                                                                               4*trigamma(2/thetaCurr.B+1)/((thetaCurr.B-2)*thetaCurr.B^2)
        }

      }
      out[["theta"]] <- out.theta
    }

    if("delta" %in% tolower(parCaller))
    {

############################################################################
      ## NOTE: This is also the unconditional gradient for Kendall's
      ## tau. The commenting out part is the way to introducing
      ## dependence between delta and theta,  i.e. theta  = ff(delta),
      ## we should prepare the following information first.
      ## tauGrad.delta <- ff'(delta)
#############################################################################
      out.delta <- delta
      out.delta[0:nObs] <- NA

      ## The gradient w.r.t. theta conditional on delta
      if(length(Idx12)>0)
      {
############################################
        ## Debugging
        ## thetaCurr <- 1.5
        ## deltaCurr <- 2.5
        ## PASSED with analytical expression in Mathematica
############################################
        thetaCurr <- theta[Idx12]
        deltaCurr <- delta[Idx12]

        out.delta[Idx12] <- -2/((thetaCurr-2)*deltaCurr^2) +
                       4*beta(2+deltaCurr, 2/thetaCurr-1)*
                       (digamma(2+deltaCurr) -
                        digamma(2/thetaCurr+deltaCurr+1)-1/deltaCurr)/
                       (thetaCurr^2*deltaCurr)

        ## B1 <- beta(2+deltaCurr, -1+1/thetaCurr)
        ## H2 <- harmonic(deltaCurr+2/thetaCurr)
        ## H4 <- harmonic(-2+2/thetaCurr)

        ## out.delta[Idx12] <- -1/(deltaCurr^2*(-2+thetaCurr)^2*thetaCurr^4)*(
        ##     2*(-2+thetaCurr)*thetaCurr^2*(-4*B1+2*B1*thetaCurr+thetaCurr^2)-
        ##     4*B1*deltaCurr*(-2+thetaCurr)^2*thetaCurr*harmonic(1+deltaCurr)+
        ##     4*B1*H2*deltaCurr*(-2+thetaCurr)^2*(thetaCurr^2-2*tauGrad.delta)+
        ##     8*B1*H4*deltaCurr*(-2+thetaCurr)^2*tauGrad.delta+
        ##     2*deltaCurr*thetaCurr*(
        ##         16*B1-16*B1*thetaCurr+4*B1*thetaCurr^2+thetaCurr^3
        ##         )*tauGrad.delta
        ##     )
      }
      if(length(IdxLarge)>0)
      {
############################################
        ## Debugging
        ## thetaCurr <- 3.5
        ## deltaCurr <- 4.5
        ## PASSED with analytical expression in Mathematica
############################################


        thetaCurr <- theta[IdxLarge]
        deltaCurr <- delta[IdxLarge]

        out.delta[IdxLarge] <- -2/((thetaCurr-2)*deltaCurr^2)-
                          4*pi*(digamma(3+deltaCurr)-
                                digamma(2/thetaCurr+deltaCurr+1)-
                                2*(1+deltaCurr)/(2*deltaCurr+deltaCurr^2))/
                          ((2+deltaCurr)*deltaCurr*thetaCurr^2*
                                                             sin(2*pi/thetaCurr)*
                                                             beta(1+deltaCurr+2/thetaCurr, 2-2/thetaCurr))

        ## B2 <- beta(1+deltaCurr+2/thetaCurr, 2-2/thetaCurr)
        ## H1 <- harmonic(1-2/thetaCurr)
        ## H2 <- harmonic(deltaCurr+2/thetaCurr)
        ## H3 <- harmonic(2+deltaCurr)

        ## out.delta[IdxLarge] <- -(
        ##     2*(B2*(2+deltaCurr)^2*thetaCurr^5+
        ##        16*pi*deltaCurr*(2+deltaCurr)*
        ##        (1-H1+H2+pi/tan(2*pi/thetaCurr))/sin(2*pi/thetaCurr)*
        ##        thetaCurr*tauGrad.delta+

        ##        4*pi/sin(2*pi/thetaCurr)*thetaCurr^3*
        ##        (4*(1+deltaCurr)-
        ##         deltaCurr*(2+deltaCurr)*
        ##         (-2*H2+2*H3+tauGrad.delta))+

        ##        thetaCurr^4*(
        ##            -4*pi*(1+deltaCurr)/sin(2*pi/thetaCurr)-
        ##            2*H2**pi*deltaCurr*(2+deltaCurr)/sin(2*pi/thetaCurr)+
        ##            2*H3*pi*deltaCurr*(2+deltaCurr)/sin(2*pi/thetaCurr)+
        ##            (2+deltaCurr)*gamma(2-2/thetaCurr)*gamma(1+deltaCurr+2/thetaCurr)*
        ##            (-2+deltaCurr*tauGrad.delta)) +

        ##        4*pi/sin(2*pi/thetaCurr)*thetaCurr^2*(
        ##            -4*(1+deltaCurr)+deltaCurr*(2+deltaCurr)*
        ##            (2*H3+H2*(-2+tauGrad.delta)+
        ##             (4-H1+pi/tan(2*pi/thetaCurr))*tauGrad.delta
        ##             )

        ##            )
        ##        )
        ##     )/(
        ##         B2*deltaCurr^2*(2+deltaCurr)^2*(-2+thetaCurr)^2*thetaCurr^4
        ##         )
      }
      if(length(Idx2)>0)
      {
############################################
        ## Debugging
        ## thetaCurr <- 2
        ## deltaCurr <- 4.5
        ## PASSED with analytical expression in Mathematica
############################################

        thetaCurr <- theta[Idx2]
        deltaCurr <- delta[Idx2]

        ## You computer will never hit this branch. But we keep it anyway.
        out.delta[Idx2] <- (digamma(2+deltaCurr)-
                            deltaCurr*trigamma(2+deltaCurr)
          -digamma(1)-1)/deltaCurr^2

        ## out[Idx2] <- (-1-digamma(1)+digamma(2+deltaCurr)-
        ##               deltaCurr*trigamma(2+deltaCurr))/deltaCurr^2
      }
      out[["delta"]] <- out.delta

    }
  }
  else if(tolower(CplNM) == "mvt")
  {
    rho <- parCpl[["rho"]]

    if("rho" %in% tolower(parCaller))
    {
      out[["rho"]] <- 2/(pi*sqrt(1-rho^2)) # n-by-lq
    }
  }
  else if(tolower(CplNM) == "gumbel")
  {
    delta <- parCpl[["delta"]]

    if("delta" %in% tolower(parCaller))
    {
      out[["delta"]] <- 1/delta^2
    }

  }
  else if(tolower(CplNM) == "clayton")
  {
    delta <- parCpl[["delta"]]

    if("delta" %in% tolower(parCaller))
    {
      out[["delta"]] <- 2/(2+delta)^2
    }
  }

  else
  {
    stop("No such copula!")
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
