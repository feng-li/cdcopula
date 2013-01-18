##' Gradient for log copula function
##'
##'
##' @title Log copula gradient
##' @param CplNM
##' @param u
##' @param parCpl
##' @param cplCaller
##' @param staticArgs
##' @return
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Fri May 11 12:42:20 CEST 2012;
##'       Current: Fri May 11 12:42:30 CEST 2012.
logCplGrad <- function(CplNM, u, parCpl, cplCaller, staticArgs, Mdl.X, Mdl.beta)
  {

###----------------------------------------------------------------------------
### Copula likelihood numerical correction if u -> 0 or 1
###----------------------------------------------------------------------------

  ## Fix u on the cliff if any u -> 0 or u -> 1.
  ## Thanks to the advice from M. Smith

  tol <- .Machine$double.eps*1e8

  u.bad1 <- (u > 1-tol)
  u.bad0 <- (u < 0+tol)
  if(any(u.bad1))
    {
      u[u.bad1] <- 1-tol
      warning("u is too close to 1. Adjusted...",
              immediate. = TRUE)

    }

  if(any(u.bad0))
    {
      u[u.bad0] <- 0 +tol
      warning("u is too close to 1. Adjusted...",
              immediate. = TRUE)
    }

###----------------------------------------------------------------------------
### Gradients for the copula
###----------------------------------------------------------------------------


    if(tolower(CplNM) == "bb7")
      {
        ## The name of marginal model
        MargisNM <- dimnames(u)[[2]]
        nObs <- dim(u)[1]

        ## Subtract the parameters list.
        ## NOTE: convert matrix into vector to match the calculation
        tau <- as.vector(parCpl[["tau"]])
        lambdaL <- as.vector(parCpl[["lambdaL"]])

        lambdaU <- as.vector(kendalltauInv(
            CplNM = CplNM,
            parRepCpl = parCpl,
            tauTabular = staticArgs[["tauTabular"]]))

        ## The standard copula parameters (recycled if necessary, should be a
        ## vector).

        delta <- -log(2)/log(lambdaL)
        theta <- log(2)/log(2-lambdaU) # ff(delta)

        if(tolower(cplCaller) == "lambdal")
          {
            ## The commenting out part of the code is for the unconditional
            ## link function which is deprecated.

            ## Gradient w.r.t delta
            ## T1 <- 1-(1-u)^theta
            ## L1 <- rowSums(T1^(-delta))-1
            ## Delta1 <- -rowSums(T1^(-delta)*log(T1))

            ## logGradCpl.delta <-  - rowSums(log(T1))-
            ##   2*(1+delta)*Delta1/(delta*L1)-
            ##     (1/theta-2)*(log(L1)-delta*Delta1/L1)/
            ##       (delta^2*(L1^(1/delta)-1))+
            ##         2*log(L1)/delta^2+
            ##           (L1^(1/delta)-(1+delta)*L1^(1/delta)*
            ##            (log(L1)-delta*Delta1/L1)/delta^2-1)/
            ##              ((1+delta)*L1^(1/delta)-delta-1/theta)

            ## tauGrad.delta <- kendalltauGrad(CplNM = CplNM,
            ##                                 theta = theta,
            ##                                 delta = delta,
            ##                                 caller = "delta")

            ###########################################################################
            ## This should be obtained through the conditional linkage
            ## TODO: This is kind of hard code, consider it in a more general way.


            ## Gradient w.r.t. tau
            gradCpl.tau.theta <- kendalltauGrad(
                CplNM = CplNM, theta = theta,
                delta = delta, caller = "theta")

            ## The gradient for the parameters in conditional link
            ## tau.b <- 1-0.1 ## NOTE: Numerical stable to not allow tau  =  1
            ## tau.a <- log(2)/(log(2)-2*log(lambdaL))
            ## if(any(tau<tau.a)) browser()
            ## linPred.tau0 <-  log(tau-tau.a) - log(tau.b-tau)

            X <- Mdl.X[[CplNM]][["tau"]]
            beta <- Mdl.beta[[CplNM]][["tau"]]

            linPred.tau <- as.vector(X%*%beta)
            grad.glogit.a <- 1/(1+exp(linPred.tau))

            grad.link.a.lambdaL <- grad.glogit.a*
              2*log(2)/(lambdaL*(log(2)-2*log(lambdaL))^2)

            ## The gradient for the reparameterized parameters
            ## lambdaL  =  2^(-1/delta)
            grad.lambdaL.delta <- 2^(-1/delta)*log(2)/delta^2
            grad.delta.theta <- (1/gradCpl.tau.theta)*
              (grad.link.a.lambdaL*grad.lambdaL.delta)
            ###########################################################################

            ub <- 1-u
            M12 <- 1-ub^theta

            ## Numeric check if M12 is too close to 1 or 0.
            tol <- .Machine$double.eps*1e8
            M12.bad1 <- (M12> (1-tol))
            if(any(M12.bad1))
              {
                M12[M12.bad1] <- 1 - tol
                warning("Numerical unstable on M12,  adjust on the cliff...",
                        immediate. = TRUE)
              }
            M12.bad0 <- (M12<tol)
            if(any(M12.bad0))
              {
                M12[M12.bad0] <- tol
                warning("Numerical unstable on M12,  adjust on the cliff...",
                        immediate. = TRUE)
              }


            M34 <- -log(M12) + ub^theta*delta*log(ub)*grad.delta.theta/M12

            M5 <- -1 + rowSums(M12^(-delta))
            M6 <- -1 + M5^(1/delta)

            M7 <- (log(M5^(1/delta))-log(M6))*grad.delta.theta/theta^2
            C1 <- rowSums((M12^(-delta)*M34))/M5
            M8 <- C1*delta + log(M5)
            M9 <- M7-(M8*(-1+1/theta))/(M6*delta^2)

            P12 <- M6*(1+delta)*theta*log(ub)*grad.delta.theta
            P34 <- -M12*log(M12)+(1+delta)*log(ub)*ub^theta*grad.delta.theta
            P56 <- (-1+theta)*log(ub)*grad.delta.theta
            P78 <- M6*P34*(1+delta)*theta/M12

            C34 <- P34*(-1+theta)/M12

            S1 <- M12^delta

            logGradCpl.delta <- (
                (S1[, 1]*S1[, 2])^(-2)*(rowSums(S1)-S1[, 1]*S1[, 2])^2*
                (
                    rowSums(C34+P12+P56+P78) - M6*theta/delta +
                    M5^(1/delta)*(1-M5^(-1/delta))*M9*(1+delta)*theta+
                    M6*(1+delta)*theta/delta-
                    (M7-M8*(-2+1/theta)/(M6*delta^2))*(-1+1/theta)*theta-
                    2*(-1+1/theta)*theta*(-C1*delta*(1+delta)+log(M5))/delta^2+
                    M5^(1/delta)*(1-M5^(-1/delta))*(1+delta)*theta*
                    (C1*(-2-1/delta)+log(M5)/delta^2)+ M6*(1+delta)*grad.delta.theta+
                    grad.delta.theta/theta+(-1+theta*grad.delta.theta)/theta
                    )
            )/(
                M5^2*(-1+(-delta+M5^(1/delta)*(1+delta))*theta)
                )


            ## Gradient w.r.t. lambdaL
            ## gradCpl.lambdaL <- 2^(-1/delta)*log(2)/delta^2

            ## The chain gradient
            out <- logGradCpl.delta/grad.lambdaL.delta
            if(any(is.na(out))) browser()
          }
        else if(tolower(cplCaller) == "tau")
          {
################################################################################
### DEBUGGING
            ## u <- matrix(c(0.6, 0.3), 1, )
            ## theta <- 3.5
            ## delta <- 2.4
### PASSED
################################################################################

            ## Gradient w.r.t theta
            T1 <- 1-(1-u)^theta

            ## Numeric check if T1 is too close to 1.
            tol <- .Machine$double.eps*1e8
            T1.bad1 <- (T1> (1-tol))
            if(any(T1.bad1))
              {
                T1[T1.bad1] <- 1 - tol
                warning("Numerical unstable on T1,  adjusted on the cliff...",
                        immediate. = TRUE)
              }

            T1.bad0 <- (T1<tol)
            if(any(T1.bad0))
              {
                T1[T1.bad0] <- tol
                warning("Numerical unstable on T1,  adjusted on the cliff...",
                        immediate. = TRUE)
              }

            L1 <- rowSums(T1^(-delta))-1
            PT1 <- matrix(T1[, 1]*T1[, 2])

            SD12 <- rowSums((T1[, 2:1])^(1+delta)*(1-u)^theta*log(1-u))
            SD34 <- rowSums(T1^(-1-delta)*(1-u)^theta*delta*log(1-u))

            C1 <- SD34*L1^(-(1+delta)/delta)/(delta-L1^(-1/delta)*delta)
            C2 <- (log(L1^(1/delta))-log(-1+L1^(1/delta)))/theta^2

            logGradCpl.theta <-
              1/(L1^3*(-1-delta*theta+L1^(1/delta)*(1+delta)*theta))*
                PT1^(-1-2*delta)*(rowSums(T1^delta)-PT1^delta)^2*
                  (1/theta*PT1^(-delta)*(-SD12*(1+delta)*theta*(
                      -2+theta-2*delta*theta+L1^(1/delta)*(1+2*delta)*theta)-L1*PT1^(1+delta)*(
                          C1*(-1+theta)*(-1+theta-theta*delta+L1^(1/delta)*(1+delta)*theta)+
                          theta*(C2+delta+C2*delta*theta-L1^(1/delta)*(1+delta)*(1+C2*theta))))+
                   L1*(-1-delta*theta+L1^(1/delta)*(1+delta)*theta)*
                   (rowSums(T1[, 2:1]*(T1+(1-u)^theta*(1+delta))*log(1-u))))


            ## Gradient w.r.t. tau
            gradCpl.tau.theta <- kendalltauGrad(
                CplNM = CplNM,
                theta = theta,
                delta = delta,
                caller = "theta")

            ## The chain gradient
            out <- logGradCpl.theta*(1/gradCpl.tau.theta)
            if(any(is.na(out))) browser()

          }
        else
          {
            ## Gradient w.r.t u. NOTE: The BB7 copula's marginal are
            ## exchangeable which means the expression for the gradient w.r.t u1
            ## and u2 are the same if swap u1 and u2
            if(tolower(cplCaller) == "u1")
              {
                u <-  u[, 1:2, drop = FALSE]
              }
            else if(tolower(cplCaller) == "u2")
              {
                u <-  u[, 2:1, drop = FALSE]
              }
            else
              {
                stop("No such copula parameter!")
              }

################################################################################
            ## DEBUGGING
            ## u <- matrix(c(0.6, 0.3), 1, )
            ## theta <- 3.5
            ## delta <- 2.4
            ## PASSED
################################################################################

            ub <- 1 - u
            ub1 <- ub[, 1, drop = FALSE]

            D12 <- 1-ub^theta
            D1 <- D12[, 1, drop = FALSE]
            D2 <- D12[, 2, drop = FALSE]

            S1 <- 1 - rowSums(D12^(-delta))
            S2 <- -1 + (-S1)^(1/delta)

            gradCpl.u <- (D1^(-1-3*delta)*D2^(-2*delta)*
                          (rowSums(D12^delta)-D1^delta*D2^delta)^2*
                          (D1^(1+delta)*S1*(-1+theta)*
                           (1+(-S1)^(1/delta)*(-1+theta)-theta+S2^2*theta+S2^2*delta*theta)+
                           D2^(-delta)*(D1^delta*(-1+D2^delta)*S2*(1+delta)*theta*
                           (-1+(1+S2+S2*delta)*theta)+D2^delta*
                            (1+theta*(-3+2*theta+S2*(1+delta)*
                                     (-2+(2+S2*delta)*theta))))*ub1^theta))/
                           (S1^3*(1+delta*theta+(-S1)^(2/delta)*(1+delta)*theta-
                                  (-S1)^(1/delta)*(1+theta+2*theta*delta))*ub1)

            ## FIXME: No idea where is the error, it is just has the opposite
            ## sign from the numerical result. Assuming the numerical result is
            ## right at the moment. Need further investigation.

            ## out <- gradCpl.u
            out <- -gradCpl.u

          }
      }
    return(matrix(out))
  }
