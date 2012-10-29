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
logCplGrad <- function(CplNM, u, parCpl, cplCaller, staticArgs)
  {
    if(tolower(CplNM) == "bb7")
      {
        ## The name of marginal model
        MargisNM <- dimnames(u)[[2]]
        nObs <- dim(u)[1]
        ## Subtract the parameters list.
        tau <- parCpl[["tau"]]
        lambdaL <- parCpl[["lambdaL"]]
        lambdaU <- kendalltauInv(CplNM = CplNM,
                                 parRepCpl = parCpl,
                                 tauTabular = staticArgs[["tauTabular"]])
        ## The standard copula parameters (recycled if necessary, should be a vector).
        delta <- as.vector(-log(2)/log(lambdaL))
        theta <- as.vector(log(2)/log(2-lambdaU)) # ff(delta)

        if(tolower(cplCaller) == "lambdal")
          {
            ## FIXME: This part needs review when the conditional linkage used.

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

            ## Gradient w.r.t. tau, i.e.  ff'(delta)
            tauGrad.delta <- kendalltauGrad(CplNM = CplNM,
                                            theta = theta,
                                            delta = delta,
                                            caller = "delta")

            ub <- 1-u
            ub1 <- ub[, 1, drop = FALSE]
            ub2 <- ub[, 2, drop = FALSE]

            M12 <- 1-ub^theta

            M34 <- -log(M12) + ub^theta*log(ub)*tauGrad.delta/M12

            M5 <- -1 + rowSums(M12^(-delta))
            M6 <- -1 + M5^(1/delta)


            M7 <- (log(M5^(1/delta))-log(M6))*tauGrad.delta/theta^2
            C1 <- rowSums((M12^(-delta)*M34))M5
            M8 <- C1*delta + log(M5)
            M9 <- M7-(M8*(-1+1/theta))/(M6*delta^2)

            P12 <- M6*(1+delta)*theta*log(ub)*tauGrad.delta
            P34 <- -M12*log(M12)+(1+delta)log(ub)*ub^theta*tauGrad.delta
            P56 <- (-1+theta)*log(ub)*tauGrad.delta
            P78 <- M6*P34*(1+delta)*theta/M12

            C34 <- P34*(-1+theta)/M12

            ############################
            S1 <- M12^delta

            gradCpl.lambdaL <- (
                (S1[, 1]*S1[, 2])^(-2)*(rowSums(M12^delta)-S1[, 1]*S1[, 2])^2*
                (
                    rowSums(C34+P12+P78) - M6*theta/delta +
                    M5^(1/delta)*(1-M5^(-1/delta))*M9*(1+delta)*theta+
                    M6*(1+delta)*theta/delta-
                    (M7-M8*(-2+1/theta)/(M6*delta^2))*(-1+1/theta)*theta-
                    2*(-1+1/theta)*theta*(-C1*delta*(1+delta)+log(M5))/delta^2+
                    M5^(1/delta)*(1-M5^(-1/delta))*(1+delta)*theta*
                    (C1*(-2-1/delta)+log(M5)/delta^2)+ M6*(1+delta)*tauGrad.delta+
                    tauGrad.delta/theta+(-1+theta*tauGrad.delta)/theta
                    )
            )/(
                M5^2*(-1+(-delta+M5^(1/delta)*(1+delta))*theta)
                )


            ## Gradient w.r.t. lambdaL
            gradCpl.lambdaL <- 2^(-1/delta)*log(2)/delta^2

            ## The chain gradient
            out <- logGradCpl.delta/gradCpl.lambdaL
          }
        else if(tolower(cplCaller) == "tau")
          {
            ## Gradient w.r.t theta
            T1 <- 1-(1-u)^theta
            T2 <- (1-u)^(theta-1)
            L1 <- rowSums(T1^(-delta))-1
            Delta2.A <- -rowSums(T1^(-1)*(1-u)^theta*log(1-u))
            Delta2.B <- -rowSums(T1^(-delta-1)*(1-u)^theta*log(1-u))
            Delta2.C <- -rowSums(T1^(1/delta-1)*(1-u)^theta*log(1-u))
            Delta2.D <- -rowSums(T1^(-1/delta-1)*(1-u)^theta*log(1-u))
            ## Delta3 <- rowSums(log(1-u))
            logGradCpl.theta <- -(1+theta)*Delta2.A + rowSums(log(1-u))+
              2*(1+delta)*Delta2.B/L1+
                (1+2*theta)*Delta2.C/((1-L1^(-1/delta))*theta*delta)-
                  log(1-L1^(-1/delta))/theta^2+
                    ((1+delta)*L1^(1/delta)+
                     theta*(1+delta)/delta*L1^(1/delta-1)*Delta2.D-delta)/
                       ((1+delta)*theta*L1^(1/delta)-theta*delta-1)

            ## Gradient w.r.t. tau
            gradCpl.tau <- kendalltauGrad(CplNM = CplNM,
                                          theta = theta,
                                          delta = delta,
                                          caller = "theta")

            ## The chain gradient
            out <- logGradCpl.theta*gradCpl.tau
          }
        else
          {
            ## Gradient w.r.t u. NOTE: The BB7 copula's marginal are
            ## exchangeable which means the expression for the gradient w.r.t u2
            ## and u2 are the same if swap u1 and u2
            if(tolower(cplCaller) == "u1")
              {
                u <-  u[, 1:2, drop = FALSE]
              }

            else if(tolower(cplCaller) == "u1")
              {
                u <-  u[, 2:1, drop = FALSE]
              }
            else
              {
                stop("No such copula parameter!")
              }

            T1 <- 1-(1-u)^theta
            T2 <- (1-u)^(theta-1)
            L1 <- rowSums(T1^(-delta))-1

            Delta4.A <- -rowSums(T1^(-1)*(1-u)^(theta-1)*theta)
            Delta4.B <- -rowSums(T1^(-delta -1)*(1-u)^(theta-1)*theta)

            gradCpl.u <- (1+delta)*theta*Delta4.A+
              (1-theta)*(colSums(1/(1-u)))-2*(1+delta)*Delta4.B/L1-
                (1/theta-2)*L1^(-1/delta-1)*Delta4.B/(1-L1^(-1/delta))-
                  (1+delta)*theta*L1^(1/theta-1)*Delta4.B/
                    ((1+delta)*theta*L1^(1/delta)-theta*delta-1)


            out <- gradCpl.u
          }
      }
    return(matrix(out))
  }
