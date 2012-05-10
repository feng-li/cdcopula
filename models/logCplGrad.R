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
##' @note Created: ; Current: .
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
        lambdaU <- kendalltauInv(CplNM = CplNM, parRepCpl = parCpl,
                                tauTabular = staticArgs[["tauTabular"]])
        ## The standard copula parameters (recycled if necessary, should be a vector).
        delta <- as.vector(-log(2)/log(lambdaL))
        theta <- as.vector(log(2)/log(2-lambdaU))

        if(tolower(cplCaller) == "lambdal")
          {
            ## Gradient w.r.t delta
            T1 <- 1-(1-u)^theta
            L1 <- rowSums(T1^(-delta))-1
            Delta1 <- -rowSums(T1^(-delta)*log(T1))

            logGradCpl.delta <-  - rowSums(log(T1))-
              2*(1+delta)*Delta1/(delta*L1)-
                (1/theta-2)*(log(L1)-delta*Delta1/L1)/
                  (delta^2*(L1^(1/delta)-1))+
                    2*log(L1)/delta^2+
                      (L1^(1/delta)-(1+delta)*L1^(1/delta)*
                       (log(L1)-delta*Delta1/L1)/delta^2-1)/
                         ((1+delta)*L1^(1/delta)-delta-1/theta)

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
            gradCpl.tau <- kendalltauGrad(copula, theta, delta, parCaller)

            ## The chain gradient
            out <- logGradCpl.theta*gradCpl.tau
          }
        else if(tolower(cplCaller) == "u")
          {
            ## Gradient w.r.t u. NOTE: The BB7 copula's marginal are
            ## exchangeable which means the expression for the gradient w.r.t u
            ## and v are the same.

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
        else
          {
            stop("No such copula parameter!")
          }
      }
    return(out)
  }
