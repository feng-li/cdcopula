##' Gradient for log copula function
##'
##'
##' @title Log copula gradient
##' @param CplNM
##' @param u
##' @param parCpl
##' @param cplCaller
##' @return
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Fri May 11 12:42:20 CEST 2012;
##'       Current: Fri May 11 12:42:30 CEST 2012.
logCplGrad <- function(CplNM, u, parCplRep, cplCaller, Mdl.X, Mdl.beta)
  {

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
        tau <- as.vector(parCplRep[["tau"]])
        lambdaL <- as.vector(parCplRep[["lambdaL"]])

        lambdaU <- as.vector(kendalltauInv(
            CplNM = CplNM,
            parRepCpl = parCpl))

        ## The standard copula parameters (recycled if necessary, should be a
        ## vector).

        delta <- -log(2)/log(lambdaL)
        theta <- log(2)/log(2-lambdaU) # ff(delta)

        if(tolower(cplCaller) == "lambdal")
            {
                T1 <- 1-(1-u)^theta
                Tu1 <- T1[, 1]
                Tv1 <- T1[, 2]

                L1 <- -1 + Tu1^(-delta)+Tv1^(-delta)

                L34 <- -1 + T1^delta
                L3 <-  L34[, 1]
                L4 <-  L34[, 2]
                D12 <- T1^(-2*delta)
                D1 <- D12[, 1]
                D2 <- D12[, 2]
                L5 <- Tu1^delta-L3*Tv1^delta

                logGradCpl.delta <- (
                    L5^2*D1*D2*(
                        (-1+L1^(1/delta))^2*delta^2*theta^2+log(L1)-
                            1/L5*(-L5*theta*(delta+L1^(2/delta)*(1+delta)*theta-
                                                 L1^(1/delta)*(3+delta+(-1+delta)*theta)
                                             )*log(L1) +
                                                 delta*(
                                  -L1^(2/delta)*(1+delta)*
                                  (L4*Tu1^delta*delta+Tv1^delta*(1+delta))*theta^2+
                                  (1+delta*theta)*(Tu1^delta*delta*theta-
                                                   Tv1^delta*(1+(1+Tu1^delta)*delta*theta))+
                                  L1^(1/delta)*theta*(
                                      L4*Tu1^delta*delta*(1+theta+2*delta*theta)
                                      + Tv1^delta*(3-theta+2*delta*(1+theta+theta*delta))
                                      ))*log(Tu1)+
                              delta*(
                                  -L1^(2/delta)*(1+delta)*
                                  (L3*Tv1^delta*delta+Tu1^delta*(1+delta))*theta^2+
                                   (1+delta*theta)*(
                                       -L3*Tv1^delta*delta*theta-Tu1^delta*(1+delta*theta))+
                                       L1^(1/delta)*theta*(
                                           L3*Tv1^delta*delta*(1+theta+2*delta*theta)+
                                           Tu1^delta*(3-theta+2*delta*(1+theta+delta*theta))
                                           ))*log(Tv1)
                              )))/(L1^2*delta^2*theta*(
                                  1+delta*theta+L1^(2/delta)*(1+delta)*theta-
                                  L1^(1/delta)*(1+theta+2*delta*theta)))


                gradCpl.tau.delta <- kendalltauGrad(
                    CplNM = CplNM, parCpl = parCpl, caller = "delta")

                out <- logGradCpl.delta*(1/gradCpl.tau.delta)
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
            ## tol <- .Machine$double.eps*1e8
            ## T1.bad1 <- (T1> (1-tol))
            ## if(any(T1.bad1))
            ##   {
            ##     T1[T1.bad1] <- 1 - tol
            ##     warning("Numerical unstable on T1,  adjusted on the cliff...",
            ##             immediate. = immediate.)
            ##   }
            ## T1.bad0 <- (T1<tol)
            ## if(any(T1.bad0))
            ##   {
            ##     T1[T1.bad0] <- tol
            ##     warning("Numerical unstable on T1,  adjusted on the cliff...",
            ##             immediate. = immediate.)
            ##   }

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
                      -2+theta-2*delta*theta+L1^(1/delta)*(1+2*delta)*theta)-
                          L1*PT1^(1+delta)*(C1*(-1+theta)*(
                              -1+theta-theta*delta
                              +L1^(1/delta)*(1+delta)*theta)+
                                  theta*(C2+delta+C2*delta*theta-
                                             L1^(1/delta)*(1+delta)*(1+C2*theta))))+
                                                 L1*(-1-delta*theta+
                                                         L1^(1/delta)*(1+delta)*theta)*
                   (rowSums(T1[, 2:1]*(T1+(1-u)^theta*(1+delta))*log(1-u))))


            ## Gradient w.r.t. tau
            gradCpl.tau.theta <- kendalltauGrad(
                CplNM = CplNM,
                parCpl = parCpl,
                caller = "theta")

            ## The chain gradient
            out <- logGradCpl.theta*(1/gradCpl.tau.theta)

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
            ## u <- matrix(c(0.2, 0.3), 1, )
            ## theta <- 1.5
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

            gradCpl.u <- -(D1^(-1-3*delta)*D2^(-2*delta)*
                          (rowSums(D12^delta)-D1^delta*D2^delta)^2*
                          (D1^(1+delta)*S1*(-1+theta)*
                               (1+(-S1)^(1/delta)*(-1+theta)-
                                    theta+S2^2*theta+S2^2*delta*theta)+
                           D2^(-delta)*(D1^delta*(-1+D2^delta)*S2*(1+delta)*theta*
                           (-1+(1+S2+S2*delta)*theta)+D2^delta*
                            (1+theta*(-3+2*theta+S2*(1+delta)*
                                     (-2+(2+S2*delta)*theta))))*ub1^theta))/
                           (S1^3*(1+delta*theta+(-S1)^(2/delta)*(1+delta)*theta-
                                  (-S1)^(1/delta)*(1+theta+2*theta*delta))*ub1)

            out <- gradCpl.u

          }
    }
    else if(tolower(CplNM) == "mvt")
        {
            ## The name of marginal model
            MargisNM <- dimnames(u)[[2]]
            nObs <- dim(u)[1]

            parCpl <- parCplRep2Std(CplNM = CplNM, parCplRep = parCplRep)
            df <- parCpl[["df"]] # n-by-1
            rho <- parCpl[["rho"]] # n-by-lq

            u.quantile <- qt(u, df)
            if(tolower(cplCaller) == "lambdal")
                {
                    ## CopulaDensity-MVT.nb

                    gradFun <- function(i, rho, df, u.quantile)
                        {
                            Sigma <- vech2m(rho[i, ])
                            v <- df[i]
                            x <- u.quantile[i, , drop = FALSE]
                            mu <- 0
                            q <- dim(Sigma)[1]

                            C2 <- t(x-mu)%*%Sigma%*%(x-mu)
                            logGradCpl.v <- (C2-q -(C2+v)*log((C2+v)/2)
                                              +(C2+v)*(-digamma(v/2)+digamma((q+v)/2)))/
                                                  (2*(C2+v))
                            retun(logGradCpl.v)
                        }
                    logGradCpl.df <- apply(matrix(1:nObs), 1,
                                           gradFun,
                                           rho = rho,
                                           df = df,
                                           u.quantile = u.quantile) # n-by-q


                    C1 <- sqrt(1-rho)/sqrt(1+rho) # n-by-lq
                    gradCpl.lambda.df <- (1/(1 + C1^2))^((1 + df)/2)*
                        (-1 + df*log(1/(1 + C1^2)) - df*digamma(df/2)
                         + df*digamma((1 + df)/2))/(2*df^(3/2)*beta(df/2, 1/2)) #n-by-lq

                    ## The chain gradient
                    out <- logGradCpl.df*(1/gradCpl.lambda.df) # n-by-lq

                }
            else if(tolower(cplCaller) == "tau")
                {
                    ## FIXME: This is only for bivariate case.
                    u.quantile <- qt(u, df)

                    gradFun <- function(i, rho, df, u.quantile)
                        {
                            Sigma <- vech2m(rho[i, ])
                            v <- df[i]
                            x <- u.quantile[i, , drop = FALSE]
                            mu <- 0

                            logGradCpl.Sigma <-
                                -(v+p)/2*(1+1/v*t(x-mu)%*%Sigma%*%(x-mu))^(-1)*
                                    1/v*(x-mu)%*%t(x-mu)-1/2*solve(Sigma)

                            logGradCpl.rho <-
                                logGradCpl.Sigma[lower.tri(logGradCpl.Sigma,
                                                           diag = FALSE)]
                            retun(logGradCpl.rho)
                        }

                    logGradCpl.rho <- apply(matrix(1:nObs), 1,
                                            gradFun,
                                            rho = rho,
                                            df = df,
                                            u.quantile = u.quantile) # n-by-lq

                    gradCpl.tau.rho <- 2/(pi*sqrt(1-rho^2)) # n-by-lq

                    out <- logGradCpl.rho*(1/gradCpl.tau.rho) # n-by-lq

                }
            else
                {
                    ## Reorder the parameters.
                    q <- dim(u)[2]
                    imar <- substr(cplCaller, 2, nchar(cplCaller))

                    uIdx <- 1:q
                    uIdx[1] <- imar
                    uIdx[imar] <- 1

                    u <- u[, uIdx]
                    x <- u.quantile[, uIdx]
                    mu <- 0

                    ## The t copula and related copulas EQ.(6)
                    fx <- dt(x = x, df = df) # n-by-1
                    f1x <- -(x-mu)*df*df*(df/(df+(x-mu)^2))^(-1+(1+df)/2)/
                        ((df+(x-mu)^2)^2*sqrt(df)*beta(df/2, 1/2)) # n-by-1

                    F1X <- matrix(NA, nObs, q)
                    I0 <- (x<mu)
                    I <- !I0
                    if(any(I))
                        {
                            F1X.I <- (2*(x-mu)^3/((x-mu)^2+df)^2-2*(x-mu)/((x-mu)^2+df))*
                                (1-(x-mu)^2/((x-mu)^2+df))^(-1+df/2)/
                                    2*sqrt((x-mu)^2/((x-mu)^2+df))/beta(1/2, df/2)

                            F1X[I] <- F1X.I
                        }

                    if(any(I0))
                        {
                            F1X.I0 <- (x-mu)*df*(df/((x-mu)^2+df))^(-1+df/2)/
                                ((x-mu)^2+df)^2/sqrt(1-df/((x-mu)^2+df))/beta(df/2, 1/2)

                            F1X[I0] <- F1X.I
                        }


                    FUN <- function(i, x, mu, df, rho, uIdx, q)
                        {
                            Sigma0 <- vech2m(rho[i, ], diag = FALSE)
                            Sigma <- Sigma0[uIdx, uIdx]

                            v <- df[i]

                            C2 <- t(x-mu)%*%Sigma%*%(x-mu)
                            gradCpl.u <- -(v+q)/2*1/C2*
                                1/v*(2*Sigma%*%(x-mu))*(1/F1x)-
                                    1/fx*f1x/F1x
                            return(gradCpl.u)
                        }

                    gradCpl.u <- apply(matrix(1:nObs), 1, FUN,
                                       df = df, x = x, mu = mu,
                                       rho = rho, uIdx = uIdx, q = q)

        }

    return(matrix(out))
  }
