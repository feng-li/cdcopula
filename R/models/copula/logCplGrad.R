##' Gradient for log copula function
##'
##'
##' @title Log copula gradient
##' @param CplNM
##' @param u
##' @param parCpl
##' @param parCaller
##' @return
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Fri May 11 12:42:20 CEST 2012; Current: Fri Mar 27 17:47:58 CST 2015.
logCplGrad <- function(CplNM, u, parCpl, parCaller)
{
    out <- list()
    q <- dim(u)[2]

    if(tolower(CplNM) == "bb7")
    {
        ## The standard copula parameters (recycled if necessary, should be a vector).
        delta <- as.numeric(parCpl[["delta"]])
        theta <- as.numeric(parCpl[["theta"]]) # ff(delta)

        if("delta" %in% tolower(parCaller))
        {
            gradFun4delta <- function(u, delta, theta)
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

                out <- ((L5^2*D1*D2*((-1+L1^(1/delta))^2*delta^2*theta^2+log(L1)- 1/L5*(
                    -L5*theta*(delta+L1^(2/delta)*(1+delta)*theta-
                                        L1^(1/delta)*(3+delta+(-1+delta)*theta))*log(L1)
                    + delta*(-L1^(2/delta)*(1+delta)*
                                 (L4*Tu1^delta*delta+Tv1^delta*(1+delta))*theta^2
                        + (1+delta*theta)*(Tu1^delta*delta*theta- Tv1^delta*(1+(1+Tu1^delta)*delta*theta))+
                        L1^(1/delta)*theta*(L4*Tu1^delta*delta*(1+theta+2*delta*theta)
                            + Tv1^delta*(3-theta+2*delta*(1+theta+theta*delta))))*log(Tu1)
                    + delta*(-L1^(2/delta)*(1+delta)*(L3*Tv1^delta*delta
                        +Tu1^delta*(1+delta))*theta^2
                        +(1+delta*theta)*(-L3*Tv1^delta*delta*theta-Tu1^delta*(1+delta*theta))
                        +L1^(1/delta)*theta*(L3*Tv1^delta*delta*(1+theta+2*delta*theta)
                            + Tu1^delta*(3-theta+2*delta*(1+theta+delta*theta))))*log(Tv1)
                )))/(L1^2*delta^2*theta*(1+delta*theta+L1^(2/delta)*(1+delta)*theta-
                                                          L1^(1/delta)*(1+theta+2*delta*theta))))
                return(out)
            }

            gradout <- gradFun4delta(u = u, theta = theta, delta = delta)
            redo.idx <- (!is.finite(gradout))
            if(any(redo.idx))
            {
                require("Rmpfr", quietly = TRUE)
                precBits <- 1024
                ## MPFR class used for u, theta,  delta
                gradout.redoMPFR <- gradFun4delta(u = mpfr(u[redo.idx, , drop = FALSE], precBits = precBits),
                                              theta = mpfr(theta[redo.idx], precBits = precBits),
                                              delta = mpfr(delta[redo.idx], precBits = precBits))
                gradout.redo <- as.numeric(gradout.redoMPFR)
                gradout[redo.idx] <- gradout.redo
                if(any(!is.finite(gradout.redo)))
                    warning("MPFR used with insufficient ", precBits,
                            " precBits in BB7 gradient for ", parCaller)
            }
            out[["delta"]] <- gradout
        }

        if( "theta" %in% tolower(parCaller))
        {
################################################################################
### DEBUGGING
            ## u <- matrix(c(0.6, 0.3), 1, )
            ## theta <- 3.5
            ## delta <- 2.4
### PASSED
################################################################################

            gradFun4theta <- function(u, theta, delta)
            {
                ## Gradient w.r.t theta
                T1 <- 1-(1-u)^theta

                L1 <- rowSums(T1^(-delta))-1
                PT1 <- matrix(T1[, 1]*T1[, 2])

                SD12 <- rowSums((T1[, 2:1])^(1+delta)*(1-u)^theta*log(1-u))
                SD34 <- rowSums(T1^(-1-delta)*(1-u)^theta*delta*log(1-u))

                C1 <- SD34*L1^(-(1+delta)/delta)/(delta-L1^(-1/delta)*delta)
                C2 <- (log(L1^(1/delta))-log(-1+L1^(1/delta)))/theta^2

                out <- (1/(L1^3*(-1-delta*theta+L1^(1/delta)*(1+delta)*theta))*
                        PT1^(-1-2*delta)*(rowSums(T1^delta)-PT1^delta)^2*(
                            1/theta*PT1^(-delta)*(-SD12*(1+delta)*theta*(
                                -2+theta-2*delta*theta+L1^(1/delta)*(1+2*delta)*theta)-
                                L1*PT1^(1+delta)*(C1*(-1+theta)*(
                                    -1+theta-theta*delta
                                    +L1^(1/delta)*(1+delta)*theta)+
                                    theta*(C2+delta+C2*delta*theta-
                                           L1^(1/delta)*(1+delta)*(1+C2*theta))))+
                                L1*(-1-delta*theta+
                                    L1^(1/delta)*(1+delta)*theta)*
                                (rowSums(T1[, 2:1]*(T1+(1-u)^theta*(1+delta))*log(1-u)))))
                return(out)
            }

            ## The chain gradient
            gradout <- gradFun4theta(u = u, theta = theta, delta = delta)
            redo.idx <- (!is.finite(gradout))
            if(any(redo.idx))
            {
                require("Rmpfr", quietly = TRUE)
                precBits <- 1024
                ## MPFR class used for u, theta,  delta
                gradout.redoMPFR <- gradFun4theta(u = mpfr(u[redo.idx, , drop = FALSE], precBits = precBits),
                                                  theta = mpfr(theta[redo.idx], precBits = precBits),
                                                  delta = mpfr(delta[redo.idx], precBits = precBits))
                gradout.redo <- as.numeric(gradout.redoMPFR)
                gradout[redo.idx] <- gradout.redo

                if(any(!is.finite(gradout.redo)))
                    warning("MPFR used with insufficient ", precBits,
                            " precBits in BB7 gradient for ", parCaller)
            }
            out[["theta"]] <- gradout
        }

        if(any(paste("u", 1:q, sep = "") %in% tolower(parCaller)))
        {
            ## Gradient w.r.t u. NOTE: The BB7 copula's marginal are exchangeable which means
            ## the expression for the gradient w.r.t u1 and u2 are the same if swap u1 and u2
            imar <- as.numeric(substr(parCaller, 2, nchar(parCaller)))
            uIdx <- 1:q
            uIdx[1] <- imar
            uIdx[imar] <- 1
            u <- u[, uIdx]

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

            gradCpl.u <- (
                -(D1^(-1-3*delta)*D2^(-2*delta)*
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
            )
            out[["u"]] <- gradCpl.u
        }
        if(any(!sapply(out, is.finite)))
        {
            warning("Bad gradients with NA/NaN or Inf.")
        }
    }
    else if(tolower(CplNM) == "mvt")
    {
        nObs <- dim(u)[1]

        df <- parCpl[["df"]] # n-by-1
        rho <- parCpl[["rho"]] # n-by-lq

        u.quantile <- qt(u, df) # the x in t(x, df)
        if("df" %in% tolower(parCaller))
        { ## CopulaDensity-MVT.nb
            gradFun4df <- function(i, rho, df, u.quantile)
            {
                Sigma <- vech2m(rho[i, ], diag = FALSE)
                v <- df[i]
                x <- matrix(u.quantile[i, ]) # col-vector
                mu <- 0
                q <- dim(Sigma)[1]

                ## C2 <- as.vector(t(x-mu)%*%solve(Sigma)%*%(x-mu))
                C2 <- as.vector(t(x-mu)%*%solve(Sigma, (x-mu)))

                out <- ((C2-q -(C2+v)*log((C2+v)/v)+
                         (C2+v)*(-digamma(v/2)+digamma((q+v)/2)))/(2*(C2+v)))
                return(out)
            }

            ## The gradient of the t copula numerator,  Demarta & Department (2006) Eq(6)
            logCplGrad.df.upper <- apply(matrix(1:nObs), 1, gradFun4df,
                                         rho = rho, df = df, u.quantile = u.quantile) # n-by-1

            ## The gradient of the t copula denominator, (split-t in MargiModelGrad())
            logCplGrad.df.lowerMat <- apply(u.quantile, 2,
                                            function(y, par, type, parCaller, densCaller)
                                            {
                                                MargiModelGrad(y = y, par = par,
                                                               type = type,
                                                               parCaller = parCaller,
                                                               densCaller = densCaller)$d
                                            },
                                            par = list(mu = 0, df = df, phi = 1, lmd = 1),
                                            type = "splitt", parCaller = "df",
                                            densCaller = "d")

            logCplGrad.df.lower <- rowSums(logCplGrad.df.lowerMat)

            out[["df"]] <- matrix(logCplGrad.df.upper -logCplGrad.df.lower) # n-by-1
        }

        if("rho" %in% tolower(parCaller))
        {

            logCplGrad.rho <- matrix(NA, nObs, ncol(rho))  # n-by-lq
            for(i in 1:nObs)
            {
                Sigma <- vech2m(rho[i, ], diag = FALSE)
                v <- df[i]
                x <- matrix(u.quantile[i, ]) # col-vector
                mu <- 0
                p <- dim(Sigma)[1]

                ## C0 <- as.vector(t(x-mu)%*%solve(Sigma)%*%(x-mu))
                C1 <- solve(Sigma, (x-mu))
                C0 <- as.vector(t(x-mu)%*%C1)
                logGradCpl.Sigma <- (-1/2*solve(Sigma)-(v+p)/2*(1+C0/v)^(-1)*(-C1%*%t(C1))/v)

                logCplGrad.rho[i, ] <- logGradCpl.Sigma[lower.tri(logGradCpl.Sigma,
                                                                  diag = FALSE)]
            }
            out[["rho"]] <- logCplGrad.rho # n-by-lq
        }

        if(any(paste("u", 1:q, sep = "") %in% tolower(parCaller)))
        {
            ## The gradient with respect to u_i Reorder the parameters.
            imar <- as.numeric(substr(parCaller, 2, nchar(parCaller)))
            uIdx <- 1:q
            uIdx[1] <- imar
            uIdx[imar] <- 1

            u <- u[, uIdx]

            x <- u.quantile[, uIdx]
            x1 <- x[, 1]
            mu <- 0

            ## The t copula and related copulas EQ.(6)
            fx1 <- dt(x = x1, df = df, log = FALSE) # n-by-1
            f1x1 <- -(x1-mu)*df*(1+df)*(df/(df+(x1-mu)^2))^(-1+(1+df)/2)/
                                                               ((df+(x1-mu)^2)^2*sqrt(df)*beta(df/2, 1/2)) # n-by-q


            ## The first marginal CDF derivative with respect to x1.
            F1x1 <- matrix(NA, nObs, 1) # n-by-1
            I0 <- (x1<mu)
            I <- !I0
            if(any(I))
            {
                df.1 <- df[I]
                x1.1 <- x1[I]

                F1x1.I <- {
                    -(2*(x1.1-mu)^3/((x1.1-mu)^2+df.1)^2-2*(x1.1-mu)/
                                                       ((x1.1-mu)^2+df.1))*
                        (1-(x1.1-mu)^2/((x1.1-mu)^2+df.1))^(-1+df.1/2)/
                                                               (2*sqrt((x1.1-mu)^2/((x1.1-mu)^2+df.1))*
                                                                beta(1/2, df.1/2))
                }

                F1x1[I] <- F1x1.I
            }

            if(any(I0))
            {
                df.0 <- df[I0]
                x1.0 <- x1[I0]

                F1X.I0 <- {
                    -(x1.0-mu)*df.0*
                        (df.0/((x1.0-mu)^2+df.0))^(-1+df.0/2)/((x1.0-mu)^2+df.0)^2/
                                                                                     sqrt(1-df.0/((x1.0-mu)^2+df.0))/beta(df.0/2, 1/2)
                }
                F1x1[I0] <- F1X.I0

            }


            ## The gradient for copula with respect to x1.
            FUN <- function(i, x, mu, df, rho, uIdx)
            {
                Sigma0 <- vech2m(rho[i, ], diag = FALSE)
                Sigma <- Sigma0[uIdx, uIdx]

                if(!is.positivedefinite(Sigma))
                {
                    out <- NA
                }

                p <- dim(Sigma)[1]
                v <- df[i]
                x_i <- matrix(x[i, ])

                C2 <- as.vector(t(x_i-mu)%*%solve(Sigma)%*%(x_i-mu))
                out <- -(v+p)/2*C2^(-1)*
                                       1/v*(2*solve(Sigma)%*%(x_i-mu))
                return(out)
            }
            gradLogCpl.x1 <- t(apply(matrix(1:nObs), 1, FUN,
                                     df = df, x = x, mu = mu,
                                     rho = rho, uIdx = uIdx))[, 1] # n-by-1

            logCplGrad.u <- gradLogCpl.x1*(1/F1x1)- 1/fx1*f1x1/F1x1
            out[["u"]] <-  logCplGrad.u
        }
    }
    else if(tolower(CplNM) == "gumbel")
    {
        delta <- parCpl[["delta"]] # n-by-1
        uvt <- -log(u)

        u1t <- uvt[, 1]
        u2t <- uvt[, 2]

        uDelta <- u1t^delta + u2t^delta

        if("delta" %in% tolower(parCaller))
        {
            u1 <- u[, 1]
            u2 <- u[, 2]

            A1 <- (-1+uDelta^(1/delta)+delta)
            A2 <- uDelta*log(uDelta)
            A3 <- u1t^delta*log(u1t)
            A4 <- u2t^delta*log(u2t)
            A6 <- uDelta^(-1+1/delta)*(A2 - (A3 + A4)*delta)/delta^2

            logCplGrad.delta <-(-log(uDelta)/delta^2+log(u1t) + log(u2t)
                + (1-2*delta)*(A3+A4)/(uDelta*delta) + A6 + (1 -A6)/A1)

            out[["delta"]] <- logCplGrad.delta
        }


        if(any(paste("u", 1:q, sep = "") %in% tolower(parCaller)))
        {
            ## The gradient with respect to u_i
            ## Reorder the parameters.
            imar <- as.numeric(substr(parCaller, 2, nchar(parCaller)))
            uIdx <- 1:q
            uIdx[1] <- imar
            uIdx[imar] <- 1
            u <- u[, uIdx]

            u1 <- u[, 1]
            logCplGrad.u <- (-1/u1 + ((-1+delta - (1+uDelta^(2/delta)
                + 3*uDelta^(1/delta)*(-1+delta)
                -3*delta+2*delta^2)*(u1t)^delta/(
                    uDelta*(-1+uDelta^(1/delta)+delta)))/(-u1t*u1)))

            out[["u"]] <- logCplGrad.u
        }
    }
    else if(tolower(CplNM) == "clayton")
    {
        delta <- parCpl[["delta"]] # n-by-1
        u1 <- u[, 1]
        u2 <- u[, 2]

        if("delta" %in% tolower(parCaller))
        {
            A <- (-1 + u1^(-delta) + u2^(-delta))
            logCplGrad.delta <- (1/(1 + delta) - log(u1) - log(u2) +
                                 (-2 - 1/delta)*(-u1^(-delta)*log(u1)
                                     - u2^(-delta)*log(u2))/A +
                                 log(A)/(delta^2))
            out[["delta"]] <- logCplGrad.delta
        }

        if(any(paste("u", 1:q, sep = "") %in% tolower(parCaller)))
        {
            ## The gradient with respect to u_i
            ## Reorder the parameters.
            imar <- as.numeric(substr(parCaller, 2, nchar(parCaller)))
            uIdx <- 1:q
            uIdx[1] <- imar
            uIdx[imar] <- 1
            u <- u[, uIdx]

            A <- u1^delta*(-1+u2^delta)

            logCplGrad.u <- (-(u2^delta*delta + A*(1+delta))/(u1*(-u2^delta+ A)))

            out[["u"]] <- logCplGrad.u
        }
    }
    else
    {
        stop("No such copula defined!")
    }


    if(length(out) == 0)
    {
        stop("No such copula parameter defined or implemented.")
    }
    return(out)
}
