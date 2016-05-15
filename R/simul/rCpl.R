##' Generate random variables using copulas
##'
##' The copula expression can be find at Trivedi and Zimmer 2005
##' @title Random *uniform* variables generator by bivariate copula
##' @param n
##' @param copula
##' @param parCpl "list" Any additional parameters input for the copula. For
##' example in the t copula, par$df: the degrees of freedom is needed.
##' @return "list" The random variables with some detailed information
##' @references Trivedi and Zimmer 2005
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     DEPENDS: mvtnorm
##'     Created: Mon Sep 26 09:43:12 CEST 2011;
##'     Current: Mon Apr 16 13:40:47 CEST 2012.
rCpl <- function(n, parCpl, CplNM, exArgs = NA)
{
    if(tolower(CplNM) == "bb7") # Joe 1997
    {
        ## Subtract the parameters list.
        delta <- parCpl[["delta"]]
        theta <- parCpl[["theta"]]


        ## Parameters healthy conditions TODO: more serious check
        ## hCond <- (theta[1] >= 1) && (theta[2]>0)
        ## if(!hCond) stop("BB7 copula should have: theta >= 1, delta >0.")

        p <- 2 #Hard code for bivariate copula
        v <- matrix(runif(n*p, 0, 1), n, p)
        u <- matrix(NA, n, p)
        u[, 1] <- v[, 1]


        ## The conditional method for sampling copula
        ## Sample v1 and v2 from U(0, 1)

        ## The conditional density equation C(u2|u1)  =  v2
        ## Not in vector form
        logcondEQ <- function(x, delta, theta, v, power)
        {
            u1 <- v[1]
            u2 <- x

            v2 <- v[2]
            L1 <- 1-(1-u1)^theta
            L2 <- 1-(1-u2)^theta
            L5 <- -1 + L1^(-delta) + L2^(-delta)

            L6 <- 1- L5^(-1/delta) # FIXME: log(L6)->Inf when u->1,  v->1.

            logcondDens <- -(1+delta)*log(L1) - (1+1/delta)*log(L5) +
                (-1+1/theta)*log(L6) + (-1+theta)*log(1-u1)

            out <- (logcondDens - log(v2))^power

            return(out)
        }


        logcondEQ.stable <- function(x, delta, theta, v, power = 1)
        {
            out <- logcondEQ(x = x, delta = delta, theta = theta, v = v,
                             power = power)
            if(!is.finite(out) || is.na(out))
            {
                require("Rmpfr")
                precBits <- 1024
                out.mpfr <- logcondEQ(x = mpfr(x, precBits = precBits),
                                      delta = mpfr(delta, precBits = precBits),
                                      theta = mpfr(theta, precBits = precBits),
                                      v = mpfr(v, precBits = precBits),
                                      power = power)
                out <- as.numeric(out.mpfr)

            }
            return(out)
        }



        ## Solve u2 from v2 = C(u2|u1) No direct solution, use nonlinear
        ## solver. TODO: Potential overflow when theta is big. Consider using
        ## Chris Sims solver.
        tol = 1e-10
        for(i in 1:n)
        {
            ## uniroot() is usually fast but instable for the boundary.
            u2.uniroot <- try(uniroot(f = logcondEQ.stable,
                                      interval = c(0+tol, 1-tol),
                                      theta = theta[i], delta = delta[i],
                                      v = v[i,], power = 1),
                              silent = TRUE)

            if(is(u2.uniroot, "try-error"))
            {
                u2.optim <- optim(par = 0.5, logcondEQ.stable,
                                  lower = 0.001, upper = 0.999,
                                  theta = theta[i], delta = delta[i],
                                  v = v[i,], power = 2,
                                  method = "L-BFGS-B")

                out.u2 <- u2.optim$par
                ## warning("NA produced in generating BB7 copula random numbers.")
            }
            else
            {
                out.u2 <- u2.uniroot[["root"]]
            }

            u[i, 2] <- out.u2
            ## Very simple weighted sampling with upper tail dependence as
            ## prob to sample extreme situations.
            ## u[i, 2] <- sample(x=c(runif(1, 0, v[i, 1]),
            ##                     runif(1, v[i, 1], 1)),
            ##                   size = 1,  prob = c(1-lambdaU, lambdaU))

        }

        ## Kendall's tau,  empirical
        emptau <- cor(u, method = "kendall")[2]

        ## Kendall's tau,  theoretical
        theotau <- kendalltau(CplNM=CplNM, parCpl = list(delta = delta, theta = theta))

        out <- list(u = u, emptau = emptau, theotau = theotau,
                    theta = theta, delta = delta)
    }
    else if(tolower(CplNM) == "gaussian")
    {
        ## Elliptical sampling,  Trivedi and Zimmer 2005 (p. 109)
        p <- 2 #Hard code for bivariate copula
        v <- matrix(rnorm(n*p, 0, 1), n, p)
        v1 <- v[, 1]
        v2 <- v[, 2]

        y <- matrix(NA, n, p)
        y[, 1] <- v1
        y[, 2] <- v1*theta + v2*sqrt(1-theta^2)
        u <- pnorm(y)

        ## Kendall's tau,  empirical
        emptau <- cor(u, method = "kendall")[2]

        ## Kendall's tau,  theoretical
        theotau <- kendalltau(CplNM=CplNM, parCpl = list(rho = theta))
        out <- list(u = u, emptau = emptau, theotau = theotau)
    }
    else if(tolower(CplNM) == "mvt") ## the multivariate t Demarta Mcneil (2005)
    {
        p <- (1+sqrt(8*length(theta)+1))/2 # p is obtained from the lower
                                        # triangular matrix.
        df <- par[["df"]]
        ## theta is the vector for the lower triangular correlation matrix. P is
        ## the full correlation matrix.  TODO: Better construction of the
        ## correlation matrix.
        P.mat <- matrix(1, p, p) # The correlation matrix
        P.mat[lower.tri(P.mat)] <- theta
        P.mat[upper.tri(P.mat)] <- t(P.mat)[upper.tri(P.mat)]
        corr <- P.mat # The covariance matrix with scale 1.

        x <- matrix(rmvt(n*p, df = df, sigma = corr, type = "shifted"), n, p)
        u <- pt(x, df = df) # The percentile

        ## Kendall's tau,  empirical
        emptau <- cor(u, method = "kendall")[2]

        ## Kendall's tau,  theoretical
        theotau <- kendalltau(CplNM=CplNM, parCpl = list(rho = corr))

        out <- list(u = u, emptau = emptau, theotau = theotau)

    }
    else if(tolower(CplNM) == "fgm")
    {
        theta <- parCpl[["theta"]]
        if(!(theta >= -1 && theta <= 1))
        {stop("FGM copula should have -1 <= theta <= 1.")}

        ## Conditional sampling Trivedi and Zimmer 2005(p. 108)
        p <- 2 #Hard code for bivariate copula
        v <- matrix(runif(n*p, 0, 1), n, p)
        v1 <- v[, 1]
        v2 <- v[, 2]

        u <- matrix(NA, n, p)
        u1 <- v1
        u[, 1] <- u1

        u2 <- (-1 - theta + 2*u1*theta +
               sqrt(4*(-1 + 2*u1)*v2*theta + (1 + theta - 2*u1*theta)^2))/
            (2*(-1 + 2*u1)*theta)

        u[, 2] <- u2

        ## Kendall's tau,  empirical
        emptau <- cor(u, method = "kendall")[2]

        ## Kendall's tau,  theoretical
        theotau <- kendalltau(CplNM=CplNM, parCpl = parCpl)

        out <- list(u = u, emptau = emptau, theotau = theotau)
    }
    else if(tolower(CplNM) == "gumbel")
    {
        ## Mixture of power simulation,  Trivedi and Zimmer 2005(p. 110)
        theta <- parCpl[["theta"]]

        alpha <- 1/theta
        gamma <- rps(n, alpha)
        v <-matrix(runif(2*n, 0, 1), n, 2)
        t <- -1/gamma*log(v)
        u <- exp(-t^(1/theta))

        emptau <- cor(u, method = "kendall")[2]

        out <- list(u = u, emptau = emptau)
    }
    else stop("Not implemented for given copula name.")
    return(out)
}
