##' Compute the copula density
##'
##' @title Calculate the density function for give bivariate copula.
##'
##' @param u "matrix"
##'     Each column of the matrix is the value in the unit interval of the real
##' line for the corresponding margin of the copula function.
##' @param theta "list";
##'     The parameter list supplied to the copula function.
##' @param copula "character";
##'     This is the option for different types of copulas used in the function.
##' Currently "gaussian" for Gaussian copula.
##' @param par "list";
##'     Any further parameters need to pass to the copula. Which is equivalent
##' to the ... argument but less possible to make error.
##' @return "list" see below.
##'     \item {copula} {"matrix"; The copula function with the same length of
##' that in each column entries in "u".}
##'     \item {tao} {The theoretical Kendall's tao.}
##'     \item {rho} {The theoretical Spearman's rho.}
##' @references
##'     Trivedi and Zimmer (2007).
##' @author
##'     Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     Created: Mon Sep 12 13:36:10 CEST 2011;
##'     Current: Tue Sep 13 11:38:16 CEST 2011.
##' TODO: the first argument maybe should be x instead of u,
##'       Return log density instead of current form.
dCpl <- function(CplNM, u, parCpl, log = TRUE)
{
    if(any(!unlist(lapply(parCpl, is.finite))))
    {
        stop("NA/NaN/Inf are not allowed in copula parameters input.")
    }


    if(tolower(CplNM) == "bb7")
    {
        theta <- as.numeric(parCpl[["theta"]])
        delta <- as.numeric(parCpl[["delta"]])

        logDensFun <- function(u, theta, delta)
        {
            ## The density function
            TC1 <- (1-(1-u)^theta) # # FIXME: Numerically instable if theta -> Inf,  then TC1-> 1
            TC2.log <- (-1+theta)*log(1-u) # TC2 = (1-u)^(-1+theta)

            L5 <- (rowSums(TC1^(-delta)) - 1)
            L6 <- 1-L5^(-1/delta) # FIXME: log(L6)->Inf when u->1,  v->1

            logCplDensObs <- ((-1-delta)*rowSums(log(TC1))+
                              rowSums(TC2.log)-
                              2*(1+delta)/delta*log(L5)+
                              (-2+1/theta)*log(L6)+
                              log(-1+theta+L5^(1/delta)*L6*(1+delta)*theta))

            out.log <- matrix(logCplDensObs)
            return(out.log)
        }
        ## The usual log density
        out.log <- logDensFun(u = u, theta = theta, delta = delta)

        ## BB7 density is very unstable numerically. Use "Multiple Precision Floating-Point
        ## Reliable" based on GNU Multiple Precision Library for "those errors only (NA, NAN,
        ## Inf)" found in the result.
        redo.idx <- (!is.finite(out.log))
        if(any(redo.idx))
        {
            require("Rmpfr")
            precBits <- 1024
            ## MPFR class used for u, theta,  delta
            out.logredoMPFR <- logDensFun(u = mpfr(u[redo.idx, , drop = FALSE], precBits = precBits),
                                          theta = mpfr(theta[redo.idx], precBits = precBits),
                                          delta = mpfr(delta[redo.idx], precBits = precBits))
            out.logredo <- as.numeric(out.logredoMPFR)
            out.log[redo.idx] <- out.logredo

            if(any(!is.finite(out.logredo)))
                warning("MPFR used with insufficient ", precBits, " precBits in BB7 density.")
        }
    }
    else if(tolower(CplNM) == "gaussian")
    {
        rho <- parCpl[["rho"]] # n-by-lq

        ## The quantile for normal CDF
        u.quantile <- qnorm(u) # n-by-q

        nObs <- nrow(u.quantile)
        ## The CplNM density function C_12(u1, u2)

        dmvNormVecFun <- function(i, x, rho)
        {
            Sigma = vech2m(rho[i, ], diag = FALSE)
            out <- dmvnorm(x = x[i, , drop = FALSE],
                           sigma = Sigma,
                           log = TRUE)
            return(out)
        }
        logDensUpper <- apply(matrix(1:nObs),1,dmvNormVecFun,
                              x = u.quantile, rho = rho)

        logDensLower <- apply(dnorm(u.quantile, log = TRUE), 1, sum)
        logDens <- logDensUpper-logDensLower

        ## The output
        out.log <- matrix(logDens)
    }
    else if(tolower(CplNM) == "mvt") # The multivariate t-copula
    {## Demarta & McNeil (2005),  The t copula and related copulas
        require("mvtnorm")

        ## df, corr
        df <- parCpl[["df"]] # n-by-1
        rho <- parCpl[["rho"]] # n-by-lq

        ## The quantile for *univariate* t, that is x in t(x, df)
        u.quantile <- qt(u, df = df)
        ## u.quantile <- X

        ## The log copula density function C_12(u1, u2)
        nObs <- length(df)

        ## the formula right before formula (1.4) in Genz and Bretz (2009), also in
        ## Wikipedia.

        ## The density of the t copula, Demarta & Department (2006) Eq(6)

        logDensUpper <- matrix(NA, nObs, 1)

        for(i in 1:nObs)
        {
            logDensUpper[i] <- dmvt(x = u.quantile[i, , drop = FALSE],
                                    sigma = vech2m(rho[i, ], diag = TRUE),
                                    type = "shifted", # wikipedia type
                                    df = df[i], log = TRUE)
        }

        logDensLower <- apply(dt(u.quantile, df = df, log = TRUE), 1, sum)
        logDens <- logDensUpper-logDensLower

        ## The output
        if(any(is.infinite(logDens))) browser()
        out.log <- matrix(logDens)
    }
    else if(tolower(CplNM) == "fgm")
    {
        u1 <- u[, 1]
        u2 <- u[, 2]
        density <- 1+theta*(1-2*u1)*(1-2*u2)
        out <- density
    }
    else if(tolower(CplNM) == "gumbel")
    {# Joe 1997. p.142 Family B6
        delta <- as.vector(parCpl[["delta"]])
        pctl.log <- pCpl(u = u, parCpl = parCpl, CplNM = "gumbel", log = TRUE)
        u.tilde <- -log(u)

        u.tildeSumdelta <- rowSums(u.tilde^delta)

        out.log <- (pctl.log+u.tilde[, 1]+u.tilde[, 2]+
                    (delta-1)*(log(u.tilde[, 1])+log(u.tilde[, 2]))-
                    (2-1/delta)*log(u.tildeSumdelta)+
                    log(u.tildeSumdelta^(1/delta)+delta-1))
    }
    else if(tolower(CplNM)  == "frank")
    {# Joe 1997. p.142 Family B3
        delta <- as.vector(parCpl[["delta"]]) # delta >= 0

        u1 <- u[, 1]
        u2 <- u[, 2]
        ## The copula function
        eta <- 1-exp(-delta)
        out.log <- (log(delta)+log(eta)-delta*(u1+u2)-
                    2*log(eta-(1-exp(-delta*u1))*-(1-exp(-delta*u1))))
    }
    else if(tolower(CplNM)  == "clayton")
    { # Joe 1997. p 141 Family B4
        delta <- as.vector(parCpl[["delta"]]) # delta >= 0
        u1 <- u[, 1]
        u2 <- u[, 2]
        out.log <- (log(1+delta) + (-delta -1)*(log(u1)+ log(u2)) +
                    (-2-1/delta)*log(u1^(-delta) + u2^(-delta) -1))
    }
    else
    {
        stop("Given copula name is not implemented.")
    }



    if(log)
    {
        out <- out.log
    }
    else
    {
        out <- exp(out.log)
    }

    return(out)
}
