#' Compute the Copula function for given margins
#'
#' The function is used for bivariate copula C(U1<u2, U2<u1)
#' @title Copula distribution
#' @param u "matrix".
#' @param CplNM  "character".
#'
#'        The copula name.
#'
#' @param parCpl "list"
#'     Any additional parameters needed in the copula. In the
#' t-copula, parCpl$df: indicates the degrees of freedom.
#' @param copula "character"
#' @return "vector"
#' @references Nelsen 2006
#' @author
#'     Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @note
#'     DEPENDS: mvtnorm
#'     Created: Mon Sep 26 13:54:13 CEST 2011; Current: Tue Jul 14 13:06:55 CST 2015.
#' @export
pCpl <- function(u, CplNM, parCpl, log = FALSE)
{
    ## Rotation check
    rotation.check <- u2rotation(u = u, CplNM = CplNM)
    u.r <- rotation.check[["u.r"]]
    CplNM0 <- rotation.check[["CplNM0"]] # The unrotated copula
    rotation <- rotation.check[["rotation"]]

    if(tolower(CplNM0) == "bb7")
    {
        p <- 2
        theta <- parCpl[["theta"]]
        delta <- parCpl["delta"]
        u1 <- u.r[, 1]
        u2 <- u.r[, 2]
        percentile <- (1 - (1 - ((1 - (1 - u1)^theta)^(-delta) +
                                 (1 - (1 -  u2)^theta)^(-delta) - 1)^(-1/delta))^(1/theta))
        out <- matrix(percentile)
    }
    else if(tolower(CplNM0) == "sjc")
    {
        ## The symmetric Joe-Clayton copula
        out <- 0.5*(pCpl(u.r, CplNM = "BB7", parCpl = parCpl, log = FALSE) +
                    pCpl(1-u.r, CplNM = "BB7", parCpl = parCpl, log = FALSE) +
                    rowSums(u.r)-1)
    }
    else if(tolower(CplNM0) == "gaussian")
    {
        theta <- parCpl[["theta"]]
        u.quantile <- qnorm(u.r)
        p <- 2 # Hard code
        corr <- matrix(c(1, theta, theta, 1), 2, 2) # The theta is the
                                        # correlation
        percentile <- apply(u.r.quantile, 1, function(x) pmvnorm(lower = c(-Inf, -Inf),
                                                               upper = x, corr = corr)[1])
        ## TODO: apply pmvnorm is very slow
        out <- matrix(percentile)
    }
    else if(tolower(CplNM0) == "mvt")
    {
        p <- ncol(u.r)
        df <- parCpl[["df"]]
        ## theta is the vector for the lower triangular correlation matrix. P is
        ## the full correlation matrix.  TODO: Better construction of the
        ## correlation matrix.
        P.mat <- matrix(1, p, p) # The correlation matrix
        P.mat[lower.tri(P.mat)] <- theta
        P.mat[upper.tri(P.mat)] <- t(P.mat)[upper.tri(P.mat)]
        corr <- P.mat # The covariance matrix with scale 1.
        u.quantile <- qt(u.r, df = df)
        percentile <- apply(u.quantile, 1,
                            function(x) pmvt(lower = rep(-Inf, p),
                                             type = "shifted",
                                             upper = x, corr = corr, df = df)[1])
        ## TODO: apply pmvt is very slow
        out <- matrix(percentile)
    }
    else if(tolower(CplNM0) == "fgm")
    {
        u1 <- u.r[, 1]
        u2 <- u.r[, 2]
        percentile.log <- log(u1)+log(u2)+log(1+theta*(1-u1)*(1-u2))
        out.log <- matrix(percentile.log)
        if(log)
        {
            out <- out.log
        }
        else
        {
            out <- exp(out.log)
        }

    }
    else if(tolower(CplNM0) == "gumbel")
    {# Joe 1997. p.142 Family B6
        delta <- as.vector(parCpl[["delta"]])
        u.tilde <- -log(u.r)
        percentile.log <- -rowSums(u.tilde^delta)^(1/delta)
        out.log <- matrix(percentile.log)

        if(log)
        {
            out <- out.log
        }
        else
        {
            out <- exp(out.log)
        }
    }
    else if(tolower(CplNM0) == "frechet-upper")
    {
        out <- apply(u.r, 1, min)
    }
    else if(tolower(CplNM0) == "frechet-lower")
    {
        d <- ncol(u.r)
        z <- cbind((rowSums(u.r)+1-d), 0)
        out <- apply(z, 1, max)
    }
    else
    {
        stop("Given copula is not implemented.")
    }

    ## Rotate copula
    out.rotated <- pCpl2rotation(percentile = out, rotation = rotation, u = u)

    return(out.rotated)
}

#' @export
pCpl2rotation <- function(percentile, rotation, u)
{
    if(rotation  == 0)
    {
        ## Brechmann & Schepsmeier (2011) Modeling dependence with C- and D-vine copulas:
        ## The R-package CDVine
        out <- percentile
    }
    else if (rotation  == 90)
    {
        out <- u[, 2] - percentile
    }
    else if (rotation  == 180)
    {
        out <- u[, 1] + u[, 2] - 1 + percentile
    }
    else if (rotation  == 270)
    {
        out <- u[, 1] - percentile
    }
    else
    {
        stop("No such copula rotation!")
    }

    return(out)
}
