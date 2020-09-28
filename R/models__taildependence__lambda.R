#' @export
lambda <- function(CplNM, parCpl)
{
    CplNM0 <- strsplit(CplNM, split = "_")[[1]][1]

    out <- list()
    if(tolower(CplNM0) == "mvt")
    {
        df = parCpl[["df"]] # n-by-1
        rho <- parCpl[["rho"]] # n-by-lq
        out[["lambdaL"]] <- 2*pt(-sqrt((df+1)*(1-rho)/(1+rho)),
                                 df = df+1, log = FALSE)
    }
    else if (tolower(CplNM0) == "bb7")
    {
        delta <- parCpl[["delta"]]
        theta <- parCpl[["theta"]]

        if(length(delta)>0)
        {
            out[["lambdaL"]] <- 2^{-1/delta}
        }

        if(length(theta)>0)
        {
            out[["lambdaU"]] <- 2-2^{1/theta}
        }
    }
    else if (tolower(CplNM0) == "clayton")
    { # Joe 1997, p.78,  Example 3.2
        delta <- parCpl[["delta"]]
        out[["lambdaL"]] <- 2^{-1/delta}

    }
    else if (tolower(CplNM0) == "gumbel")
    {
        delta <- parCpl[["delta"]]
        out[["lambdaU"]] <- 2-2^{1/delta}
    }

    else
    {
        stop("No such copula!")
    }

    return(out)
}
