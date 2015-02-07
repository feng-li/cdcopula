lambda <- function(CplNM, parCpl)
{

    out <- list()
    if(tolower(CplNM) == "mvt")
        {
            df = parCpl[["df"]]
            rho <- parCpl[["rho"]]
            out[["lambdaL"]] <- 2*dt(-sqrt(df+1)*sqrt(1-rho)/sqrt(1+rho), df = df+1)
        }
    else if (tolower(CplNM) == "bb7")
        {
            if(exists("delta"))
                {
                    out[["lambdaL"]] <- 2^{-1/delta}
                }

            if(exists("theta"))
                {
                    out[["lambdaU"]] <- 2-2^{1/theta}
                }


        }
    else
        {
            stop("No such copula!")
        }

    return(out)
}
