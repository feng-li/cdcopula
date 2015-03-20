lambda <- function(CplNM, parCpl)
{
  out <- list()
  if(tolower(CplNM) == "mvt")
    {
      df = parCpl[["df"]] # n-by-1
      rho <- parCpl[["rho"]] # n-by-lq
      out[["lambdaL"]] <- 2*dt(-sqrt((df+1)*(1-rho)/(1+rho)),
                               df = df+1, log = FALSE)
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
