lambdaGrad <- function(CplNM, parCpl, caller)
  {
    out <- list()

    if((tolower(CplNM) == "bb7"))
      {
        delta <- parCpl[["delta"]] # n-by-lq

        if( "delta" %in% tolower(caller))
          {
            out[["delta"]] <- 2^(-1/delta)*log(2)/delta^2
          }
      }
    else if(tolower(CplNM) == "mvt")
      {
        df = parCpl[["df"]] # n-by-1
        rho <- parCpl[["rho"]] # n-by-lq

        if( "df" %in% tolower(caller))
          {
            grad4df.bvt <- function(df1, rho1)
              {
                out <- 1/2*(
                        -1/beta((1+df1)/2, 1/2)*2^(1/2*(-1-df1))*(1+rho1)^((1+df1)/2)*
                        gamma((1+df1)/2)^2*
                        regghypergeo(cbind((1+df1)/2, (1+df1)/2, 1/2),
                                     cbind((3+df1)/2, (3+df1)/2),
                                     (1+rho1)/2)+
                        ibeta((1+rho1)/2, (1+rho1)/2, 1/2, reg = TRUE)*
                        (-harmonic((-1+df1)/2)+ harmonic(df1/2) + log((1+rho1)/2)))
                return(out)
              }
            out[["df"]] <- apply(rho, 2, grad4df.bvt, df1 = df)
          }

        if( "rho" %in% tolower(caller))
          {
            grad4rho.bvt <- function(df1, rho1)
              {
                out <- 2^(-df1/2)*(1+rho1)^((-1+df1)/2)/
                (sqrt(1 - rho1)*beta((1+df1)/2, 1/2))
                return(out)
              }
            out[["rho"]] <- apply(rho, 2, grad4rho.bvt, df1 = df)
          }

      }
    else
      {
        stop("No such copula.")
      }

    return(out)
  }
