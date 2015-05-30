lambdaGrad <- function(CplNM, parCpl, caller)
  {
    if(tolower(CplNM) == "mvt")
      {
        df = parCpl[["df"]] # n-by-1
        rho <- parCpl[["rho"]] # n-by-lq

        if(tolower(caller) == "df")
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

            out <- apply(rho, 2, grad4df.bvt, df1 = df)

          }
        else
          {
            stop("No such caller.")
          }
      }
    else
      {
        stop("No such copula")
      }

    return(out)
  }
