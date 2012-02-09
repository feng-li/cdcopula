chainFracCplGrad <- function(copula, parCpl, chainCaller)
  {
    if(tolower(copula) == "bb7")
      {
        delta <- parCpl[["delta"]]
        theta <- parCpl[["theta"]]

        if(all.equal(tolower(chainCaller), c("delta", "lambdal")))
          {
            out <-  1/(2^(-1/delta)*log(2)/delta^2)
          }
        else if(all.equal(tolower(chainCaller), c("delta", "tau")))
          {
            out <- 1/(kendalltauGrad(copula = copula, theta = theta,
                                          parCaller = "delta"))
          }
        else if(all.equal(tolower(chainCaller), c("theta", "tau")))
          {
            out <- 1/(kendalltauGrad(copula = copula, theta = theta,
                                     parCaller = "theta"))
          }
        else 
          {
            ## There is no such chain (reparametrization), i.e. estimate delta
            ## and theta directly. 
            out <- 1
          }
      }
    return(out)
  }
