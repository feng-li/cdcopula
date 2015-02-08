chainFracCplGrad <- function(CplNM, parCpl, chainCaller)
  {
    if(tolower(CplNM) == "bb7")
      {
        delta <- parCpl[["delta"]]
        theta <- parCpl[["theta"]]

        if(all.equal(tolower(chainCaller), c("delta", "lambdal")))
          {
            out <-  1/(2^(-1/delta)*log(2)/delta^2)
          }
        else if(all.equal(tolower(chainCaller), c("delta", "tau")))
          {
            out <- 1/(kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
                                          parCaller = "delta"))
          }
        else if(all.equal(tolower(chainCaller), c("theta", "tau")))
          {
            out <- 1/(kendalltauGrad(CplNM = CplNM, parCpl = parCpl,
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
