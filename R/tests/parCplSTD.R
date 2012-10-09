##' <description>
##'
##' <details>
##' @title <short tile>
##' @param copula 
##' @param parRepCpl 
##' @return 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: ; Current: .
parCplSTD <- function(CplNM, parRepCpl)
  {
    if(tolower(CplNM) == "bb7")
      {
        ## The name for the reparametrized parameters
        CplParNM <- names(parRepCpl)
        ## The parameter names in the standard parametrization
        CplParSTDNM <- c("theta", "delta")

        ## Initialize output Storage and name it
        out <- parRepCpl
        names(out) <- CplParSTDNM
        
        lambdaL <- parRepCpl[["lambdaL"]]
        tau <- parRepCpl[["tau"]]


        ## The first parameter
        out[["delta"]] <- -log(2)/log(lambdaL)
        
        ## The second parameter
        
        out[["theta"]] <- kendalltauInv(CplNM = CplNM, parCpl = out,
                                                tau = tau, parCaller = "theta")
      }
    else
      {
        stop("This copula is not implemented!")
      }
    return(out)
  }
