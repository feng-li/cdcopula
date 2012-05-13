CplLinkConstrainGrad <- function(CplNM, Mdl.par, Mdl.parLink, chainCaller)
  {
    if(tolower(CplNM) == "bb7")
      {
        CompCaller <- chainCaller[1]
        parCaller <- chainCaller[2]
        if(tolower(parCaller) == "lambdal")
          {
            ## The conditional links
            tau <- Mdl.par[[CplNM]][["tau"]]
            a <- 0 ## The lower bound of generalized logit link
            b <- 2^(1/2-1/(2*tau)) ## the upper bound
            extArgs <- list(a = a, b = b)

            LinkGradRaw <-  parMeanFunGrad(
                              par = Mdl.par[[CompCaller]][[parCaller]],
                              link = Mdl.parLink[[CompCaller]][[parCaller]],
                              extArgs = extArgs)
            LinkGradObs <-
          }
        else
          {
            ## The unconditional gradients for the linkage
            LinkGradObs <-  parMeanFunGrad(
                              par = Mdl.par[[CompCaller]][[parCaller]],
                              link = Mdl.parLink[[CompCaller]][[parCaller]]
                              extArgs = NA)
          }

        out <- LinkGradObs
      }
    return(out)
  }
