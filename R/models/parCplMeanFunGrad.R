parCplMeanFunGrad <- function(CplNM, Mdl.par, Mdl.parLink, chainCaller)
  {
    if(tolower(CplNM) == "bb7")
      {
        CompCaller <- chainCaller[1]
        parCaller <- chainCaller[2]



        ## The linkage for all parameters

        if(tolower(parCaller) == "tau")
          {
            ## The conditional links
            ## TODO: link function for tau are actually connected with
            ## lambdaL. Therefore the whole link gradient should be very
            ## complicated. Here we assume the condition is known a priori
            ## which is a bit sloppy.

            linkCurr <- Mdl.parLink[[CompCaller]][[parCaller]]

            ## Construct the extArgs
            if(tolower(linkCurr) == "glogit")
              {
                ## tau <- as.numeric(Mdl.par[[CplNM]][["tau"]])

                ## ## The upper and lower bounds are dynamic
                ## a <- 0 ## The lower bound
                ## b <- 2^(1/2-1/(2*tau)) ## the upper bound

                lambdaL <- Mdl.par[[CplNM]][["lambdaL"]]

                a <- log(2)/(log(2)-log(lambdaL))
                b <- 1 ## the upper bound


                extArgs <- list(a = a, b = b)
              }
            else
              {
                extArgs <- NA
              }

            ## The gradient
            LinkGradRaw <-  parMeanFunGrad(
                par = Mdl.par[[CompCaller]][[parCaller]],
                link = linkCurr,
                extArgs = extArgs)
            LinkGradObs <- LinkGradRaw
          }
        else
          {
            ## The unconditional gradients for the linkage
            LinkGradObs <-  parMeanFunGrad(
                par = Mdl.par[[CompCaller]][[parCaller]],
                link = Mdl.parLink[[CompCaller]][[parCaller]],
                extArgs = NA)
          }
        out <- LinkGradObs
      }
    return(out)
  }
