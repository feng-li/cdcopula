parCplMeanFunGrad <- function(CplNM, Mdl.par, Mdl.parLink, chainCaller)
  {
    if(tolower(CplNM) == "bb7")
      {
        CompCaller <- chainCaller[1]
        parCaller <- chainCaller[2]

        linkCurr <- Mdl.parLink[[CompCaller]][[parCaller]]

        if(tolower(parCaller) == "tau")
          {
            ## The conditional links
            ## TODO: link function for tau are actually connected with
            ## lambdaL. Therefore the whole link gradient should be very
            ## complicated. Here we assume the condition is known a priori
            ## which is a bit sloppy.

            ## Construct the extArgs
            if(tolower(linkCurr[["type"]]) == "glogit")
              {
                ## tau <- as.numeric(Mdl.par[[CplNM]][["tau"]])

                ## ## The upper and lower bounds are dynamic
                ## a <- 0 ## The lower bound
                ## b <- 2^(1/2-1/(2*tau)) ## the upper bound

                ## The lower and upper bounds
                lambdaL <- Mdl.par[[CplNM]][["lambdaL"]]
                tau.a <- log(2)/(log(2)-log(lambdaL))
                linkCurr$a <- tau.a

                ## tau.b <- linkCurr$b  ## NOTE: Numerical stable. keep it slightly away from 1.
                ## linkCurr$b <- tau.b
              }
            else
              {
                stop("No such link available for tau")
              }

          }

        ## The unconditional gradients for the linkage
        LinkGradObs <-  parMeanFunGrad(
            par = Mdl.par[[CompCaller]][[parCaller]],
            linkArgs = linkCurr)
      }
    out <- LinkGradObs

    return(out)
  }
