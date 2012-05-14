##' Link function constrains for copula models.
##'
##'
##' @title Constrains of link functions
##' @param CplNM
##' @return
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Fri May 11 09:34:09 CEST 2012;
##'       Current: Fri May 11 09:34:16 CEST 2012.
##' TODO: Write it in a more elegant way
CplLinkConstrain <- function(CplNM, Mdl.X,  Mdl.parLink, Mdl.beta,
                             parUpdate, Mdl.par)
  {
    if(tolower(CplNM) == "bb7")
      {
        ## BB7 copula (tau, lambdaL) representation requires that
        ## 0 < lambdaL < 2^(1/2-1/(2tau)) See Li 2012 for the proof.
        ## The condition should be calculated in the end for safety.

        CompNM <- names(Mdl.beta)
        ## First update all the independent linkages
        for(CompCurr in CompNM)
          {
            parUpdateNM <- names(parUpdate[[CompCurr]] == TRUE)
            for(ParCurr in parUpdateNM)
              {
                ## Check if particular constrain is needed.
                if(!(tolower(ParCurr) %in% c("lambdal")))
                  {
                    extArgs <- NA
                    ## Update the parameters for the updated part
                    Mdl.par[[CompCurr]][[ParCurr]] <-
                      parMeanFun(X = Mdl.X[[CompCurr]][[ParCurr]],
                                 beta = Mdl.beta[[CompCurr]][[ParCurr]],
                                 link = Mdl.parLink[[CompCurr]][[ParCurr]],
                                 extArgs = extArgs)
                  }
              }
          }

        ## Check if update the conditional parameters
        ## Special case for "lambdaL" as it is conditional on "tau"
        if(parUpdate[[CplNM]][["tau"]] == TRUE ||
           parUpdate[[CplNM]][["lambdaL"]] == TRUE)
          {
            linkCurr <- Mdl.parLink[[CplNM]][["lambdaL"]]
            XCurr <- Mdl.X[[CplNM]][["lambdaL"]]
            betaCurr <- Mdl.beta[[CplNM]][["lambdaL"]]

            if(tolower(linkCurr) == "glogit")
              {
                tau <- Mdl.par[[CplNM]][["tau"]]
                a <- 0 ## The lower bound of generalized logit link
                b <- 2^(1/2-1/(2*tau)) ## the upper bound

                extArgs <- list(a = a, b = b)
              }
            else
              {
                stop("Such conditional not implemented!")
              }

            Mdl.par[[CplNM]][["lambdaL"]] <-
              parMeanFun(X = XCurr,
                         beta = betaCurr,
                         link = linkCurr,
                         extArgs = extArgs)
          }
        out <- Mdl.par
      }
    return(out)
  }
