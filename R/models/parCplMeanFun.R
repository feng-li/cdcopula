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
parCplMeanFun <- function(CplNM, Mdl.X,  Mdl.parLink, Mdl.beta,
                          parUpdate, Mdl.par)
  {
    if(tolower(CplNM) == "bb7")
      {
        ## BB7 copula (tau, lambdaL) representation requires that
        ## 0 < lambdaL < 2^(1/2-1/(2*tau)) See Li 2012 for the proof.
        ## The condition should be calculated in the end for safety.

        CompNM <- names(Mdl.beta)

        ## Which parameter are conditionally considered
        ## Hard coded, maybe should treat it as an input
        ## condPar <- c("lambdaL")
        condPar <- c("tau")

###----------------------------------------------------------------------------
### (1) update all the independent linkages
###----------------------------------------------------------------------------
        for(CompCurr in CompNM)
          {
            parUpdateNM <- names(parUpdate[[CompCurr]] == TRUE)
            for(ParCurr in parUpdateNM)
              {
                ## Check if particular constrain is needed.
                if(!(tolower(ParCurr) %in% tolower(condPar)))
                  {
                    ## Update the parameters for the updated part
                    Mdl.par[[CompCurr]][[ParCurr]] <- parMeanFun(
                        X = Mdl.X[[CompCurr]][[ParCurr]],
                        beta = Mdl.beta[[CompCurr]][[ParCurr]],
                        linkArgs = Mdl.parLink[[CompCurr]][[ParCurr]])
                  }
              }
          }

###----------------------------------------------------------------------------
### (2) Special case for conditional linkage
###----------------------------------------------------------------------------

        ## The parameter tau is updated individually. NOTE that the parameter
        ## tau depends on lambdaL. So when lambdal is updated, the information
        ## of tau should also be updated.

        if(parUpdate[[CplNM]][["tau"]] == TRUE |
           parUpdate[[CplNM]][["lambdaL"]] == TRUE)
          {
            ## linkCurr <- Mdl.parLink[[CplNM]][["lambdaL"]]
            ## XCurr <- Mdl.X[[CplNM]][["lambdaL"]]
            ## betaCurr <- Mdl.beta[[CplNM]][["lambdaL"]]


            linkCurr <- Mdl.parLink[[CplNM]][["tau"]]
            XCurr <- Mdl.X[[CplNM]][["tau"]]
            betaCurr <- Mdl.beta[[CplNM]][["tau"]]

            if(tolower(linkCurr[["type"]]) == "glogit")
              {
                ## tau <- Mdl.par[[CplNM]][["tau"]]
                ## tau.a <- 0 ## The lower bound of generalized logit link
                ## tau.b <- 2^(1/2-1/(2*tau)) ## the upper bound

                ## The lower and upper bounds of generalized logit link
                lambdaL <- Mdl.par[[CplNM]][["lambdaL"]]

                tau.a <- log(2)/(log(2)-log(lambdaL))
                tau.b <- linkCurr[["b"]] ## NOTE: Numerical stable. keep it slightly away
                ## from 1.

                linkCurr$a <- tau.a
                linkCurr$b <- tau.b

              }
            else
              {
                stop("Such conditional linkage is not implemented!")
              }


            Mdl.par[[CplNM]][["tau"]] <- parMeanFun(
                X = XCurr, beta = betaCurr, linkArgs = linkCurr)
          }
        out <- Mdl.par
      }

    return(out)
  }
