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
CplLinkConstrain <- function(CplNM, Mdl.X,  Mdl.parLink, Mdl.beta, parUpdate,
                             Mdl.par)
  {
    if(tolower(CplNM) == "bb7")
      {
        ## BB7 copula (tau, lambdaL) representation requires that
        ## 0 < lambdaL < 2^(1/2-1/(2tau)) See Li 2012 for the proof.

        CompNM <- names(Mdl.beta)
        for(CompCurr in CompNM)
          {
            parUpdateNM <- names(parUpdate[[CompCurr]] == TRUE)
            for(ParCurr in parUpdateNM)
              {
                ## Check if particular constrain is needed.
                if(tolower(ParCurr) == "lambdaL")
                  {
                    tau <- Mdl.par[[CplNM]][["tau"]]
                    a <- 0
                    b <- 2^(1/2-1/(2*tau))
                    extArgs <- list(a = a, b = b)
                  }
                else
                  {
                    extArgs <- NA
                  }

                ## Update the parameters for the updated part
                Mdl.par[[CompCurr]][[ParCurr]] <-
                  parMeanFun(X = Mdl.X[[CompCurr]][[ParCurr]],
                             beta = Mdl.beta[[CompCurr]][[ParCurr]],
                             link = Mdl.parLink[[CompCurr]][[ParCurr]],
                             extArgs = extArgs)
              }
          }

        out <- Mdl.par
      }
    return(out)
  }
