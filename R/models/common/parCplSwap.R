##' Swap the parameters between model structure and vector
##'
##' This function is used for numerical optimization in which it usually need
##' vector arguments.
##' @param betaInput "list or vector" depends on the contents.
##' @param Mdl.beta "list"
##' @param Mdl.betaIdx "list"
##' @param parUpdate "list"
##' @return NA
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Mar 14 13:36:29 CET 2013;
##'       Current: Thu Mar 14 13:36:35 CET 2013.
parCplSwap <- function(betaInput, Mdl.beta = NA, Mdl.betaIdx = NA, parUpdate = NA)
  {
    ## Convert Mdl.beta -> betaVec
    if(class(betaInput)  == "list")
      {
        for(CompCaller in names(Mdl.betaIdx))
          {
            for(parCaller in names(Mdl.betaIdx[[CompCaller]]))
              {
                ## delete current entry if not updated
                if(parUpdate[[CompCaller]][[parCaller]] == FALSE)
                  {
                    Mdl.beta[[CompCaller]][[parCaller]] <- NULL
                    Mdl.betaIdx[[CompCaller]][[parCaller]] <- NULL
                  }
              }
          }

        betaVec <- unlist(Mdl.beta)
        betaIdxVec <- unlist(Mdl.betaIdx)
        betaVec <- betaVec[betaIdxVec]

        ## names(betaVec) <- NULL
        out <- betaVec
      }
    else ## Convert betaVec -> Mdl.beta
      {
        ## Convert the parameter vector into the model structure
        Idx0 <- NA # the initial parameter index
        Idx1 <- 0
        for(CompCaller in names(Mdl.betaIdx))
          {
            for(parCaller in names(Mdl.betaIdx[[CompCaller]]))
              {
                if(parUpdate[[CompCaller]][[parCaller]] == TRUE)
                  {
                    betaIdxCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]]
                    betaLen <- length(betaIdxCurr)
                    betaLenNoZero <-sum(betaIdxCurr)

                    Idx0 <- Idx1+1
                    Idx1 <- Idx0 + betaLenNoZero -1
                    IdxCurr <-Idx0:Idx1

                    betaCurr <- array(0, dim(betaIdxCurr))
                    betaCurr[betaIdxCurr] <- betaInput[IdxCurr]

                    if(!is.matrix(betaCurr)) browser()
                    Mdl.beta[[CompCaller]][[parCaller]] <- betaCurr
                  }
              }

          }
        out <- Mdl.beta
      }
    return(out)
  }
