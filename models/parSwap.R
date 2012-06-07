parSwap <- function(betaInput, Mdl.beta = NA, Mdl.betaIdx = NA, parUpdate = NA)
  {
    ## Convert Mdl.beta -> betaVec
    if(class(betaInput)  == "list")
      {
        for(CompCurr in names(Mdl.betaIdx))
          {
            for(parCurr in names(Mdl.betaIdx[[CompCurr]]))
              {
                if(parUpdate[[CompCurr]][[parCurr]] == FALSE)
                  {
                    Mdl.beta[[CompCurr]][[parCurr]] <- NULL
                    Mdl.betaIdx[[CompCurr]][[parCurr]] <- NULL
                  }
              }
          }
        betaVec <- unlist(Mdl.beta)
        betaIdxVec <- unlist(Mdl.betaIdx)
        betaVec <- betaVec[betaIdxVec]
        names(betaVec) <- NULL
        out <- betaVec
      }
    else ## Convert betaVec -> Mdl.beta
      {
        ## Convert the parameter vector into the model structure
        Idx0 <- NA # the initial parameter index
        Idx1 <- 0
        for(CompCurr in names(Mdl.betaIdx))
          {
            for(parCurr in names(Mdl.betaIdx[[CompCurr]]))
              {
                if(parUpdate[[CompCurr]][[parCurr]] == TRUE)
                  {
                    betaIdxCurr <- Mdl.betaIdx[[CompCurr]][[parCurr]]
                    betaLen <- length(betaIdxCurr)
                    betaLenNoZero <-sum(betaIdxCurr)

                    Idx0 <- Idx1+1
                    Idx1 <- Idx0 + betaLenNoZero -1
                    IdxCurr <-Idx0:Idx1

                    betaCurr <- matrix(0, betaLen, 1)
                    betaCurr[betaIdxCurr] <- betaInput[IdxCurr]
                    Mdl.beta[[CompCurr]][[parCurr]] <- betaCurr
                  }
              }

          }
        out <- Mdl.beta
      }
    return(out)
  }
