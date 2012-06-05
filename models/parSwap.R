parSwap <- function(betaInput, Mdl.betaIdx = NA, parUpdate = NA)
  {
    ## Convert betaVec to Mdl.beta
    if(is.vector(betaInput))
      {
        ## Convert the parameter vector into the model structure
        Mdl.beta <- Mdl.betaIdx # reserve the structure
        Idx0 <- NA # the initial parameter index
        Idx1 <- 0
        for(CompCurr in names(Mdl.betaIdx))
          {
            for(parCurr in names(Mdl.betaIdx[[CompCurr]]))
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

        out <- Mdl.beta
      }

    if(is.list(betaInput))
      {
        betaVec <- unlist(Mdl.beta)
        betaVec <- betaVec[betaVec != 0]
        names(betaVec) <- NULL


        out <- betaVec
      }
    return(out)
  }
