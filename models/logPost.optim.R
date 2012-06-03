logPost.optim <- function(betaVec, CplNM, Mdl.Y, Mdl.X, Mdl.betaIdx, Mdl.parLink,
                          varSelArgs, MargisTypes, priArgs, parUpdate,
                          staticArgs)
  {
    ## a wrapper of the log posterior function that can be used for directly
    ## optimization via Newton's method

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
            betaCurr[betaIdxCurr] <- betaVec[IdxCurr]
            Mdl.beta[[CompCurr]][[parCurr]] <- betaCurr
          }
      }
    ## update staticArgs

    ## Update logPost
    out <- logPost(CplNM = CplNM,
                   Mdl.Y = Mdl.Y,
                   Mdl.X = Mdl.X,
                   Mdl.beta = Mdl.beta,
                   Mdl.betaIdx = Mdl.betaIdx,
                   Mdl.parLink = Mdl.parLink,
                   varSelArgs = varSelArgs,
                   MargisTypes = MargisTypes,
                   priArgs = priArgs,
                   parUpdate = parUpdate,
                   staticArgs = staticArgs)
    return(out)
  }
