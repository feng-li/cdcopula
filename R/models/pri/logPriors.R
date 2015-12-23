##' This is the prior settings for the Covariates dependent copula model
##'
##' The detailed documentation is in the main setting file for each parameter.
##' @title Copula model prior settings.
##' @param Mdl.X "list"
##' @param Mdl.parLink "list"
##' @param Mdl.beta "list"
##' @param Mdl.betaIdx "list"
##' @param varSelArgs "list"
##' @param priArgs "list"
##' @param parUpdate "list"
##' @param priCurr "list"
##' @return "list" synced
##' @references "list"
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Dec 15 10:45:56 CET 2011;
logPriors <- function(Mdl.X, Mdl.parLink, Mdl.beta, Mdl.betaIdx,
                      varSelArgs, priArgs, parUpdate, Mdl.logPri)
{
  ## Allocate the structure.
    if(missing(Mdl.logPri))
    {
        Mdl.logPri <- rapply(Mdl.beta, function(x) NA, how = "replace")
    }

    if(missing(parUpdate))
    {
        parUpdate <- rapply(Mdl.beta, function(x) TRUE, how = "replace")
    }


    ## Loop over all updated parameter candidates

  CompNM <- names(priArgs)
  for(CompCaller in CompNM)
  {
###----------------------------------------------------------------------------
### Only update priors for parameters that need to update.
###----------------------------------------------------------------------------
    parUpdateIdx <- names(parUpdate[[CompCaller]])[unlist(parUpdate[[CompCaller]])]
    for(parCaller in parUpdateIdx)
    {   ## Initial the storage structure for current log prior
      outCurr <-  priArgs[[CompCaller]][[parCaller]]

###----------------------------------------------------------------------------
### Prior for variable selection indicators
###----------------------------------------------------------------------------
      priArgsCurr <- priArgs[[CompCaller]][[parCaller]][["indicators"]]
      betaIdxCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]] # p-yb-lq
      nPar <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]

      if(tolower(priArgsCurr[["type"]]) == "bern") # Bernoulli prior
      {
        ## The probability is set for each variable involved in the
        ## variable selection procedure via "valSelArgs"
        ## The probability is recycled if necessary
        prob <- priArgsCurr[["prob"]]

        ## Variable section candidates.
        if(dim(betaIdxCurr)[1] == 1)
        {
          ## Only intercept is in, do not allow for variable
          ## selection.
          candIdx <- NULL
        }
        else
        {
          varSelCandConfigCurr <- varSelArgs[[CompCaller]][[parCaller]][["cand"]]

          if(class(varSelCandConfigCurr) == "character" &&
             tolower(varSelCandConfigCurr) == "2:end")
          {
            ncolX.ij <- nrow(Mdl.beta[[CompCaller]][[parCaller]])
            candIdx <- 2:ncolX.ij
          }
          else
          {
            candIdx <- varSelCandConfigCurr
          }

        }


        if(length(candIdx)>0)
        {
          probMat <- matrix(prob, length(candIdx), nPar)
          ## TRUE or FALSE of variable selection candidates

          varSelCandTF <- betaIdxCurr[candIdx, , drop = FALSE]

          logDens <- sum(dbinom(x = varSelCandTF, size = 1,
                                prob = probMat, log = TRUE))
        }
        else
        {
          ## No variable selection is used, thus we don't need
          ## prior.

          logDens <- NULL
        }
        outCurr[["indicators"]] <- logDens
      }

###----------------------------------------------------------------------------
### Prior for coefficients
###----------------------------------------------------------------------------

### intercept as special case. The intercept should alway be included.
      priArgsCurr <- priArgs[[CompCaller]][[parCaller]][["beta"]][["intercept"]]

      betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][1,,drop = FALSE]#intercepts

      ## if(is(betaCurr, "try-error")) browser()

      linkCurr <- Mdl.parLink[[CompCaller]][[parCaller]]

      if(tolower(priArgsCurr[["type"]]) == "custom")
      {
        ## Call the any2any() function
        densOutput <- any2any(densArgs = priArgsCurr, linkArgs = linkCurr)

        mean <- densOutput$mean
        variance <- densOutput$variance
        shrinkage <- priArgsCurr[["output"]][["shrinkage"]]

        logDens <- dnorm(x = betaCurr, mean = mean,
                         sd = sqrt(variance*shrinkage), log = TRUE)
        outCurr[["beta"]][["intercept"]] <- logDens
      }

### coefficients (conditional on variable selection indicators)
      priArgsCurr <- priArgs[[CompCaller]][[parCaller]][["beta"]][["slopes"]]

      ## Slopes and variable selection indicators(taking away intercept)
      betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][-1, , drop = FALSE]
      betaIdxNoIntCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]][-1,,drop = FALSE]


      if(length(betaIdxNoIntCurr) == 0L)
      {
        ## No covariates at all (only intercept in the model)
        outCurr[["beta"]][["slopes"]] <- NULL
      }
      else
      {
        if(tolower(priArgsCurr[["type"]]) == "cond-mvnorm")
        {
          ## Normal distribution condition normal The full beta
          ## vector is assumed as normal. Since variable
          ## selection is included in the MCMC, The final
          ## proposed beta are those non zeros. We need to using
          ## the conditional normal density See Mardia p. 63

          ## Subtract the prior information for the full beta
          mean <- priArgsCurr[["mean"]] # mean of density
          covariance <- priArgsCurr[["covariance"]] # Covariates
          shrinkage <- priArgsCurr[["shrinkage"]] # Shrinkage

          ## Split the beta vector by selected and non-selected.
          Idx1 <- which(betaIdxNoIntCurr == TRUE)
          Idx0 <- which(betaIdxNoIntCurr == FALSE)

          Idx0Len <- length(Idx0)
          Idx1Len <- length(Idx1)

          betaLen <- length(betaIdxNoIntCurr)

          ## The mean vector (recycled if necessary)
          meanVec <- matrix(mean, 1, betaLen)

          ## The covariance matrix for the whole beta vector
          if(tolower(covariance) == "g-prior")
          {
            ## The covariance matrix for the whole beta vector
            X <- Mdl.X[[CompCaller]][[parCaller]][, -1, drop = FALSE]
            coVar0 <- qr.solve(crossprod(X))
            coVar0Lst <- lapply(as.vector(rep(NA, nPar),"list"),
                                function(x) coVar0)

            coVar <- block.diag(coVar0Lst)

          }
          else if(tolower(covariance) == "identity")
          {
            coVar <- diag1(betaLen)
          }

          ## Calculate the log density
          if(Idx0Len == 0L)
          {   ## 1. all are selected or
            ## Switch to unconditional prior.
            outCurr[["beta"]][["slopes"]] <- dmvnorm(matrix(betaCurr, 1, betaLen,byrow = TRUE),
                                                     meanVec, coVar*shrinkage, log = TRUE)
          }
          else if( Idx0Len == betaLen)
          {   ## 2. non are selected:
            ## Switch to unconditional prior.
            outCurr[["beta"]][["slopes"]] <- dmvnorm(matrix(betaCurr, 1, betaLen, byrow = TRUE),
                                                     meanVec, coVar*shrinkage, log = TRUE)
          }
          else if(Idx0Len > 0 & Idx0Len < betaLen)
          {
            ## 3. some are selected
            ## Using the conditional prior
            A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
            condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
            condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]

            outCurr[["beta"]][["slopes"]] <- dmvnorm(matrix(betaCurr[Idx1], 1),
                                                     condMean, condCovar*shrinkage, log = TRUE)
          }
          else
          {
            stop("Debug me: Unknown situation for conditional priors.")
          }
        }
      }

###----------------------------------------------------------------------------
### Update the output for prior
###----------------------------------------------------------------------------
      Mdl.logPri[[CompCaller]][[parCaller]] <- outCurr
    }
  }
  return(Mdl.logPri)
}
