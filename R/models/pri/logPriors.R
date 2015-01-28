##' This is the prior settings for the Covariates dependent copula model
##'
##' The detailed documentation is in the main setting file for each parameter.
##' @title Copula model prior settings.
##' @param Mdl.X
##' @param Mdl.parLink
##' @param Mdl.beta
##' @param Mdl.betaIdx
##' @param varSelArgs
##' @param priArgs "list"
##' @param parUpdate
##' @param priCurr "list"
##' @return "list" synced
##' @references
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
    ## Loop over all updated parameter candidates

    CompNM <- names(priArgs)
    for(iComp in CompNM)
        {
###----------------------------------------------------------------------------
### Only update priors for parameters that need to update.
###----------------------------------------------------------------------------
            parUpdateIdx <- which(parUpdate[[iComp]] == TRUE)

            for(iPar in parUpdateIdx)
                {
                    ## Initial the storage structure for current log prior
                    outCurr <-  priArgs[[iComp]][[iPar]]

###----------------------------------------------------------------------------
### Prior for variable selection indicators
###----------------------------------------------------------------------------
                    priArgsCurr <- priArgs[[iComp]][[iPar]][["indicators"]]
                    betaIdxCurr <- Mdl.betaIdx[[iComp]][[iPar]]

                    if(tolower(priArgsCurr[["type"]]) == "bern") # Bernoulli prior
                        {
                            ## The probability is set for each variable involved in the
                            ## variable selection procedure via "valSelArgs"
                            ## The probability is recycled if necessary
                            prob <- priArgsCurr[["prob"]]

                            ## Variable section candidates
                            candIdx <- varSelArgs[[iComp]][[iPar]][["cand"]]
                            candLen <- length(candIdx)
                            probVec <- matrix(prob, candLen, 1)

                            ## TRUE or FALSE of variable selection candidates
                            varSelCandTF <- betaIdxCurr[candIdx]

                            logDens <- sum(log(probVec[varSelCandTF])) +
                                sum(log(1-probVec[!varSelCandTF]))

                            outCurr[["indicators"]] <- logDens
                        }

###----------------------------------------------------------------------------
### Prior for coefficients
###----------------------------------------------------------------------------

### intercept as special case
### The intercept should alway be included.
                    priArgsCurr <- priArgs[[iComp]][[iPar]][["beta"]][["intercept"]]
                    betaCurr <- Mdl.beta[[iComp]][[iPar]][1] # the intercept
                    linkCurr <- Mdl.parLink[[iComp]][[iPar]]

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

### coefficients conditional on variable selection indicators
                    priArgsCurr <- priArgs[[iComp]][[iPar]][["beta"]][["slopes"]]
                    betaCurr <- Mdl.beta[[iComp]][[iPar]][-1] # Slopes(taking away intercept)
                    betaIdxNoIntCurr <- Mdl.betaIdx[[iComp]][[iPar]][-1] # Variable section
                                        # indicator without intercept

                    ## The covariance matrix for the whole beta vector
                    X <- Mdl.X[[iComp]][[iPar]][, -1, drop = FALSE]

                    if(length(X) == 0L)
                        {
                            ## No covariates at all (only intercept in the model)
                            outCurr[["beta"]][["slopes"]] <- NULL
                        }
                    else
                        {

                            if(tolower(priArgsCurr[["type"]]) == "cond-mvnorm")
                                {
                                    ## Normal distribution condition normal The full beta
                                    ## vector is assumed as normal. Since variable selection is
                                    ## included in the MCMC, The final proposed beta are those
                                    ## non zeros. We need to using the conditional normal
                                    ## density See Mardia p. 63

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
                                            coVar <- qr.solve(crossprod(X))
                                        }
                                    else if(tolower(covariance) == "identity")
                                        {
                                            coVar <- diag1(betaLen)
                                        }

                                    ## Calculate the log density
                                    if(Idx0Len == 0L )
                                        {
                                            ## 1. all are selected or
                                            ## Switch to unconditional prior.

                                            outCurr[["beta"]][["slopes"]] <- dmvnorm(
                                                matrix(betaCurr, 1, betaLen,byrow = TRUE),
                                                meanVec, coVar*shrinkage, log = TRUE)
                                        }
                                    else if( Idx0Len == betaLen)
                                        {
                                            ## 2. non are selected:
                                            ## Switch to unconditional prior.
                                            outCurr[["beta"]][["slopes"]] <- dmvnorm(
                                                matrix(betaCurr, 1, betaLen, byrow = TRUE),
                                                meanVec, coVar*shrinkage, log = TRUE)
                                        }
                                    else if(Idx0Len > 0 & Idx0Len < betaLen)
                                        {
                                            ## 3. some are selected
                                            ## Using the conditional prior
                                            A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
                                            condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
                                            condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]

                                            outCurr[["beta"]][["slopes"]] <- dmvnorm(
                                                matrix(betaCurr[Idx1], 1),
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
                    Mdl.logPri[[iComp]][[iPar]] <- outCurr
                }
        }
    return(Mdl.logPri)
}
