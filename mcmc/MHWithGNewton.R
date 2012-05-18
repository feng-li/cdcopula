##' Metropolis–Hastings algorithm with K-step Newton method with variable
##' selection.
##'
##' Metropolis-Hastings with generalized Newton method
##' @param CplNM
##' @param Mdl.Y
##' @param Mdl.X
##' @param Mdl.beta
##' @param Mdl.betaIdx
##' @param Mdl.parLink
##' @param parUpdate
##' @param priArgs
##' @param varSelArgs
##' @param propArgs
##' @param MargisTypes
##' @param staticArgs
##' @return "list"
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Thu Feb 17 14:03:14 CET 2011;
##'       Current: Fri May 18 11:05:06 CEST 2012.
##' DEPENDS: mvtnorm
MHWithGNewton <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                          Mdl.parLink, parUpdate, priArgs, varSelArgs,
                          propArgs, MargisTypes, staticArgs)
{

  ## The updating component parameter chain
  chainCaller <- parCaller(parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  ## Set a reject flag to handle unexpected situations. If TRUE, the proposal
  ## is rejected anyway regardless of other situations.
  rejectFlag <- FALSE

###----------------------------------------------------------------------------
### VARIABLE SELECTION PROPOSAL
###----------------------------------------------------------------------------
  ## Randomly propose a subset for covariates to change
  varSelCand <- varSelArgs[[CompCaller]][[parCaller]]$cand
  betaIdxArgs <- propArgs[[CompCaller]][[parCaller]][["indicators"]]
  if(betaIdxArgs$type == "binom")
    {
      ## Binomial proposal a small subset
      betaIdx.propCand <- which(rbinom(n = length(varSelCand) , size = 1,
                                       prob = betaIdxArgs$prob) == 1)
      ## Special case to make sure at least one variable is proposed a change
      if(length(betaIdx.propCand) == 0)
        {
          betaIdx.propCand <- sample(varSelCand, 1)
        }
    }

  ## The current and proposed variable selection indicators
  betaIdx.curr <- Mdl.betaIdx[[CompCaller]][[parCaller]]

  betaIdx.prop <- betaIdx.curr
  betaIdx.prop[betaIdx.propCand] <- !betaIdx.curr[betaIdx.propCand]

  ## No. of covariates
  nCovs <- length(betaIdx.curr)

  ## The jump density for the variable selection indicators
  ## TODO: Add adaptive scheme
  logJump.betaIdx.currATprop <- 1
  logJump.betaIdx.propATcurr <- 1

###----------------------------------------------------------------------------
### Make good proposal via K-steps Newton's method
###----------------------------------------------------------------------------
  beta.curr <- Mdl.betaIdx[[CompCaller]][[parCaller]][betaIdx.curr]
  Mdl.beta.curr <- Mdl.beta
  Mdl.betaIdx.curr <- Mdl.betaIdx

  ## Newton method to approach the posterior based on the current draw
  beta.curr2mode <- GNewtonMove(propArgs = propArgs,
                                varSelArgs = varSelArgs,
                                priArgs = priArgs,
                                betaIdxProp = betaIdx.prop,
                                parUpdate = parUpdate,
                                CplNM = CplNM,
                                Mdl.Y = Mdl.Y,
                                Mdl.X = Mdl.X,
                                Mdl.parLink = Mdl.parLink,
                                Mdl.beta = Mdl.beta.curr,
                                Mdl.betaIdx = Mdl.betaIdx.curr,
                                MargisTypes = MargisTypes,
                                staticArgs = staticArgs)

  ## The information for proposed density via K-step Newton's method
  beta.curr2mode.mean <- matrix(beta.curr2mode$param, 1) # 1-by-p
  beta.curr2mode.sigma <- -beta.curr2mode$HessObsInv

  ## Check if it is a good proposal
  if(any(is.na(beta.curr2mode.sigma)) ||
     is(try(chol(beta.curr2mode.sigma), silent=TRUE), "try-error"))
    {## Something is wrong.
      rejectFlag <- TRUE
    }
  else # Continues with the Metropolis-Hastings
    { ## Propose a draw from multivariate t-distribution based on the proposed
      ## An idea (out of loud) : Generate a matrix of param. Select the one
      ## that give max acceptance probability (but need correct the acceptance
      ## probability.)

      beta.propArgs <- propArgs[[CompCaller]][[parCaller]][["beta"]]
      beta.prop.type <- beta.propArgs[["type"]]
      if(tolower(beta.prop.type) == "mvt")
        {
          ## The proposal parameters block
          beta.prop.df <- beta.propArgs[["df"]]
          beta.prop <- beta.curr2mode.mean + rmvt(n = 1,
                                                  sigma = beta.curr2mode.sigma,
                                                  df = beta.prop.df)
        }
    }

###----------------------------------------------------------------------------
### Jump via K-steps Newton's method
###----------------------------------------------------------------------------

  ## Newton method to approach the posterior for the proposed draw
  if(rejectFlag == FALSE)
    {
      beta.prop.full <- matrix(0, nCovs, 1)
      beta.prop.full[betaIdx.prop] <- beta.prop
      Mdl.beta.prop <- Mdl.beta
      Mdl.beta.prop[[CompCaller]][[parCaller]] <- beta.prop.full

      Mdl.betaIdx.prop <- Mdl.betaIdx
      Mdl.betaIdx.prop[[CompCaller]][[parCaller]] <- betaIdx.prop

      beta.prop2mode <- GNewtonMove(propArgs = propArgs,
                                    varSelArgs = varSelArgs,
                                    priArgs = priArgs,
                                    betaIdxProp = betaIdx.curr,
                                    parUpdate = parUpdate,
                                    CplNM = CplNM,
                                    Mdl.Y = Mdl.Y,
                                    Mdl.X = Mdl.X,
                                    Mdl.parLink = Mdl.parLink,
                                    Mdl.beta = Mdl.beta.prop,
                                    Mdl.betaIdx = Mdl.betaIdx.prop,
                                    MargisTypes = MargisTypes,
                                    staticArgs = staticArgs)

      ## The information for proposed density via K-step Newton's method
      beta.prop2mode.mean <- matrix(beta.prop2mode$param, 1) # 1-by-p
      beta.prop2mode.sigma <- - beta.prop2mode$HessObsInv

      if(any(is.na(beta.prop2mode.sigma)))
        { ## Something is wrong at GNewton2,  reject it soon.
          rejectFlag <- TRUE
        }
      else # all are right
        {
          ## The jump density for current point at proposal mode
          logJump.beta.currATprop2mode <- dmvt(x = beta.curr - beta.prop2mode.mean,
                                               sigma = beta.prop2mode.sigma,
                                               df = beta.prop.df, log = TRUE)

          ## The jump density for propose draw to current draw.
          logJump.beta.propATcurr2mode<- dmvt(x = beta.prop - beta.curr2mode.mean,
                                              sigma = beta.curr2mode.sigma,
                                              df = beta.prop.df, log = TRUE)

          ## The log posterior for the proposed draw
          logPost.prop <- logPost(CplNM = CplNM,
                                  Mdl.Y = Mdl.Y,
                                  Mdl.X = Mdl.X,
                                  Mdl.beta = Mdl.beta.prop,
                                  Mdl.betaIdx = Mdl.betaIdx.prop,
                                  Mdl.parLink = Mdl.parLink,
                                  varSelArgs = varSelArgs,
                                  MargisTypes = MargisTypes,
                                  priArgs = priArgs,
                                  parUpdate = parUpdate,
                                  staticArgs = staticArgs)[["Mdl.logPost"]]

          ## The log posterior for the current draw
          logPost.curr <- logPost(CplNM = CplNM,
                                  Mdl.Y = Mdl.Y,
                                  Mdl.X = Mdl.X,
                                  Mdl.beta = Mdl.beta.curr,
                                  Mdl.betaIdx = Mdl.betaIdx.curr,
                                  Mdl.parLink = Mdl.parLink,
                                  varSelArgs = varSelArgs,
                                  MargisTypes = MargisTypes,
                                  priArgs = priArgs,
                                  parUpdate = parUpdate,
                                  staticArgs = staticArgs)[["Mdl.logPost"]]

          ## compute the MH ratio.
          logMHRatio <- logPost.prop - logPost.curr +
            logJump.beta.currATprop2mode - logJump.beta.propATcurr2mode +
              logJump.betaIdx.currATprop - logJump.betaIdx.propATcurr

          MHRatio <- exp(logMHRatio)
        }
    }
###----------------------------------------------------------------------------
### The MH acceptance probability and keep/update the proposed draw.
###----------------------------------------------------------------------------

  ## the acceptance probability
  accept.prob <- min(1, MHRatio)
  if(!is.na(accept.prob) && runif(1) < accept.prob) # keep update
    {
      out <- list(betaIdx = betaIdx.prop,
                  beta = beta.prop,
                  accept.prob = accept.prob)
    }
  else # keep current
    {
      out <- list(betaIdx = betaIdx.curr,
                  beta = beta.curr,
                  accept.prob = 0)
    }
  return(out)
}
