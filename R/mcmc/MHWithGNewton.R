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
  chainCaller <- parCplCaller(parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  ## Set a reject flag to handle unexpected situations. If TRUE, the proposal
  ## is rejected anyway regardless of other situations.
  rejectFlag <- FALSE
  accept.prob <- NA

###----------------------------------------------------------------------------
### VARIABLE SELECTION PROPOSAL
###----------------------------------------------------------------------------
  ## Randomly propose a subset for covariates to change

  varSelCand <- varSelArgs[[CompCaller]][[parCaller]]$cand
  betaIdxArgs <- propArgs[[CompCaller]][[parCaller]][["indicators"]]

  ## The current and proposed variable selection indicators are the same by
  ## assuming there is no variable selection
  betaIdx.curr <- Mdl.betaIdx[[CompCaller]][[parCaller]]
  betaIdx.prop <- betaIdx.curr

  ## The jump density for the variable selection indicators
  ## TODO: Add adaptive scheme
  logJump.Idx.currATprop <- 1
  logJump.Idx.propATcurr <- 1

  ## No. of covariates
  nCovs <- length(betaIdx.curr)

  ## If variable selection is available, make a change proposal
  ## Otherwise do nothing.
  if(length(varSelCand) > 0)
    {
      if(betaIdxArgs$type == "binom")
        {
          ## Binomial proposal a small subset
          betaIdx.propCandIdx <- which(rbinom(n = length(varSelCand) , size = 1,
                                              prob = betaIdxArgs$prob) == 1)

          ## Special case to make sure at least one variable is proposed a
          ## change. FIXME: This may break the binomial settings
          if(length(betaIdx.propCandIdx) == 0)
            {
              betaIdx.propCand <- sample(varSelCand, 1)
            }
          else
            {
              betaIdx.propCand <- varSelCand[betaIdx.propCandIdx]
            }
        }
      ## Propose a change
      betaIdx.prop[betaIdx.propCand] <- !betaIdx.curr[betaIdx.propCand]
    }

###----------------------------------------------------------------------------
### PROPOSAL A DRAW VIA K-STEP NEWTON'S METHOD
###----------------------------------------------------------------------------
  beta.curr.full <- Mdl.beta[[CompCaller]][[parCaller]]
  beta.curr <- beta.curr.full[betaIdx.curr]

  Mdl.beta.curr <- Mdl.beta
  Mdl.betaIdx.curr <- Mdl.betaIdx
  staticArgs.curr <- staticArgs

  ## Newton method to approach the posterior based on the current draw
  beta.NTProp <- GNewtonMove(
      propArgs = propArgs,
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
      staticArgs = staticArgs.curr)

  ## The information for proposed density via K-step Newton's method
  beta.propArgs <- propArgs[[CompCaller]][[parCaller]][["beta"]]
  beta.prop.type <- beta.propArgs[["type"]]
  beta.prop.mean <- matrix(beta.NTProp[["param"]], 1) # 1-by-p
  beta.prop.sigma <- -beta.NTProp[["HessObsInv"]]
  staticArgs.prop <- beta.NTProp[["staticArgs"]]

  ## Check if it is a good proposal
  if(any(is.na(beta.prop.sigma)) ||
     is(try(chol(beta.prop.sigma), silent=TRUE), "try-error"))
    {## Something is wrong.
      rejectFlag <- TRUE
    }
  else # Continues with the Metropolis-Hastings
    { ## Propose a draw from multivariate t-distribution based on the proposed
      ## An idea (out of loud) : Generate a matrix of param. Select the one
      ## that give max acceptance probability (but need correct the acceptance
      ## probability.)
      if(tolower(beta.prop.type) == "mvt")
        {
          ## The proposal parameters block
          beta.prop.df <- beta.propArgs[["df"]]
          beta.prop <- beta.prop.mean + rmvt(
              sigma = beta.prop.sigma,
              n = 1, df = beta.prop.df)
        }
    }

###----------------------------------------------------------------------------
### REVERSE OF THE NEWTON METHOD
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

      beta.NTPropRev <- GNewtonMove(
          propArgs = propArgs,
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
          staticArgs = staticArgs.prop)

      ## The information for proposed density via K-step Newton's method
      beta.propRev.mean <- matrix(beta.NTPropRev[["param"]], 1) # 1-by-p
      beta.propRev.sigma <- -beta.NTPropRev[["HessObsInv"]] # p-by-p

###----------------------------------------------------------------------------
### COMPUTING THE METROPOLIS-HASTINGS RATIO
###----------------------------------------------------------------------------
      if(all(!is.na(beta.propRev.sigma)))
        {
          ## The jump density for proposed point at proposed mode
          logJump.propATprop <- dmvt(
              x = beta.prop - beta.prop.mean,
              sigma = beta.prop.sigma,
              df = beta.prop.df, log = TRUE)

          ## The jump density for curr draw at reverse proposed mode.
          logJump.currATpropRev<- dmvt(
              x = beta.curr - beta.propRev.mean,
              sigma = beta.propRev.sigma,
              df = beta.prop.df, log = TRUE)

          ## The log posterior for the proposed draw
          logPost.propOut <- logPost(
              CplNM = CplNM,
              Mdl.Y = Mdl.Y,
              Mdl.X = Mdl.X,
              Mdl.beta = Mdl.beta.prop,
              Mdl.betaIdx = Mdl.betaIdx.prop,
              Mdl.parLink = Mdl.parLink,
              varSelArgs = varSelArgs,
              MargisTypes = MargisTypes,
              priArgs = priArgs,
              parUpdate = parUpdate,
              staticArgs = staticArgs.prop)

          logPost.prop <- logPost.propOut[["Mdl.logPost"]]
          staticArgs.prop <- logPost.propOut[["staticArgs"]]

          ## The log posterior for the current draw
          logPost.curr <- logPost(
              CplNM = CplNM,
              Mdl.Y = Mdl.Y,
              Mdl.X = Mdl.X,
              Mdl.beta = Mdl.beta.curr,
              Mdl.betaIdx = Mdl.betaIdx.curr,
              Mdl.parLink = Mdl.parLink,
              varSelArgs = varSelArgs,
              MargisTypes = MargisTypes,
              priArgs = priArgs,
              parUpdate = parUpdate,
              staticArgs = staticArgs.curr)[["Mdl.logPost"]]

          ## compute the MH ratio.
          MHRatio <- exp(logPost.prop - logPost.curr +
                         logJump.currATpropRev - logJump.propATprop +
                         logJump.Idx.currATprop - logJump.Idx.propATcurr)

          ## cat(logPost.prop, logPost.curr, "\n")

          ## the acceptance probability
          accept.prob <- min(1, MHRatio)
        }
    }
###----------------------------------------------------------------------------
### THE MH ACCEPTANCE PROBABILITY AND ACCEPT/REJECT THE PROPOSED DRAW.
###----------------------------------------------------------------------------

  if(is.na(accept.prob)) browser()

  if(rejectFlag == FALSE &&
     !is.na(accept.prob) &&
     runif(1) < accept.prob) # keep update
    {
      out <- list(betaIdx = betaIdx.prop,
                  beta = beta.prop.full,
                  accept.prob = accept.prob,
                  staticArgs = staticArgs.prop)
    }
  else # keep current
    {
      out <- list(betaIdx = betaIdx.curr,
                  beta = beta.curr.full,
                  accept.prob = accept.prob,
                  staticArgs = staticArgs.curr)
    }
  return(out)
}
