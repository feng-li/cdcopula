##' Metropolis–Hastings algorithm with K-step Newton method with variable selection.
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
##' @param staticCache
##' @return "list"
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Thu Feb 17 14:03:14 CET 2011; Current: Fri Mar 27 11:29:18 CST 2015.
MetropolisHastings <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                              Mdl.parLink, parUpdate, priArgs, varSelArgs,
                              propArgs, MargisTypes, staticCache, MCMCUpdateStrategy)
{
  ## The updating component parameter chain
  chainCaller <- parCplRepCaller(CplNM = CplNM, parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  ## The proposal methods
  algmArgs <- propArgs[[CompCaller]][[parCaller]][["algorithm"]]
  beta.propArgs <- propArgs[[CompCaller]][[parCaller]][["beta"]]
  betaIdx.propArgs <- propArgs[[CompCaller]][[parCaller]][["indicators"]]

###----------------------------------------------------------------------------
### INITIAL COPY OF CURRENT VALUES
###----------------------------------------------------------------------------

  Mdl.beta.curr <- Mdl.beta
  Mdl.betaIdx.curr <- Mdl.betaIdx
  staticCache.curr <- staticCache

  beta.curr.full <- Mdl.beta[[CompCaller]][[parCaller]] # p-by-lq
  betaIdx.curr <- Mdl.betaIdx[[CompCaller]][[parCaller]] # p-by-lq
  beta.curr <- beta.curr.full[betaIdx.curr] # p*lq

  ## Assume initial no variable selection
  betaIdx.prop <- betaIdx.curr

###----------------------------------------------------------------------------
### VARIABLE SELECTION PROPOSAL
###----------------------------------------------------------------------------

  ## No. of covariates
  nCovs <- dim(betaIdx.curr)[1]
  nPar <- dim(betaIdx.curr)[2]

  ## Randomly propose a subset for covariates to change
  beta01Mat <- matrix(1:(nCovs*nPar), nCovs, nPar)
  varSelCandRow <- varSelArgs[[CompCaller]][[parCaller]][["cand"]] # sub.q-by-1
  varSelCand <- beta01Mat[varSelCandRow, ] # sub.p-by-lq

  ## If variable selection is available, make a change proposal. Otherwise no variable
  ## selection.
  if(length(varSelCand) > 0)
    {
      if(betaIdx.propArgs[["type"]] == "binom")
        {
          ## Binomial proposal a small subset
          betaIdx.propCandIdx <- which(rbinom(
                  n = length(varSelCand),
                  size = 1L,
                  prob = betaIdx.propArgs[["prob"]]) == 1L)

          if(length(betaIdx.propCandIdx) != 0)
            {
              ## Situation when propose a change with prob = betaIdxArgs$prob
              betaIdx.propCand <- varSelCand[betaIdx.propCandIdx]
              betaIdx.prop[betaIdx.propCand] <- !betaIdx.curr[betaIdx.propCand]

              VSProp <- "VSProp"
            }
          else
            {
              ## Situation when propose no changes  =  no variable selection
              VSProp <- NULL
            }
        }
      else
        {
          stop("No such variable selection proposals!")
        }
    }
  else
    {   ## No variable selection
      VSProp <- NULL
    }

  ## The jump density for the variable selection indicators. TODO: Add adaptive scheme

  logJump.Idx.currATprop <- 1
  logJump.Idx.propATcurr <- 1

###----------------------------------------------------------------------------
### INITIALIZING THE MH SCHEMES.
###----------------------------------------------------------------------------

  ## Initially we assume there is only mh step.
  MHUpdate <- "MH"

  ## If variable selection (VS) is proposed, add an extra MH step with VS at the
  ## beginning. otherwise we do the usual update without variable selection,
  ## i.e. betaidx.prop = betaidx.curr.
  MHUpdate <- c(VSProp, MHUpdate)

  ## For each mh step, we obtain an acceptance probability
  nMH <- length(MHUpdate)
  accept.probs <- rep(NA, nMH)

  ## Set a reject flag to handle unexpected situations. If TRUE, the proposal is rejected
  ## anyway regardless of other situations.
  errorFlags <- rep(FALSE, nMH)

###----------------------------------------------------------------------------
### METROPOLIS-HASTINGS ALGORITHM
###----------------------------------------------------------------------------
  for(iMH in 1:nMH)
    {
      if(tolower(algmArgs[["type"]]) == "gnewtonmove")
        { ## Newton method to approach the posterior based on the current draw
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
                  staticCache = staticCache.curr,
                  MCMCUpdateStrategy = MCMCUpdateStrategy)
        }
      else if(tolower(algmArgs[["type"]])  == "randomwalk")
        { ## Random walk metropolis (with/without variable selection)
          stop("Not implement yet!")

          ## beta.NTProp <- list(errorFlag = FALSE,
          ##                     param = ??,
          ##                     HessObsInv = ??)
        }
      else
        {
          stop("Not implement yet!")
        }

      ## Check if it is a good proposal
      if(beta.NTProp$errorFlag)
        { ## Something is wrong. Skip this MH update(if is VS, jump to pure update without
          ## variable selection) or terminate the algorithm.

          errorFlags[iMH] <- TRUE
          betaIdx.prop <- betaIdx.curr
          next
        }

      ## else Continues with the Metropolis-Hastings Propose a draw from multivariate
      ## t-distribution based on the proposed An idea (out of loud) : Generate a matrix of
      ## param. Select the one that give max acceptance probability (but need correct the
      ## acceptance probability.)

      ## The information for proposed density via K-step Newton's method
      beta.prop.mean <- matrix(beta.NTProp[["param"]], 1) # 1-by-p
      beta.prop.sigma <- -beta.NTProp[["HessObsInv"]]


      staticCache.prop <- beta.NTProp[["staticCache"]]
      if(tolower(beta.propArgs[["type"]]) == "mvt")
        {
          ## require("mvtnorm")
          ## The proposal parameters block
          beta.prop <- beta.prop.mean + rmvt(
                  sigma = (beta.prop.sigma+t(beta.prop.sigma))/2,
                  n = 1, df = beta.propArgs[["df"]])
        }

###----------------------------------------------------------------------------
### REVERSE PROBABILITY OF THE PROPOSALS
###----------------------------------------------------------------------------

      ## Newton method to approach the posterior for the proposed draw
      beta.prop.full <- matrix(0, nCovs, nPar)
      beta.prop.full[betaIdx.prop] <- beta.prop
      Mdl.beta.prop <- Mdl.beta.curr
      Mdl.beta.prop[[CompCaller]][[parCaller]] <- beta.prop.full

      Mdl.betaIdx.prop <- Mdl.betaIdx.curr
      Mdl.betaIdx.prop[[CompCaller]][[parCaller]] <- betaIdx.prop

      if(tolower(algmArgs[["type"]]) == "gnewtonmove")
        {
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
                  staticCache = staticCache,
                  MCMCUpdateStrategy = MCMCUpdateStrategy)
        }
      else if(tolower(algmArgs[["type"]])  == "randomwalk")
        {
          stop("Not implement yet!")

          ## beta.NTPropRev <- list(errorFlag = FALSE,
          ##                        param = ??,
          ##                        HessObsInv = ??)
        }
      else
        {
          stop("Not implement yet!")
        }

      if(beta.NTPropRev$errorFlag)
        { ## Something is wrong.
          errorFlags[iMH] <- TRUE
          betaIdx.prop <- betaIdx.curr
          next
        }

      ## The information for proposed density via K-step Newton's method
      beta.propRev.mean <- matrix(beta.NTPropRev[["param"]], 1) # 1-by-p
      beta.propRev.sigma <- -beta.NTPropRev[["HessObsInv"]] # p-by-p

      ## browser()
###----------------------------------------------------------------------------
### COMPUTING THE METROPOLIS-HASTINGS RATIO
###----------------------------------------------------------------------------

      ## The jump density for proposed point at proposed mode and the jump density for
      ## current draw at reverse proposed mode.
      if(tolower(beta.propArgs[["type"]]) == "mvt")
        {
          logJump.propATprop <- dmvt(
                  x = beta.prop - beta.prop.mean,
                  sigma = (beta.prop.sigma+t(beta.prop.sigma))/2,
                  df = beta.propArgs[["df"]], log = TRUE)

          logJump.currATpropRev<- dmvt(
                  x = beta.curr - beta.propRev.mean,
                  sigma = (beta.propRev.sigma+t(beta.propRev.sigma))/2,
                  df = beta.propArgs[["df"]], log = TRUE)
        }

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
              staticCache = staticCache,
              MCMCUpdateStrategy = MCMCUpdateStrategy)

      logPost.prop <- logPost.propOut[["Mdl.logPost"]]
      staticCache.prop <- logPost.propOut[["staticCache"]]

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
              staticCache = staticCache,
              MCMCUpdateStrategy = MCMCUpdateStrategy)[["Mdl.logPost"]]

      ## Compute the (log) MH ratio.
      logMHRatio <- logPost.prop - logPost.curr +
                     logJump.currATpropRev - logJump.propATprop +
                     logJump.Idx.currATprop - logJump.Idx.propATcurr

      ## The acceptance probability
      accept.prob.curr <- exp(min(0, logMHRatio))

###----------------------------------------------------------------------------
### THE MH ACCEPTANCE PROBABILITY AND ACCEPT/REJECT THE PROPOSED DRAW.
###----------------------------------------------------------------------------

      if(!is.na(accept.prob.curr) &&
         runif(1) < accept.prob.curr)
        {## keep the proposal
          betaIdx.curr <- betaIdx.prop
          beta.curr <- beta.prop

          Mdl.beta.curr <- Mdl.beta.prop
          Mdl.betaIdx.curr <- Mdl.betaIdx.prop

          staticCache.curr <- staticCache.prop
          ## browser()
        }
      else
        {## keep the current
          betaIdx.prop <- betaIdx.curr
        }
      accept.probs[iMH] <- accept.prob.curr
    }

###----------------------------------------------------------------------------
### THE FINAL OUTPUT
###----------------------------------------------------------------------------
  ## The acceptance prob are from the last MH update or keep current draw.
  if(errorFlags[nMH] == TRUE)
    { ## epic failure
      out <- list(errorFlag = TRUE)
    }
  else
    {
      out <- list(beta =  Mdl.beta.curr[[CompCaller]][[parCaller]],
                  betaIdx = Mdl.betaIdx.curr[[CompCaller]][[parCaller]],
                  accept.prob = accept.probs[nMH],
                  staticCache = staticCache.curr,
                  errorFlag = FALSE)
    }
  return(out)
}
