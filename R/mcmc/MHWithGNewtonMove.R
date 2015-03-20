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
##' @param staticCache
##' @return "list"
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Thu Feb 17 14:03:14 CET 2011;
##'       Current: Fri May 18 11:05:06 CEST 2012.
##' DEPENDS: mvtnorm
MHWithGNewtonMove <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                              Mdl.parLink, parUpdate, priArgs, varSelArgs,
                              propArgs, MargisTypes, staticCache, MCMCUpdateStrategy)
{
  ## The updating component parameter chain
  chainCaller <- parCplRepCaller(CplNM = CplNM, parUpdate)

  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  ## Set a reject flag to handle unexpected situations. If TRUE, the proposal
  ## is rejected anyway regardless of other situations.
  MHUpdate <- "MH"

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
  varSelCandRow <- varSelArgs[[CompCaller]][[parCaller]]$cand # sub.q-by-1
  varSelCand <- beta01Mat[varSelCandRow, ] # sub.p-by-lq

  betaIdxArgs <- propArgs[[CompCaller]][[parCaller]][["indicators"]]


  ## If variable selection is available, make a change proposal Otherwise no variable
  ## selection.
  if(length(varSelCand) > 0)
    {
      if(betaIdxArgs$type == "binom")
        {
          ## Binomial proposal a small subset
          betaIdx.propCandIdx <- which(rbinom(
              n = length(varSelCand),
              size = 1L,
              prob = betaIdxArgs$prob) == 1L)

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



  ## The jump density for the variable selection indicators
  ## TODO: Add adaptive scheme
  logJump.Idx.currATprop <- 1
  logJump.Idx.propATcurr <- 1


###----------------------------------------------------------------------------
### LOOP OVER THE MH SCHEMES.
###----------------------------------------------------------------------------
  ## If variable selection is proposed,  add an extra MHNewton at the
  ## beginning,  else the usual update without variable selection,
  ## i.e. betaIdx.prop = betaIdx.curr.

  MHUpdate <- c(VSProp, MHUpdate)
  nMH <- length(MHUpdate)
  errorFlags <- rep(FALSE, nMH)
  accept.probs <- rep(NA, nMH)

  for(iMH in 1:nMH)
    {
###----------------------------------------------------------------------------
### PROPOSAL A DRAW VIA K-STEP NEWTON'S METHOD
###----------------------------------------------------------------------------

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
          staticCache = staticCache.curr,
          MCMCUpdateStrategy = MCMCUpdateStrategy)

      ## Check if it is a good proposal
      if(beta.NTProp$errorFlag ## ||
         ## any(is.na(beta.prop.sigma)) ||
         ## is(try(chol(beta.prop.sigma), silent=TRUE), "try-error")
         )
        {## Something is wrong. Skip this MH update(if is VS, jump to pure
          ## update without variable selection) or terminate the algorithm.
          errorFlags[iMH] <- TRUE
          betaIdx.prop <- betaIdx.curr
          next
        }

      ## else # Continues with the Metropolis-Hastings
      ##   { ## Propose a draw from multivariate t-distribution based on the proposed
      ## An idea (out of loud) : Generate a matrix of param. Select the one
      ## that give max acceptance probability (but need correct the acceptance
      ## probability.)

      ## The information for proposed density via K-step Newton's method
      beta.propArgs <- propArgs[[CompCaller]][[parCaller]][["beta"]]
      beta.prop.type <- beta.propArgs[["type"]]
      beta.prop.mean <- matrix(beta.NTProp[["param"]], 1) # 1-by-p
      beta.prop.sigma <- -beta.NTProp[["HessObsInv"]]
      ## beta.prop.sigma <- diag1(length(beta.prop.mean))*0.1
      ## if(all(beta.prop.sigma<0.01)) cat("beta.prop.sigma too small")


      staticCache.prop <- beta.NTProp[["staticCache"]]


      if(tolower(beta.prop.type) == "mvt")
        {
          ## The proposal parameters block
          beta.prop.df <- beta.propArgs[["df"]]
          beta.prop <- beta.prop.mean + rmvt(
              sigma = (beta.prop.sigma+t(beta.prop.sigma))/2,
              n = 1, df = beta.prop.df)
        }
      ## }

###----------------------------------------------------------------------------
### REVERSE OF THE NEWTON METHOD
###----------------------------------------------------------------------------

      ## Newton method to approach the posterior for the proposed draw

      beta.prop.full <- matrix(0, nCovs, nPar)
      beta.prop.full[betaIdx.prop] <- beta.prop
      Mdl.beta.prop <- Mdl.beta.curr
      Mdl.beta.prop[[CompCaller]][[parCaller]] <- beta.prop.full

      Mdl.betaIdx.prop <- Mdl.betaIdx.curr
      Mdl.betaIdx.prop[[CompCaller]][[parCaller]] <- betaIdx.prop

      browser()

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

      if(beta.NTPropRev$errorFlag ## ||
         ## any(is.na(beta.prop.sigma)) ||
         ## is(try(chol(beta.prop.sigma), silent=TRUE), "try-error")
         )
        {## Something is wrong.
          errorFlags[iMH] <- TRUE
          betaIdx.prop <- betaIdx.curr
          next
        }

      ## The information for proposed density via K-step Newton's method
      beta.propRev.mean <- matrix(beta.NTPropRev[["param"]], 1) # 1-by-p
      beta.propRev.sigma <- -beta.NTPropRev[["HessObsInv"]] # p-by-p
      ## beta.propRev.sigma <- diag1(length(beta.propRev.mean))*0.1
      ## browser()
###----------------------------------------------------------------------------
### COMPUTING THE METROPOLIS-HASTINGS RATIO
###----------------------------------------------------------------------------
      ## if(all(!is.na(beta.propRev.sigma)))
      ##   {
      ## The jump density for proposed point at proposed mode
      logJump.propATprop <- dmvt(
          x = beta.prop - beta.prop.mean,
          sigma = (beta.prop.sigma+t(beta.prop.sigma))/2,
          df = beta.prop.df, log = TRUE)

      ## The jump density for curr draw at reverse proposed mode.
      logJump.currATpropRev<- dmvt(
          x = beta.curr - beta.propRev.mean,
          sigma = (beta.propRev.sigma+t(beta.propRev.sigma))/2,
          df = beta.prop.df, log = TRUE)

      ## if(is(logJump.currATpropRev, "try-error")) browser()

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

      ## compute the MH ratio.
      MHRatio <- exp(logPost.prop - logPost.curr +
                       logJump.currATpropRev - logJump.propATprop +
                         logJump.Idx.currATprop - logJump.Idx.propATcurr)

      ##  cat(logPost.prop, logPost.curr, logJump.currATpropRev, logJump.propATprop, "\n")


      ## the acceptance probability
      accept.prob.curr <- min(1, MHRatio)
      ## print(round(accept.prob*100))
      ## }
      ## }

###----------------------------------------------------------------------------
### THE MH ACCEPTANCE PROBABILITY AND ACCEPT/REJECT THE PROPOSED DRAW.
###----------------------------------------------------------------------------
      ## print(accept.prob.curr)
      ## if(is.na(accept.prob.curr)) browser()

      if(## rejectFlag == FALSE &&
          !is.na(accept.prob.curr) &&
          runif(1) < accept.prob.curr) # keep update
        {
          betaIdx.curr <- betaIdx.prop
          beta.curr <- beta.prop

          Mdl.beta.curr <- Mdl.beta.prop
          Mdl.betaIdx.curr <- Mdl.betaIdx.prop

          staticCache.curr <- staticCache.prop
          ## browser()
        }
      else # keep current
        {
          betaIdx.prop <- betaIdx.curr
        }


      accept.probs[iMH] <- accept.prob.curr
    }

  ## cat("accept.probs", accept.probs, "\n")
  ## cat("errorFlags:", errorFlags, "\n")
  ## cat("logPost.prop", logPost.prop, "\n")
  ## cat("logPost.curr", logPost.curr, "\n")
  ## if(is.na(accept.probs[nMH])) browser()

  ## The final update, the acceptance prob are from the last MH update or keep
  ## current draw.
  if(errorFlags[nMH] == TRUE)
    {
      ## epic failure
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
