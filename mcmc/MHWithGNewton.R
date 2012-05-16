##' Metropolis–Hastings algorithm with K-step Newton method with variable
##' selection.
##'
##' @title
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
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Thu Feb 17 14:03:14 CET 2011;
##'       Current: Wed Feb 01 16:09:04 CET 2012.
##' DEPENDS: mvtnorm
MHWithGNewton <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                          Mdl.parLink, parUpdate, priArgs, varSelArgs,
                          propArgs, MargisTypes, staticArgs)
{

###----------------------------------------------------------------------------
### VARIABLE SELECTION PROPOSAL
### TODO: Write a general variable selection scheme
###----------------------------------------------------------------------------

  ## The updating component parameter chain
  chainCaller <- parCaller(parUpdate)
  CompCurr <- chainCaller[1]
  parCurr <- chainCaller[2]

  ## The current variable selection indicators
  betaIdxCurr <- Mdl.betaIdx[[CompCurr]][[parCurr]]

  ## Randomly propose a subset
  ## TODO: A general proposal function
  varSelCand <- varSelArgs[[CompCurr]][[parCurr]]$cand
  IdxArgs <- propArgs[[CompCurr]][[parCurr]][["indicators"]]
  betaIdxProp <- betaIdxCurr # The proposal is base on current status
  if(IdxArgs$type == "binom")
    {
      ## Binomial proposal a small subset
      betaIdxPropCand <- which(rbinom(n = length(varSelCand) , size = 1,
                                      prob = IdxArgs$prob) == 1)
      ## Special case to make sure at least one variable is proposed a change
      if(length(betaIdxPropCand) == 0)
        {
          betaIdxPropCand <- sample(varSelCand, 1)
        }
    }

  ## The proposed variable selection indicators
  betaIdxProp[betaIdxPropCand] <- 1 - betaIdxCurr[betaIdxPropCand]

###----------------------------------------------------------------------------
### Make good proposal via K-steps Newton's method
###----------------------------------------------------------------------------

  ## Newton method to approach the posterior for the current draw
  GNewton1 <- GNewtonMove(propArgs = propArgs,
                          varSelArgs = varSelArgs,
                          priArgs = priArgs,
                          betaIdxProp = betaIdxProp,
                          parUpdate = parUpdate,
                          CplNM = CplNM,
                          Mdl.Y = Mdl.Y,
                          Mdl.X = Mdl.X,
                          Mdl.parLink = Mdl.parLink,
                          Mdl.beta = Mdl.beta,
                          Mdl.betaIdx = Mdl.betaIdx,
                          MargisTypes = MargisTypes,
                          staticArgs = staticArgs)

  ## The information for proposed density via K-step Newton's method
  param.cur.prop <- GNewton1$param # mean information
  HessObs.cur.prop <- GNewton1$HessObs # variance information
  invHessObs.cur.prop <- GNewton1$HessObsInv

  browser()
  ## Check if it is a good proposal
  if(any(is.na(HessObs.cur.prop)) ||
     is(try(chol(-invHessObs.cur.prop), silent=T), "try-error"))
    {
      ## Something is wrong.
      logjump.cur2prop <- NaN
      param.prop <- NaN
    }
  else # Continues with the Metropolis-Hastings
    {
      ## Propose a draw from multivariate t-distribution based on the proposed
      ## An idea (out of loud) : Generate a ## matrix of param. Select the one
      ## that give max acceptance probability.TODO: more general proposal type.

      propArgs.beta <- propArgs[[CompCurr]][[parCurr]][["beta"]]
      prop.type <- propArgs.beta[["type"]]
      prop.df <- propArgs.beta[["df"]]

      param.prop <- rmvt(n = 1,
                         sigma = -HessObs.cur.prop,
                         df = prop.df) + matrix(param.cur.prop, 1)

      ## The jump density from the proposed draw
      logjump.cur2prop <- dmvt(x = param.prop - matrix(param.cur.prop, 1),
                               sigma = -HessObs.cur.prop,
                               df = prop.df, log = TRUE)
    }

###----------------------------------------------------------------------------
### Jump back via K-steps Newton's method
###----------------------------------------------------------------------------

  if(any(is.na(param.prop))) # Previous proposal unsuccessful
    {
      HessObs.prop.prop <- NaN
    }
  else # all are right
    {
      ## Newton method to approach the posterior for the proposed draw
      GNewton2 <- GNewtonMove(propArgs = propArgs,
                              varSelArgs = varSelArgs,
                              priArgs = priArgs,
                              betaIdxProp = betaIdxProp,
                              parUpdate = parUpdate,
                              CplNM = CplNM,
                              Mdl.Y = Mdl.Y,
                              Mdl.X = Mdl.X,
                              Mdl.parLink = Mdl.parLink,
                              Mdl.beta = Mdl.beta,
                              Mdl.betaIdx = Mdl.betaIdx,
                              MargisTypes = MargisTypes,
                              staticArgs = staticArgs)

      ## The information for proposed density via K-step Newton's method
      param.prop.prop <- GNewton2$param
      HessObs.prop.prop <- GNewton2$hessObs
      invHessObs.prop.prop <- GNewton2$invHessObs
    }

  if(any(is.na(HessObs.prop.prop))) # Something is wrong at GNewton2,  reject it.
    {
      logpost.cur <- NaN
      logpost.prop <- NaN
      logjump.prop2cur <- NaN
    }
  else # all are right
    {
      ## The jump density from propose draw to current draw.
      logjump.prop2cur <- dmvt(x = param.cur - mattrix(param.prop.prop, 1),
                               sigma = -invHessObs.prop.prop,
                               df = prop.df, log = TRUE)

      ## Params.prop <- Params
      ## Params.prop[[callParam$id]][callParam$subset] <- param.prop
      ## Params.cur <- Params
      ## Params.cur[[callParam$id]][callParam$subset] <- param.cur

      ## The log posterior for the proposed draw
      logpost.prop <- logPost(CplNM = CplNM,
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

      ## The log posterior for the current draw
      logpost.cur <- logPost(CplNM = CplNM,
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
    }

###----------------------------------------------------------------------------
### Compute the MH ratio and make decision to update or keep current draw.
###----------------------------------------------------------------------------

  ## compute the MH ratio.
  log.r <- logpost.prop - logpost.cur + logjump.prop2cur - logjump.cur2prop
  r <- exp(log.r)

  ## the acceptance probability
  accept.prob <- min(1, r)

  if(!is.na(accept.prob) && runif(1) < accept.prob)
    {
      param.out <- param.prop  # keep update
    }
  else
    {
      param.out <- param.cur # keep current
      accept.prob <- 0 # set acceptance probs to zero.
    }

###----------------------------------------------------------------------------
### The output
###----------------------------------------------------------------------------
  out <- list(param.out = param.out, accept.prob = accept.prob)

  browser()

  return(out)
}
