##' This is the Newton move with dimension changes for the copula model.
##'
##' @title Generalized Newton method for the copula model
##' @param propArgs
##' @param varSelArgs
##' @param priArgs
##' @param betaIdxProp
##' @param parUpdate
##' @param CplNM
##' @param Mdl.Y
##' @param Mdl.X
##' @param Mdl.beta
##' @param Mdl.betaIdx
##' @param Mdl.parLink
##' @param MargisTypes
##' @param staticArgs
##' @param param.cur "matrix".
##'         The initial values for the Newton update.
##' @return "list". See bellow.
##' \item   {gradObs}
##'         {"matrix". The gradient}
##' \item   {HessObsInv}
##'         {"matrix". The inverse Hessian matrix.}
##' \item   {param}
##'         {"matrix". The updated parameters after K step Newton integrations.}
##' @references Li Villani Kohn 2010
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Wed Sep 29 17:18:22 CEST 2010;
##'       Current: Mon Mar 05 10:33:29 CET 2012.
GNewtonMove <- function(propArgs,
                        varSelArgs,
                        priArgs,
                        betaIdxProp,
                        parUpdate,
                        CplNM,
                        Mdl.Y,
                        Mdl.X,
                        Mdl.beta,
                        Mdl.betaIdx,
                        Mdl.parLink,
                        MargisTypes,
                        staticArgs)
{

  ## The updating component parameter chain
  chainCaller <- parCplCaller(parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]


  ## if(parCaller == "tau") browser()

  ## The current parameters
  X <- Mdl.X[[CompCaller]][[parCaller]]
  betaCurr <- Mdl.beta[[CompCaller]][[parCaller]] # p-by-1
  betaIdxCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]] # p-by-1

  ## Finite Newton move. K steps to approach the mode plus one more step to
  ## update the gradient and Hessian at i:th step.
  kSteps <- propArgs[[CompCaller]][[parCaller]][["algorithm"]][["ksteps"]]
  hessMethod <- propArgs[[CompCaller]][[parCaller]][["algorithm"]][["hess"]]

  ## Initialize the Newton move with the proposed variable selection indicator
  Mdl.betaIdx[[CompCaller]][[parCaller]] <- betaIdxProp
  param <- betaCurr[betaIdxCurr, , drop = FALSE]

  ## Initial update staticArgs for current Newton move

  staticArgs.curr <- logPost(
      CplNM = CplNM,
      Mdl.Y = Mdl.Y,
      Mdl.X = Mdl.X,
      Mdl.beta = Mdl.beta,
      Mdl.betaIdx = Mdl.betaIdx,
      Mdl.parLink = Mdl.parLink,
      varSelArgs = varSelArgs,
      MargisTypes = MargisTypes,
      priArgs = priArgs,
      parUpdate = parUpdate,
      staticArgs = staticArgs,
      staticArgsOnly = TRUE)[["staticArgs"]]

###----------------------------------------------------------------------------
### The k-step Generalized Newton Move
###----------------------------------------------------------------------------

  for(iStep in 1:(kSteps+1))
    {
      ## The gradient and Hessian in the likelihood
      logLikGradHess.prop <- logLikelihoodGradHess(
          CplNM = CplNM,
          Mdl.Y = Mdl.Y,
          Mdl.X = Mdl.X,
          Mdl.parLink = Mdl.parLink,
          Mdl.beta = Mdl.beta,
          MargisTypes = MargisTypes,
          Mdl.betaIdx = Mdl.betaIdx,
          parUpdate = parUpdate,
          varSelArgs = varSelArgs,
          staticArgs = staticArgs.curr)

      ## Gradient Hessian for the prior *including non selected covariates*
      ## NOTE: The Hessian matrix of the prior is also approximated, we should
      ## use the explicit Hessian whenever possible.
      logPriGradHess.prop <- logPriorsGradHess(
          Mdl.X = Mdl.X,
          Mdl.beta = Mdl.beta,
          Mdl.betaIdx = Mdl.betaIdx,
          Mdl.parLink = Mdl.parLink,
          varSelArgs = varSelArgs,
          priArgs = priArgs,
          chainCaller = chainCaller)

      ## Gradient and Hessian for the likelihood
      logLikGrad.prop <- logLikGradHess.prop[["logLikGradObs"]] # n-by-pp
      logLikHess.prop <- hessApprox(logLikGrad.prop, hessMethod)

      logPriGrad.prop <- logPriGradHess.prop[["gradObs"]] # pp-by-1
      logPriHess.prop <- logPriGradHess.prop[["HessObs"]] # pp-by-pp


      ## The gradient and Hessian subsets in the priors due to variable selection
      logPriGrad.pp <- logPriGrad.prop[betaIdxProp, , drop = FALSE] # pp-by-1
      logPriHess.pp <- logPriHess.prop[betaIdxProp, betaIdxProp, drop = FALSE] # pp-by-pp
      logPriHess.pc <- logPriHess.prop[betaIdxProp, betaIdxCurr, drop = FALSE] # pp-by-pc

      ## The selected covariates in the proposed and current draw
      X.prop <- X[ , betaIdxProp, drop = FALSE] # n-by-pp
      X.curr <- X[ , betaIdxCurr, drop = FALSE] # n-by-pc

      ## The gradient and Hessian in the general Newton's update
      gradObs.pp <- matrix(rowSums(Md(t(X.prop), logLikGrad.prop)) + logPriGrad.pp) # pp-by-1
      HessObs.pp <- tMdN(X.prop, logLikHess.prop, X.prop) + logPriHess.pp # pp-by-pp
      HessObs.pc <- tMdN(X.prop, logLikHess.prop, X.curr) + logPriHess.pc # pp-by-pc

      HessObsInv.pp <- solve(HessObs.pp) # pp-by-pp

      ## The general Newton's Update
      if((iStep <= kSteps))
        {
          ## if(iStep == 2) browser()
          ## update the proposed parameters via the general Newton formula
          param <- HessObsInv.pp%*%(HessObs.pc%*%param - gradObs.pp)

          ## Update the parameter with current updated results.
          ## If variable selection did not chose pth covariate,  then the pth
          ## coefficient is zero naturally. After the first
          ## variable-dimensional move, the algorithm switches to usual
          ## Newton's move.
          betaIdxCurr <- betaIdxProp

          ## the full parameters including zeros.
          param.full <- matrix(0, length(betaIdxCurr), 1)
          param.full[betaIdxProp] <- param
          Mdl.beta[[CompCaller]][[parCaller]] <- param.full

          ## Update the staticArgs
          ## Initial update staticArgs for current Newton move

          ## staticArgs.curr <- logPost(
          ##     CplNM = CplNM,
          ##     Mdl.Y = Mdl.Y,
          ##     Mdl.X = Mdl.X,
          ##     Mdl.beta = Mdl.beta,
          ##     Mdl.betaIdx = Mdl.betaIdx,
          ##     Mdl.parLink = Mdl.parLink,
          ##     varSelArgs = varSelArgs,
          ##     MargisTypes = MargisTypes,
          ##     priArgs = priArgs,
          ##     parUpdate = parUpdate,
          ##     staticArgs = staticArgs.curr,
          ##     staticArgsOnly = TRUE)[["staticArgs"]]

          staticArgs.curr <- try(logPost(
              CplNM = CplNM,
              Mdl.Y = Mdl.Y,
              Mdl.X = Mdl.X,
              Mdl.beta = Mdl.beta,
              Mdl.betaIdx = Mdl.betaIdx,
              Mdl.parLink = Mdl.parLink,
              varSelArgs = varSelArgs,
              MargisTypes = MargisTypes,
              priArgs = priArgs,
              parUpdate = parUpdate,
              staticArgs = staticArgs.curr,
              staticArgsOnly = TRUE)[["staticArgs"]], silent = FALSE)


          ## if(is(staticArgs.curr, "try-error")) browser()

        }
      else # (k+1):th step.  Make a output
        {
          out <- list(param = param,
                      gradObs = gradObs.pp,
                      HessObs = HessObs.pp,
                      HessObsInv = HessObsInv.pp,
                      staticArgs = staticArgs.curr)
          ## print(gradObs.pp)
        }
    }

  return(out)
}
