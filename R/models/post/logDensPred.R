##' The log predictive likelihood density.
##'
##' This is used for prediction and LPDS.
##' @param CplOut
##' @param Mdl.Idx.testing
##' @param Mdl.X.testing
##' @param Mdl.Y.testing
##' @param pred "character" The predictive likelihood for marginal or copula
##' likelihood.
##' @return "matrix" "mcmc sample by Lik.len"
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Feb 25 19:20:57 CET 2013; Current: Sat Jul 18 09:30:58 CST 2015.
logDensPred <- function(CplOut, Mdl.Idx.testing, Mdl.X.testing, Mdl.Y.testing,
                        MCMC.beta, MCMC.betaIdx, pred)
{
###----------------------------------------------------------------------------
### Extract the MMCMC output list
###----------------------------------------------------------------------------
  list2env(CplOut, envir = environment())

###----------------------------------------------------------------------------
### The testing covariates
###----------------------------------------------------------------------------
  ## Unless user specify the predict covariates, use the default in the
  ## configure files.
  if(missing(Mdl.X.testing) || missing(Mdl.Y.testing))
  {
    subsetFun <- function(x, idx)x[idx, , drop = FALSE]
    ## Mdl.Idx.testing <- crossValidIdx[["testing"]][[iCross]]

    Mdl.Y.testing <- rapply(object=Mdl.Y,
                            f = subsetFun,
                            idx = Mdl.Idx.testing,
                            how = "replace")

    if(any(rapply(Mdl.X, class) != "matrix"))
    {## Foreign marginal models are used.
      Mdl.X.Pred <- MargiModelForeignPred(MargisNM = MargisNM,
                                          MargisType = MargisType,
                                          Mdl.ForeignFit =Mdl.ForeignFit,
                                          Mdl.Y = Mdl.Y.testing)

      Mdl.X.testing <- c(Mdl.X.Pred[["Mdl.X"]],
                          rapply(object=Mdl.X[MargisNM[length(MargisNM)]],
                                 f = subsetFun,
                                 idx = Mdl.Idx.testing,
                                 how = "replace"))
    }
    else
    {## The native model structure
      Mdl.X.testing <- rapply(object=Mdl.X,
                              f = subsetFun,
                              idx = Mdl.Idx.testing,
                              how = "replace")
    }

  }

###----------------------------------------------------------------------------
### The MCMC burnin
###----------------------------------------------------------------------------

  n.burnin <- floor(nIter*burnin)
  nUsed <- nIter - n.burnin

  partiMethod <- crossValidArgs$partiMethod

  ## The sample indices for LPDS after burn-in Make a decreasing sequence and
  ## sort it by increasing order. Just to make sure last draw is allays used
  ## The exact size may not same as the result from sample proportion.
  MCMC.sampleIdxInv <- seq(from = nIter,
                           to = (nIter-nUsed+1),
                           by = -round(1/sampleProp))
  MCMC.sample.len <- length(MCMC.sampleIdxInv)
  MCMC.sampleIdx <- MCMC.sampleIdxInv[MCMC.sample.len:1]

  nPred <- length(Mdl.Y.testing[[1]])
  if(partiMethod == "time-series")
  {
    ## Special case for time series where the dependence are used The LPDS
    ## is approximated by computing each term p(y_{t+1}|y_{1:t}) using the
    ## same posterior sample base on data update to time t.  See Villani et
    ## al 2009 or Li et al 2010

    LikLst.Idx <- as.list(1:nPred) # length of nPred
  }
  else
  {
    ## The independent likelihood
    LikLst.Idx <- list(1:nPred) # length of 1
  }

  ## Allocate the log predictive matrix
  LikLst.len <- length(LikLst.Idx)

  out.pred <- matrix(NA, MCMC.sample.len, LikLst.len)
###----------------------------------------------------------------------------
### Calculate the predictive densities in all likelihood segments
###----------------------------------------------------------------------------
  ## The predictive likelihood
  if(tolower(pred) == "joint")
  {
    MCMCUpdateStrategy <- "joint"
    ## Use the whole copula likelihood
    parUpdate <- rapply(parUpdate,  function(x) TRUE,
                        how  = "replace")

  }
  else
  {
    MCMCUpdateStrategy <- "margin"

    ## The marginal likelihood
    parUpdate <- rapply(parUpdate,
                        function(x) FALSE,
                        how  = "replace")
    parUpdate[[pred]] <- rapply(parUpdate[[pred]],
                                function(x) TRUE,
                                how  = "replace")
  }


  for(i in 1:LikLst.len)
  {
    idx.curr <- LikLst.Idx[[i]]
    Mdl.X.testing.curr <- rapply(object=Mdl.X.testing,
                                 f = subsetFun,
                                 idx = idx.curr,
                                 how = "replace")
    Mdl.Y.testing.curr <- rapply(object=Mdl.Y.testing,
                                 f = subsetFun,
                                 idx = idx.curr,
                                 how = "replace")

    which.j <- 0
    for(j in MCMC.sampleIdx) ## Just the likelihood function with posterior samples
    {

      Mdl.beta.curr <- rapply(object=MCMC.beta,
                              f = subsetFun,
                              idx = j,
                              how = "replace")

      ## Mdl.betaIdx.curr <- rapply(object=MCMC.betaIdx,
      ##                            f = subsetFun,
      ##                            idx = j,
      ##                            how = "replace")

      ## The log predictive likelihood.  Note that the corresponding
      ## updating flags should be switched on
      Mdl.par <- parCplMeanFun(Mdl.X = Mdl.X.testing.curr,
                               Mdl.parLink = Mdl.parLink,
                               Mdl.beta = Mdl.beta.curr,
                               parUpdate = parUpdate,
                               Mdl.par = parUpdate)

      Mdl.ud <- logDens(MargisType = MargisType,
                        Mdl.Y = Mdl.Y.testing.curr,
                        Mdl.par = Mdl.par,
                        Mdl.u = staticCache$Mdl.u,
                        Mdl.d = staticCache$Mdl.d,
                        parUpdate = parUpdate,
                        MCMCUpdateStrategy = MCMCUpdateStrategy)
      Mdl.d <- Mdl.ud[["Mdl.d"]]
      ## Mdl.u <- Mdl.ud[["Mdl.u"]]
      ## Mdl.PostComp <- Mdl.ud[["Mdl.PostComp"]]

      ## logPred <- logPost(MargisType = MargisType,
      ##                    Mdl.Y = Mdl.Y.testing.curr,
      ##                    Mdl.X = Mdl.X.testing.curr,
      ##                    Mdl.beta = Mdl.beta.curr,
      ##                    Mdl.betaIdx = Mdl.betaIdx.curr,
      ##                    Mdl.parLink = Mdl.parLink,
      ##                    varSelArgs = varSelArgs,
      ##                    priArgs = priArgs,
      ##                    parUpdate = parUpdate,
      ##                    MCMCUpdateStrategy = MCMCUpdateStrategy
      ##                    )[["Mdl.logLik"]]

      which.j <- which.j + 1

      out.pred[which.j, i] <- sum(Mdl.d)
    }


    ## Simple progress bar
    ## progressbar(((iCross-1)*nSample + which.j), nFold*nSample)
  }

  if(partiMethod == "time-series")
  {
    ## Sum to get the loglikelihood when the likelihood are calculated
    ## conditionally.
    out <- matrix(apply(out.pred, 1, sum))
  }
  else
  {
    out <- out.pred
  }
  return(out)
}
