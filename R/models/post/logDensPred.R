##' The MCMC samples for log predictive likelihood/density.
##'
##' This is used for calculating mean prediction, posterior creditable interval and LPDS.
##' @param CplOut
##' @param Mdl.Idx.testing
##' @param Mdl.X.testing
##' @param Mdl.Y.testing
##' @param pred "character" The predictive likelihood for marginal or copula likelihood.
##' @return "matrix" "No. of MCMC samples-by- length of predictive likelihood/density"
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

  ## MCMC.sampleProp
###----------------------------------------------------------------------------
### The testing covariates
###----------------------------------------------------------------------------
  ## Unless user specify the predict covariates, use the default in the
  ## configure files.
  if(missing(Mdl.X.testing) || missing(Mdl.Y.testing))
  {
    subsetFun <- function(x, idx)
    {
      x[idx, , drop = FALSE]
    }
    ## Mdl.Idx.testing <- crossValidIdx[["testing"]][[iCross]]

    Mdl.Y.testing <- rapply(object=Mdl.Y, f = subsetFun,
                            idx = Mdl.Idx.testing, how = "replace")

    if(any(rapply(Mdl.X, class) != "matrix"))
    {## Foreign marginal models are used.
      Mdl.X.Pred <- MargiModelForeignPred(MargisNM = MargisNM,
                                          MargisType = MargisType,
                                          Mdl.ForeignFit =Mdl.ForeignFit,
                                          Mdl.Y = Mdl.Y.testing)

      Mdl.X.testing <- c(Mdl.X.Pred[["Mdl.X"]],
                         rapply(object=Mdl.X[MargisNM[length(MargisNM)]],
                                f = subsetFun, idx = Mdl.Idx.testing,
                                how = "replace"))
    }
    else
    {## The native model structure
      Mdl.X.testing <- rapply(object=Mdl.X, f = subsetFun,
                              idx = Mdl.Idx.testing, how = "replace")
    }

  }

###----------------------------------------------------------------------------
### The MCMC burnin
###----------------------------------------------------------------------------

  n.burnin <- floor(MCMC.nIter*MCMC.burninProp)
  nUsed <- MCMC.nIter - n.burnin

  ## The sample indices for LPDS after burn-in Make a decreasing sequence and
  ## sort it by increasing order. Just to make sure last draw is always used
  ## The exact size may not same as the result from sample proportion.
  MCMC.sampleIdxRev <- seq(from = MCMC.nIter,
                           to = (MCMC.nIter-nUsed+1),
                           by = -round(1/MCMC.sampleProp))
  MCMC.sample.len <- length(MCMC.sampleIdxRev)
  MCMC.sampleIdx <- MCMC.sampleIdxRev[MCMC.sample.len:1]

  nPred <- length(Mdl.Y.testing[[1]])

  partiMethod <- crossValidArgs[["partiMethod"]]

  ## The independent likelihood. TODO: Special case to be considered for time series where
  ## the dependence are used The LPDS is approximated by computing each term
  ## p(y_{t+1}|y_{1:t}) using the same posterior sample base on data update to time t.
  ## See Villani et al 2009 or Li et al 2010
  LikLst.Idx <- 1:nPred # length of 1

  ## Allocate the log MCMC predictive matrix
  MCMC.logPredDens <- matrix(NA, MCMC.sample.len, 1)

  MCMC.Y.Pred <- list()
  for(iComp in (names(MargisType)[-length(MargisType)]))
  {
    MCMC.Y.Pred[[iComp]] <- array(NA, c(MCMC.sample.len, nPred, 4),
                                  dimnames = list(NULL, rownames(Mdl.X.testing[[1]][[1]]),
                                                  c("mean", "variance", "skewness",  "kurtosis")))
  }

###----------------------------------------------------------------------------
### Calculate the predictive densities in all likelihood segments
###----------------------------------------------------------------------------
  ## Specify the predictive likelihood components
  if(tolower(pred) == "joint")
  { ## Use the whole copula likelihood/density
    MCMCUpdateStrategy4LPDS <- "joint"
    parUpdate <- rapply(parUpdate,  function(x) TRUE, how  = "replace")

  }
  else
  {## Use the marginal likelihood/density
    MCMCUpdateStrategy4LPDS <- "margin"
    parUpdate <- rapply(parUpdate, function(x) FALSE, how  = "replace")
    parUpdate[[pred]] <- rapply(parUpdate[[pred]], function(x) TRUE, how  = "replace")
  }

  ## Calculate the predictive likelihood
  Mdl.X.testing.curr <- rapply(object=Mdl.X.testing, f = subsetFun,
                               idx = LikLst.Idx, how = "replace")
  Mdl.Y.testing.curr <- rapply(object=Mdl.Y.testing, f = subsetFun,
                               idx = LikLst.Idx, how = "replace")

  j <- 1
  for(iMCMC.sampleIdx in MCMC.sampleIdx)
  {## Just the likelihood function with posterior samples

    subsetFun4beta <- function(x, idx)
    {
      if((dim(x)[1] == 1 && length(idx)>1) ||
         (dim(x)[1] == 1 && length(idx) == 1 && idx != 1))
      {# check whether some parameters are not updated
        out <- x
      }
      else
      {
        out <- x[idx, , drop = FALSE]
      }
      return(out)
    }

    Mdl.beta.curr <- rapply(object=MCMC.beta,
                            f = subsetFun4beta,
                            idx = iMCMC.sampleIdx,
                            how = "replace")

    ## The log predictive likelihood.  Note that the corresponding
    ## updating flags should be switched on
    Mdl.par.curr <- parCplMeanFun(Mdl.X = Mdl.X.testing.curr,
                                  Mdl.parLink = Mdl.parLink,
                                  Mdl.beta = Mdl.beta.curr,
                                  parUpdate = parUpdate,
                                  Mdl.par = parUpdate)
    Mdl.Y.pred.curr <- logCplPredict(MargisType = MargisType,
                                     Mdl.par = Mdl.par.curr)

    for(iComp in names(Mdl.Y.pred.curr))
    {
      MCMC.Y.Pred[[iComp]][j, , ] <-  Mdl.Y.pred.curr[[iComp]]
    }

    Mdl.ud <- logDens(MargisType = MargisType,
                      Mdl.Y = Mdl.Y.testing.curr,
                      Mdl.par = Mdl.par.curr,
                      parUpdate = parUpdate,
                      MCMCUpdateStrategy = MCMCUpdateStrategy4LPDS)
    Mdl.d <- Mdl.ud[["Mdl.d"]]
    Mdl.PostComp <- Mdl.ud[["Mdl.PostComp"]]

    MCMC.logPredDens[j, ] <- sum(Mdl.d[, unlist(Mdl.PostComp)])
    j <- j + 1

    ## Simple progress bar
    ## progressbar(((iCross-1)*nSample + which.j), nFold*nSample)
  }

  ## "mean", "variance", "skewness",  "kurtosis"
  applyFun <- function(x, mar, fun,...) {apply(x, mar, fun, ...)}
  MVSK <- list()
  MVSK[["mean"]] <- lapply(MCMC.Y.Pred, applyFun, mar = c(2, 3), fun = mean)
  MVSK[["var"]] <- lapply(MCMC.Y.Pred, applyFun, mar = c(2, 3), fun = var)
  MVSK[["HPDL"]] <- lapply(MCMC.Y.Pred, applyFun, mar = c(2, 3),
                           fun = quantile, probs = 0.025)
  MVSK[["HPDU"]] <- lapply(MCMC.Y.Pred, applyFun, mar = c(2, 3),
                           fun = quantile, probs = 0.975)

  ## Residuals and the uncertainty
  MCMC.residual <- mapply(function(x, y) matrix(x, nrow(y), ncol(y), byrow = TRUE)-y,
                          Mdl.Y.testing, lapply(MCMC.Y.Pred, function(x) x[, ,"mean"]),
                          SIMPLIFY = FALSE)
  RESID <- list()
  RESID[["mean"]] <- lapply(MCMC.residual, applyFun, mar = 2, fun = mean)
  RESID[["var"]] <- lapply(MCMC.residual, applyFun, mar = 2, fun = var)
  RESID[["HPDL"]] <- lapply(MCMC.residual, applyFun, mar  = 2,
                            fun = quantile, probs = 0.025)
  RESID[["HPDU"]] <- lapply(MCMC.residual, applyFun, mar = 2,
                            fun = quantile, probs = 0.975)

  out <- list(MCMC.logPredDens = MCMC.logPredDens,
              MVSK = MVSK,
              RESID = RESID,
              MCMC.Y.Pred = MCMC.Y.Pred)

  return(out)
}
