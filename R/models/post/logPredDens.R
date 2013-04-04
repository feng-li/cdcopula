##' The log predictive likelihood density.
##'
##' This is used for prediction and LPDS.
##' @param CplOut
##' @param Testing.Idx
##' @param MdlTesting.X
##' @param MdlTesting.Y
##' @param pred "character" The predictive likelihood for marginal or copula
##' likelihood.
##' @return "matrix" "mcmc sample by Lik.len"
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Feb 25 19:20:57 CET 2013;
##'       Current: Mon Feb 25 19:21:03 CET 2013.
logPredDens <- function(CplOut, Testing.Idx, MdlTesting.X, MdlTesting.Y, pred = CplNM)
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
    if(missing(MdlTesting.X) || missing(MdlTesting.Y))
      {
        subsetFun <- function(x, idx)x[idx, , drop = FALSE]
        ## Testing.Idx <- crossValidIdx[["testing"]][[iCross]]

        MdlTesting.X <- rapply(object=Mdl.X,
                               f = subsetFun,
                               idx = Testing.Idx,
                               how = "replace")
        MdlTesting.Y <- rapply(object=Mdl.Y,
                               f = subsetFun,
                               idx = Testing.Idx,
                               how = "replace")
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

    nPred <- length(MdlTesting.Y[[1]])
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
    for(i in 1:LikLst.len)
      {
        idx.curr <- LikLst.Idx[[i]]
        MdlTesting.X.curr <- rapply(object=MdlTesting.X,
                                    f = subsetFun,
                                    idx = idx.curr,
                                    how = "replace")
        MdlTesting.Y.curr <- rapply(object=MdlTesting.Y,
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

            Mdl.betaIdx.curr <- rapply(object=MCMC.betaIdx,
                                       f = subsetFun,
                                       idx = j,
                                       how = "replace")

            ## The log predictive likelihood.
            ## Note that all the updating flags should be switched on
            logPredLik <- logPost(
                CplNM = CplNM,
                Mdl.Y = MdlTesting.Y.curr,
                Mdl.X = MdlTesting.X.curr,
                Mdl.beta = Mdl.beta.curr,
                Mdl.betaIdx = Mdl.betaIdx.curr,
                Mdl.parLink = Mdl.parLink,
                varSelArgs = varSelArgs,
                MargisTypes = MargisTypes,
                priArgs = priArgs,
                parUpdate = rapply(parUpdate, function(x) TRUE, how = "replace"),
                call.out = c("likelihood"))[["Mdl.logLik"]]

            ## The predictive likelihood
            if(pred == CplNM)
              {
                ## The whole copula likelihood
                logPred <-  sum(logPredLik)
              }
            else
              {
                ## Particular margin
                logPred <- sum(logPredLik[, pred])
              }

            which.j <- which.j + 1

            out.pred[which.j, i] <- logPred
          }


        ## Simple progress bar
        ## progressbar(((iCross-1)*nSample + which.j), nFold*nSample)
      }

    if(partiMethod == "time-series")
      {
        ## Sum to get the loglikelihood when the likelihood are calumniated
        ## conditionally.
        out <- matrix(apply(out.pred, 1, sum))
      }
    else
      {
        out <- out.pred
      }
    return(out)
  }
