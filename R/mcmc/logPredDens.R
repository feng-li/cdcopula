##' <title>
##'
##' <description>
##' @param CplMain.out
##' @param Testing.Idx
##' @param MdlTesting.X
##' @param MdlTesting.Y
##' @return "matrix" "mcmc sample by Lik.len"
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Feb 25 19:20:57 CET 2013;
##'       Current: Mon Feb 25 19:21:03 CET 2013.
logPredDens <- function(CplMain.out,Testing.Idx = NA,
                        MdlTesting.X = NA, MdlTesting.Y = NA)
  {

###----------------------------------------------------------------------------
### Extract the MMCMC output list
###----------------------------------------------------------------------------
    list2env(CplMain.out)

###----------------------------------------------------------------------------
### The testing covariates
###----------------------------------------------------------------------------
    ## Unless user specify the predict covariates, use the default in the
    ## configure files.
    subsetFun <- function(x, idx)x[idx, , drop = FALSE]
    if(is.na(MdlTesting.X) || is.na(MdlTesting.Y))
      {
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
    MCMC.sampleIdx <- MCMC.sampleIdxInv[length(MCMC.sampleIdxInv):1]

    nPred <- length(MdlTesting.Y[[1]])
    if(partiMethod == "time-series")
      {
        ## Special case for time series where the dependence are used
        ## The LPDS is approximated by computing each term p(y_{t+1}|y_{1:t})
        ## using the same posterior sample base on data upda to time t.
        ## See Villani et al 2009 or Li et al 2010
        LikLst.Idx <- as.list(1:nPred)
      }
    else
      {
        LikLst.Idx <- list(1:nPred)
      }

    ## Allocate the log predictive matrix
    LikLst.len <- length(LikLst.Idx)
    MCMC.sample.len <- length(MCMC.sample.len)
    out <- matrix(NA, MCMC.sample.len, LikLst.len)

###----------------------------------------------------------------------------
### Calculate the predictive densities for all folds
###----------------------------------------------------------------------------

    for(i in 1:LikLst.len)
      {
        idx.curr <- LikLst.Idx[[i]]
        MdlTesting.X.curr <- rapply(object=MdlTesting.X,
                                    f = subsetFun,
                                    idx = idx.curr,
                                    how = "replace")
        MdlTesting.Y.curr <- rapply(object=MdlTesting.Y.curr,
                                    f = subsetFun,
                                    idx = idx.curr,
                                    how = "replace")


        which.j <- 0
        for(j in MCMC.sampleIdx) ## Just the likelihood function with posterior samples
          {

            ## Params.j <- lapply(OUT.Params, function(x) apply(x[, , j, iCross, drop =
            ##                                                    FALSE], c(1, 2), "["))
            ## caller.log.like <- call(logpost.fun.name,Y = Y.iTesting, x = x.iTesting,
            ##                         Params = Params.j, callParam = list(id =
            ##                                              "likelihood"), priorArgs =
            ##                         priorArgs, splineArgs = splineArgs, Params_Transform
            ##                         = Params_Transform)

            ## log.like <- eval(caller.log.like)

            ## logPost(CplNM, Mdl.Y, Mdl.X,Mdl.beta,Mdl.betaIdx,Mdl.parLink,
            ##         varSelArgs,MargisTypes,priArgs,parUpdate,
            ##         staticCache, staticCacheOnly = FALSE, parUpdate4Pri = parUpdate)


            Mdl.beta.curr <- rapply(object=MCMC.beta,
                                    f = subsetFun,
                                    idx = j,
                                    how = "replace")

            Mdl.betaIdx.curr <- rapply(object=MCMC.betaIdx,
                                       f = subsetFun,
                                       idx = j,
                                       how = "replace")

            logPred <- logPost(
                CplNM = CplNM,
                Mdl.Y = MdlTesting.X.curr,
                Mdl.X = MdlTesting.Y.curr,
                Mdl.beta = Mdl.beta.curr,
                Mdl.betaIdx = Mdl.betaIdx.curr,
                Mdl.parLink = Mdl.parLink,
                varSelArgs = varSelArgs,
                MargisTypes = MargisTypes,
                priArgs = priArgs,
                parUpdate = parUpdate,
                staticCache = staticCache,
                call.out = "likelihood")

            which.j <- which.j + 1
            out[which.j, i] <- logPred
          }


        ## Simple progress bar
        ## progressbar(((iCross-1)*nSample + which.j), nFold*nSample)
      }
    return(out)
  }
