##' The MCMC samples for log predictive likelihood/density.
##'
##' This is used for calculating mean prediction, posterior creditable interval and LPDS.
##' @param CplFitted
##' @param Mdl.Idx.testing
##' @param Mdl.X.testing
##' @param Mdl.Y.testing
##' @param MCMC.beta
##' @param PredDens "character" The predictive likelihood for marginal or copula
##'     likelihood.
##' @return "matrix" "No. of MCMC samples-by- length of predictive likelihood/density"
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Feb 25 19:20:57 CET 2013; Current: Sat Jul 18 09:30:58 CST 2015.
logDensPred <- function(CplFitted, Mdl.Idx.testing, Mdl.X.testing, Mdl.Y.testing)
{
###----------------------------------------------------------------------------
### Extract the MMCMC output list
###----------------------------------------------------------------------------
    MCMC.nIter <- CplFitted[["MCMC.nIter"]]
    MCMC.burninProp <- CplFitted[["MCMC.burninProp"]]
    MCMC.sampleProp <- CplFitted[["MCMC.sampleProp"]]
    Mdl.crossValidArgs <- CplFitted[["Mdl.crossValidArgs"]]
    MCMC.beta <- CplFitted[["MCMC.beta"]]
    Mdl.parLink <- CplFitted[["Mdl.parLink"]]
    Mdl.Y <- CplFitted[["Mdl.Y"]]
    Mdl.X <- CplFitted[["Mdl.X"]]
    Mdl.ForeignFitted <- CplFitted[["Mdl.ForeignFitted"]]
    Mdl.MargisNM <- CplFitted[["Mdl.MargisNM"]]
    Mdl.MargisType <- CplFitted[["Mdl.MargisType"]]
    PredDens <- CplFitted[["LPDS"]]

    ## list2env(CplFitted, envir = environment())

    ## Check if testing data Y is there
    subsetFun <- function(x, idx) {x[idx, , drop = FALSE]}
    if(missing(Mdl.Y.testing))
    {
        Mdl.Y.testing <- rapply(object = Mdl.Y, f = subsetFun,
                                idx = Mdl.Idx.testing, how = "replace")
    }

###----------------------------------------------------------------------------
### Foreign multivariate model/methods used in fitting/prediction
###----------------------------------------------------------------------------
    if(length(Mdl.ForeignFitted)>0)
    {
        ForeignModelType <- Mdl.MargisType[length(Mdl.MargisType)]

        out <- ModelForeignPred(model = ForeignModelType,
                                fitted.model = Mdl.ForeignFitted,
                                data.pred = Mdl.Y.testing)
        return(out)
    }

###----------------------------------------------------------------------------
### The testing covariates
###----------------------------------------------------------------------------
    ## Unless user specify the predict covariates, use the default in the
    ## configure files.

    if(missing(Mdl.X.testing))
    {
        if(any(rapply(Mdl.X, class) != "matrix"))
        {
            ## Foreign marginal models are used.  We only allow either all margins are
            ## foreign or all should be native. mixing is not implemented yet!
            Mdl.MargisForeignFitted <- CplFitted[["Mdl.MargisForeignFitted"]]
            Mdl.X.Pred <- MargiModelForeignPred(Mdl.MargisNM = Mdl.MargisNM,
                                                Mdl.MargisType = Mdl.MargisType,
                                                Mdl.MargisForeignFitted = Mdl.MargisForeignFitted,
                                                Mdl.Y = Mdl.Y.testing)

            Mdl.X.testing <- c(Mdl.X.Pred[["Mdl.X"]],
                               rapply(object=Mdl.X[Mdl.MargisNM[length(Mdl.MargisNM)]],
                                      f = subsetFun, idx = Mdl.Idx.testing,
                                      how = "replace"))
        }
        else
        {   ## The native model structure
            Mdl.X.testing <- rapply(object=Mdl.X, f = subsetFun,
                                    idx = Mdl.Idx.testing, how = "replace")
        }
    }
###----------------------------------------------------------------------------
### The MCMC burnin
###----------------------------------------------------------------------------
    n.burnin <- floor(MCMC.nIter*MCMC.burninProp)
    nUsed <- MCMC.nIter - n.burnin

    ## The sample indices for LPDS after burn-in Make a decreasing sequence and sort it by
    ## increasing order. Just to make sure last draw is always used The exact size may not
    ## same as the result from sample proportion.
    MCMC.sampleIdxRev <- seq(from = MCMC.nIter,
                             to = (MCMC.nIter-nUsed+1),
                             by = -round(1/MCMC.sampleProp))
    MCMC.sample.len <- length(MCMC.sampleIdxRev)
    MCMC.sampleIdx <- MCMC.sampleIdxRev[MCMC.sample.len:1]

    nPred <- length(Mdl.Y.testing[[1]])

    partiMethod <- Mdl.crossValidArgs[["partiMethod"]]

    ## The independent likelihood. TODO: Special case to be considered for time series
    ## where the dependence are used The LPDS is approximated by computing each term
    ## p(y_{t+1}|y_{1:t}) using the same posterior sample base on data update to time t.
    ## See Villani et al 2009 or Li et al 2010
    LikLst.Idx <- 1:nPred # length of nPred

###----------------------------------------------------------------------------
### Calculate the predictive densities in all likelihood segments
###----------------------------------------------------------------------------
    Mdl.X.testing.curr <- rapply(object=Mdl.X.testing, f = subsetFun,
                                 idx = LikLst.Idx, how = "replace")
    ## Define parUpdate
    if(!is.null(PredDens))
    {

        ## Allocate the log MCMC predictive matrix
        Mdl.logPredDens <- matrix(NA, MCMC.sample.len, length(PredDens),
                                  dimnames = list(NULL, PredDens))

        Mdl.Y.testing.curr <- rapply(object=Mdl.Y.testing, f = subsetFun,
                                     idx = LikLst.Idx, how = "replace")

        ## Specify the predictive likelihood components
        if("joint" %in% tolower(PredDens) |
           tolower(names(MCMC.beta)[length(MCMC.beta)]) %in% tolower(PredDens))
        { ## Use the whole copula likelihood/density
            MCMC.UpdateStrategy4LPDS <- "joint"
            parUpdate <- rapply(MCMC.beta,  function(x) TRUE, how  = "replace")
        }
        else
        {## Use the marginal likelihood/density
            MCMC.UpdateStrategy4LPDS <- "margin"
            parUpdate <- rapply(MCMC.beta, function(x) FALSE, how  = "replace")
            parUpdate[PredDens] <- rapply(parUpdate[PredDens],
                                          function(x) TRUE, how  = "replace")
        }
    }

    Mdl.PredY <- list()
    CompUpdate <- lapply(parUpdate, function(x) any(unlist(x) == TRUE))
    for(iComp in (names(Mdl.MargisType)[-length(Mdl.MargisType)]))
    {
        if(CompUpdate[[iComp]])
        {
            Mdl.PredY[[iComp]] <- array(NA, c(MCMC.sample.len, nPred, 4),
                                        dimnames = list(NULL, rownames(Mdl.X.testing[[1]][[1]]),
                                                        c("mean", "variance",
                                                          "skewness",  "kurtosis")))
        }
    }


    ## Loop over all MCMC.sampleIdx to obtain Mdl.logPredDens
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
    j <- 1
    for(iMCMC.sampleIdx in MCMC.sampleIdx)
    {   ## Just the likelihood function with posterior samples
        Mdl.beta.curr <- rapply(object=MCMC.beta, f = subsetFun4beta,
                                idx = iMCMC.sampleIdx, how = "replace")

        ## The log predictive likelihood.  Note that the corresponding updating flags
        ## should be switched on
        Mdl.par.curr <- parCplMeanFun(Mdl.X = Mdl.X.testing.curr,
                                      Mdl.parLink = Mdl.parLink,
                                      Mdl.beta = Mdl.beta.curr,
                                      parUpdate = parUpdate)
        Mdl.Y.pred.curr <- logCplPredict(Mdl.MargisType = Mdl.MargisType,
                                         Mdl.par = Mdl.par.curr)

        for(iComp in names(Mdl.PredY))
        {
            Mdl.PredY[[iComp]][j, , ] <-  Mdl.Y.pred.curr[[iComp]]
        }

        if(!is.null(PredDens))
        {
            Mdl.ud <- logDens(Mdl.MargisType = Mdl.MargisType,
                              Mdl.Y = Mdl.Y.testing.curr,
                              Mdl.par = Mdl.par.curr,
                              parUpdate = parUpdate,
                              MCMC.UpdateStrategy = MCMC.UpdateStrategy4LPDS)
            Mdl.d <- Mdl.ud[["Mdl.d"]]
            Mdl.PostComp <- Mdl.ud[["Mdl.PostComp"]]

            Mdl.Lik <- apply(Mdl.d[, unlist(Mdl.PostComp), drop = FALSE], 2, sum)
            if(length(PredDens) == 1)
            {
                if(PredDens  == names(MCMC.beta)[length(MCMC.beta)])
                {
                    ## PredDens is the copula density only
                    Mdl.logPredDens[j, PredDens] <- Mdl.Lik[PredDens]
                }
                else if(tolower(PredDens)  == "joint")
                {
                    ## PredDens is the joint density only
                    Mdl.logPredDens[j, PredDens] <- sum(Mdl.Lik)
                }
                else
                {
                    Mdl.logPredDens[j, PredDens] <- Mdl.Lik
                }
            }
            else
            {
                ## First update margins
                CommMargin <- intersect(PredDens, names(Mdl.Lik))
                Mdl.logPredDens[j, CommMargin] <- Mdl.Lik[CommMargin]

                ## Then update joint
                if("joint" %in% tolower(PredDens))
                {
                    Mdl.logPredDens[j, "joint"] <- sum(Mdl.Lik)
                }
            }

        }

        j <- j + 1

        ## Simple progress bar
        ## progressbar(((iCross-1)*nSample + which.j), nFold*nSample)
    }
    ## print(PredDens)

    ## "mean", "variance", "skewness",  "kurtosis" p(Y_p|Y_old)
    MVSK <- list()
    applyFun <- function(x, mar, fun,...){apply(x, mar, fun, ...)}
    MVSK[["mean"]] <- lapply(Mdl.PredY, applyFun, mar = c(2, 3), fun = mean)
    MVSK[["var"]] <- lapply(Mdl.PredY, applyFun, mar = c(2, 3), fun = var)
    MVSK[["HPDL"]] <- lapply(Mdl.PredY, applyFun, mar = c(2, 3),
                             fun = quantile, probs = 0.025)
    MVSK[["HPDU"]] <- lapply(Mdl.PredY, applyFun, mar = c(2, 3),
                             fun = quantile, probs = 0.975)

    if(!is.null(PredDens))
    {
        ## Residuals and the uncertainty
        MCMC.residual <- mapply(FUN = function(x, y) {matrix(x, nrow(y), ncol(y), byrow = TRUE)-y},
                                x = Mdl.Y.testing[names(Mdl.PredY)],
                                y = lapply(Mdl.PredY, function(x) x[, ,"mean"]),
                                SIMPLIFY = FALSE)

        RESID <- list()
        RESID[["mean"]] <- lapply(MCMC.residual, applyFun, mar = 2, fun = mean)
        RESID[["var"]] <- lapply(MCMC.residual, applyFun, mar = 2, fun = var)
        RESID[["HPDL"]] <- lapply(MCMC.residual, applyFun, mar  = 2,
                                  fun = quantile, probs = 0.025)
        RESID[["HPDU"]] <- lapply(MCMC.residual, applyFun, mar = 2,
                                  fun = quantile, probs = 0.975)
    }
    else
    {
        Mdl.logPredDens <- NA
        RESID <- NA
    }

    out <- list(Mdl.logPredDens = Mdl.logPredDens,
                Mdl.PredMVSK = MVSK,
                Mdl.RredRESID = RESID
                # Mdl.PredY = Mdl.PredY # too big
                )
    return(out)
}
