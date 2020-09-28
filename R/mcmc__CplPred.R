#' Simulate random variables for the predictive density based on fitted MCMC results
#'
#' This can be use for computing the VaR and Expected Shortfall.
#' @param CplFitted  Fitted results
#' @param Mdl.Idx.testing  Testing index
#' @return Array nMCMC-by-2-nTesting
#' @references NA
#' @author Feng Li, Central University of Finance and Economics.
#' @export
rCplPred <- function(CplFitted, Mdl.Idx.testing)
{
    Mdl.MargisType <- CplFitted[["Mdl.ConfigEnv"]][["Mdl.MargisType"]]
    Mdl.MargisNM <- names(Mdl.MargisType)
    CplNM <- Mdl.MargisType[length(Mdl.MargisType)]

    ForeignModelSpec <- CplFitted[["ForeignModelSpec"]]
    if(length(ForeignModelSpec)>0 &&  class(ForeignModelSpec) != "logical")
    {   ## Special case when a foreign multivariate model is introduced. Use the model's
        ## predict method.

        ForeignModelType <- Mdl.MargisType[length(Mdl.MargisType)]
        Mdl.ForeignFitted <- CplFitted[["Mdl.ForeignFitted"]]
        Mdl.Y <- CplFitted[["Mdl.ConfigEnv"]][["Mdl.Y"]]

        ## Check if Mdl.Y.testing is in-sample or out-of-sample
        Mdl.Idx.training <- CplFitted[["Mdl.Idx.training"]]
        Split.Idx <- (Mdl.Idx.testing %in% Mdl.Idx.training)

        Mdl.Idx.inSample <- Mdl.Idx.testing[Split.Idx]
        Mdl.Idx.outSample <- Mdl.Idx.testing[!Split.Idx]


        if(length(Mdl.Idx.outSample)>0)
        {
            Mdl.Y.outSample <- rapply(object=Mdl.Y,
                                      f = function(x, idx){x[idx,,drop = FALSE]},
                                      idx = Mdl.Idx.outSample, how = "replace")
        }
        else
        {
            Mdl.Y.outSample <- NULL
        }

        Mdl.Pred <- ModelForeignPred(model = ForeignModelType,
                                     fitted.model = Mdl.ForeignFitted,
                                     data.pred = Mdl.Y.outSample)
        ## Simulate p(y_new|y_old)
        mean.out <- NULL
        var.out <- NULL
        if(length(Mdl.Idx.outSample)>0)
        {
            mean.out <- Mdl.Pred[["MVSK"]][["mean"]] # n-by-2 matrix
            var.out <- Mdl.Pred[["MVSK"]][["var"]] # 2-by-2-n array
        }

        ## Simulate p(y|Model)
        mean.in <- NULL
        var.in <- NULL
        if(length(Mdl.Idx.inSample)>0)
        {
            mean.in <- Mdl.Pred[["MVSK.fitted"]][["mean"]]
            var.in <- Mdl.Pred[["MVSK.fitted"]][["var"]]
        }

        ## Merge all parameters
        means <- rbind(mean.in, mean.out)
        vars <- array(c(var.in, var.out), dim = c(2, 2, length(Mdl.Idx.testing)))

        MCMC.sampleSize <- 1000 # arbitrary setting
        Mdl.Ysim <- array(NA, dim = c(MCMC.sampleSize,
                                      length(Mdl.MargisType)-1,
                                      length(Mdl.Idx.testing)),
                          dimnames = list(NULL, Mdl.MargisNM[-length(Mdl.MargisNM)], NULL))

        for(iTesting in 1:length(Mdl.Idx.testing))
        {
            Mdl.Ysim[, , iTesting] <- rmvnorm(MCMC.sampleSize, means[i, ], vars[, , i])
        }
    }
    else
    {
        ## The native model
        MCMC.beta <- CplFitted[["MCMC.beta"]]

        Mdl.parLink <- CplFitted[["Mdl.ConfigEnv"]][["Mdl.parLink"]]
        Mdl.X <- CplFitted[["Mdl.ConfigEnv"]][["Mdl.X"]]
        MCMC.Update <- CplFitted[["Mdl.ConfigEnv"]][["MCMC.Update"]]
        MCMC.nIter <- CplFitted[["Mdl.ConfigEnv"]][["MCMC.nIter"]]
        MCMC.burninProp <- CplFitted[["Mdl.ConfigEnv"]][["MCMC.burninProp"]]
        MCMC.sampleProp <- CplFitted[["Mdl.ConfigEnv"]][["MCMC.sampleProp"]]
        Mdl.MargisForeignFitted = CplFitted[["Mdl.MargisForeignFitted"]]

        ## MCMC.sampleProp <- 0.01

        n.burn <- round(MCMC.nIter*MCMC.burninProp)
        MCMC.sampleIdx <- round(seq(n.burn+1, MCMC.nIter,
                                    length.out = round((MCMC.nIter-n.burn)*MCMC.sampleProp)))

        ## Mdl.Ysim <- array(NA, dim = c(length(MCMC.sampleIdx),
        ##                               length(Mdl.MargisType)-1,
        ##                               length(Mdl.Idx.testing)),
        ##                   dimnames = list(NULL, Mdl.MargisNM[-length(Mdl.MargisNM)], NULL))

        ## TODO: Add foreign marginal model

        Mdl.X.testing = make_X.testing(Mdl.Idx.testing = Mdl.Idx.testing,
                                       Mdl.X = Mdl.X,
                                       Mdl.ConfigEnv = Mdl.ConfigEnv,
                                       Mdl.MargisForeignFitted =  Mdl.MargisForeignFitted)
        MCMC.par <- parCplMCMC(MCMC.beta = MCMC.beta,
                               Mdl.X = Mdl.X.testing,
                               Mdl.parLink = Mdl.parLink,
                               MCMC.Update = MCMC.Update,
                               MCMC.sampleIdx = MCMC.sampleIdx)

        CompCpl <- Mdl.MargisNM[length(Mdl.MargisType)] # only copula component

        ## Marginal models are also simulated
        ## CompUpLst <- unlist(lapply(MCMC.Update, function(x) any(unlist(x) == TRUE)))

        Mdl.YsimFun = function(iTesting)
        {
            Mdl.Ysim.iTesting <- array(
                NA, dim = c(length(MCMC.sampleIdx),
                            length(Mdl.MargisType)-1),
                dimnames = list(NULL, Mdl.MargisNM[-length(Mdl.MargisNM)]))

            for(iMCMC in 1:length(MCMC.sampleIdx))
            {
                subfun <- function(x, iMCMC, iTesting)
                {
                    if((dim(x)[1] == 1 && length(iMCMC)>1) ||
                       (dim(x)[1] == 1 && length(iMCMC) == 1 && iMCMC != 1))
                    {
                        sampleIdx <- 1
                    }
                    else
                    {
                        sampleIdx <- iMCMC
                    }
                    x.use <- x[sampleIdx, iTesting, , drop = FALSE]
                    return(x.use)
                }

                Mdl.par <- rapply(MCMC.par, subfun, iMCMC = iMCMC,
                                  iTesting = iTesting, how = "replace")

                parCplStd <- parCplRep2Std(CplNM = CplNM, parCplRep = Mdl.par[[CompCpl]])

                Mdl.u <- rCpl(n = 1, CplNM = CplNM, parCpl = parCplStd)[["u"]]

                colnames(Mdl.u) <- Mdl.MargisNM[-length(Mdl.MargisNM)]

                ## cat("iTesting, iMCMC:", iTesting, iMCMC,
                ##     "nTesting, nMCMC:", length(Mdl.Idx.testing), length(MCMC.sampleIdx),
                ##     "\n")
                for(iComp in setdiff(Mdl.MargisNM, CplNM))
                {
                    ## if(CompUpLst[iComp])
                    ## {
                    Mdl.Ysim.iTesting[iMCMC, iComp] <- MargiModelInv(
                        u = Mdl.u[, iComp],
                        par = lapply(Mdl.par[[iComp]], rbind),
                        type = Mdl.MargisType[[iComp]])
                    ## }
                }
            }
            return(Mdl.Ysim.iTesting)
        }

        Mdl.Idx.testingMat = matrix(1:length(Mdl.Idx.testing)) # nPred-by-1

        cl <- getDefaultCluster()
        ## ce <- clusterEvalQ(cl,{
        ##     source("~/code/flutils/R/systools/sourceDir.R")
        ##     sourceDir("~/code/cdcopula/R", "~/code/flutils/R", recursive = TRUE)
        ## })
        Mdl.YsimMat = parApply(cl = cl, X = Mdl.Idx.testingMat, MARGIN = 1,
                               FUN = Mdl.YsimFun) # nMCMC*nComp-by-nPred
        Mdl.Ysim = array(Mdl.YsimMat, dim = c(length(MCMC.sampleIdx),
                                              length(Mdl.MargisType)-1,
                                              length(Mdl.Idx.testing)),
                          dimnames = list(NULL, Mdl.MargisNM[-length(Mdl.MargisNM)], NULL))

    }
    return(Mdl.Ysim)
}



###----------------------------------------------------------------------------
### Calculate the predictive density p(y_pred|y_old)
###----------------------------------------------------------------------------
#' @export
dCplPred <- function(Out.Fitted)
{

    ## Generate mesh grid for fake Y.
    Y.Fake = do.call(cbind, lapply(Mdl.Y, quantile, probs = seq(0, 1, 0.01))) # nGrid-by-2
    Y.FakeMesh = mesh.grid(Y.Fake[, 1], Y.Fake[, 2]) # nGrid^2-by-2
    colnames(Y.FakeMesh) = names(Mdl.Y)


    iCross <- 1
    iLPDS = "joint"
    Mdl.Idx.testing <- 1

    CplFitted[["MCMC.sampleProp"]] <- 0.005
    CplFitted[["LPDS"]] <- "joint"

    nGrid <- nrow(Y.FakeMesh)
    out <- matrix(NA, nGrid, 1)
    for(i in 1:nGrid)
    {
        cat(i, ", ")
        Y.FakeMeshLst = array2list(Y.FakeMesh[i, , drop = FALSE], 2)
        Y.Fake.Pred = logDensPred(CplFitted = CplFitted,
                                  Mdl.Idx.testing = Mdl.Idx.testing,
                                  Mdl.Y.testing = Y.FakeMeshLst)
        out[i] <- logPredDensScore(
            list(Y.Fake.Pred[["Mdl.logPredDens"]][, iLPDS, drop = FALSE]))[1]
    }
}
