parCplMCMC <- function(MCMC.beta, Mdl.X, Mdl.parLink, MCMC.Update, MCMC.sampleIdx)
{
    MCMC.par <- list()
    Mdl.X.training <- Mdl.X

    nMCMCSample <- length(MCMC.sampleIdx)
    nTraining <- nrow(Mdl.X.training[[1]][[1]])
    nTrainingNames <- rownames(Mdl.X.training[[1]][[1]])

    Mdl.MargisNM <- names(MCMC.Update)


    ## Regenerate MCMC.par
    for(CompCaller in names(MCMC.Update))
    {
        for(parCaller in names(MCMC.Update[[CompCaller]]))
        {
            ncolX.ij <- ncol(Mdl.X.training[[CompCaller]][[parCaller]])
            nPar.ij <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]
            namesX.ij <- rep(colnames(Mdl.X.training[[CompCaller]][[parCaller]]), nPar.ij)


            if((CompCaller %in% Mdl.MargisNM[length(Mdl.MargisNM)]) & nPar.ij != 1)
            {
                nDim <- length(Mdl.MargisNM)-1
                namesParFull.ij <- matrix(paste(matrix(1:nDim, nDim, nDim),
                                                matrix(1:nDim, nDim, nDim,
                                                       byrow = TRUE), sep = "."), nDim)
                namesPar.ij <- namesParFull.ij[lower.tri(namesParFull.ij,
                                                         diag = FALSE)]
            }
            else
            {
                namesPar.ij <- "1.1"
            }


            MCMC.par[[CompCaller]][[parCaller]] <- array(NA, c(nMCMCSample, nTraining, nPar.ij),
                                                         dimnames = list(NULL, nTrainingNames, namesPar.ij))
        }
    }

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

    for(iMCMC.sampleIdx in 1:nMCMCSample)
    {
        Mdl.beta.curr <- rapply(object=MCMC.beta, f = subsetFun4beta,
                                idx = MCMC.sampleIdx[iMCMC.sampleIdx],
                                how = "replace")

        Mdl.par.curr <- parCplMeanFun(Mdl.X = Mdl.X.training,
                                      Mdl.parLink = Mdl.parLink,
                                      Mdl.beta = Mdl.beta.curr,
                                      parUpdate = MCMC.Update)

        for(CompCaller in names(MCMC.Update))
        {
            for(parCaller in names(MCMC.Update[[CompCaller]]))
            {
                MCMC.par[[CompCaller]][[parCaller]][iMCMC.sampleIdx, ,] <- Mdl.par.curr[[CompCaller]][[parCaller]]
            }
        }
    }


    ## Summary of MCMC.par
    subFun3 <- function(obj, fun, dim, ...)
    {
        if(any(is.na(obj)))
        {
            out <- NA
        }
        else
        {
            out <- apply(obj, dim, fun, ...)
        }
        return(out)
    }

    out <- list()
    out[["mean"]] <- rapply(MCMC.par, subFun3, how = "replace", fun = mean, dim = 3)
    out[["median"]] <- rapply(MCMC.par, subFun3, how = "replace", fun = median, dim = 3,)
    out[["sd"]] <- rapply(MCMC.par, subFun3, how = "replace",   fun = sd, dim = 3)
    out[["ts.mean"]] <- rapply(MCMC.par, subFun3, how = "replace", fun = mean, dim = c(2, 3))
    out[["ts.median"]] <- rapply(MCMC.par, subFun3, how = "replace", fun = median,
                                 dim = c(2, 3))
    out[["ts.sd"]] <- rapply(MCMC.par, subFun3, how = "replace", fun = sd, dim = c(2, 3))
    out[["ts.hpd95"]] <- rapply(MCMC.par, subFun3, how = "replace", fun = quantile,
                                dim = c(2, 3), probs = c(0.025, 0.975))


    return(out)
}
