#' @export
parCplMCMC <- function(MCMC.beta, Mdl.X, Mdl.parLink, MCMC.Update, MCMC.sampleIdx)
{
    MCMC.nIter <- length(MCMC.sampleIdx)
    nObs <- nrow(Mdl.X[[1]][[1]])
    nObsNames <- rownames(Mdl.X[[1]][[1]])

    Mdl.MargisNM <- names(MCMC.Update)

    ## Regenerate MCMC.par
    MCMC.par <- list()
    for(CompCaller in names(MCMC.Update))
    {
        for(parCaller in names(MCMC.Update[[CompCaller]]))
        {
            ncolX.ij <- ncol(Mdl.X[[CompCaller]][[parCaller]])
            nPar.ij <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]
            namesX.ij <- rep(colnames(Mdl.X[[CompCaller]][[parCaller]]), nPar.ij)

            nMCMCSample <- ifelse(MCMC.Update[[CompCaller]][[parCaller]], MCMC.nIter, 1)

            if((CompCaller %in% Mdl.MargisNM[length(Mdl.MargisNM)]) & nPar.ij != 1)
            {
                nDim <- length(Mdl.MargisNM)-1
                namesParFull.ij <- matrix(paste(matrix(1:nDim, nDim, nDim),
                                                matrix(1:nDim, nDim, nDim,
                                                       byrow = TRUE), sep = "."), nDim)
                namesPar.ij <- namesParFull.ij[lower.tri(namesParFull.ij, diag = FALSE)]
            }
            else if((CompCaller %in% Mdl.MargisNM[-length(Mdl.MargisNM)]) & nPar.ij != 1)
            {
                namesPar.ij <- 1:nPar.ij
            }
            else
            {
                namesPar.ij <- "1.1"
            }

            MCMC.par[[CompCaller]][[parCaller]] <- array(NA, c(nMCMCSample, nObs, nPar.ij),
                                                         dimnames = list(NULL, nObsNames, namesPar.ij))
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



    for(iMCMC.sampleIdx in 1:MCMC.nIter)
    {
        Mdl.beta.curr <- rapply(object=MCMC.beta, f = subsetFun4beta,
                                idx = MCMC.sampleIdx[iMCMC.sampleIdx],
                                how = "replace")

        Mdl.par.curr <- parCplMeanFun(Mdl.X = Mdl.X,
                                      Mdl.parLink = Mdl.parLink,
                                      Mdl.beta = Mdl.beta.curr,
                                      parUpdate = MCMC.Update)

        for(CompCaller in names(MCMC.Update))
        {
            for(parCaller in names(MCMC.Update[[CompCaller]]))
            {
                if(iMCMC.sampleIdx == 1 ||
                   MCMC.Update[[CompCaller]][[parCaller]] == TRUE)
                {  ## MCMC.par is updated at each MCMC draw if MCMC.Update is TRUE for
                    ## that parameter. Otherwise, only use the first MCMC draw.
                    MCMC.par[[CompCaller]][[parCaller]][iMCMC.sampleIdx, ,] <-
                        Mdl.par.curr[[CompCaller]][[parCaller]]
                }

            }
        }
    }
    return(MCMC.par)
}

#' @export
parCplMCMCSummary <- function(MCMC.par, MCMC.sampleIdx)
{
   subFun3 <- function(obj, fun, dim, idx, ...)
   {
       if((dim(obj)[1] == 1 && length(idx)>1) ||
          (dim(obj)[1] == 1 && length(idx) == 1 && idx != 1))
       {
           sampleIdx <- 1
       }
       else
       {
           sampleIdx <- idx
       }
       obj.use <- obj[sampleIdx, , , drop = FALSE]

       ## Check missing values
       if(any(is.na(obj.use)))
       {
           out <- NA
       }
       else
       {
           out <- apply(obj.use, dim, fun, ...)
       }



       return(out)
    }
   out <- list()
   out[["mean"]] <- rapply(MCMC.par, subFun3, how = "replace",
                           fun = mean, dim = 3, idx = MCMC.sampleIdx)
   out[["median"]] <- rapply(MCMC.par, subFun3, how = "replace",
                             fun = median, dim = 3,idx = MCMC.sampleIdx)
   out[["sd"]] <- rapply(MCMC.par, subFun3, how = "replace",
                         fun = sd, dim = 3, idx = MCMC.sampleIdx)
   out[["ts.mean"]] <- rapply(MCMC.par, subFun3, how = "replace",
                              fun = mean, dim = c(2, 3), idx = MCMC.sampleIdx)
   out[["ts.median"]] <- rapply(MCMC.par, subFun3, how = "replace",
                                fun = median, dim = c(2, 3), idx = MCMC.sampleIdx)
   out[["ts.sd"]] <- rapply(MCMC.par, subFun3, how = "replace",
                            fun = sd, dim = c(2, 3), idx = MCMC.sampleIdx)
   out[["ts.hpd95"]] <- rapply(MCMC.par, subFun3, how = "replace",
                               fun = quantile,idx = MCMC.sampleIdx,
                               dim = c(2, 3), probs = c(0.025, 0.975))
   return(out)
}


#' @export
parCplMCMCSummary4Tau <- function(MCMC.par, MCMC.sampleIdx)
{
    MargisNM <- names(MCMC.par)
    CompCpl <- MargisNM[length(MargisNM)]
    CplNM0 <- strsplit(CplNM, split = "_")[[1]][1]

    parCpl <- MCMC.par[[CompCpl]]

    ## TODO: change the below to parCplRep2Std()

    if(tolower(CplNM0)  == "bb7")
    {
        lambdaL <- MCMC.par[[CompCpl]][["lambdaL"]]
        lambdaU <- MCMC.par[[CompCpl]][["lambdaU"]]

        delta <- (-log(2)/log(lambdaL))
        theta <- log(2)/log(2-lambdaU)

        tau <- kendalltau(CplNM = CplNM, parCpl = list(delta = delta, theta = theta))
    }


    parCpl4Tau <- list()
    parCpl4Tau[[CompCpl]] <- list("tau" = tau)
    ## browser()
    out <- parCplMCMCSummary(MCMC.par = parCpl4Tau,
                             MCMC.sampleIdx = MCMC.sampleIdx)

    return(out)
}
