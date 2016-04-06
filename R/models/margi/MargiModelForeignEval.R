MargiModelForeignEval <- function(MargisNM, MargisType, MargisForeignConfig, Mdl.Y)
{
    Mdl.X <- list()
    Mdl.ForeignFit <- list()

    for(iComp in 1:(length(MargisNM)-1)) ## TODO: Parallelize marginal models.
    {
        cat("Evaluating foreign marginal model: ", MargisNM[iComp], "...\n")
        if(tolower(MargisType[iComp])  == "garch")
        {
            ## GARCH Model
            require("fGarch", quietly = TRUE)

            parArgs <- MargisForeignConfig[[iComp]]
            parArgs[["data"]] <- Mdl.Y[[MargisNM[iComp]]]
            MargiModel.Fit.caller <- as.call(c(garchFit, parArgs))
            MargiModel.Fit <- eval(MargiModel.Fit.caller)

            Mdl.X[[MargisNM[iComp]]] <- list()
            if(parArgs[["cond.dist"]] == "norm")
            { # Normal innovation
                Mdl.X[[MargisNM[iComp]]][["mu"]] <- matrix(MargiModel.Fit@fitted)
                Mdl.X[[MargisNM[iComp]]][["phi"]] <- matrix(MargiModel.Fit@sigma.t)
            }
            Mdl.ForeignFit[[MargisNM[iComp]]] <- MargiModel.Fit
        }
        else if(tolower(MargisType[iComp])  == "stochvol")
        {
            require("stochvol", quietly = TRUE)

            ## The Stochastic Volatility Model
            parArgs <- MargisForeignConfig[[iComp]]
            parArgs[["y"]] <- Mdl.Y[[MargisNM[iComp]]]
            MargiModel.Fit.caller <- as.call(c(svsample, parArgs))
            MargiModel.Fit <- eval(MargiModel.Fit.caller)

            Mdl.X[[MargisNM[iComp]]] <- list()
            Mdl.X[[MargisNM[iComp]]][["mu"]] <- matrix(0, dim(Mdl.Y[[MargisNM[iComp]]])[1], 1)
            Mdl.X[[MargisNM[iComp]]][["phi"]] <- matrix(sqrt(apply(exp(MargiModel.Fit$latent), 2, mean)))

            Mdl.ForeignFit[[MargisNM[iComp]]] <- MargiModel.Fit

        }
        else
        {
            stop("This foreign marginal model is not implemented!")
        }
    }

    out <- list(Mdl.X = Mdl.X, Mdl.ForeignFit = Mdl.ForeignFit)
    return(out)
}
