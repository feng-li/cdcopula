MargiModelForeignEval <- function(Mdl.MargisNM, Mdl.MargisType, MargisForeignConfig, Mdl.Y)
{
    Mdl.X <- list()
    Mdl.ForeignFit <- list()

    for(iComp in 1:(length(Mdl.MargisNM)-1)) ## TODO: Parallelize marginal models.
    {
        cat("Evaluating foreign marginal model: ", Mdl.MargisNM[iComp], "...\n")
        if(tolower(Mdl.MargisType[iComp])  == "garch")
        {
            ## GARCH Model
            require("fGarch", quietly = TRUE)

            parArgs <- MargisForeignConfig[[iComp]]
            parArgs[["data"]] <- Mdl.Y[[Mdl.MargisNM[iComp]]]
            MargiModel.Fit.caller <- as.call(c(garchFit, parArgs))
            MargiModel.Fit <- eval(MargiModel.Fit.caller)

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()
            if(parArgs[["cond.dist"]] == "norm")
            { # Normal innovation
                Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(MargiModel.Fit@fitted)
                Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(MargiModel.Fit@sigma.t)
            }
            Mdl.ForeignFit[[Mdl.MargisNM[iComp]]] <- MargiModel.Fit
        }
        else if(tolower(Mdl.MargisType[iComp])  == "stochvol")
        {
            require("stochvol", quietly = TRUE)

            ## The Stochastic Volatility Model
            parArgs <- MargisForeignConfig[[iComp]]
            parArgs[["y"]] <- Mdl.Y[[Mdl.MargisNM[iComp]]]
            MargiModel.Fit.caller <- as.call(c(svsample, parArgs))
            MargiModel.Fit <- eval(MargiModel.Fit.caller)

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()
            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(0, dim(Mdl.Y[[Mdl.MargisNM[iComp]]])[1], 1)
            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(sqrt(apply(exp(MargiModel.Fit$latent), 2, mean)))

            Mdl.ForeignFit[[Mdl.MargisNM[iComp]]] <- MargiModel.Fit

        }
        else
        {
            stop("This foreign marginal model is not implemented!")
        }
    }

    out <- list(Mdl.X = Mdl.X, Mdl.ForeignFit = Mdl.ForeignFit)
    return(out)
}
