#' @export
MargiModelForeignEval <- function(Mdl.MargisNM, Mdl.MargisType, MargisForeignConfig, Mdl.Y)
{
    Mdl.X <- list()
    Mdl.MargisForeignFitted <- list()

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
            MargiModel.Fitted <- eval(MargiModel.Fit.caller)

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()
            if(parArgs[["cond.dist"]] == "norm")
            { # Normal innovation
                Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(MargiModel.Fitted@fitted)
                Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(MargiModel.Fitted@sigma.t)
            }
            Mdl.MargisForeignFitted[[Mdl.MargisNM[iComp]]] <- MargiModel.Fitted
        }
        else if(tolower(Mdl.MargisType[iComp])  == "stochvol")
        {
            require("stochvol", quietly = TRUE)

            ## The Stochastic Volatility Model
            parArgs <- MargisForeignConfig[[iComp]]
            parArgs[["y"]] <- Mdl.Y[[Mdl.MargisNM[iComp]]]
            MargiModel.Fit.caller <- as.call(c(svsample, parArgs))
            MargiModel.Fitted <- eval(MargiModel.Fit.caller)

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()
            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(0, dim(Mdl.Y[[Mdl.MargisNM[iComp]]])[1], 1)
            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(sqrt(apply(exp(MargiModel.Fitted$latent), 2, mean)))

            Mdl.MargisForeignFitted[[Mdl.MargisNM[iComp]]] <- MargiModel.Fitted

        }
        else if(tolower(Mdl.MargisType[iComp]) ==  "teigen")
        {
            require("teigen", quietly = TRUE)

            parArgs <- MargisForeignConfig[[iComp]]
            parArgs[["x"]] <- Mdl.Y[[Mdl.MargisNM[iComp]]]

            MargiModel.Fitted <- do.call(teigen, parArgs)

            nObs <- length(Mdl.Y[[Mdl.MargisNM[iComp]]])
            nGroups <- MargiModel.Fitted[["G"]]

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()
            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <-  matrix(
                MargiModel.Fitted[["parameters"]][["mean"]], nObs, nGroups, byrow = TRUE)
            Mdl.X[[Mdl.MargisNM[iComp]]][["df"]] <-  matrix(
                MargiModel.Fitted[["parameters"]][["df"]], nObs, nGroups, byrow = TRUE)
            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(
                MargiModel.Fitted[["parameters"]][["sigma"]], nObs, nGroups, byrow = TRUE)
            Mdl.X[[Mdl.MargisNM[iComp]]][["weights"]] <-  matrix(
                MargiModel.Fitted[["parameters"]][["pig"]], nObs, nGroups, byrow = TRUE)
            Mdl.MargisForeignFitted[[Mdl.MargisNM[iComp]]] <- MargiModel.Fitted

        }
        else
        {
            stop("This foreign marginal model is not implemented!")
        }
    }

    out <- list(Mdl.X = Mdl.X, Mdl.MargisForeignFitted = Mdl.MargisForeignFitted)
    return(out)
}
