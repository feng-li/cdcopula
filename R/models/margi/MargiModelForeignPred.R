MargiModelForeignPred <- function(Mdl.MargisNM, Mdl.MargisType, Mdl.MargisForeignFitted, Mdl.Y)
{
    Mdl.X <- list()
    Mdl.ForeignPred <- list()
    for(iComp in 1:length(Mdl.MargisForeignFitted)) ## TODO: Parallelize marginal models.
    {
        if(tolower(Mdl.MargisType[iComp]) ==  "garch")
        {
            require("fGarch", quietly = TRUE)

            MargiModel.Pred.caller <- as.call(c(predict, Mdl.MargisForeignFitted[[iComp]],
                                                n.ahead = length(Mdl.Y[[iComp]])))

            MargiModel.Pred <- eval(MargiModel.Pred.caller)

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()

            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(MargiModel.Pred[["meanForecast"]])
            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(MargiModel.Pred[["standardDeviation"]])

            Mdl.ForeignPred[[Mdl.MargisNM[iComp]]] <- MargiModel.Pred
        }
        else if(tolower(Mdl.MargisType[iComp]) ==  "stochvol")
        {
            ## MargiModel.Pred.caller <- as.call(c(predict.svdraws,
            ##                                     Mdl.MargisForeignFitted[[iComp]],
            ##                                     steps = length(Mdl.Y[[iComp]])))
            ## MargiModel.Pred <- eval(MargiModel.Pred.caller)
            require("stochvol", quietly = TRUE)

            MargiModel.Pred <- predict.svdraws(object = Mdl.MargisForeignFitted[[iComp]],
                                               steps = length(Mdl.Y[[iComp]]))

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()

            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(apply(MargiModel.Pred, 2, mean))
            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(apply(MargiModel.Pred, 2, sd))

            Mdl.ForeignPred[[Mdl.MargisNM[iComp]]] <- MargiModel.Pred

        }
        else
        {
            stop("No such foreign marginal model!")
        }
    }

    out <- list(Mdl.X = Mdl.X, Mdl.ForeignPred = Mdl.ForeignPred)
    return(out)
}
