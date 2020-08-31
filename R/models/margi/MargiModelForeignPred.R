#' @export
MargiModelForeignPred <- function(Mdl.MargisNM, Mdl.MargisType, Mdl.MargisForeignFitted, n.ahead)
{
    Mdl.X <- list()
    Mdl.ForeignPred <- list()
    for(iComp in 1:length(Mdl.MargisForeignFitted)) ## TODO: Parallelize marginal models.
    {
        if(tolower(Mdl.MargisType[iComp]) ==  "garch")
        {
            require("fGarch", quietly = TRUE)
            MargiModel.Pred.caller <- as.call(c(predict, Mdl.MargisForeignFitted[[iComp]],
                                                n.ahead = n.ahead))

            MargiModel.Pred <- eval(MargiModel.Pred.caller)

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()

            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(MargiModel.Pred[["meanForecast"]])
            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(MargiModel.Pred[["standardDeviation"]])

            Mdl.ForeignPred[[Mdl.MargisNM[iComp]]] <- MargiModel.Pred
        }
        else if(tolower(Mdl.MargisType[iComp]) ==  "stochvol")
        {
            ## MargiModel.Pred.caller <- as.call(c(predict.svdraws,
            ##                                     Mdl.MargisForeignFitted[[iComp]][[iComp]],
            ##                                     steps = length(Mdl.Y[[iComp]])))
            ## MargiModel.Pred <- eval(MargiModel.Pred.caller)
            require("stochvol", quietly = TRUE)

            MargiModel.Pred <- predict.svdraws(object = Mdl.MargisForeignFitted[[iComp]], steps = n.ahead)

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()

            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <- matrix(apply(MargiModel.Pred, 2, mean))
            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(apply(MargiModel.Pred, 2, sd))

            Mdl.ForeignPred[[Mdl.MargisNM[iComp]]] <- MargiModel.Pred

        }
        else if(tolower(Mdl.MargisType[iComp]) ==  "admit")
        {
            require("AdMit", quietly = TRUE)
            stop("Not yet ready!")
        }
        else if(tolower(Mdl.MargisType[iComp]) ==  "teigen")
        {
            require("teigen", quietly = TRUE)

            ## n.ahead <- length(Mdl.Y[[Mdl.MargisNM[iComp]]])
            nGroups <-  Mdl.MargisForeignFitted[[iComp]][["G"]]

            Mdl.X[[Mdl.MargisNM[iComp]]] <- list()
            Mdl.X[[Mdl.MargisNM[iComp]]][["mu"]] <-  matrix(
                Mdl.MargisForeignFitted[[iComp]][["parameters"]][["mean"]],
                n.ahead, nGroups, byrow = TRUE)

            Mdl.X[[Mdl.MargisNM[iComp]]][["df"]] <-  matrix(
                Mdl.MargisForeignFitted[[iComp]][["parameters"]][["df"]],
                n.ahead, nGroups, byrow = TRUE)

            Mdl.X[[Mdl.MargisNM[iComp]]][["phi"]] <- matrix(
                Mdl.MargisForeignFitted[[iComp]][["parameters"]][["sigma"]],
                n.ahead, nGroups, byrow = TRUE)

            Mdl.X[[Mdl.MargisNM[iComp]]][["weights"]] <-  matrix(
                Mdl.MargisForeignFitted[[iComp]][["parameters"]][["pig"]],
                n.ahead, nGroups, byrow = TRUE)

            Mdl.ForeignPred[[Mdl.MargisNM[iComp]]] <- NA

        }

        else
        {
            stop("No such foreign marginal model!")
        }
    }

    out <- list(Mdl.X = Mdl.X, Mdl.ForeignPred = Mdl.ForeignPred)
    return(out)
}
