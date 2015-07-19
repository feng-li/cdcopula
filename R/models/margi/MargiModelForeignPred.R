MargiModelForeignPred <- function(MargisNM, MargisType, Mdl.ForeignFit, Mdl.Y)
{
  Mdl.X <- list()
  Mdl.ForeignPred <- list()
  for(iComp in 1:length(Mdl.ForeignFit)) ## TODO: Parallelize marginal models.
  {
    if(tolower(MargisType[iComp]) %in% c("garch-normal", "garch-t"))
    {
      require("fGarch")

      MargiModel.Pred.caller <- as.call(c(predict, Mdl.ForeignFit[[iComp]],
                                          n.ahead = length(Mdl.Y[[iComp]])))

      MargiModel.Pred <- eval(MargiModel.Pred.caller)

      Mdl.X[[MargisNM[iComp]]] <- list()
      Mdl.X[[MargisNM[iComp]]][["mu"]] <- matrix(MargiModel.Pred$meanForecast)
      Mdl.X[[MargisNM[iComp]]][["phi"]] <- matrix(MargiModel.Pred$standardDeviation)

      Mdl.ForeignPred[[MargisNM[iComp]]] <- MargiModel.Pred
    }
  }

  out <- list(Mdl.X = Mdl.X, Mdl.ForeignPred = Mdl.ForeignPred)
  return(out)
}
