MargiModelForeignEval <- function(MargisNM, MargisType, MargiModelForeignConfig, Mdl.Y)
{
  out <- list()
  for(iComp in 1:(length(MargisNM)-1))
  {
    if(tolower(MargisType[iComp]) %in% c("garch-normal", "garch-t"))
    {
      require("fGarch")

      parArgs <- MargiModelForeignConfig[[iComp]]
      parArgs[["data"]] <- Mdl.Y[[MargisNM[iComp]]]
      MargiModel.Fit.caller <- as.call(c(garchFit, parArgs))
      MargiModel.Fit <- eval(MargiModel.Fit.caller)

      out[[MargisNM[iComp]]][["mu"]] <- matrix(MargiModel.Fit@fitted)
      out[[MargisNM[iComp]]][["phi"]] <- matrix(MargiModel.Fit@sigma.t)
    }
  }

  return(out)
}
