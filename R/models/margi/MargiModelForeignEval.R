MargiModelForeignEval <- function(MargisNM, MargisType, MargisForeignConfig, Mdl.Y)
{
  Mdl.X <- list()
  Mdl.ForeignFit <- list()

  for(iComp in 1:(length(MargisNM)-1)) ## TODO: Parallelize marginal models.
  {
    cat("Evaluating foreign marginal model: ", MargisNM[iComp], "...\n")
    if(tolower(MargisType[iComp])  == "garch")
    {
      require("fGarch")

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
    else
    {
      stop("This foreign model is not implemented!")
    }
  }

  out <- list(Mdl.X = Mdl.X, Mdl.ForeignFit = Mdl.ForeignFit)
  return(out)
}
