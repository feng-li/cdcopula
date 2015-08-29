logCplPredict <- function(MargisType, Mdl.par)
{
  MargisNM <- names(MargisType[-length(MargisType)])

  out <- list()
  for(iComp in MargisNM)
  {
    if("splitt" %in% tolower(MargisType[iComp]))
    {

      mu <- Mdl.par[[iComp]][["mu"]]
      df <- Mdl.par[[iComp]][["df"]]
      phi <- Mdl.par[[iComp]][["phi"]]
      lmd <- Mdl.par[[iComp]][["lmd"]]

      mean <- splitt.mean(mu, df, phi, lmd)
      var <- splitt.var(df, phi, lmd)
      skewness <- splitt.skewness(df, phi, lmd)
      kurtosis <- splitt.kurtosis(df, phi, lmd)
    }
    else
    {
      stop("No such marginal component!")
    }

    out.iComp <- cbind(mean, var, skewness, kurtosis)
    colnames(out.iComp) <- c("mean", "variance", "skewness", "kurtosis")

    out[[iComp]] <- out.iComp
  }

  return(out)
}
