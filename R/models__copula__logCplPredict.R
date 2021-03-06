#' @export
logCplPredict <- function(Mdl.MargisType, Mdl.par)
{
    Mdl.MargisNM <- names(Mdl.MargisType[-length(Mdl.MargisType)])

    out <- list()
    for(iComp in Mdl.MargisNM)
    {
        type <- tolower(Mdl.MargisType[iComp])

        out.iComp <- list()

        ## Convert model types into standard version
        if(type  %in% c("garch", "stochvol"))
        {
            typeStd <- "gaussian"
        }
        else
        {
            typeStd <- type
        }
###----------------------------------------------------------------------------
### Standard model type
###----------------------------------------------------------------------------
        if(typeStd == "gaussian")
        {
            mu <- Mdl.par[[iComp]][["mu"]]
            phi <- Mdl.par[[iComp]][["phi"]]

            out.iComp[["mean"]] <- mu
            out.iComp[["var"]] <- phi^2

        }
        else if(typeStd == "splitt")
        {
            mu <- Mdl.par[[iComp]][["mu"]]
            df <- Mdl.par[[iComp]][["df"]]
            phi <- Mdl.par[[iComp]][["phi"]]
            lmd <- Mdl.par[[iComp]][["lmd"]]

            out.iComp[["mean"]] <- splitt.mean(mu, df, phi, lmd)
            out.iComp[["var"]] <- splitt.var(df, phi, lmd)
            out.iComp[["skewness"]] <- splitt.skewness(df, phi, lmd)
            out.iComp[["kurtosis"]] <- splitt.kurtosis(df, phi, lmd)
        }
        else if(typeStd == "teigen")
        {

            mu <- Mdl.par[[iComp]][["mu"]]
            df <- Mdl.par[[iComp]][["df"]]
            phi <- Mdl.par[[iComp]][["phi"]]
            lmd <- array(1, dim(mu))
            weights <- Mdl.par[[iComp]][["weights"]]

            out.iComp[["mean"]] <- NA
            out.iComp[["var"]] <- NA
            ## out.iComp[["phi"]] <- NA
            ## out.iComp[["lmd"]] <- NA
            ## out.iComp[["weights"]] <- weights


        }
        else if(typeStd == "poisson")
        {
            mu <- Mdl.par[[iComp]][["mu"]]
            out.iComp[["mean"]] <- mu
        }
        else
        {
            stop("No such marginal component!")
        }

        out[[iComp]] <- do.call(cbind, out.iComp)
    }

    return(out)
}
