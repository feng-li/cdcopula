logDens <- function(Mdl.MargisType, Mdl.Y,  Mdl.u, Mdl.d, Mdl.par,
                    parUpdate, MCMC.UpdateStrategy)
{
###----------------------------------------------------------------------------
### UPDATE MARGINAL PDF AND/OR CDF THE MARGINAL U AND ONLY UPDATED IF THE CORRESPONDING
### PARAMETERS ARE UPDATED.
###----------------------------------------------------------------------------
    if(missing(Mdl.u))
    {
        Mdl.u <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.Y),
                        dimnames = list(NULL, names(Mdl.Y)))
    }

    if(missing(Mdl.d))
    {
        Mdl.d <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.par),
                        dimnames = list(NULL, names(Mdl.par)))
    }

    CompNM <- names(Mdl.par)
    ## Allocate the output structure: Margins + Copula (NA)
    ## Mdl.logLik <- cbind(Mdl.d, NA)
    ## colnames(Mdl.logLik) <- CompNM

    CplNM <- Mdl.MargisType[length(Mdl.MargisType)]
    CompUpNM <- unlist(lapply(parUpdate, function(x) any(unlist(x) == TRUE)))
    Mdl.MargisNM <- CompNM[CompNM  != CplNM]
    MargisUpNM <- CompNM[(CompNM  != CplNM) & CompUpNM]

    ## THE COPULA LOG LIKELIHOOD
    if(tolower(MCMC.UpdateStrategy) == "joint")
    {
        evalCpl <- TRUE
        Mdl.PostComp <- lapply(parUpdate, function(x) TRUE)
    }
    else if(tolower(MCMC.UpdateStrategy) == "twostage" ||
            tolower(MCMC.UpdateStrategy) == "margin")
    {
        if(length(MargisUpNM) == 0)
        {
            ## Stage two of the two stage approach, only updating copula density
            evalCpl <- TRUE
            Mdl.PostComp <- lapply(parUpdate, function(x) FALSE)
            Mdl.PostComp[[CplNM]] <- TRUE
        }
        else
        {
            ## Only updating marginal density
            evalCpl <- FALSE
            Mdl.PostComp <- lapply(parUpdate, function(x) any(unlist(x) == TRUE))
            Mdl.PostComp[[CplNM]] <- FALSE
        }
    }
    else
    {
        stop(paste("MCMC update strategy:", MCMC.UpdateStrategy,
                   "not implemented!"))
    }

    ## THE MARGINAL LIKELIHOODS
    densCaller <- list()
    for(CompCaller in Mdl.MargisNM)
    {
        if(CompCaller %in% MargisUpNM)
        {
            if(tolower(MCMC.UpdateStrategy) == "joint")
            {
                densCaller[[CompCaller]] <- c("u", "d")
            }
            else if(tolower(MCMC.UpdateStrategy) == "twostage")
            {
                ## Stage two of the two stage approach
                densCaller[[CompCaller]] <- c("u", "d")
            }
            else if(tolower(MCMC.UpdateStrategy) == "margin")
            {
                densCaller[[CompCaller]] <- c(NA, "d")
            }
            else
            {
                stop(paste("MCMC update strategy:", MCMC.UpdateStrategy,
                           "not implemented!"))
            }
        }
        else
        {
            if(evalCpl == TRUE && any(is.na(Mdl.u)))
            {
                densCaller[[CompCaller]] <- c("u", "d")
            }
        }
    }

###----------------------------------------------------------------------------
### UPDATING THE LIKELIHOOD
###----------------------------------------------------------------------------
    ## Updating the marginal models TODO: parallel version
    for(CompCaller in names(densCaller))
    {
        ## if(CompCaller == "BABA") browser()

        Mdl.ud <- MargiModel(y = Mdl.Y[[CompCaller]],
                             type = Mdl.MargisType[which(Mdl.MargisNM == CompCaller)],
                             par = Mdl.par[[CompCaller]],
                             densCaller = densCaller[[CompCaller]])

        if("u" %in% densCaller[[CompCaller]])
        {
            Mdl.u[,  CompCaller] <- Mdl.ud[["u"]]
        }
        if("d" %in% densCaller[[CompCaller]])
        {
            Mdl.d[,  CompCaller] <- Mdl.ud[["d"]]
        }
    }

    ## Updating the copula model if required
    R_CPL_NPARALLEL <- as.numeric(Sys.getenv("R_CPL_NPARALLEL"))
    if(evalCpl == TRUE)
    {

        if(!is.na(R_CPL_NPARALLEL) && R_CPL_NPARALLEL > 1)
        {
            logCplLikFUN.NM <- "logCplLikParallel"
        }
        else
        {
            logCplLikFUN.NM <- "logCplLik"
        }


        Mdl.logLikCpl.caller <- call(logCplLikFUN.NM,
                                     u = Mdl.u,
                                     CplNM = CplNM,
                                     parCplRep = Mdl.par[[CplNM]],
                                     sum = FALSE)
        Mdl.logLikCpl <- eval(Mdl.logLikCpl.caller) # n-by-1

        Mdl.d[, CplNM] <- Mdl.logLikCpl
    }
    ## browser()
    ## if(any(!is.finite(Mdl.u))) browser()
    out <- list(Mdl.d = Mdl.d, Mdl.u = Mdl.u,
                Mdl.PostComp = Mdl.PostComp)
    return(out)
}
