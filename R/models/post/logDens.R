#' @export
logDens <- function(Mdl.MargisType, Mdl.Y,  Mdl.u, Mdl.d, Mdl.par,
                    parUpdate, MCMC.UpdateStrategy)
{
###----------------------------------------------------------------------------
### UPDATE MARGINAL PDF AND/OR CDF THE MARGINAL U AND ONLY UPDATED IF THE CORRESPONDING
### PARAMETERS ARE UPDATED.
###----------------------------------------------------------------------------
    if(missing(Mdl.u))
    {
        ## Mdl.u <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.Y),
        ##                 dimnames = list(NULL, names(Mdl.Y)))
        Mdl.u <- list()
    }

    if(missing(Mdl.d))
    {
        ## Mdl.d <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.par),
        ##                 dimnames = list(NULL, names(Mdl.par)))
        Mdl.d <- list()
    }

    # CompNM <- names(Mdl.par)
    ## Allocate the output structure: Margins + Copula (NA)
    ## Mdl.logLik <- cbind(Mdl.d, NA)
    ## colnames(Mdl.logLik) <- CompNM

    CplNM <- Mdl.MargisType[length(Mdl.MargisType)]
    Mdl.MargisNM <- names(Mdl.MargisType)

    CompCpl <- Mdl.MargisNM[length(Mdl.MargisType)] # only copula component
    CompMargis <- Mdl.MargisNM[-length(Mdl.MargisType)] # only marginal components

    CompUp <- unlist(lapply(parUpdate, function(x) any(unlist(x) == TRUE)))
    CompUpMargis <- CompUp[CompMargis]
    ## Mdl.MargisNM <- CompNM[-length(Mdl.MargisType)]
    ## MargisUpNM <- Mdl.MargisNM[CompUpNM[-length(Mdl.MargisType)]]

    ## THE COPULA LOG LIKELIHOOD
    if(tolower(MCMC.UpdateStrategy) == "joint")
    {
        ## Joint updating should update all marginal components and copula component.
        evalCpl <- TRUE
        Mdl.PostComp <- lapply(parUpdate, function(x) TRUE)
    }
    else if(tolower(MCMC.UpdateStrategy) == "twostage" ||
            tolower(MCMC.UpdateStrategy) == "margin")
    {
        if(length(CompUpMargis) == 0)
        {
            ## Stage two of the two-stage approach: only updating copula density
            evalCpl <- TRUE
            Mdl.PostComp <- lapply(parUpdate, function(x) FALSE)
            Mdl.PostComp[[CplComp]] <- TRUE
        }
        else
        {
            ## Stage one of the two-stage approach: Only updating marginal density
            evalCpl <- FALSE
            Mdl.PostComp <- lapply(parUpdate, function(x) any(unlist(x) == TRUE))
            Mdl.PostComp[[CompCpl]] <- FALSE
        }
    }
    else
    {
        stop(paste("MCMC update strategy:", MCMC.UpdateStrategy,
                   "not implemented!"))
    }

    ## CHECK THE DENSITY CALLER (PDF and/or CDF) FOR THE MARGINAL LIKELIHOODS
    densCaller <- list()
    for(CompCaller in names(CompUpMargis))
    {
        if(tolower(MCMC.UpdateStrategy) == "joint")
        {
            densCaller[[CompCaller]] <- c("u", "d")
        }
        else if(tolower(MCMC.UpdateStrategy) == "twostage")
        {
            ## Stage one of the two-stage approach
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

        if(evalCpl == TRUE && length(Mdl.u) == 0)
        {
            ## Special case when "u" is not supplied when evalCpl = TRUE
            densCaller[[CompCaller]] <- c("u", "d")
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
                             type = Mdl.MargisType[CompCaller],
                             par = Mdl.par[[CompCaller]],
                             densCaller = densCaller[[CompCaller]])

        if("u" %in% densCaller[[CompCaller]])
        {
            ## Mdl.u[,  CompCaller] <- Mdl.ud[["u"]]
            Mdl.u[[CompCaller]] <- Mdl.ud[["u"]]

        }
        if("d" %in% densCaller[[CompCaller]])
        {
            Mdl.d[[CompCaller]] <- Mdl.ud[["d"]]
        }
    }

    ## print(Mdl.u)
    ## print(Mdl.d)
    ## Updating the copula model if required
    CDCOPULA_NPARALLEL <- as.numeric(Sys.getenv("CDCOPULA_NPARALLEL"))
    if(evalCpl == TRUE)
    {

        if(!is.na(CDCOPULA_NPARALLEL) && CDCOPULA_NPARALLEL > 1)
        {
            logCplLikFUN.NM <- "logCplLikParallel"
        }
        else
        {
            logCplLikFUN.NM <- "logCplLik"
        }

        Mdl.logLikCpl.caller <- call(logCplLikFUN.NM,
                                     Mdl.u = Mdl.u,
                                     CplNM = CplNM,
                                     parCplRep = Mdl.par[[CompCpl]],
                                     sum = FALSE)
        Mdl.logLikCpl <- eval(Mdl.logLikCpl.caller) # n-by-1

        Mdl.d[[CompCpl]] <- Mdl.logLikCpl
    }
    ## browser()
    ## if(any(!is.finite(Mdl.u))) browser()
    out <- list(Mdl.d = Mdl.d, Mdl.u = Mdl.u,
                Mdl.PostComp = Mdl.PostComp)
    return(out)
}
