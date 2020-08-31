#' @export
logPostVine <- function(Vine.Cpl.beta,
                        Vine.Cpl.betaIdx,
                        Vine.Margis.beta,
                        Vine.Margis.betaIdx,
                        Vine.RVM,
                        Vine.X,
                        Vine.Y,
                        Vine.InitTree,
                        Mdl.varSelArgs,
                        Mdl.priArgs,
                        Mdl.MargisType,
                        Mdl.parLink,
                        Mdl.parUpdate,
                        MCMC.UpdateStrategy)

{

###----------------------------------------------------------------------------
### Check the updating strategy
###----------------------------------------------------------------------------
    ## THE COPULA LOG LIKELIHOOD
    if(tolower(MCMC.UpdateStrategy) == "joint")
    {
        ## Joint updating should update all marginal components and copula component.
        evalCpl <- TRUE

        ## Mdl.PostComp <- lapply(parUpdate, function(x) TRUE)
    }
    else if(tolower(MCMC.UpdateStrategy) == "twostage" ||
            tolower(MCMC.UpdateStrategy) == "margin")
    {
        ## if(length(CompUpMargis) == 0)
        ## {
            ## Stage two of the two-stage approach: only updating copula density
            evalCpl <- TRUE
            ## Mdl.PostComp <- lapply(parUpdate, function(x) FALSE)
            ## Mdl.PostComp[[CplComp]] <- TRUE
        ## }
        ## else
        ## {
            ## Stage one of the two-stage approach: Only updating marginal density
            ## evalCpl <- FALSE
            ## Mdl.PostComp <- lapply(parUpdate, function(x) any(unlist(x) == TRUE))
            ## Mdl.PostComp[[CompCpl]] <- FALSE
        ## }
    }
    else
    {
        stop(paste("MCMC update strategy:", MCMC.UpdateStrategy,
                   "not implemented!"))
    }


###----------------------------------------------------------------------------
### Calculating conditional marginal models (parallel)
###----------------------------------------------------------------------------

    ## Marginal parameters
    MargisType.LLst <- unlist(Vine.InitTree[["MargisType"]], recursive = FALSE)
    Margis.beta.LLst <- unlist(Vine.Margis.beta, recursive = FALSE)
    MargisName.LLst <- names(MargisType.LLst)

    Margis.X.LLst <- unlist(Vine.X[["Margis.X"]], recursive = FALSE)
    names(Margis.X.LLst) <- MargisName.LLst # Margis.X.LLst did not have a name.

    Margis.parLink <- lapply(MargisType.LLst, function(x) Mdl.parLink[[x]])
    Margis.parUpdate <- lapply(MargisType.LLst, function(x) Mdl.parUpdate[[x]])


    Margis.par <-  parCplMeanFun(Mdl.X = Margis.X.LLst,
                                 Mdl.parLink = Margis.parLink,
                                 Mdl.beta = Margis.beta.LLst,
                                 parUpdate = Margis.parUpdate)


    ## Calculate the marginal densities
    Margis.Y.LLst <- unlist(Vine.Y, recursive = FALSE) # long list
    names(Margis.Y.LLst) <- MargisName.LLst


    ## TODO: If the updating strategy is two-stage, the `u.LLst` is
    ## static all the time. Not need to recalculate. See `logDens` for
    ## the solution in bivariate cases.

    u.LLst <- mapply(MargiModel, y = Margis.Y.LLst, type = MargisType.LLst,
                     par = Margis.par, SIMPLIFY = FALSE)

    logMargiLik <- 0


    ## Copula parameters
    BiCplInfo <- RVineCondBiCpl(Vine.RVM)
    BiCplFamilyNM.Lst <- as.list(BiCplInfo[["familyname"]])
    Cpl.Name.LLst <- paste("pair", as.character(1:length(BiCplFamilyNM.Lst)), sep = "")

    Cpl.X <- Vine.X[["Cpl.X"]]
    Cpl.parLink <- lapply(BiCplFamilyNM.Lst, function(x) Mdl.parLink[[x]])
    Cpl.parUpdate <- lapply(BiCplFamilyNM.Lst, function(x) Mdl.parUpdate[[x]])

    names(Cpl.X) <- Cpl.Name.LLst
    names(Cpl.parLink) <- Cpl.Name.LLst
    names(Cpl.parUpdate) <- Cpl.Name.LLst
    names(Vine.Cpl.beta) <- Cpl.Name.LLst


    Vine.Cpl.par <- parCplMeanFun(Mdl.X = Cpl.X,
                                  Mdl.parLink = Cpl.parLink,
                                  Mdl.beta = Vine.Cpl.beta,
                                  parUpdate = Cpl.parUpdate)



###----------------------------------------------------------------------------
### Log Priors
###----------------------------------------------------------------------------

    warning("Prior not set yet.")
    ## Vine.Margis.logPri <- logPriors(Mdl.X = Mdl.X,
    ##                                 Mdl.parLink = Mdl.parLink,
    ##                                 Mdl.beta = Mdl.beta,
    ##                                 Mdl.betaIdx = Mdl.betaIdx,
    ##                                 Mdl.varSelArgs = Mdl.varSelArgs,
    ##                                 Mdl.priArgs = Mdl.priArgs,
    ##                                 parUpdate = parUpdate)

    ## Vine.Cpl.logPri <- logPriors(Mdl.X = Mdl.X,
    ##                              Mdl.parLink = Mdl.parLink,
    ##                              Mdl.beta = Mdl.beta,
    ##                              Mdl.betaIdx = Mdl.betaIdx,
    ##                              Mdl.varSelArgs = Mdl.varSelArgs,
    ##                              Mdl.priArgs = Mdl.priArgs,
    ##                              parUpdate = parUpdate)

    Vine.logCplPriOut <- 0
    logMargiPriOut <- 0
###----------------------------------------------------------------------------
### Vine Copula Log  Likelihood
###----------------------------------------------------------------------------
    Vine.parAry <- par2RVMpar(Vine.Cpl.par, BiCplFamilyNM.Lst)

    par.Ary <- Vine.parAry[["par.Ary"]]
    par2.Ary <- Vine.parAry[["par2.Ary"]]

    Vine.data <- do.call(cbind, lapply(u.LLst, function(x) x[["u"]][, "u"]))

    Vine.dataminus <- do.call(cbind, lapply(u.LLst, function(x)
    {
        u <- x[["u"]]
        dim <- dim(u)
        if(dim[2] == 2)  u[, "um"] else matrix(0, dim[1])
    }))

    Vine.isdiscrete <- unlist(Vine.InitTree[["isdiscrete"]])


    ## Generate RVineMatrix structure. In this setp we just cheat RVineMatrix() with fake
    ## par and par2 just to generate the correct RVM because we have arrays. We my need to
    ## write a new RVineMatrix to reflect the array parameters we have.

    Vine.RVM[["par"]] <- par.Ary[, , 1]
    Vine.RVM[["par2"]] <- par2.Ary[, , 1]
    Vine.RVM[["names"]] <- names(unlist(Vine.InitTree[["MargisNM"]]))

    RVM = RVineMatrix(Vine.RVM$Matrix, Vine.RVM$family,
                      Vine.RVM$par, Vine.RVM$par2, Vine.RVM$names)

    Vine.logCplLik <- RVineLogLik_discretecontinuous(
        data = Vine.data,
        dataminus = Vine.dataminus,
        RVM = RVM,
        isdiscrete = Vine.isdiscrete,
        par = par.Ary,
        par2 = par2.Ary,
        separate = TRUE)


    out <- sum(Vine.logCplLik[["loglik"]]) + logMargiLik  + Vine.logCplPriOut + logMargiPriOut

    return(out)
}
