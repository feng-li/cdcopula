##' DGP for the copula model.
##'
##'
##' @title Copula model DGP
##' @param configfile "character"
##'        The configuration file for the Copula DGP
##' @param export
##' \item {character string "list"}{Return a list containing the DGP results}
##' \item {character string "parent.env"}{The DGP ouput are written to the parent
##' environment directly}
##' @return See "export" argument.
##' @references Li, F., 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Wed Mar 07 17:33:13 CET 2012; Current: Tue Apr 12 19:09:49 CST 2016
DGPCpl <- function(DGPconfigfile, export = "list")
{
    CplNM <- NA
    MCMCUpdate <- NA
    MdlDGP.par <- NA
    MdlDGP.intercept <- NA
    MdlDGP.parLink <- NA
    MdlDGP.nCovs <- NA
    nObs <- NA

    ## source the configure file
    source(file = DGPconfigfile, local = TRUE)

    nMargis <- length(MargisType)
    out <- vector("list", nMargis)
    names(out) <- names(MargisType)
    MargisNM <- names(MargisType)

    ## THE RANDOM CDF VARIABLE IN THE COPULA
    CplNM <- MargisType[length(MargisType)]

    parCplStd <- parCplRep2Std(CplNM = CplNM, parCplRep = MdlDGP.par[[CplNM]])

    MdlDGP.u <- rCpl(n = nObs, parCpl = parCplStd, CplNM = CplNM)
    colnames(MdlDGP.u[["u"]]) <- MargisNM[-length(MargisNM)]

    ## Marginal models are also simulated
    CompUpLst <- unlist(lapply(MCMCUpdate, function(x) any(unlist(x) == TRUE)))
    Mdl.Y <- list()
    for(iComp in setdiff(MargisNM, CplNM))
    {
        if(CompUpLst[iComp])
        {
            Mdl.Y[[iComp]] <- MargiModelInv(u = MdlDGP.u$u[, iComp],
                                            par = MdlDGP.par[[iComp]],
                                            type = MargisType[[iComp]])
        }
        else
        {
            Mdl.Y[[iComp]] <- NA
        }
    }

    ## The base covariates
    Mdl.X <- MCMCUpdate
    for(i in names(MCMCUpdate))
    {
        for(j in names(MCMCUpdate[[i]]))
        {
            parLin <- parLinkFun(MdlDGP.par[[i]][[j]],
                                 linkArgs = Mdl.parLink[[i]][[j]])
            if(is.na(MdlDGP.beta[[i]][[j]]))
            {
                Mdl.X[[i]][[j]] <- matrix(1, nObs, 1)
                MdlDGP.beta[[i]][[j]] <- matrix(parLin[1])

            }
            else
            {
                Mdl.X[[i]][[j]] <- DGPlm(Y = parLin, beta = MdlDGP.beta[[i]][[j]],
                                         Xlim = c(0, 1),
                                         intercept = MdlDGP.intercept[[i]][[j]])
            }
        }
    }



    dev.width <- getOption("width")
    cat("DGP DATA SUMMARY...")
    cat("\n", rep("-", dev.width-1), "\n", sep = "")
    MdlDGP.Summary <- rbind(unlist(rapply(MdlDGP.par, mean, how = "replace")),
                            unlist((rapply(MdlDGP.par, sd, how = "replace"))),
                            unlist(MdlDGP.intercept))
    rownames(MdlDGP.Summary) <- c("par.mean", "par.sd", "intercept(Y/N)")
    print(MdlDGP.Summary)

    cat("DGP True Coefficients:\n")
    print(unlist(rapply(MdlDGP.beta, t, how = "replace"), recursive = T))
    cat("\n", rep("-", dev.width-1), "\n", sep = "")
    ## print(Mdl.parLink)


    ## The output
    out <- list(Mdl.Y = Mdl.Y, Mdl.X = Mdl.X,
                MdlDGP.u = MdlDGP.u, MdlDGP.beta = MdlDGP.beta,
                Mdl.parLink = Mdl.parLink)

    if(tolower(export)  == "list")
    {
        return(out)
    }
    else if(tolower(export)  == "parent.env")
    {
        list2env(x = out, envir = sys.frame(sys.parent(1)))
    }
}
