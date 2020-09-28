#' DGP for the copula model.
#'
#'
#' @title Copula model DGP
#' @param configfile "character"
#'        The configuration file for the Copula DGP
#' @param export
#'  {character string "list"}{Return a list containing the DGP results}
#'  {character string "parent.env"}{The DGP ouput are written to the parent
#' environment directly}
#' @return See "export" argument.
#' @references Li, F., 2012
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @note Created: Wed Mar 07 17:33:13 CET 2012; Current: Tue Apr 12 19:09:49 CST 2016
#' @export
DGPCpl <- function(DGPconfigfile, export = "list")
{
    CplNM <- NA
    MdlDGP.Update <- NA
    MdlDGP.par <- NA
    MdlDGP.intercept <- NA
    nObs <- NA
    Mdl.MargisType <- NA
    Mdl.parLink <- NA

    ## source the configure file
    source(file = DGPconfigfile, local = TRUE)

    nMargis <- length(Mdl.MargisType)
    out <- vector("list", nMargis)
    names(out) <- names(Mdl.MargisType)
    Mdl.MargisNM <- names(Mdl.MargisType)

    ## THE RANDOM CDF VARIABLE IN THE COPULA
    CplNM <- Mdl.MargisType[length(Mdl.MargisType)]
    Mdl.MargisNM <- names(Mdl.MargisType)
    CompCpl <- Mdl.MargisNM[length(Mdl.MargisType)] # only copula component

    parCplStd <- parCplRep2Std(CplNM = CplNM, parCplRep = MdlDGP.par[[CompCpl]])

    MdlDGP.u <- rCpl(n = nObs, parCpl = parCplStd, CplNM = CplNM)
    colnames(MdlDGP.u[["u"]]) <- Mdl.MargisNM[-length(Mdl.MargisNM)]

    ## Marginal models are also simulated
    CompUpLst <- unlist(lapply(MdlDGP.Update, function(x) any(unlist(x) == TRUE)))
    Mdl.Y <- list()
    for(iComp in setdiff(Mdl.MargisNM, CplNM))
    {
        if(CompUpLst[iComp])
        {
            Mdl.Y[[iComp]] <- MargiModelInv(u = MdlDGP.u$u[, iComp],
                                            par = MdlDGP.par[[iComp]],
                                            type = Mdl.MargisType[[iComp]])
        }
        else
        {
            Mdl.Y[[iComp]] <- NA
        }
    }

    ## The base covariates
    Mdl.X <- MdlDGP.Update
    for(i in names(MdlDGP.Update))
    {
        for(j in names(MdlDGP.Update[[i]]))
        {
            parLin <- parLinkFun(MdlDGP.par[[i]][[j]],
                                 linkArgs = Mdl.parLink[[i]][[j]])
            if(all(is.na(MdlDGP.beta[[i]][[j]])))
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

    cat("\n", rep("-", dev.width-1), "\n", sep = "")
    cat("DGP DATA SUMMARY ...")
    cat("\n", rep("-", dev.width-1), "\n", sep = "")
    MdlDGP.Summary <- rbind(unlist(rapply(MdlDGP.par, mean, how = "replace")),
                            unlist((rapply(MdlDGP.par, sd, how = "replace"))),
                            unlist(MdlDGP.intercept))
    rownames(MdlDGP.Summary) <- c("par.mean", "par.sd", "intercept(Y/N)")
    print(MdlDGP.Summary)

    cat("\nDGP True Coefficients:\n")
    cat(rep("-", dev.width-1), "\n", sep = "")

    for(i in names(MdlDGP.beta))
    {
        for(j in names(MdlDGP.beta[[i]]))
        {
            cat(i,"." , j, ": ", sep = "")
            cat(MdlDGP.beta[[i]][[j]], "\n")
        }
    }

    cat(rep("-", dev.width-1), "\n", sep = "")

    ## The output
    out <- list(Mdl.Y = Mdl.Y, Mdl.X = Mdl.X,
                MdlDGP.u = MdlDGP.u, MdlDGP.beta = MdlDGP.beta,
                MdlDGP.par = MdlDGP.par,
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


## Slightly modify the DGP par to avoid under/over flow in parLinkFun()
#' @export
DGPCplRestrictPar <- function(MdlDGP.par, Mdl.parLink)
{
    out <- Map(function(i, j) Map(parRestricFun, i, j), MdlDGP.par, Mdl.parLink)

    return(out)
}
