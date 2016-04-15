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


    ## THE RANDOM CDF VARIABLE IN THE COPULA
    CplNM <- MargisType[length(MargisType)]

    parCplStd <- parCplRep2Std(CplNM = CplNM, parCplRep = MdlDGP.par[[CplNM]])
    Mdl.u <- rCpl(n = nObs, parCpl = parCplStd, CplNM = CplNM)

    if(length(MargisType)>1)
    {
        Mdl.Y <- vector("list", nMargis-1)
        for(iComp in 1:(length(MargisType)-1))
        {
            Mdl.Y[[iComp]] <- MargiModelInv(u = Mdl.u$u[, iComp], par = MdlDGP.par[[iComp]],
                                            type = MargisType[[iComp]])
        }
    }

    ## The base covariates
    Mdl.X <- MCMCUpdate
    Mdl.XFixed <- MCMCUpdate

    for(i in names(MCMCUpdate))
    {
        for(j in names(MCMCUpdate[[i]]))
        {
            browser()

            ParResp <- parLinkFun(MdlDGP.par[[i]][[j]],
                                  linkArgs = Mdl.parLink[[i]][[j]])
            Mdl.XFixed[[i]][[j]] <- DGPlm(Y = ParResp, beta = MdlDGP.beta[[i]][[j]],
                                          Xlim = c(0, 1),
                                          intercept = MdlDGP.intercept[[i]][[j]])
        }
    }

    ## The output
    out <- list(Mdl.Y = Mdl.Y, Mdl.X = Mdl.X, Mdl.u = Mdl.u, MdlDGP.beta = MdlDGP.beta)
    if(tolower(export)  == "list")
    {
        return(out)
    }
    else if(tolower(export)  == "parent.env")
    {
        list2env(x = out, envir = sys.frame(sys.parent(1)))
    }
}
