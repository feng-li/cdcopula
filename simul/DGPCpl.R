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
##' @note Created: Wed Mar 07 17:33:13 CET 2012;
##'       Current: Wed Mar 07 17:33:20 CET 2012.
DGPCpl <- function(configfile, export = "list")
  {
    ## source the configure file
    source(file = configfile, local = TRUE)
    
    ## COVARIATES USED IN THE MODEL
    X <- list()
    for(j in 1:length(MargisNM))
      {
        browser()
        nCovs <- 
        X[[MargisNM[j]]] <- matrix(runif(nObs*nCovs), nObs, nCovs)
      }
        
    ## COVARIATES USED FOR THE MARGINAL AND COPULA PARAMETERS
    Mdl.X <- MdlDataStruc
    Mdl.X[[1]][[1]] <- cbind(1, X[[1]])
    Mdl.X[[1]][[2]] <- cbind(1, X[[1]])
    Mdl.X[[2]][[1]] <- cbind(1, X[[2]])
    Mdl.X[[2]][[2]] <- cbind(1, X[[2]])
    Mdl.X[[3]][[1]] <- cbind(1, X[[1]], X[[2]])
    Mdl.X[[3]][[2]] <- cbind(1, X[[1]], X[[2]])

    ## PARAMETERS IN COPULA FUNCTION
    DGP.par <- MdlDataStruc
    for(i in 1:length(MdlDataStruc))
      {
        
        for(j in 1:length(MdlDataStruc[[i]]))
          {
            DGP.par[[i]][[j]] <- parMeanFun(X = Mdl.X[[i]][[j]],
                                            beta = DGP.betaTRUE[[i]][[j]],
                                            link = Mdl.parLink[[i]][[j]])
          } 
      }

    
    ## TRANSFORM COPULA PARAMETERS INTO THE STANDARD PARAMETRIZATION.
    DGP.parCpl <- kendalltauInv(CplNM = CplNM, parRepCpl = DGP.par[[CplNM]],
                                tauTabular = tauTabular)

    ## THE RANDOM CDF VARIABLE IN THE COPULA
    uOut <- ruCpl(n = nObs, parCpl = DGP.parCpl, copula = CplNM)

    ## THE RESPONSE VARIABLE
    Mdl.Y <- u2qtl(u = uOut$u, parMargis = DGP.par[MargisNM],
                   MargisTypes = MargisTypes)

    out <- list(Mdl.Y = Mdl.Y, X = X)
    if(tolower(export)  == "list")
      {
        return(out)
      }
    else if(tolower(export)  == "parent.env")
      {
        envirOut <- sys.fram(sys.parent(1))
        list2env(x = out, envir = exvirOut)        
      }
  }
