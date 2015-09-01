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
##' @note
##'       DEPENDS: flutils
##'       Created: Wed Mar 07 17:33:13 CET 2012;
##'       Current: Wed Mar 07 17:33:20 CET 2012.
DGPCpl <- function(DGPconfigfile, export = "list")
{
  ## TODO: check no visible bindings for DGPCpl,  DGP.par, MdlDGP.*

MdlDGP.par <- NA
MdlDGP.intercept <- NA
MdlDGP.parLink <- NA
MdlDGP.nCovs <- NA

  ## source the configure file
  source(file = DGPconfigfile, local = TRUE)


  ## THE RANDOM CDF VARIABLE IN THE COPULA
  uOut <- rCpl(n = nObs, parCpl = MdlDGP.par[[CplNM]], CplNM = CplNM)

  ## Generate the response variables
  Mdl.Y <- qCpl(u = uOut$u, parMargis = MdlDGP.par[MargisNM],
                MargisType = MargisType)

  ## The base covariates
  MdlDGP.beta <- MCMCUpdate
  Mdl.X <- MCMCUpdate
  Mdl.XFixed <- MCMCUpdate

  for(i in names(MCMCUpdate))
  {
    for(j in names(MCMCUpdate[[i]]))
    {
      Intercept <- ifelse(MdlDGP.intercept[[i]][[j]], TRUE, FALSE)

      linkCurr <- MdlDGP.parLink[[i]][[j]]
      if(tolower(linkCurr)  == "glogit")
      {
        warning("DGPCpl function needs review with conditional link,  Feng")
        tau <- MdlDGP.par[[CplNM]][["tau"]]
        a <- 0 ## The lower bound of generalized logit link
        b <- 2^(1/2-1/(2*tau)) ## the upper bound
        linkArgs <- list(a = a, b = b, type = linkCurr)
      }
      else
      {
        linkArgs <- list(type = linkCurr)
      }

      ## FIXME: Conditional link function
      ParResp <- parLinkFun(MdlDGP.par[[i]][[j]],
                            linkArgs = linkArgs)

      nCovsTol <- MdlDGP.nCovs[[i]][[j]]$total
      nCovsFixed <- MdlDGP.nCovs[[i]][[j]]$fixed

      betaFixed <- runif(n = nCovsFixed, min = 0, max = 1)

      Mdl.XFixed[[i]][[j]] <- DGPlm(Y = ParResp, beta = betaFixed,
                                    Xlim = c(0, 1),
                                    intercept = Intercept)
      MdlDGP.beta[[i]][[j]] <- matrix(c(betaFixed,
                                        rep(0, nCovsTol-nCovsFixed)))
    }
  }

  ## The extended covariates that are from the combination of the base
  ## covariates FIXME: it is better to select the non-fixed covariates from
  ## the known fixed covariates.
  for(i in names(MCMCUpdate))
  {
    for(j in names(MCMCUpdate[[i]]))
    {
      nCovsTol <- MdlDGP.nCovs[[i]][[j]]$total
      nCovsFixed <- MdlDGP.nCovs[[i]][[j]]$fixed

      XFinal1 <- Mdl.XFixed[[i]][[j]]
      XFinal2 <- matrix(runif(nObs*(nCovsTol-nCovsFixed)), nObs)

      XFinal <- cbind(XFinal1, XFinal2)

      Mdl.X[[i]][[j]] <- XFinal
    }
  }

  ## The output
  out <- list(Mdl.Y = Mdl.Y, Mdl.X = Mdl.X, MdlDGP.beta = MdlDGP.beta)
  if(tolower(export)  == "list")
  {
    return(out)
  }
  else if(tolower(export)  == "parent.env")
  {
    list2env(x = out, envir = sys.frame(sys.parent(1)))
  }
}
