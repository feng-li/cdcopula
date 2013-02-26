##' The log posterior of the copula model
##'
##' The structure of the input are constructed via the design of variable
##' "MdlDataStuc" in the main setting file. See the individual description for
##' each variable in the setting files.
##' @param CplNM "character".
##'        The copula name.
##'
##' @param Mdl.Y "list"
##'        The responses of each marginal model.
##'
##' @param Mdl.X "list".
##'        The covariate used in each parameter components. The structure is
##'        designed by "MdlDataStruc" variable in the main file. The intercept
##'        is included if called in the data construction procedure.
##'
##' @param Mdl.beta "list".
##'
##' @param Mdl.betaIdx "list".
##'
##' @param Mdl.parLink "list".
##'        The link function used in the MCMC procedure. See the main setting
##'        file for details.
##'
##' @param varSelArgs "list"
##'
##' @param MargisTypes "list".
##'        The model type in each marginal distribution.
##'
##' @param priArgs "list".
##'        The prior settings for each parameter components.
##'
##' @param parUpdate "list".
##'        The parameters list to be updated. In the MCMC draw. Most time we
##'        are doing conditional posterior which means some components are kept
##'        uncaged. This can reduce computing time.
##'
##' @param staticCache "list"
##'        Miscellaneous arguments that are needed in the model.
##'
##' @param call.out
##' @param staticCacheOnly "logical"
##'        If TRUE,  only update the staticCache,  otherwise, do a full log
##'        posterior updating.
##'
##' @param Mdl.beta "list".
##'        The beta coefficients for each parameter. Note that the length
##'        should be same as the length of covariates in the covariates.
##'
##' @param Mdl.betaIdx "list".
##'        The variable selection index. For each parameter. The second column
##'        shows the proposal results, i.e. 1 for selected and 0 for not
##'        selection in current MCMC draw.
##'
##' @param tauTabular "list".
##'        The extra input for numerical dictionary looking up in the inverse
##'        calculation for the inverse Kendall's tau. See the documentation for
##'        "kendalltauTabular()" for details.
##'
##' @return "list".
##'        The list should contain the updated components.
##'
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Oct 24 15:07:01 CEST 2011;
##'       Current: Thu May 10 20:17:09 CEST 2012.
logPost <- function(CplNM, Mdl.Y, Mdl.X,Mdl.beta,Mdl.betaIdx,Mdl.parLink,
                    varSelArgs,MargisTypes,priArgs,parUpdate,staticCache = NA,
                    call.out = c("prior", "likelihood", "posterior", "staticCache")[3])
{
  ## Debugging symbol: if the warning should be printed out immediately.
  immediate. <- FALSE

  ## The pre-saved information. The idea is to make even staticCache is not
  ## available, the log posterior is still working.

  if(!is.na(staticCache))
    {
      Mdl.par <- staticCache[["Mdl.par"]]
      Mdl.u <- staticCache[["Mdl.u"]]
      Mdl.d <- staticCache[["Mdl.d"]]
      Mdl.logPri <- staticCache[["Mdl.logPri"]]
    }
  else
    {
      Mdl.par <- NA
      Mdl.u <- NA
      Mdl.d <- NA
      Mdl.logPri <- NA
    }

  ## Allocate the output structure
  Mdl.logLik <- NA
  Mdl.logPost <- NA

###----------------------------------------------------------------------------
### UPDATE THE LOG PRIORS
###----------------------------------------------------------------------------
  if(any(c("prior", "posterior", "staticCache") %in% call.out))
    {

      Mdl.logPri <- logPriors(
          Mdl.X = Mdl.X,
          Mdl.parLink = Mdl.parLink,
          Mdl.beta = Mdl.beta,
          Mdl.betaIdx = Mdl.betaIdx,
          varSelArgs = varSelArgs,
          priArgs = priArgs,
          Mdl.logPri = Mdl.logPri,
          parUpdate = parUpdate)
      ## Mdl.logPri <- unlist(Mdl.logPri, recursive = FALSE)[unlist(parUpdate)]

    }

###----------------------------------------------------------------------------
### THE MARGINAL LIKELIHOOD
###----------------------------------------------------------------------------

  if(any(c("likelihood", "posterior", "staticCache") %in% call.out))
    {
### Update Mdl.par
      Mdl.par <- parCplMeanFun(
          CplNM = CplNM,
          Mdl.X = Mdl.X,
          Mdl.parLink = Mdl.parLink,
          Mdl.beta = Mdl.beta,
          parUpdate = parUpdate,
          Mdl.par = Mdl.par)

      ## print(Mdl.par$BB7$lambdaL)
      ## print(Mdl.beta$BB7$lambdaL)

      if(any(is.na(unlist(Mdl.par))))
        {
          warning("DEBUGGING: NA happens when updating ``Mdl.par''...",
                  immediate. = immediat.)
        }

### Update marginal pdf and cdf
### the marginal u and only updated if the corresponding parameters are updated.
      CompNM <- names(Mdl.beta)
      MargisNM <- CompNM[CompNM != CplNM]
      for(CompCaller in MargisNM)
        {
          CompUpdate <- any(parUpdate[[CompCaller]] == TRUE)
          ## Marginal Update available
          if(CompUpdate == TRUE)
            {
              Margi.ud <- MargiModel(y = Mdl.Y[[CompCaller]],
                                     type = MargisTypes[CompCaller],
                                     par = Mdl.par[[CompCaller]])
              Mdl.u[, CompCaller] <- Margi.ud[["u"]] # the marginal cdf
              Mdl.d[, CompCaller] <- Margi.ud[["d"]] # the marginal pdf

              ## plot(Mdl.u, xlim = c(0, 1), ylim = c(0, 1))
            }

        }

    }


###----------------------------------------------------------------------------
### THE LOG LIKELIHOOD
###----------------------------------------------------------------------------

  if(any(c("likelihood", "posterior") %in% call.out))
    {
      Mdl.logLikCpl.sum <- logCplLik(
          u = Mdl.u,
          CplNM = CplNM,
          parCpl = Mdl.par[[CplNM]],
          staticCache = staticCache, logLik = TRUE) # n-by-1

      Mdl.logLikMargis.sum <- sum(Mdl.d)
      Mdl.logLik <- Mdl.logLikMargis.sum  + Mdl.logLikCpl.sum

      ## cat("Prior:    ", Mdl.logPri.sum, "\n")
      ## cat("CplLik:   ", Mdl.logLikCpl.sum, "\n")
      ## cat("MargisLik:", Mdl.logLikMargis.sum, "\n")

    }

###----------------------------------------------------------------------------
### THE STATIC ARGUMENT UPDATE
###----------------------------------------------------------------------------

  if("staticCache" %in% call.out)
    {
      staticCache[["Mdl.logPri"]] <- Mdl.logPri
      staticCache[["Mdl.par"]] <- Mdl.par
      staticCache[["Mdl.u"]] <- Mdl.u
      staticCache[["Mdl.d"]] <- Mdl.d
    }


###----------------------------------------------------------------------------
### THE LOG POSTERIOR
###----------------------------------------------------------------------------

  if("posterior" %in% call.out)
    {
      Mdl.logPri.sum <- sum(unlist(Mdl.logPri), na.rm = TRUE)
      Mdl.logPost <- Mdl.logLik + Mdl.logPri.sum
    }


  out <- list(Mdl.logPost = Mdl.logPost,
              Mdl.logLik = Mdl.logLik,
              Mdl.logPri = Mdl.logPri,
              staticCache = staticCache)
  return(out)
}
