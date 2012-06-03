##' The log posterior of the copula model.
##'
##' The structure of the input are constructed via the design of variable
##' "MdlDataStuc" in the main setting file. See the individual description for
##' each variable in the setting files.
##' @title The log posterior function of the copula model
##'
##' @param CplNM "character".
##'        The copula name.
##'
##' @param Mdl.Y ""
##'
##' @param Mdl.X "list".
##'        The covariate used in each parameter components. The structure is
##'        designed by "MdlDataStruc" variable in the main file. The intercept
##'        is included if called in the data construction procedure.
##'
##' @param Mdl.beta
##' @param Mdl.betaIdx
##' @param Mdl.parLink "list".
##'        The link function used in the MCMC procedure. See the main setting
##'        file for details.
##'
##' @param varSelArgs
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
##' @param staticArgs
##' @param MdlCurr.beta "list".
##'        The beta coefficients for each parameter. Note that the length
##'        should be same as the length of covariates in the covariates.
##'
##' @param MdlCurr.betaIdx "list".
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
logPost <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx, Mdl.parLink,
                    varSelArgs, MargisTypes, priArgs, parUpdate, staticArgs)
{
  ## The pre-saved information
  Mdl.par <- staticArgs[["Mdl.par"]]
  Mdl.u <- staticArgs[["Mdl.u"]]
  Mdl.d <- staticArgs[["Mdl.d"]]
  Mdl.logPri <- staticArgs[["Mdl.logPri"]]

###----------------------------------------------------------------------------
### THE MARGINAL LIKELIHOOD
### The idea is to make even staticArgs is not available,  the log posterior is
### still working.
###----------------------------------------------------------------------------

### Update Mdl.par
  Mdl.par <- CplLinkConstrain(CplNM = CplNM,
                              Mdl.X = Mdl.X,
                              Mdl.parLink = Mdl.parLink,
                              Mdl.beta = Mdl.beta,
                              parUpdate = parUpdate,
                              Mdl.par = Mdl.par)

### Update marginal pdf and cdf
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
        }
    }

###----------------------------------------------------------------------------
### THE COPULA LIKELIHOOD
###----------------------------------------------------------------------------
  Mdl.logLikCpl <- logCplLik(u = Mdl.u,
                             CplNM = CplNM,
                             parCpl = Mdl.par[[CplNM]],
                             staticArgs = staticArgs) # n-by-1

###----------------------------------------------------------------------------
### THE LOG PRIORS
###----------------------------------------------------------------------------
  Mdl.logPri <- logPriors(Mdl.X = Mdl.X,
                          Mdl.parLink = Mdl.parLink,
                          Mdl.beta = Mdl.beta,
                          Mdl.betaIdx = Mdl.betaIdx,
                          varSelArgs = varSelArgs,
                          priArgs = priArgs,
                          Mdl.logPri = Mdl.logPri,
                          parUpdate = parUpdate)

###----------------------------------------------------------------------------
### THE FINAL LOG POSTERIOR AND STATIC ARGUMENT UPDATE
###----------------------------------------------------------------------------

  Mdl.logPri <- unlist(Mdl.logPri, recursive = FALSE)[unlist(parUpdate)]
  Mdl.logPost <- sum(unlist(Mdl.logPri)) + Mdl.logLikCpl + sum(Mdl.d)

  staticArgs[["Mdl.logPri"]] <- Mdl.logPri
  staticArgs[["Mdl.par"]] <- Mdl.par
  staticArgs[["Mdl.u"]] <- Mdl.u
  staticArgs[["Mdl.d"]] <- Mdl.d

  out <- list(Mdl.logPost = Mdl.logPost,
              staticArgs = staticArgs)

  return(out)
}
