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
##' @param staticArgs "list"
##'        Miscellaneous arguments that are needed in the model.
##'
##' @param staticArgsOnly "logical"
##'        If TRUE,  only update the staticArgs,  otherwise, do a full log
##'        posterior updating.
##'
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
logPost <- function(CplNM, Mdl.Y, Mdl.X,Mdl.beta,Mdl.betaIdx,Mdl.parLink,
                    varSelArgs,MargisTypes,priArgs,parUpdate,
                    staticArgs, staticArgsOnly = FALSE, parUpdate4Pri = parUpdate)
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
              immediate. = TRUE)
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

          plot(Mdl.u, xlim = c(0, 1), ylim = c(0, 1))
        }

    }

###----------------------------------------------------------------------------
### UPDATE THE LOG PRIORS
###----------------------------------------------------------------------------
  Mdl.logPri <- logPriors(
      Mdl.X = Mdl.X,
      Mdl.parLink = Mdl.parLink,
      Mdl.beta = Mdl.beta,
      Mdl.betaIdx = Mdl.betaIdx,
      varSelArgs = varSelArgs,
      priArgs = priArgs,
      Mdl.logPri = Mdl.logPri,
      parUpdate = parUpdate4Pri)

  ## Mdl.logPri <- unlist(Mdl.logPri, recursive = FALSE)[unlist(parUpdate)]

###----------------------------------------------------------------------------
### THE STATIC ARGUMENT UPDATE
###----------------------------------------------------------------------------

  staticArgs[["Mdl.logPri"]] <- Mdl.logPri
  staticArgs[["Mdl.par"]] <- Mdl.par
  staticArgs[["Mdl.u"]] <- Mdl.u
  staticArgs[["Mdl.d"]] <- Mdl.d

###----------------------------------------------------------------------------
### THE LOG POSTERIOR
###----------------------------------------------------------------------------

  if(staticArgsOnly == FALSE)
    {
      Mdl.logLikCpl.sum <- logCplLik(
          u = Mdl.u,
          CplNM = CplNM,
          parCpl = Mdl.par[[CplNM]],
          staticArgs = staticArgs, logLik = TRUE) # n-by-1

      Mdl.logPri.sum <- sum(unlist(Mdl.logPri), na.rm = TRUE)
      Mdl.logLikMargis.sum <- sum(Mdl.d)


      Mdl.logPost <-  Mdl.logLikMargis.sum  + Mdl.logLikCpl.sum + Mdl.logPri.sum

      ## cat("Prior:    ", Mdl.logPri.sum, "\n")
      ## cat("CplLik:   ", Mdl.logLikCpl.sum, "\n")
      ## cat("MargisLik:", Mdl.logLikMargis.sum, "\n")

    }
  else
    {
      Mdl.logPost <- NA
    }

  out <- list(Mdl.logPost = Mdl.logPost,
              staticArgs = staticArgs)
  return(out)
}
