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
##' @param priArgs "list".
##'        The prior settings for each parameter components.
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
##' @param Mdl.parLink "list".
##'        The link function used in the MCMC procedure. See the main setting
##'        file for details.
##' 
##' @param MargisTypes "list".
##'        The model type in each marginal distribution.
##' 
##' @param parUpdate "list".
##'        The parameters list to be updated. In the MCMC draw. Most time we
##'        are doing conditional posterior which means some components are kept
##'        uncaged. This can reduce computing time.
##' 
##' @param tauTabular "list".
##'        The extra input for numerical dictionary looking up in the inverse
##'        calculation for the inverse Kendall's tau. See the documentation for
##'        "kendalltauTabular()" for details.
##' 
##' @return "list".
##'        The list should contain the updated components.
##'        
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Oct 24 15:07:01 CEST 2011;
##'       Current: Tue Jan 10 17:10:10 CET 2012.
logPost <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx, Mdl.parLink,
                    varSelArgs, MargisTypes, priArgs, parUpdate, staticArgs)  
{


  Mdl.par <- staticArgs[["Mdl.par"]]
  Mdl.u <- staticArgs[["Mdl.u"]]
###----------------------------------------------------------------------------
### THE MARGINAL LIKELIHOOD
### The idea is to make even staticArgs is not available,  the log posterior is
### still working. 
###----------------------------------------------------------------------------

  CompNM <- names(Mdl.beta)
  MargisNM <- CompNM[CompNM != CplNM]
  for(i in CompNM)
    {
      parUpdateNM <- names(parUpdate[[i]] == TRUE) 
      for(j in parUpdateNM)
        {
          ## Update the parameters for the updated part
          Mdl.par[[i]][[j]] <- parMeanFun(X = Mdl.X[[i]][[j]],
                                          beta = Mdl.beta[[i]][[j]],
                                          link = Mdl.parLink[[i]][[j]])
        }

      ## Marginal Update available 
      if(length(parUpdateNM)>0 &&
         tolower(i) != tolower(CompNM)) # updating marginal parameters
        {
          MargisOut <- MargisModles(Mdl.Y = Mdl.Y,
                                    MargisTypes = MargisTypes,
                                    parMargis = Mdl.par[MargisNM],
                                    whichMargis = i)

          Mdl.u[, i] <- MargisOut[["Mdl.u"]] # the marginal cdf
          Mdl.d[, i] <- MargisOut[["Mdl.d"]] # the marginal pdf
        }
    }

###----------------------------------------------------------------------------
### THE COPULA LIKELIHOOD 
###----------------------------------------------------------------------------
  logLikCplOut <- logLikCpl(u = Mdl.u, CplNM = CplNM, parCpl = parCpl)

###----------------------------------------------------------------------------
### THE LOG PRIORS
###----------------------------------------------------------------------------
  Mdl.logPriOut <- logPriors(Mdl.X = Mdl.X,
                         Mdl.parLink = Mdl.parLink,
                         Mdl.beta = Mdl.beta,
                         Mdl.betaIdx = Mdl.betaIdx,
                         varSelArgs = varSelArgs,
                         priArgs = priArgs,
                         Mdl.logPri = Mdl.logPri,
                         parUpdate = parUpdate)
  

###----------------------------------------------------------------------------
### THE FINAL LOG POSTERIOR
###----------------------------------------------------------------------------
  
  logPostOut <- sum(unlist(logPriOut)) + logLikCplOut + sum(Mdl.d)


  staticArgs[["Mdl.logPri"]] <- Mdl.logPriOut
  staticArgs[["Mdl.par"]] <- Mdl.par
  staticArgs[["Mdl.u"]] <- Mdl.u
  staticArgs[["Mdl.d"]] <- Mdl.d

  
  out <- list(logPostOut = logPostOut,
              staticArgs = staticArgs)
  return(out)
}
