##' The log posterior of the copula model
##'
##' The structure of the input are constructed via the design of variable "MdlDataStuc" in
##' the main setting file. See the individual description for each variable in the setting
##' files. This is used to calculate the conditional log posterior for the full copula
##' model.
##' @param CplNM "character".  The copula name.
##'
##' @param Mdl.Y "list" The responses of each marginal model.
##'
##' @param Mdl.X "list".  The covariate used in each parameter components. The structure
##'   is designed by "MCMCUpdate" variable in the main file. The intercept is included if
##'   called in the data construction procedure.
##'
##' @param Mdl.beta "list".
##'
##' @param Mdl.betaIdx "list".
##'
##' @param Mdl.parLink "list".  The link function used in the MCMC procedure. See the main
##'   setting file for details.
##'
##' @param varSelArgs "list"
##'
##' @param MargisType "list".  The model type in each marginal distribution.
##'
##' @param priArgs "list".  The prior settings for each parameter components.
##'
##' @param parUpdate "list".  The parameters list to be updated. In the MCMC draw. Most
##'   time we are doing conditional posterior which means some components are kept
##'   uncaged. This can reduce computing time.
##'
##' @param staticCache "list" Arguments that are cached in the model.
##'
##' @param call.out "character vector"
##'
##' @param split "logical"
##'
##'        If TRUE, the marginal model and copula model are split. This can be used in the
##' two stage method.
##'
##' @return "list".  The list should contain the updated components.
##'
##' @references Li 2012
##' @author Feng Li, Central University of Finance and Economics.
##' @note Created: Mon Oct 24 15:07:01 CEST 2011; Current: Sat Jul 18 10:47:12 CST 2015
logPost <- function(MargisType, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx, Mdl.parLink,
                    varSelArgs,priArgs,parUpdate,staticCache, MCMCUpdateStrategy)
{
    ## Assume no error a priori
    errorFlag <- FALSE
    ## Debugging symbol: if the warning should be printed out immediately.
    ## use options(warn = 1)
    ## immediate. <- FALSE

    ## The cached (pre-saved) information. The idea is to make even staticCache
    ## is not available, the log posterior is still working.
    ## TODO: Change staticCache.

    if(missing(staticCache))
    {
        ## Initialize "staticCache" structure.
        Mdl.u <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.Y),
                        dimnames = list(NULL, names(Mdl.Y)))

        Mdl.d <- cbind(Mdl.u, NA)
        colnames(Mdl.d) <- names(Mdl.beta)

        Mdl.par <- parCplMeanFun(Mdl.X = Mdl.X,
                                 Mdl.parLink = Mdl.parLink,
                                 Mdl.beta = Mdl.beta)

        Mdl.logPri <- logPriors(Mdl.X = Mdl.X,
                                Mdl.parLink = Mdl.parLink,
                                Mdl.beta = Mdl.beta,
                                Mdl.betaIdx = Mdl.betaIdx,
                                varSelArgs = varSelArgs,
                                priArgs = priArgs)

        staticCache <- list(Mdl.logPri =  parUpdate,
                            Mdl.par = Mdl.par,
                            Mdl.d = Mdl.d,
                            Mdl.u = Mdl.u)
    }
    else
    {
        Mdl.u <- staticCache[["Mdl.u"]]
        Mdl.d <- staticCache[["Mdl.d"]]
        Mdl.par <- staticCache[["Mdl.par"]]
        Mdl.logPri <- staticCache[["Mdl.logPri"]]
    }
###----------------------------------------------------------------------------
### UPDATE THE LOG PRIORS
###----------------------------------------------------------------------------
    Mdl.logPri <- logPriors(Mdl.X = Mdl.X,
                            Mdl.parLink = Mdl.parLink,
                            Mdl.beta = Mdl.beta,
                            Mdl.betaIdx = Mdl.betaIdx,
                            varSelArgs = varSelArgs,
                            priArgs = priArgs,
                            parUpdate = parUpdate,
                            Mdl.logPri = Mdl.logPri)

    Mdl.logPri.SubSum <- sum(unlist(mapply(function(x, y) ifelse(x == TRUE, y, 0),
                                           x = parUpdate, y = Mdl.logPri)))

###----------------------------------------------------------------------------
### UPDATING THE MODEL LIKELIHOOD PARAMETERS
###----------------------------------------------------------------------------
    ## browser()
    Mdl.par <- parCplMeanFun(Mdl.X = Mdl.X,
                             Mdl.parLink = Mdl.parLink,
                             Mdl.beta = Mdl.beta,
                             parUpdate = parUpdate,
                             Mdl.par = Mdl.par)

    ## if(any(is.nan(unlist(Mdl.beta)))) browser()

    Mdl.ud <- logDens(MargisType = MargisType,
                      Mdl.Y = Mdl.Y,
                      Mdl.par = Mdl.par,
                      Mdl.u = Mdl.u,
                      Mdl.d = Mdl.d,
                      parUpdate = parUpdate,
                      MCMCUpdateStrategy = MCMCUpdateStrategy)
    Mdl.d <- Mdl.ud[["Mdl.d"]]
    Mdl.u <- Mdl.ud[["Mdl.u"]]
    Mdl.PostComp <- Mdl.ud[["Mdl.PostComp"]]

    Mdl.logLik.SubSum <- sum(Mdl.d[, unlist(Mdl.PostComp)])
###----------------------------------------------------------------------------
### THE STATIC ARGUMENT UPDATE
###----------------------------------------------------------------------------
    staticCache[["Mdl.logPri"]] <- Mdl.logPri
    staticCache[["Mdl.par"]] <- Mdl.par
    staticCache[["Mdl.u"]] <- Mdl.u
    staticCache[["Mdl.d"]] <- Mdl.d

###----------------------------------------------------------------------------
### THE LOG POSTERIOR
###----------------------------------------------------------------------------
    Mdl.logPost <-  Mdl.logLik.SubSum+Mdl.logPri.SubSum

    out <- list(Mdl.logPost = Mdl.logPost,
                Mdl.logLik = Mdl.logLik.SubSum,
                Mdl.logPri = Mdl.logPri.SubSum,
                staticCache = staticCache,
                errorFlag = errorFlag)


    return(out)
}
