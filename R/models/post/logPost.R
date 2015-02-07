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
##' is designed by "MdlDataStruc" variable in the main file. The intercept is included if
##' called in the data construction procedure.
##'
##' @param Mdl.beta "list".
##'
##' @param Mdl.betaIdx "list".
##'
##' @param Mdl.parLink "list".  The link function used in the MCMC procedure. See the main
##' setting file for details.
##'
##' @param varSelArgs "list"
##'
##' @param MargisTypes "list".  The model type in each marginal distribution.
##'
##' @param priArgs "list".  The prior settings for each parameter components.
##'
##' @param parUpdate "list".  The parameters list to be updated. In the MCMC draw. Most
##' time we are doing conditional posterior which means some components are kept
##' uncaged. This can reduce computing time.
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
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Mon Oct 24 15:07:01 CEST 2011; Current: Mon Jan 05 22:58:03 CST 2015
logPost <- function(CplNM, Mdl.Y, Mdl.X,Mdl.beta,Mdl.betaIdx,Mdl.parLink,
                    varSelArgs,MargisTypes,priArgs,parUpdate,staticCache,
                    MCMCUpdateStrategy)
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
            ## Initialize "staticCache" structure
            u <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.Y),
                        dimnames = list(NULL, names(Mdl.Y)))
            MdlDataStruc <- rapply(object = Mdl.beta,
                                   f = function(x)NA,
                                   how = "replace")
            staticCache <- list(Mdl.logPri =  MdlDataStruc,
                                Mdl.par = MdlDataStruc,
                                Mdl.d = u,
                                Mdl.u = u)
        }

    Mdl.par <- staticCache[["Mdl.par"]]
    Mdl.u <- staticCache[["Mdl.u"]]
    Mdl.d <- staticCache[["Mdl.d"]]
    Mdl.logPri <- staticCache[["Mdl.logPri"]]

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

###----------------------------------------------------------------------------
### UPDATING THE LIKELIHOOD
###----------------------------------------------------------------------------
    Mdl.par <- parCplMeanFun(CplNM = CplNM,
                             Mdl.X = Mdl.X,
                             Mdl.parLink = Mdl.parLink,
                             Mdl.beta = Mdl.beta,
                             parUpdate = parUpdate,
                             Mdl.par = Mdl.par)

### Update marginal pdf and/or cdf the marginal u and only updated if the corresponding
### parameters are updated.

    CompNM <- names(Mdl.beta)

    ## Allocate the output structure
    Mdl.logLik <- cbind(Mdl.d, NA)
    colnames(Mdl.logLik) <- CompNM

    CompUpNM <- unlist(lapply(parUpdate, function(x) any(unlist(x) == TRUE)))
    MargisNM <- CompNM[CompNM  != CplNM]
    MargisUpNM <- CompNM[(CompNM  != CplNM) & CompUpNM]

    ## The Marginal likelihoods
    for(iComp in MargisNM)
        {
            if(iComp %in% MargisUpNM)
                {
                    if(tolower(MCMCUpdateStrategy) == "joint")
                        {
                            densCaller <- c("u", "d")
                        }
                    else if(tolower(MCMCUpdateStrategy) == "twostage")
                        {
                            ## Stage two of the two stage approach
                            densCaller <- c("u", "d")
                        }
                    else if(tolower(MCMCUpdateStrategy) == "margin")
                        {
                            densCaller <- c(NA, "d")
                        }
                    else
                        {
                            stop(paste("MCMC update strategy:", MCMCUpdateStrategy,
                                       "not implemented!"))
                        }
                }
            else
                {
                    if(tolower(MCMCUpdateStrategy) == "joint" && any(is.na(Mdl.u)))
                        {
                            densCaller <- c("u", "d")
                        }
                }

            ## MARGINAL LIKELIHOOD
            Mdl.ud <- MargiModel(
                y = Mdl.Y[[iComp]],
                type = MargisTypes[which(MargisNM == iComp)],
                par = Mdl.par[[iComp]],
                densCaller = densCaller)

            if("u" %in% densCaller)
                {
                    Mdl.u[, iComp] <- Mdl.ud[["u"]]
                }
            if("d" %in% densCaller)
                {
                    Mdl.d[, iComp] <- Mdl.ud[["d"]]
                    Mdl.logLik[, iComp] <- Mdl.ud[["d"]]
                }

        }

    ## THE COPULA LOG LIKELIHOOD
    if(tolower(MCMCUpdateStrategy) == "joint")
        {
            evalCpl <- TRUE
            PostComp <- lapply(parUpdate, function(x) TRUE)
        }

    else if(tolower(MCMCUpdateStrategy) == "twostage" |
            tolower(MCMCUpdateStrategy) == "margin")
        {
            if(length(MargisUpNM) == 0)
                {
                    ## Stage two of the two stage approach
                    evalCpl <- TRUE
                    PostComp <- lapply(parUpdate, function(x) FALSE)
                    PostComp[[CplNM]] <- TRUE

                    if(any(is.na(Mdl.u)))
                        {
                            ## In stage two, Mdl.u is required to compute the copula
                            ## density.
                            for(iComp in MargisNM)
                                {
                                    Mdl.ud <- MargiModel(
                                        y = Mdl.Y[[iComp]],
                                        type = MargisTypes[which(MargisNM == iComp)],
                                        par = Mdl.par[[iComp]],
                                        densCaller = "u")
                                    Mdl.u[, iComp] <- Mdl.ud[["u"]]
                                }


                        }

                }
            else
                {
                    evalCpl <- FALSE
                    PostComp <- lapply(parUpdate, function(x) any(unlist(x) == TRUE))
                }
        }
    else
        {
            stop(paste("MCMC update strategy:", MCMCUpdateStrategy,
                       "not implemented!"))
        }


    if(evalCpl == TRUE)
        {
            Mdl.logLikCpl <- logCplLik(
                u = Mdl.u,
                CplNM = CplNM,
                parCplRep = Mdl.par[[CplNM]],
                sum = FALSE) # n-by-1

            Mdl.logLik[, CplNM] <- Mdl.logLikCpl
        }
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
    Mdl.logPri.SubSum <- sum(unlist(Mdl.logPri[unlist(PostComp)]))
    Mdl.logLik.SubSum <- sum(Mdl.logLik[, unlist(PostComp)])
    Mdl.logPost <-  Mdl.logLik.SubSum+Mdl.logPri.SubSum

    out <- list(Mdl.logPost = Mdl.logPost,
                Mdl.logLik = Mdl.logLik.SubSum,
                Mdl.logPri = Mdl.logPri.SubSum,
                staticCache = staticCache,
                errorFlag = errorFlag)

    return(out)
}
