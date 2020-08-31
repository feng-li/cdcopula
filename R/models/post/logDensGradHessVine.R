#' The gradient for the log posterior in the copula model
#'
#' This function usually does not need to modify. Note that this function does not give
#' the complete posteriors but its components, i.e. gradient for the likelihood function,
#' gradient for the prior, gradient for the linkage.
#' @param CplNM "character".
#'
#' @param Mdl.MargisType "list".
#'
#' @param Mdl.Y "list".
#'
#' @param Mdl.parLink "list".
#'
#' @param parUpdate "list".
#'
#' @param Mdl.varSelArgs "list".
#'
#' @param staticCache "list".
#'
#' @param gradMethods "character" If is "numeric", the numeric gradient is returned; else
#' return the analytical gradient.
#'
#' @return "list".
#'
#'        errorFlag: If anything went wrong, quit with TRUE.
#'
#'        logLilGradObs: Gradient for the log likelihood
#'
#' @references Li 2012
#' @author Feng Li, Central University of Finance and Economics.
#' @note Created: Thu Feb 02 22:45:42 CET 2012; Current: Mon Dec 22 20:25:44 CST 2014
#' @export
logDensGradHessVine <- function(Vine.Cpl.parUpdate,
                                Vine.Cpl.par,
                                Vine.Y,
                                Vine.u,
                                Vine.RVM,
                                Vine.InitTree,
                                gradMethods = c("analytic", "numeric")[1],
                                MCMC.UpdateStrategy)
{

    ## The updating component parameter chain
    chainCaller <- parCplRepCallerVine(Vine.Cpl.parUpdate)
    CompCaller <- chainCaller[1]
    parCaller <- chainCaller[2]


    BiCplInfo <- RVineCondBiCpl(Vine.RVM)
    BiCplCaller <- BiCplInfo[["family"]][CompCaller]
    BiCplCallerNM <- BiCplInfo[["familyname"]][CompCaller]
    MargisCaller <- BiCplInfo[["conditioned"]][CompCaller, ]


    CDCOPULA_NPARALLEL <- as.numeric(Sys.getenv("CDCOPULA_NPARALLEL"))
###----------------------------------------------------------------------------
### SPLIT THE GRADIENT INTO COPULA AND MARGINAL ACCORDING TO MCMC STRATEGY
###----------------------------------------------------------------------------
    if(tolower(MCMC.UpdateStrategy) == "joint")
    {
        stop("`MCMC.UpdateStrategy = " , MCMC.UpdateStrategy,  "` is not allowed in Vine!")
        evalCpl <- TRUE
        cplCaller <- paste("u", which(Mdl.MargisNM%in%CompCaller), sep = "")

        evalMargi <- TRUE
        densCaller <- c("u", "d")
    }
    else if(tolower(MCMC.UpdateStrategy) == "twostage")
    {
        ## Only update the vine copula parameters. This is different from the bivariate
        ## case in `logDens()`
        evalCpl <- TRUE
        cplCaller <- NA

        evalMargi <- FALSE
        densCaller <- c("d")
    }
    else if(tolower(MCMC.UpdateStrategy) == "margin")
    {
        ## Only evaluate the marginal models
        stop("`MCMC.UpdateStrategy = " , MCMC.UpdateStrategy,  "` is not allowed in Vine!")

        evalCpl <- FALSE
        cplCaller <- NA

        evalMargi <- TRUE
        densCaller <- c("d")
    }
    else
    {
        stop(paste("MCMC update strategy:", MCMC.UpdateStrategy,
                   "not implemented!"))
    }

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE MARGINAL LIKELIHOOD
###----------------------------------------------------------------------------
    if(evalMargi == TRUE)
    {
        yCurr <- Mdl.Y[[CompCaller]]
        parCurr <- Mdl.par[[CompCaller]]
        typeCurr <- Mdl.MargisType[CompCaller]

        ## Gradient Fraction in the marginal component. n-by-1
        if("analytic" %in% tolower(gradMethods))
        {
            if(!is.na(CDCOPULA_NPARALLEL) && CDCOPULA_NPARALLEL > 1)
            {
                MargiGradFUN.NM <- "MargiModelGradParallel"
            }
            else
            {
                MargiGradFUN.NM <- "MargiModelGrad"
            }
            MargiGradObs.ana.caller <- call(MargiGradFUN.NM, par = parCurr, y = yCurr,
                                            type = typeCurr, parCaller = parCaller,
                                            densCaller = densCaller)
            MargiGradObs.ana <- eval(MargiGradObs.ana.caller)
            MargiGradObs.u <- MargiGradObs.ana[["u"]]
            MargiGradObs.d <- MargiGradObs.ana[["d"]]
        }

        ## Numerical Gradients, check with analytical gradients.
        if("numeric" %in% tolower(gradMethods))
        {
            MargiGradObs.num <- MargiModelGradNum(y = yCurr, par = parCurr, type = typeCurr,
                                                  parCaller = parCaller,
                                                  densCaller = densCaller)
            MargiGradObs.u <- MargiGradObs.num[["u"]]
            MargiGradObs.d <- MargiGradObs.num[["d"]]

            ## DEBUG: Check if any gradient component is not correctly computed.  To check
            ## the overall gradient chain, look at the "PropGNewtonMove()" function. Below
            ## evaluates if the numeric and analytic gradients are consistent
            ## try(plot(sort(MargiGradObs.ana), MargiGradObs.num[order(MargiGradObs.ana)],
            ## type = "l", pch = 20, main = chainCaller), silent = TRUE)
        }
    }
    else
    {
        ## Only update the gradient for copula parameters
        ## Gradient Fraction in the copula component.
        MargiGradObs.u <- 1
        MargiGradObs.d <- 0
    }

    ## Error checking
    if(length(MargiGradObs.u) != 0 && (any(!is.finite(MargiGradObs.u))))
    {
        return(list(errorFlag = TRUE))
    }

    if(any(!is.finite(MargiGradObs.d)))
    {
        return(list(errorFlag = TRUE))
    }

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE COPULA LIKELIHOOD
###----------------------------------------------------------------------------
    if(evalCpl == TRUE)
    {
        if("analytic" %in% tolower(gradMethods))
        {
            if(!is.na(CDCOPULA_NPARALLEL) && CDCOPULA_NPARALLEL > 1)
            {
                ## parallel version
                logCplGrad.caller <- logCplRepGradVineParallel(
                    CplNM = CplNM,
                    Mdl.u = staticCache$Mdl.u,
                    parCplRep = Mdl.par[[CompCpl]],
                    parCaller = cplCaller)
            }
            else
            {
                BiCplInfo <- RVineCondBiCpl(Vine.RVM)
                BiCplFamilyNM.Lst <- as.list(BiCplInfo[["familyname"]])

                Vine.parAry <- par2RVMpar(Vine.Cpl.par, BiCplFamilyNM.Lst)

                par.Ary <- Vine.parAry[["par.Ary"]]
                par2.Ary <- Vine.parAry[["par2.Ary"]]

                Vine.data <- do.call(cbind, lapply(Vine.u, function(x) x[["u"]][, "u"]))

                Vine.dataminus <- do.call(cbind, lapply(Vine.u, function(x)
                {
                    u <- x[["u"]]
                    dim <- dim(u)
                    if(dim[2] == 2)  u[, "um"] else matrix(0, dim[1])
                }))

                Vine.isdiscrete <- unlist(Vine.InitTree[["isdiscrete"]])

                ## Fake the RVM matrix
                Vine.RVM[["par"]] <- par.Ary[, , 1]
                Vine.RVM[["par2"]] <- par2.Ary[, , 1]
                Vine.RVM[["names"]] <- names(unlist(Vine.InitTree[["MargisNM"]]))
                RVM = RVineMatrix(Vine.RVM$Matrix, Vine.RVM$family,
                                  Vine.RVM$par, Vine.RVM$par2, Vine.RVM$names)


                logCplGrad <- RVineLogLikGradient_discretecontinuous(
                    data = Vine.data,
                    dataminus = Vine.dataminus,
                    RVM = RVM,
                    isdiscrete = Vine.isdiscrete,
                    par = par.Ary,
                    par2 = par2.Ary,
                    start.V=NA,
                    posParams=(RVM$family > 0))

                ## logCplGrad.caller <- CplParRepGrad(
                ##     CplNM = BiCplCallerNM,
                ##     parCplRep = ???,
                ##     parCaller = ???)

            }


            ## The gradient for the copula function. n-by-1
            logCplGradObs.ana <- eval(logCplGrad.caller)
            logCplGradObs <- logCplGradObs.ana
        }
        if("numeric" %in% tolower(gradMethods))
        {
            logCplGradObs.num <- logCplRepGradNum(CplNM = CplNM,
                                                  Mdl.u = staticCache$Mdl.u,
                                                  parCplRep = Mdl.par[[CompCpl]],
                                                  parCaller = cplCaller)
            logCplGradObs <- logCplGradObs.num
        }
    }
    else
    {
        logCplGradObs <- 1
        MargiGradObs.u <- 0
    }

    ## Error checking
    if(length(logCplGradObs) == 0 || ## NULL happens
       any(!is.finite(logCplGradObs))) ## NA/NaN, Inf
    {
        return(list(errorFlag = TRUE))
    }


###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------
    ## The gradient for the full likelihood,  n-by-1
    logDensGradObs <- (logCplGradObs*MargiGradObs.u) + MargiGradObs.d #*LinkGradObs

    ## The output
    out <- list(logGradObs = logDensGradObs, # n-by-1
                logHessObs = NA,
                errorFlag = FALSE)
    return(out)
}
