##' Generalized Newton move with dimension changes for the copula model.
##'
##' @param MCMC.propArgs "list".
##' @param MCMC.varSelArgs "list".
##' @param Mdl.priArgs "list".
##' @param betaIdxProp "list".
##' @param parUpdate "list".
##' @param CplNM "list".
##' @param Mdl.Y "list".
##' @param Mdl.X "list".
##' @param Mdl.beta "list".
##' @param Mdl.betaIdx "list".
##' @param Mdl.parLink "list".
##' @param Mdl.MargisType "list".
##' @param staticCache "list".
##' @param param.cur "matrix".
##'         The initial values for the Newton update.
##' @return "list".
##' \item   {gradObs}
##'         {"matrix". The gradient}
##' \item   {HessObsInv}
##'         {"matrix". The inverse Hessian matrix.}
##' \item   {param}
##'         {"matrix". The updated parameters after K step Newton integrations.}
##' @references Li Villani Kohn 2010
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Wed Sep 29 17:18:22 CEST 2010; Current: Mon Mar 05 10:33:29 CET 2012.
PropGNewtonMove <- function(MCMC.propArgs, MCMC.varSelArgs, Mdl.priArgs, betaIdxProp, parUpdate,
                            Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx, Mdl.parLink,
                            Mdl.MargisType, staticCache, MCMC.UpdateStrategy)
{
    require("MASS")

    ## browser()
    ## The updating component parameter chain
    chainCaller <- parCplRepCaller(parUpdate)
    CompCaller <- chainCaller[1]
    parCaller <- chainCaller[2]

    errorFlag <- FALSE

    ## The current parameters
    X <- Mdl.X[[CompCaller]][[parCaller]]
    betaCurr <- Mdl.beta[[CompCaller]][[parCaller]] # p-by-lq
    betaIdxCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]] # p-by-lq

    ## if(any(is.nan(betaCurr))) browser()


    ## Finite Newton move. K steps to approach the mode plus one more step to update the
    ## gradient and Hessian at i:th step.
    kSteps <- MCMC.propArgs[[CompCaller]][[parCaller]][["algorithm"]][["ksteps"]]
    hessMethod <- MCMC.propArgs[[CompCaller]][[parCaller]][["algorithm"]][["hess"]]

    ## Initialize the Newton move with the proposed variable selection indicator
    Mdl.betaIdx[[CompCaller]][[parCaller]] <- betaIdxProp
    param <- matrix(betaCurr[betaIdxCurr]) # col-vector

    ## Initial update staticCache for current Newton move
    staticCache.curr <- logPost(Mdl.MargisType = Mdl.MargisType,
                                Mdl.Y = Mdl.Y,
                                Mdl.X = Mdl.X,
                                Mdl.beta = Mdl.beta,
                                Mdl.betaIdx = Mdl.betaIdx,
                                Mdl.parLink = Mdl.parLink,
                                MCMC.varSelArgs = MCMC.varSelArgs,
                                Mdl.priArgs = Mdl.priArgs,
                                parUpdate = parUpdate,
                                staticCache = staticCache,
                                MCMC.UpdateStrategy = MCMC.UpdateStrategy)[["staticCache"]]
###----------------------------------------------------------------------------
### The K-step Generalized Newton Move
###----------------------------------------------------------------------------
    for(iStep in 1:(kSteps+1))
    {

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE LINK FUNCTION
###----------------------------------------------------------------------------
        ## The gradient for the link function n-by-1
        LinkGradObs <- parCplMeanFunGrad(CplNM = Mdl.MargisType[length(Mdl.MargisType)],
                                         Mdl.par = staticCache[["Mdl.par"]],
                                         Mdl.parLink = Mdl.parLink,
                                         chainCaller = chainCaller)

        ## Error checking
        ## if(any(!is.finite(LinkGradObs)))
        ##   {
        ##     return(list(errorFlag = TRUE))
        ##   }


        ## The gradient and Hessian in the likelihood
        ## a <- proc.time()
        logDensGradHess.prop <- logDensGradHess(Mdl.MargisType = Mdl.MargisType,
                                                Mdl.Y = Mdl.Y,
                                                Mdl.parLink = Mdl.parLink,
                                                parUpdate = parUpdate,
                                                staticCache = staticCache.curr,
                                                gradMethods = c("analytic", "numeric")[1],
                                                MCMC.UpdateStrategy = MCMC.UpdateStrategy)

        ## cat("Analytical gradient:\n")
        ## print(proc.time()-a)

        ## DEBUGING FIXME: DEBUGING code
        DEBUGGING <- FALSE
        if(DEBUGGING == TRUE)
        {
            ## APPROACH ONE: This version calculates the numerical gradient with respect to
            ## the model parameters (not the covariate dependent beta parameters) via the
            ## chain rule.
            ## a <- proc.time()
            logDensGradHess.prop.num.split <- logDensGradHess(Mdl.MargisType = Mdl.MargisType,
                                                              Mdl.Y = Mdl.Y,
                                                              Mdl.parLink = Mdl.parLink,
                                                              parUpdate = parUpdate,
                                                              staticCache = staticCache.curr,
                                                              gradMethods = "numeric",
                                                              MCMC.UpdateStrategy = MCMC.UpdateStrategy)

            ## cat("Numerical gradient (split):\n")
            ## print(proc.time()-a)

            ## a <- proc.time()
            logDensGradHess.prop.num.joint <- logDensGradHessNum(Mdl.MargisType = Mdl.MargisType,
                                                                 Mdl.Y = Mdl.Y,
                                                                 Mdl.parLink = Mdl.parLink,
                                                                 parUpdate = parUpdate,
                                                                 staticCache = staticCache.curr,
                                                                 MCMC.UpdateStrategy = MCMC.UpdateStrategy) # n-by-pp

            ## cat("Numerical gradient (joint):\n")
            ## print(proc.time()-a)


            ## Define gradient accuracy coefficient. The TRUE coefficient should be one if
            ## analytical and numerical methods are of the same.
            if(!logDensGradHess.prop[["errorFlag"]])
            {

                g.ana <- logDensGradHess.prop[["logGradObs"]]
                g.num.split <- logDensGradHess.prop.num.split[["logGradObs"]]
                g.num.joint <- logDensGradHess.prop.num.joint[["logGradObs"]]

                nplot = ncol(g.ana)
                ncol = min(3, nplot)
                par(mfrow = c(ceiling(nplot/ncol), ncol))
                for(i in 1:nplot)
                {
                    try(plot(g.num.joint[, i], g.ana[, i], main = as.character(chainCaller),
                             pch = 20, col = "blue"), silent = TRUE)

                    g.lm <- try(lm(g.ana[, i]~0+g.num.joint[, i]), silent = TRUE)
                    if(is(g.lm, "try-error") || is.na(g.lm$coef) ||abs(g.lm$coef-1)>0.1)
                    {
                        ## Sys.sleep(1)
                        warning("Analytical and numerical gradients fail to match.")
                    }

                }
            }
            else
            {
                stop(paste("Analytical gradient error occurs in:", chainCaller[1], chainCaller[2]))
            }

        }

        ## Break the loop if something went wrong in the gradient
        if(logDensGradHess.prop[["errorFlag"]])
        {
            errorFlag <- TRUE
            ## browser()
            break
        }

        ## Gradient Hessian for the prior *including non selected covariates*
        ## NOTE: The Hessian matrix of the prior is also approximated, we should
        ## use the explicit Hessian whenever possible.
        logPriGradHess.prop <- logPriorsGradHess(Mdl.X = Mdl.X,
                                                 Mdl.beta = Mdl.beta,
                                                 Mdl.betaIdx = Mdl.betaIdx,
                                                 Mdl.parLink = Mdl.parLink,
                                                 MCMC.varSelArgs = MCMC.varSelArgs,
                                                 Mdl.priArgs = Mdl.priArgs,
                                                 chainCaller = chainCaller)

        ## Gradient and Hessian for the likelihood with linkage
        logDensGrad.prop <- logDensGradHess.prop[["logGradObs"]]*LinkGradObs # n-by-pp
        logDensHess.prop <- hessApprox(logDensGrad.prop, hessMethod) # diag elements only

        logPriGrad.prop <- logPriGradHess.prop[["gradObs"]] # pp-by-1
        logPriHess.prop <- diag(as.vector(hessApprox(logPriGrad.prop, hessMethod)),
                                nrow = length(logPriGrad.prop))
        ## logPriHess.prop <- logPriGradHess.prop[["HessObs"]] # pp-by-pp


        ## The gradient and Hessian subsets in the priors due to variable selection
        logPriGrad.pp <- matrix(logPriGrad.prop[betaIdxProp]) # pp-by-1
        logPriHess.pp <- logPriHess.prop[betaIdxProp, betaIdxProp, drop = FALSE] # pp-by-pp
        logPriHess.pc <- logPriHess.prop[betaIdxProp, betaIdxCurr, drop = FALSE] # pp-by-pc

        ## The selected covariates in the proposed and current draw

        ## Let mu = X%*%beta where X is n-by-q, beta is q-by-lq. Then we have "vec(mu) = (I
        ## %k% X) vec(beta)" where I is lq-by-lq. This generalized Villani (2009)'s
        ## generalized Newton algorithm that allows mu to be a matrix and

        nPar <- dim(betaCurr)[2]
        X.KR <- kronecker(diag1(nPar), X) # (n*lq)-by-(q*lq)

        X.prop <- X.KR[ , betaIdxProp, drop = FALSE]# n-by-pp
        X.curr <- X.KR[ , betaIdxCurr, drop = FALSE]# n-by-pp


        ## The gradient and Hessian in the general Newton's update
        gradObs.pp <- t(X.prop)%*%matrix(logDensGrad.prop) + logPriGrad.pp # pp-by-1
        HessObs.pp <- tMdN(X.prop, logDensHess.prop, X.prop) + logPriHess.pp # pp-by-pp
        HessObs.pc <- tMdN(X.prop, logDensHess.prop, X.curr) + logPriHess.pc # pp-by-pc
        HessObs.pp.Inv <- ginv(HessObs.pp) # pp-by-pp

        ## TODO: Testing if a subset of gradients works, seems not
        ## nObs <- dim(X)[1]
        ## ratio <- 0.3
        ## idx <- sample(1:nObs, round(ratio*nObs))

        ## gradObs.pp.sample <- matrix(1/ratio*rowSums(Md(t(X.prop[idx,, drop = FALSE ]),
        ##                                        logDensGrad.prop[idx]))
        ##                             + logPriGrad.pp) #

        ## plot(sort(gradObs.pp.sample), gradObs.pp[order(gradObs.pp.sample)], pch = 20,
        ##      type = "b")


        ## if(is(HessObs.pp.Inv, "try-error"))
        ##   {
        ##     errorFlag <- TRUE
        ##     break
        ##   }

        ## The general Newton's Update
        if((iStep <= kSteps))
        {
            ## update the proposed parameters via the general Newton formula
            ## param <- HessObs.pp.Inv%*%(HessObs.pc%*%param - gradObs.pp)
            param <- solve(HessObs.pp, (HessObs.pc%*%param - gradObs.pp)) ## robust solution

            ## Update the parameter with current updated results. If variable selection did
            ## not choose pth covariate, then the pth coefficient is zero naturally. After
            ## the first variable-dimensional move, the algorithm switches to usual Newton's
            ## move.
            betaIdxCurr <- betaIdxProp

            ## the full parameters including zeros.
            param.full <- array(0, dim(betaCurr))
            param.full[betaIdxProp] <- param
            Mdl.beta[[CompCaller]][[parCaller]] <- param.full

            ## Update the staticCache
            ## Initial update staticCache for current Newton move
            staticCache.curr <- logPost(Mdl.MargisType = Mdl.MargisType,
                                        Mdl.Y = Mdl.Y,
                                        Mdl.X = Mdl.X,
                                        Mdl.beta = Mdl.beta,
                                        Mdl.betaIdx = Mdl.betaIdx,
                                        Mdl.parLink = Mdl.parLink,
                                        MCMC.varSelArgs = MCMC.varSelArgs,
                                        Mdl.priArgs = Mdl.priArgs,
                                        parUpdate = parUpdate,
                                        staticCache = staticCache.curr,
                                        MCMC.UpdateStrategy = MCMC.UpdateStrategy)[["staticCache"]]

            ## Mdl.par0 <- unlist(staticCache.curr[["Mdl.par"]])
            ## if(any(is.na(Mdl.par0))) browser()

        }
        else # (k+1):th step.  Make a output
        {
            out <- list(param = param,
                        gradObs = gradObs.pp,
                        HessObs = HessObs.pp,
                        HessObsInv = HessObs.pp.Inv,
                        staticCache = staticCache.curr,
                        errorFlag = errorFlag)
            ## print(gradObs.pp)
        }
    }

    ## if(any(is.nan(param))) browser()

    if(errorFlag)
    {
        out <- list(errorFlag = errorFlag)
    }

    return(out)
}
