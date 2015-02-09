##' Generalized Newton move with dimension changes for the copula model.
##'
##' @param propArgs "list".
##' @param varSelArgs "list".
##' @param priArgs "list".
##' @param betaIdxProp "list".
##' @param parUpdate "list".
##' @param CplNM "list".
##' @param Mdl.Y "list".
##' @param Mdl.X "list".
##' @param Mdl.beta "list".
##' @param Mdl.betaIdx "list".
##' @param Mdl.parLink "list".
##' @param MargisTypes "list".
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
##' @note Created: Wed Sep 29 17:18:22 CEST 2010;
##'       Current: Mon Mar 05 10:33:29 CET 2012.
GNewtonMove <- function( propArgs, varSelArgs, priArgs, betaIdxProp, parUpdate,
                        CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                        Mdl.parLink, MargisTypes, staticCache, MCMCUpdateStrategy)
{
    ## The updating component parameter chain
    chainCaller <- parCplRepCaller(CplNM = CplNM, parUpdate)
    CompCaller <- chainCaller[1]
    parCaller <- chainCaller[2]

    errorFlag <- FALSE

    ## The current parameters
    X <- Mdl.X[[CompCaller]][[parCaller]]
    betaCurr <- Mdl.beta[[CompCaller]][[parCaller]] # p-by-lq
    betaIdxCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]] # p-by-lq

    ## Finite Newton move. K steps to approach the mode plus one more step to update the
    ## gradient and Hessian at i:th step.
    kSteps <- propArgs[[CompCaller]][[parCaller]][["algorithm"]][["ksteps"]]
    hessMethod <- propArgs[[CompCaller]][[parCaller]][["algorithm"]][["hess"]]

    ## Initialize the Newton move with the proposed variable selection indicator
    Mdl.betaIdx[[CompCaller]][[parCaller]] <- betaIdxProp
    param <- matrix(betaCurr[betaIdxCurr]) # col-vector

    ## Initial update staticCache for current Newton move
    staticCache.curr <- logPost(
        CplNM = CplNM,
        Mdl.Y = Mdl.Y,
        Mdl.X = Mdl.X,
        Mdl.beta = Mdl.beta,
        Mdl.betaIdx = Mdl.betaIdx,
        Mdl.parLink = Mdl.parLink,
        varSelArgs = varSelArgs,
        MargisTypes = MargisTypes,
        priArgs = priArgs,
        parUpdate = parUpdate,
        staticCache = staticCache,
        MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

###----------------------------------------------------------------------------
### The K-step Generalized Newton Move
###----------------------------------------------------------------------------

    for(iStep in 1:(kSteps+1))
        {
            browser()

            ## The gradient and Hessian in the likelihood
            logLikGradHess.prop <- logLikelihoodGradHess(
                CplNM = CplNM,
                Mdl.Y = Mdl.Y,
                Mdl.X = Mdl.X,
                Mdl.parLink = Mdl.parLink,
                Mdl.beta = Mdl.beta,
                MargisTypes = MargisTypes,
                Mdl.betaIdx = Mdl.betaIdx,
                parUpdate = parUpdate,
                varSelArgs = varSelArgs,
                staticCache = staticCache.curr,
                MCMCUpdateStrategy = MCMCUpdateStrategy)


            ## DEBUGING
            ## FIXME: DEBUGING code
            ## browser()

            DEBUGGING <- FALSE
            if(DEBUGGING == TRUE)
                {
                    logLikGradHess.prop.num <- logLikelihoodGradHess(
                        CplNM = CplNM,
                        Mdl.Y = Mdl.Y,
                        Mdl.X = Mdl.X,
                        Mdl.parLink = Mdl.parLink,
                        Mdl.beta = Mdl.beta,
                        MargisTypes = MargisTypes,
                        Mdl.betaIdx = Mdl.betaIdx,
                        parUpdate = parUpdate,
                        varSelArgs = varSelArgs,
                        staticCache = staticCache.curr,
                        gradMethods = "numeric",
                        MCMCUpdateStrategy = MCMCUpdateStrategy)

                    g.math <- logLikGradHess.prop[["logLikGradObs"]]
                    g.num <- logLikGradHess.prop.num[["logLikGradObs"]]
                    try(plot(g.num, g.math, main = as.character(chainCaller),
                             pch = 20, col = "blue"), silent = TRUE)

                    ## print(g.num)

                    ## Define gradient accuracy coefficient. The TRUE coefficient should be one
                    ## if analytical and numerical methods are of the same.
                    g.lm <- try(lm(g.math~0+g.num), silent = TRUE)
                    if(is(g.lm, "try-error") || abs(g.lm$coef-1)>0.1)
                        {
                            ## Sys.sleep(1)
                            ## browser(text = "Something Wrong!")
                        }
                }

            ## Break the loop if something went wrong in the gradient
            if(logLikGradHess.prop$errorFlag)
                {
                    errorFlag <- TRUE
                    break
                }

            ## Gradient Hessian for the prior *including non selected covariates*
            ## NOTE: The Hessian matrix of the prior is also approximated, we should
            ## use the explicit Hessian whenever possible.
            logPriGradHess.prop <- logPriorsGradHess(
                Mdl.X = Mdl.X,
                Mdl.beta = Mdl.beta,
                Mdl.betaIdx = Mdl.betaIdx,
                Mdl.parLink = Mdl.parLink,
                varSelArgs = varSelArgs,
                priArgs = priArgs,
                chainCaller = chainCaller)

            ## Gradient and Hessian for the likelihood
            logLikGrad.prop <- logLikGradHess.prop[["logLikGradObs"]] # n-by-pp


            ## cbind(logLikGrad.prop, logLikGrad.prop.num)
            ## browser()

            logLikHess.prop <- hessApprox(logLikGrad.prop, hessMethod)

            logPriGrad.prop <- logPriGradHess.prop[["gradObs"]] # pp-by-1
            logPriHess.prop <- logPriGradHess.prop[["HessObs"]] # pp-by-pp


            ## The gradient and Hessian subsets in the priors due to variable selection
            logPriGrad.pp <- logPriGrad.prop[betaIdxProp, , drop = FALSE] # pp-by-1
            logPriHess.pp <- logPriHess.prop[betaIdxProp, betaIdxProp, drop = FALSE] # pp-by-pp
            logPriHess.pc <- logPriHess.prop[betaIdxProp, betaIdxCurr, drop = FALSE] # pp-by-pc

            ## The selected covariates in the proposed and current draw
            X.prop <- X[ , betaIdxProp, drop = FALSE] # n-by-pp
            X.curr <- X[ , betaIdxCurr, drop = FALSE] # n-by-pc


            ## The gradient and Hessian in the general Newton's update
            gradObs.pp <- matrix(rowSums(Md(t(X.prop), logLikGrad.prop)) + logPriGrad.pp) # pp-by-1
            HessObs.pp <- tMdN(X.prop, logLikHess.prop, X.prop) + logPriHess.pp # pp-by-pp
            HessObs.pc <- tMdN(X.prop, logLikHess.prop, X.curr) + logPriHess.pc # pp-by-pc


            ## TODO: Testing if a subset of gradients works, seems not
            ## nObs <- dim(X)[1]
            ## ratio <- 0.3
            ## idx <- sample(1:nObs, round(ratio*nObs))

            ## gradObs.pp.sample <- matrix(1/ratio*rowSums(Md(t(X.prop[idx,, drop = FALSE ]),
            ##                                        logLikGrad.prop[idx]))
            ##                             + logPriGrad.pp) #

            ## plot(sort(gradObs.pp.sample), gradObs.pp[order(gradObs.pp.sample)], pch = 20,
            ##      type = "b")

            ## HessObsInv.pp <- qr.solve(HessObs.pp)
            HessObsInv.pp <- try(qr.solve(HessObs.pp), silent = TRUE) # pp-by-pp

            if(is(HessObsInv.pp, "try-error"))
                {
                    errorFlag <- TRUE
                    break
                }

            ## The general Newton's Update
            if((iStep <= kSteps))
                {
                    ## if(iStep == 2) browser()
                    ## update the proposed parameters via the general Newton formula
                    ## if(any(is.na(gradObs.pp))) browser()

                    param <- HessObsInv.pp%*%(HessObs.pc%*%param - gradObs.pp)

                    ## param <- param - diag(length(gradObs.pp))%*%gradObs.pp

                    ## print(HessObsInv.pp%*%gradObs.pp)

                    ## Update the parameter with current updated results.
                    ## If variable selection did not chose pth covariate,  then the pth
                    ## coefficient is zero naturally. After the first
                    ## variable-dimensional move, the algorithm switches to usual
                    ## Newton's move.
                    betaIdxCurr <- betaIdxProp

                    ## the full parameters including zeros.
                    param.full <- matrix(0, length(betaIdxCurr), 1)
                    param.full[betaIdxProp] <- param
                    Mdl.beta[[CompCaller]][[parCaller]] <- param.full

                    ## Update the staticCache
                    ## Initial update staticCache for current Newton move


                    staticCache.curr <- logPost(
                        CplNM = CplNM,
                        Mdl.Y = Mdl.Y,
                        Mdl.X = Mdl.X,
                        Mdl.beta = Mdl.beta,
                        Mdl.betaIdx = Mdl.betaIdx,
                        Mdl.parLink = Mdl.parLink,
                        varSelArgs = varSelArgs,
                        MargisTypes = MargisTypes,
                        priArgs = priArgs,
                        parUpdate = parUpdate,
                        staticCache = staticCache.curr,
                        MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

                    ## Mdl.par0 <- unlist(staticCache.curr[["Mdl.par"]])
                    ## if(any(is.na(Mdl.par0))) browser()

                }
            else # (k+1):th step.  Make a output
                {
                    out <- list(param = param,
                                gradObs = gradObs.pp,
                                HessObs = HessObs.pp,
                                HessObsInv = HessObsInv.pp,
                                staticCache = staticCache.curr,
                                errorFlag = errorFlag)
                    ## print(gradObs.pp)
                }
        }

    if(errorFlag)
        {
            out <- list(errorFlag = errorFlag)
        }

    return(out)
}
