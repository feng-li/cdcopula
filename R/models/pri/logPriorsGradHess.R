##' A collection of gradient for the copula model
##'
##' @title Gradient for priors
##'
##' @param Mdl.X
##' @param Mdl.beta
##' @param Mdl.betaIdx
##' @param Mdl.parLink
##' @param Mdl.varSelArgs
##' @param Mdl.priArgs
##' @param chainCaller
##' @return "list". The gradient and Hessian matrix, see below.
##' \item   {gradObsPri}
##'         {"matrix". One-column.}
##' \item   {hessObsPri}
##'         {"matrix". A squre matrix. Dimension same as prior_type$Sigma0.}
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Mar 30 09:33:23 CEST 2010;
##'       Current: Mon Feb 27 15:22:27 CET 2012.
logPriorsGradHess <- function(Mdl.X, Mdl.beta, Mdl.betaIdx,Mdl.parLink,
                              Mdl.varSelArgs, Mdl.priArgs, chainCaller)
    {
        ## Only update priors for parameters that need to update.
        ## Initial the storage structure for current log prior
        CompCaller <- chainCaller[[1]]
        parCaller <- chainCaller[[2]]

        ## Reserve list structure of the gradient and Hessian
        gradObsLst <- list()
        HessObsLst <- list()

###----------------------------------------------------------------------------
### Gradient and Hessian for the intercept as a special case
###----------------------------------------------------------------------------
        Mdl.priArgsCurr <- Mdl.priArgs[[CompCaller]][[parCaller]][["beta"]][["intercept"]]
        betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][1, , drop = FALSE] # the intercepts
        linkCurr <- Mdl.parLink[[CompCaller]][[parCaller]]

        if(tolower(Mdl.priArgsCurr[["type"]]) == "custom")
            {
                ## Call the any2any() function
                densOutput <- any2any(densArgs = Mdl.priArgsCurr, linkArgs = linkCurr)

                mean <- densOutput$mean
                variance <- diag(densOutput$variance, length(betaCurr))
                shrinkage <- Mdl.priArgsCurr[["output"]][["shrinkage"]]

                ## Gradient and Hessian for the intercept
                GradHessInt <- DensGradHess(B = matrix(betaCurr),
                                            mean = mean,
                                            covariance = variance*shrinkage,
                                            grad = TRUE, Hess = FALSE)

                ## The output
                gradObsLst[["Int"]] <- t(GradHessInt[["grad"]]) # row-vector
                # HessObsLst[["Int"]] <- GradHessInt[["Hess"]]
            }
###----------------------------------------------------------------------------
### Gradient for beta|I and Hessian for beta (unconditional)
###----------------------------------------------------------------------------

        Mdl.priArgsCurr <- Mdl.priArgs[[CompCaller]][[parCaller]][["beta"]][["slopes"]]
        nPar <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]

        ## Slopes and variable selection indicators(taking away intercept)
        betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][-1, , drop = FALSE]
        betaIdxNoIntCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]][-1, , drop = FALSE]

        X <- Mdl.X[[CompCaller]][[parCaller]][, -1, drop = FALSE]
        if(length(X) == 0L)
            {
                ## No covariates at all (only intercept in the model)
                gradObsLst[["Slop"]] <- NULL
                HessObsLst[["SlopFull"]] <- NULL
            }
        else
            {
                if(tolower(Mdl.priArgsCurr[["type"]]) == "cond-mvnorm")
                    {
                        ## Normal distribution condition normal The full beta vector is assumed
                        ## as normal. Since variable selection is included in the MCMC, The
                        ## final proposed beta are those non zeros. We need to using the
                        ## gradient for the conditional normal density See Mardia p. 63.

                        ## Subtract the prior information for the full beta
                        mean <- Mdl.priArgsCurr[["mean"]] # mean of density
                        covariance <- Mdl.priArgsCurr[["covariance"]] # Covariates
                        shrinkage <- Mdl.priArgsCurr[["shrinkage"]] # Shrinkage

                        ## Split the beta vector by zero and nonzero.
                        Idx1 <- which(betaIdxNoIntCurr == TRUE)
                        Idx0 <- which(betaIdxNoIntCurr == FALSE)
                        Idx0Len <- length(Idx0)
                        Idx1Len <- length(Idx1)
                        betaLen <- length(betaIdxNoIntCurr)

                        SlopCondGrad <- array(NA, dim(betaIdxNoIntCurr))

                        ## The mean vector for the whole beta vector (recycled if necessary)
                        meanVec <- array(mean, dim(betaIdxNoIntCurr))

                        ## The covariance matrix for the whole beta vector
                        if(tolower(covariance) == "g-prior")
                            {
                                coVar0 <- qr.solve(crossprod(X))
                                coVar0Lst <- lapply(
                                    as.vector(rep(NA, nPar),"list"),
                                    function(x) coVar0)
                                coVar <- block.diag(coVar0Lst)
                            }
                        else if(tolower(covariance) == "identity")
                            {
                                coVar <- diag(betaLen)
                            }

                        ## The conditional gradient. Consider three situations for
                        ## gradient.
                        if(Idx0Len == 0)
                            {
                                ## 1. all are selected. Switch to unconditional prior.
                                SlopCondGrad[Idx1] <- DensGradHess(
                                    B = matrix(betaCurr),
                                    mean = matrix(meanVec),
                                    covariance = coVar*shrinkage,
                                    grad = TRUE, Hess = FALSE)[["grad"]]
                            }

                        else if(Idx0Len > 0 && Idx0Len < betaLen)
                            {
                                ## 2. some are selected (the most common situation)
                                ## The conditional prior
                                A <- coVar[Idx1, Idx0, drop = FALSE]%*%
                                    solve(coVar[Idx0, Idx0, drop = FALSE])

                                condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]

                                condCovar <- coVar[Idx1, Idx1, drop = FALSE] -
                                    A%*%coVar[Idx0, Idx1, drop = FALSE]

                                ## The conditional gradient
                                SlopCondGrad[Idx1] <- DensGradHess(
                                    B = matrix(betaCurr[Idx1]),
                                    mean = condMean,
                                    covariance = condCovar*shrinkage,
                                    grad = TRUE, Hess = FALSE)[["grad"]]
                            }
                        else
                            {
                                ## 3. non are selected, do nothing
                                ##  Switch to unconditional prior.
                                ## browser()
                                ## SlopCondGrad[Idx1] <- DensGradHess(
                                ##     B = betaCurr,
                                ##     mean = meanVec,
                                ##     covariance = coVar*shrinkage,
                                ##     grad = TRUE, Hess = FALSE)[["grad"]]
                                ## SlopCondGrad[Idx1] <- NA
                            }

                        gradObsLst[["Slop"]] <- SlopCondGrad

                        ## The unconditional full Hessian matrix
                        ## HessObsLst[["SlopFull"]] <- DensGradHess(
                        ##     B = matrix(betaCurr),
                        ##     mean = matrix(meanVec),
                        ##     covariance = coVar*shrinkage,
                        ##     grad = FALSE, Hess = TRUE)[["Hess"]]
                    }
            }
###----------------------------------------------------------------------------
### The output
###----------------------------------------------------------------------------
        ## The final gradient output.
        ## The intercept and the conditional gradient; The unconditional Hessian

        gradObs <- rbind(gradObsLst[["Int"]], gradObsLst[["Slop"]])

        ## TESTING code to reorganize the hesssian matrix, failed.
        HessObs <- NA
        ## HessObs <- block.diag(HessObsLst)

        ## IdxMat <- matrix(1:15, 5, 3)

        ## IdxSlop <- IdxMat[-1, , drop = FALSE]
        ## IdxInt <- IdxMat[1, ]
        ## HessTest <- matrix(0, 15, 15)

        ## HessTest[IdxSlop, IdxSlop] <- HessObsLst[["SlopFull"]]
        ## HessTest[IdxInt, IdxInt] <- HessObsLst[["Int"]]

        out <- list(gradObs = gradObs, HessObs = HessObs)
        return(out)
    }
