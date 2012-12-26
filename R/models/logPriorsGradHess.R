##' A collection of gradient for the copula model
##'
##' @title Gradient for priors
##'
##' @param Mdl.X
##' @param Mdl.beta
##' @param Mdl.betaIdx
##' @param Mdl.parLink
##' @param varSelArgs
##' @param priArgs
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
logPriorsGradHess <- function(
    Mdl.X,
    Mdl.beta,
    Mdl.betaIdx,
    Mdl.parLink,
    varSelArgs,
    priArgs,
    chainCaller)
  {
    ## Only update priors for parameters that need to update.
    ## Initial the storage structure for current log prior
    CompCurr <- chainCaller[[1]]
    parCurr <- chainCaller[[2]]

    ## Reserve list structure of the gradient and Hessian
    gradObsLst <- list()
    HessObsLst <- list()

    ## if(parCurr == "tau") browser()

###----------------------------------------------------------------------------
### Gradient and Hessian for the intercept as a special case
###----------------------------------------------------------------------------
    priArgsCurr <- priArgs[[CompCurr]][[parCurr]][["beta"]][["intercept"]]
    betaCurr <- Mdl.beta[[CompCurr]][[parCurr]][1] # the intercept
    linkCurr <- Mdl.parLink[[CompCurr]][[parCurr]]

    if(tolower(priArgsCurr[["type"]]) == "custom")
      {
        ## Call the any2any() function
        densOutput <- any2any(densArgs = priArgsCurr, linkArgs = linkCurr)
        mean <- densOutput$mean
        variance <- densOutput$variance
        shrinkage <- priArgsCurr[["output"]][["shrinkage"]]

        ## Gradient and Hessian for the intercept
        GradHessInt <- DensGradHess(B = betaCurr,
                                    mean = mean,
                                    covariance = variance*shrinkage,
                                    grad = TRUE, Hess = TRUE)

        ## The output
        gradObsLst[["Int"]] <- GradHessInt[["grad"]]
        HessObsLst[["Int"]] <- GradHessInt[["Hess"]]
      }
###----------------------------------------------------------------------------
### Gradient for beta|I and Hessian for beta (unconditional)
###----------------------------------------------------------------------------

    priArgsCurr <- priArgs[[CompCurr]][[parCurr]][["beta"]][["slopes"]]
    betaCurr <- Mdl.beta[[CompCurr]][[parCurr]][-1] # Slopes(taking away intercept)
    betaIdxNoIntCurr <- Mdl.betaIdx[[CompCurr]][[parCurr]][-1] # Variable section
                                        # indicator without intercept

    X <- Mdl.X[[CompCurr]][[parCurr]][, -1, drop = FALSE]
    if(length(X) == 0L)
      {
        ## No covariates at all (only intercept in the model)
        gradObsLst[["Slop"]] <- NULL
        HessObsLst[["SlopFull"]] <- NULL
      }
    else
      {
        if(tolower(priArgsCurr[["type"]]) == "cond-mvnorm")
          {
            ## Normal distribution condition normal The full beta vector is assumed
            ## as normal. Since variable selection is included in the MCMC, The
            ## final proposed beta are those non zeros. We need to using the
            ## gradient for the conditional normal density See Mardia p. 63.

            ## Subtract the prior information for the full beta
            mean <- priArgsCurr[["mean"]] # mean of density
            covariance <- priArgsCurr[["covariance"]] # Covariates
            shrinkage <- priArgsCurr[["shrinkage"]] # Shrinkage

            ## Split the beta vector by zero and nonzero.
            Idx1 <- which(betaIdxNoIntCurr == TRUE)
            Idx0 <- which(betaIdxNoIntCurr == FALSE)
            Idx0Len <- length(Idx0)
            Idx1Len <- length(Idx1)
            betaLen <- length(betaIdxNoIntCurr)

            SlopCondGrad <- matrix(NA, betaLen, 1)

            ## The mean vector for the whole beta vector (recycled if necessary)
            meanVec <- matrix(mean, betaLen, 1)

            ## The covariance matrix for the whole beta vector
            if(tolower(covariance) == "g-prior")
              {
                coVar <- qr.solve(crossprod(X))
              }
            else if(tolower(covariance) == "identity")
              {
                coVar <- diag(length(betaLen))
              }

            ##--------------------The conditional gradient--------------------------

            ## Consider three situations for gradient:
            if(Idx0Len == 0)
              {
                ## 1. all are selected. Switch to unconditional prior.
                ## The conditional gradient
                SlopCondGrad[Idx1] <- DensGradHess(
                    B = betaCurr,
                    mean = meanVec,
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
                    B = betaCurr[Idx1],
                    mean = condMean,
                    covariance = condCovar*shrinkage,
                    grad = TRUE, Hess = FALSE)[["grad"]]
              }
            else
              {
                ## 3. non are selected
                SlopCondGrad[Idx] <- NA
              }

            gradObsLst[["Slop"]] <- SlopCondGrad

            ## -------------The unconditional full Hessian matrix-------------------
            HessObsLst[["SlopFull"]] <- DensGradHess(
                B = betaCurr,
                mean = meanVec,
                covariance = coVar*shrinkage,
                grad = FALSE, Hess = TRUE)[["Hess"]]
          }
      }
###----------------------------------------------------------------------------
### The output
###----------------------------------------------------------------------------
    ## The final gradient output.
    ## The intercept and the conditional gradient; The unconditional Hessian

    gradObs <- matrix(unlist(gradObsLst))
    HessObs <- block.diag(HessObsLst)

    out <- list(gradObs = gradObs, HessObs = HessObs)
    return(out)
  }
