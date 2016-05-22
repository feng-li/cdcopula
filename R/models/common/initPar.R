##' Setup the initial values for the copula model
##'
##' This function is only once for the first iteration. The values are updated
##' via MCMC scheme.
##' @param Mdl.varSelArgs "list".
##'
##'        Variable selection argument
##'
##' @param Mdl.betaInit "character or numeric".
##'
##'        If is "character", the corresponding method are used to generate the
##' initial value.
##'
##'        If is "numeric", The initial value are taken as is.
##'
##' @param Mdl.X "list".
##'
##'        The covariates list
##'
##' @return "list"
##'
##'       Both the initial value for the beta coefficients ("Mdl.beta") and the
##' initial value for the variable selection indicator ("Mdl.betaIdx") are
##' returned.
##'
##' @references Li 2012
##' @author Feng Li, Central University of Finance and Economics.
##' @note Created: Thu Dec 22 15:57:14 CET 2011; Current: Tue Sep 29 23:50:07 CST 2015.
initPar <- function(Mdl.varSelArgs, Mdl.priArgs, Mdl.betaInit, Mdl.X, Mdl.Y, Mdl.MargisType,
                    Mdl.parLink, MCMC.Update, MCMC.optimInit)
{
    initParOut <- initPar0(Mdl.varSelArgs = Mdl.varSelArgs,
                           Mdl.betaInit = Mdl.betaInit,
                           Mdl.X = Mdl.X,
                           Mdl.Y = Mdl.Y,
                           Mdl.parLink = Mdl.parLink,
                           parUpdate = MCMC.Update,
                           MCMC.optimInit = MCMC.optimInit)
    Mdl.beta = initParOut[["Mdl.beta"]]
    Mdl.betaIdx = initParOut[["Mdl.betaIdx"]]

###----------------------------------------------------------------------------
### STABILIZE THE INITIAL VALUES VIA NEWTON ITERATIONS
###----------------------------------------------------------------------------
    ## Generate initial values that does not let log posterior be -Inf.
    ## Loop and count how many times tried for generating initial values
    if(MCMC.optimInit == TRUE &&
       any(tolower(unlist(Mdl.betaInit)) == "random"))
    {
        require("optimx", quietly = TRUE)

        ## Sample size to be used. Usually set it as 10% of the full data but no more than 500
        ## observations. Use the full data unless the it is really small.

        nObs <- length(Mdl.Y[[1]])

        if(nObs<50)
        {
            nOptim.Sample <- nObs
        }
        else
        {
            nOptim.Sample <- max(min(floor(nObs*.10), 500), 50)
        }

        subsetFun <- function(x, idx)x[idx, , drop = FALSE]

        Mdl.Idx.optim.sample <- floor(seq(1, nObs, length.out = nOptim.Sample))
        Mdl.X.optim.sample <- rapply(object=Mdl.X, f = subsetFun,
                                     idx = Mdl.Idx.optim.sample, how = "replace")

        if(exists("Mdl.Y"))
        {
            Mdl.Y.optim.sample <- rapply(object=Mdl.Y, f = subsetFun,
                                         idx = Mdl.Idx.optim.sample,
                                         how = "replace")
        }
        else
        {
            stop("Mdl.u should be supplied.")
            ## Mdl.u.optim.sample <- Mdl.u.training[Mdl.Idx.optim.sample,, drop = FALSE]
        }

        ## browser()
        cat("Optimizing initial values, may take a few minutes...\n\n")

        for(CompCaller in names(Mdl.beta))
        {
            ## If nothing to update, skip optimizing this component.
            if(all(unlist(MCMC.Update[[CompCaller]]) == FALSE)) next

            cat("\nInitializing model component:", CompCaller, "...\n")
            InitGoodCompCurr <- FALSE
            nLoopInit <- 0
            maxLoopInit <- 1

            ## Only current component to be updated.
            parUpdate <- rapply(MCMC.Update, function(x) FALSE, how = "replace")
            parUpdate[[CompCaller]] <- MCMC.Update[[CompCaller]]

            while(InitGoodCompCurr == FALSE){

                Mdl.betaIdx[[CompCaller]] <- initParOut[["Mdl.betaIdx"]][[CompCaller]]
                Mdl.beta[[CompCaller]] <- initParOut[["Mdl.beta"]][[CompCaller]]


                ## Optimize the initial values via BFGS. NOTE: The variable selection
                ## indicators are fixed (not optimized) loop over all the marginal
                ## models and copula via TWO STAGE OPTIMIZATION

                betaVecInitComp <- parCplSwap(betaInput = Mdl.beta,
                                              Mdl.beta = Mdl.beta,
                                              Mdl.betaIdx = Mdl.betaIdx,
                                              parUpdate = parUpdate)
                ## if(CompCaller == "MVT") browser()
                ## Optimize the initial values
                betaVecOptimComp <-  optimx(par = betaVecInitComp,
                                            fn = logPostOptim,
                                            control = list(maximize = TRUE,
                                                           all.methods = FALSE,
                                                           kkt = FALSE,
                                                           maxit = 100),
                                            method = "BFGS",
                                            hessian = FALSE,
                                            Mdl.MargisType = Mdl.MargisType,
                                            Mdl.Y = Mdl.Y.optim.sample,
                                            Mdl.X = Mdl.X.optim.sample,
                                            Mdl.beta = Mdl.beta,
                                            Mdl.betaIdx = Mdl.betaIdx,
                                            Mdl.parLink = Mdl.parLink,
                                            Mdl.varSelArgs = Mdl.varSelArgs,
                                            Mdl.priArgs = Mdl.priArgs,
                                            ## staticCache = staticCache.sample,
                                            parUpdate = parUpdate,
                                            MCMC.UpdateStrategy = "twostage")
                if(any(is.na(as.numeric(betaVecOptimComp[1, 1:length(betaVecInitComp)]))))
                {# It does not have to be converged.
                    cat("Initializing algorithm failed,  retry...\n")

                    InitGoodCompCurr <- FALSE
                    ## break

                    ## Reassign the initial values
                    initParOut <- initPar0(Mdl.varSelArgs = Mdl.varSelArgs,
                                           Mdl.betaInit = Mdl.betaInit,
                                           Mdl.X = Mdl.X.optim.sample,
                                           Mdl.Y = Mdl.Y.optim.sample,
                                           Mdl.parLink = Mdl.parLink,
                                           parUpdate = parUpdate,
                                           MCMC.optimInit = MCMC.optimInit)
                }
                else
                {
                    InitGoodCompCurr <- TRUE
                    Mdl.beta <- parCplSwap(betaInput = as.numeric(
                                               betaVecOptimComp[1, 1:length(betaVecInitComp)]),
                                           Mdl.beta = Mdl.beta,
                                           Mdl.betaIdx = Mdl.betaIdx,
                                           parUpdate = parUpdate)
                }
                nLoopInit <- nLoopInit +1
                if((nLoopInit >= maxLoopInit) & InitGoodCompCurr  == FALSE)
                {## Too many failures,  abort!
                    cat("The initializing algorithm failed more that", nLoopInit, "times.\n")
                    cat("Trying to continue without initial value optimization in this component.\n\n")
                    break
                }
            }
        }
    }

    cat("\nINITIAL VALUES FOR BETA COEFFICIENTS:\n",
        "(conditional on variable selection indicators)\n")
    for(iComp in names(Mdl.beta))
    {
        for(jPar in names(Mdl.beta[[iComp]]))
        {
            prtval <- t(Mdl.beta[[iComp]][[jPar]])
            rownames(prtval) <- paste(iComp, jPar, sep = ".")
            print(prtval)
        }
    }
    ## browser()
    out <- list(Mdl.beta = Mdl.beta,
                Mdl.betaIdx = Mdl.betaIdx)
    return(out)
}

initPar0 <- function(Mdl.varSelArgs, Mdl.betaInit, Mdl.X, Mdl.Y, Mdl.parLink,
                     parUpdate, MCMC.optimInit)
{
    ## The output structure.
    Mdl.betaIdx <- Mdl.betaInit
    Mdl.beta <- Mdl.betaInit
    nObs <- length(Mdl.Y[[1]])


    ## Loop to assign the initial values
    CompNM <- names(Mdl.betaInit)
    if(length(CompNM) == 0)
    {
        stop("Model component name must be specified first!")
    }


    for(i in CompNM)
    {
        CompParNM <- names(Mdl.parLink[[i]])
        for(j in CompParNM)
        {
            ## No. of col for covariates, including intercept if applicable.
            ncolX.ij <- ncol(Mdl.X[[i]][[j]])
            nPar.ij <- Mdl.parLink[[i]][[j]][["nPar"]]
###----------------------------------------------------------------------------
### Initialize the variable section indicator
###----------------------------------------------------------------------------
            ## Initial value for variable selection indicators which can be
            ## character or vector
            varSelCandConfigCurr <- Mdl.varSelArgs[[i]][[j]][["cand"]]

            if(class(varSelCandConfigCurr) == "character")
            {
                cand <- strsplit(varSelCandConfigCurr, ":")[[1]]
                cand1 <- as.numeric(cand[1])
                if(cand[2] == "end")
                    {
                        cand2 <- ncolX.ij
                    }
                else
                {
                    cand2 <- as.numeric(cand[2])
                }
                varSelCandCurr <- 2:ncolX.ij
            }
            else
            {
                varSelCandCurr <- varSelCandConfigCurr
            }

            varSelInitCurr <- Mdl.varSelArgs[[i]][[j]][["init"]]

            ## Check if variable selection candidates are all subsets of X.
            if(!all(varSelCandCurr %in% 1:ncolX.ij))
            {
                stop("Variable selection candidates are not subset of covariates in component: ",  i, "-", j, "!")
            }

            if(class(varSelInitCurr) == "character" &&
               tolower(varSelInitCurr) == "all-in")
            {
                Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)
            }
            else if(class(varSelInitCurr) == "character" &&
                    tolower(varSelInitCurr) == "all-out")
            {
                Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)

                Mdl.betaIdx[[i]][[j]][varSelCandCurr, ] <- FALSE
            }
            else if(class(varSelInitCurr) == "character" &&
                    tolower(varSelInitCurr) == "random")
            {

                Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)

                ## Randomly pick up half in
                betaIdxCurrSubOut <- sample(varSelCandCurr,
                                            round(length(varSelCandCurr)/2))
                Mdl.betaIdx[[i]][[j]][betaIdxCurrSubOut, ] <- FALSE
            }
            else # Do nothing, use user input
            {
                Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)
                Mdl.betaIdx[[i]][[j]][varSelCandCurr, ] <- FALSE
                Mdl.betaIdx[[i]][[j]][varSelInitCurr, ] <- TRUE
            }

###----------------------------------------------------------------------------
### Initial value for coefficients
###----------------------------------------------------------------------------
            Mdl.betaInitCurr <- Mdl.betaInit[[i]][[j]]
            if(class(Mdl.betaInitCurr) == "character" &&
               tolower(Mdl.betaInitCurr) == "random")
            {

                if(parUpdate[[i]][[j]] == FALSE && MCMC.optimInit  == FALSE)
                {
                    stop("Using random initial values but don't update (",  i, ", " ,  j ,
                         ") component is neither logical nor allowed.")
                }

                ## A simple version of random initial value.
                Mdl.beta[[i]][[j]] <- array(runif(ncolX.ij*nPar.ij, -1, 1) , c(ncolX.ij, nPar.ij))

            }
            else if (class(Mdl.betaInitCurr) == "character" &&
                     tolower(Mdl.betaInitCurr) == "ols")
            {
                Y <- Mdl.Y[[i]]
                X <- Mdl.X[[i]][[j]]

                lmcoef <- lm(Y~0+X)$coef
                Mdl.beta[[i]][[j]] <- array(lmcoef, c(ncolX.ij, nPar.ij))
            }
            else # Do nothing and use user input
            {
                Mdl.beta[[i]][[j]] <- array(Mdl.betaInitCurr, c(ncolX.ij, nPar.ij))
            }

        }
    }

    out <- list(Mdl.beta = Mdl.beta,
                Mdl.betaIdx = Mdl.betaIdx)
    return(out)
}
