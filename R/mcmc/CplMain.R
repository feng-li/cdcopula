##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param CplConfigFile "character".
##'        The path where the setup script is located.
##'
##' @param Mdl.Idx.training "vector"
##'
##' @return "MCMC-details"
##'
##'        This will return the MCMC details into the global environment.
##'
##' @references Li, F., 2014 Covariate-dependent copula
##'
##' @author Feng Li, Central University of Finance and Economics.
##'
##' @note Created: Thu Feb 02 13:33:06 CET 2012;
##'       Current: Sun Dec 21 22:58:16 CST 2014.
CplMain <- function(Mdl.Idx.training, CplConfigFile)
{
###----------------------------------------------------------------------------
### Load user setup file
###----------------------------------------------------------------------------
    ## Load dependences
    require("mvtnorm", quietly = TRUE)
    require("optimx", quietly = TRUE)


    ## Load the sourceDir tool
    R_CPL_LIB_ROOT_DIR <- Sys.getenv("R_CPL_LIB_ROOT_DIR")
    if(length(R_CPL_LIB_ROOT_DIR) == 0L)
        {
            stop("R_CPL_LIB_ROOT_DIR is not set properly!")
        }

    sys.source(file.path(R_CPL_LIB_ROOT_DIR, "R/flutils/sourceDir.R"),
               envir = .GlobalEnv)

    ## Load the whole library
    Cpl.source <- sourceDir(file.path(R_CPL_LIB_ROOT_DIR, "R"),
                            byte.compile = 0,
                            recursive = TRUE,
                            ignore.error = TRUE)

    ## Source the configuration file for the model
    source(CplConfigFile, local = TRUE)

###----------------------------------------------------------------------------
### INITIALIZE THE STORAGE AND DATA STRUCTURE
###----------------------------------------------------------------------------

###----------------------------------------------------------------------------
### Initialize the MCMC
###----------------------------------------------------------------------------

    ## Extract the training and testing data according to cross-validation
    subsetFun <- function(x, idx)x[idx, , drop = FALSE]

    ## Mdl.Idx.training <- crossValidIdx[["training"]][[iCross]] # obtained from inputs
    nTraining <- length(Mdl.Idx.training)
    Mdl.X.training <- rapply(object=Mdl.X, f = subsetFun,
                             idx = Mdl.Idx.training, how = "replace")
    Mdl.Y.training <- rapply(object=Mdl.Y, f = subsetFun,
                             idx = Mdl.Idx.training, how = "replace")

    ## Switch all the updating indicators ON
    parUpdate <- MCMCUpdate

    ## Assign the initial values
    initParOut <- initPar(varSelArgs = varSelArgs, betaInit = betaInit,
                          Mdl.X = Mdl.X.training, Mdl.Y = Mdl.Y.training,
                          Mdl.parLink = Mdl.parLink)

    Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
    Mdl.beta <- initParOut[["Mdl.beta"]]

###----------------------------------------------------------------------------
### STABILIZE THE INITIAL VALUES VIA NEWTON ITERATIONS
###----------------------------------------------------------------------------

    ## Generate initial values that does not let log posterior be -Inf.
    ## Loop and count how many times tried for generating initial values
    optimInit <- TRUE

    if(optimInit == TRUE &&
       any(tolower(unlist(betaInit)) == "random"))
        {
            require(optimx)

            InitGood <- FALSE
            nLoopInit <- 0
            maxLoopInit <- 3

            Mdl.Idx.training.sample <-
                Mdl.Idx.training[seq(1, nTraining, length.out = 30)]

            Mdl.X.training.sample <- rapply(object=Mdl.X, f = subsetFun,
                                            idx = Mdl.Idx.training.sample, how = "replace")
            Mdl.Y.training.sample <- rapply(object=Mdl.Y, f = subsetFun,
                                            idx = Mdl.Idx.training.sample, how = "replace")

            cat("Optimizing initial values, may take a few minutes...\n\n")
            while(InitGood == FALSE)
                {
                    ## Reassign the initial values
                    initParOut <- initPar(
                        varSelArgs = varSelArgs,
                        betaInit = betaInit,
                        Mdl.X = Mdl.X.training.sample,
                        Mdl.Y = Mdl.Y.training.sample,
                        Mdl.parLink = Mdl.parLink)
                    Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
                    Mdl.beta <- initParOut[["Mdl.beta"]]

                    ## Dry run to obtain the initial "staticCache" NOTE: As this is the
                    ## very first run, the "parUpdate" should be all on for this time.


                    staticCache.sample <- logPost(
                        CplNM = CplNM,
                        Mdl.Y = Mdl.Y.training.sample,
                        Mdl.X = Mdl.X.training.sample,
                        Mdl.beta = Mdl.beta,
                        Mdl.betaIdx = Mdl.betaIdx,
                        Mdl.parLink = Mdl.parLink,
                        varSelArgs = varSelArgs,
                        MargisTypes = MargisTypes,
                        priArgs = priArgs,
                        parUpdate = MCMCUpdate,
                        MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

                    ## Optimize the initial values via BFGS. NOTE: The variable selection
                    ## indicators are fixed (not optimized) loop over all the marginal
                    ## models and copula via TWO STAGE OPTIMIZATION

                    for(iComp in names(Mdl.beta))
                        {

                            ## If nothing to update, optimization inside this components
                            ## skipped.
                            if(all(unlist(parUpdate[[iComp]]) == FALSE)) next

                            cat("Initializing model component:", iComp, "...\n")


                            ## Only current component is updated.
                            parUpdate <- rapply(MCMCUpdate, function(x) FALSE,
                                                    how = "replace")
                            parUpdate[[iComp]] <- MCMCUpdate[[iComp]]

                            betaVecInitComp <- parCplSwap(
                                betaInput = Mdl.beta,
                                Mdl.beta = Mdl.beta,
                                Mdl.betaIdx = Mdl.betaIdx,
                                parUpdate = parUpdate)

                            ## Optimize the initial values
                            betaVecOptimComp <- try(optimx(
                                par = betaVecInitComp,
                                fn = logPostOptim,
                                control = list(maximize = TRUE, maxit = 100,
                                    all.methods = FALSE),
                                method = "BFGS",
                                CplNM = CplNM,
                                Mdl.Y = Mdl.Y.training.sample,
                                Mdl.X = Mdl.X.training.sample,
                                Mdl.beta = Mdl.beta,
                                Mdl.betaIdx = Mdl.betaIdx,
                                Mdl.parLink = Mdl.parLink,
                                varSelArgs = varSelArgs,
                                MargisTypes = MargisTypes,
                                priArgs = priArgs,
                                staticCache = staticCache.sample,
                                parUpdate = parUpdate,
                                MCMCUpdateStrategy = "twostage"
                                ), silent = FALSE)

                            if(is(betaVecOptimComp, "try-error") == TRUE)
                                {# It does not have to be converged.
                                    cat("Initializing algorithm failed,  retry...\n")
                                    InitGood <- FALSE
                                    break
                                }
                            else
                                {
                                    InitGood <- TRUE
                                    ## print(betaVecOptimComp)
                                    Mdl.beta <- parCplSwap(
                                        betaInput = as.numeric(betaVecOptimComp[1,
                                            1:length(betaVecOptimComp)]),
                                        Mdl.beta = Mdl.beta,
                                        Mdl.betaIdx = Mdl.betaIdx,
                                        parUpdate = parUpdate)
                                }
                        }

                    nLoopInit <- nLoopInit +1

                    cat("The initial values for beta (conditional on variable selection indicators) are:\n")
                    print(Mdl.beta)

                    ## Too many failures,  abort!
                    if(nLoopInit >= maxLoopInit)
                        {
                            ## InitGood <- TRUE
                            cat("The initializing algorithm failed more that", nLoopInit, "times.\n")
                            cat("Trying to continue without initial value optimization.\n\n")
                            break
                        }

                }
        }

    ## Dry run to obtain staticCache for the initial values
    ## Again this time all the parameters should be updated.
    ## browser()
    staticCache <- logPost(
        CplNM = CplNM,
        Mdl.Y = Mdl.Y.training,
        Mdl.X = Mdl.X.training,
        Mdl.beta = Mdl.beta,
        Mdl.betaIdx = Mdl.betaIdx,
        Mdl.parLink = Mdl.parLink,
        varSelArgs = varSelArgs,
        MargisTypes = MargisTypes,
        priArgs = priArgs,
        parUpdate = MCMCUpdate,
        MCMCUpdateStrategy = MCMCUpdateStrategy)[["staticCache"]]

    ## Clear all warnings during initial value
    ## optimization. NOTE: not working
    ## warningsClear(envir = environment())

###----------------------------------------------------------------------------
### HARD DEBUGGING CODE， should be removed
###----------------------------------------------------------------------------
    Mdl.par <- staticCache$Mdl.par
    plot <- FALSE
    if(plot == TRUE)
        {
            ## nTraining <- length(Mdl.Y[[1]])
            X.ID0 <- as.Date(ID[1:nTraining])

            ## par(mfcol = c(5, 2), mar = c(2.5, 4, 2, 0))
            plot(X.ID0, Mdl.par[[1]][[1]], type = "l", col = "blue", xlab = "",
                 ylab = expression(mu), main = names(Mdl.Y)[1], ylim = c(-1, 2))
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")
            plot(X.ID0, Mdl.par[[1]][[2]], type = "l", col = "blue", xlab = "",
                 ylab = expression(phi), , log = "y")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[1]][[3]], type = "l", col = "blue", xlab = "",
                 ylab = expression(kappa),  log = "y")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[1]][[4]], type = "l", col = "blue", xlab = "",
                 ylab = expression(lambda), log = "y")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[3]][[1]], type = "l", col = "blue", xlab = "Time",
                 ylab = expression(tau), main = "Copula")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[2]][[1]], type = "l", col = "blue", xlab = "",
                 ylab = "", main = names(Mdl.Y)[2], , ylim = c(-1, 2))
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[2]][[2]], type = "l", col = "blue", xlab = "",
                 ylab = "", log = "y")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[2]][[3]], type = "l", col = "blue", xlab = "",
                 ylab = "", log = "y")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[2]][[4]], type = "l", col = "blue", xlab = "",
                 ylab = "", log = "y")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

            plot(X.ID0, Mdl.par[[3]][[2]], type = "l", col = "blue", xlab = "Time",
                 ylab = expression(lambda[L]), main = "Copula")
            abline(v = X.ID0[3000], col = "red", lwd = 2, lty = "dashed")

        }

###----------------------------------------------------------------------------
###  ALLOCATE THE STORAGE
###----------------------------------------------------------------------------

    ## The final parameters in current fold
    MCMC.beta <- MdlDataStruc
    MCMC.betaIdx <- MdlDataStruc
    MCMC.par <- MdlDataStruc
    MCMC.AccProb <- MdlDataStruc

    if(!exists("MCMC.density"))
        {
            MCMC.density <- list()
            MCMC.density[["d"]] <- array(NA, c(nTraining, length(MargisNM),  nIter))
            MCMC.density[["u"]] <- array(NA, c(nTraining, length(MargisNM),  nIter))
            ## FIXME: This is really big ~ 1G
        }

    for(i in names(MdlDataStruc))
        {
            for(j in names(MdlDataStruc[[i]]))
                {
                    ncolX.ij <- ncol(Mdl.X.training[[i]][[j]])
                    nPar.ij <- Mdl.parLink[[i]][[j]][["nPar"]]
                    namesX.ij <- rep(colnames(Mdl.X.training[[i]][[j]]), nPar.ij)

                    ## The MCMC storage
                    MCMC.betaIdx[[i]][[j]] <- matrix(
                        Mdl.betaIdx[[i]][[j]], nIter, ncolX.ij*nPar.ij, byrow = TRUE,
                        dimnames = list(NULL, namesX.ij))
                    MCMC.beta[[i]][[j]] <- matrix(
                        Mdl.beta[[i]][[j]], nIter, ncolX.ij*nPar.ij, byrow = TRUE,
                        dimnames = list(NULL, namesX.ij))
                    MCMC.par[[i]][[j]] <- array(
                        staticCache[["Mdl.par"]][[i]][[j]], c(nIter, nTraining, nPar.ij))

                    ## The Metropolis-Hasting acceptance rate
                    MCMC.AccProb[[i]][[j]] <- matrix(NA, nIter, 1)
                }
        }

###----------------------------------------------------------------------------
### THE METROPOLIS-HASTINGS WITHIN GIBBS
###----------------------------------------------------------------------------

    cat("Posterior sampling using Metropolis-Hastings within Gibbs\n")
    ## The updating matrix

    Starting.time <- Sys.time()
    UpdateMat <- parCplRepCaller(CplNM = CplNM,
                              parUpdate = MCMCUpdate,
                              parUpdateOrder = MCMCUpdateOrder)
    nInner <- nrow(UpdateMat)

    for(iUpdate in 1:(nInner*nIter))
        {
            iInner <- ifelse((iUpdate%%nInner) == 0, nInner, iUpdate%%nInner)
            iIter <- floor((iUpdate-1)/nInner)+1

            ##print(c(iUpdate, iInner, iIter))

            chainCaller <- UpdateMat[iInner, ]
            CompCaller <- UpdateMat[iInner, 1] ## SP100, SP600,  BB7
            parCaller <- UpdateMat[iInner, 2] ## mean,  df...

            ## Switch all the updating indicators OFF Switch current updating parameter
            ## indicator ON
            parUpdate <- rapply(MCMCUpdate, function(x) FALSE, how = "replace")
            parUpdate[[CompCaller]][[parCaller]] <- TRUE

            ## Call the proposal algorithm
            algmArgs <- propArgs[[CompCaller]][[parCaller]][["algorithm"]]

            if(tolower(algmArgs[["type"]]) == "gnewtonmove")
                {
                    ## staticCache <- list(u = u, Mdl.par = Mdl.par)
                    MHOut <- MHWithGNewtonMove(
                        CplNM = CplNM,
                        propArgs = propArgs,
                        varSelArgs = varSelArgs,
                        priArgs = priArgs,
                        parUpdate = parUpdate,
                        Mdl.Y = Mdl.Y.training,
                        Mdl.X = Mdl.X.training,
                        Mdl.beta = Mdl.beta,
                        Mdl.betaIdx = Mdl.betaIdx,
                        MargisTypes = MargisTypes,
                        Mdl.parLink = Mdl.parLink,
                        staticCache = staticCache,
                        MCMCUpdateStrategy = MCMCUpdateStrategy)

                    if(MHOut$errorFlag == FALSE)
                        {
                            ## Update the MH results to the current parameter structure
                            staticCache <- MHOut[["staticCache"]]
                            Mdl.beta[[CompCaller]][[parCaller]] <- MHOut[["beta"]]
                            Mdl.betaIdx[[CompCaller]][[parCaller]] <- MHOut[["betaIdx"]]

                            ## Export the parameters in each cross-validation fold
                            MCMC.beta[[CompCaller]][[parCaller]][iIter, ] <-
                                MHOut[["beta"]]
                            MCMC.betaIdx[[CompCaller]][[parCaller]][iIter, ] <-
                                MHOut[["betaIdx"]]
                            MCMC.AccProb[[CompCaller]][[parCaller]][iIter,] <-
                                MHOut[["accept.prob"]]
                            MCMC.par[[CompCaller]][[parCaller]][iIter, ,] <-
                                staticCache[["Mdl.par"]][[CompCaller]][[parCaller]]

                            ## Save the marginal densities.
                            ## FIXME: This is not updated for every iteration.
                            MCMC.density[["u"]][, , iIter] <-
                                MHOut[["staticCache"]][["Mdl.u"]]
                            MCMC.density[["d"]][, , iIter] <-
                                MHOut[["staticCache"]][["Mdl.d"]]

                        }
                    else
                        {
                            ## Set acceptance probability to zero if this
                            ## iteration fails.
                            MCMC.AccProb[[CompCaller]][[parCaller]][iIter,] <- 0
                        }

                    ## print(MHOut[["accept.prob"]])

                }
            else
                {
                    stop("Unknown proposal algorithm!")
                }

            ## MCMC trajectory
            if(track.MCMC == TRUE && iInner == nInner)
                {
                    CplMCMC.summary(iIter = iIter, nIter = nIter,
                                    interval = 0.01, burnin = burnin,
                                    MCMC.beta = MCMC.beta,
                                    MCMC.betaIdx = MCMC.betaIdx,
                                    MCMC.par = MCMC.par,
                                    MCMC.AccProb = MCMC.AccProb,
                                    Starting.time = Starting.time,
                                    MCMCUpdate = MCMCUpdate)
                }

        }


    ## Fetch everything at current environment to a list
    ## list2env(out, envir = .GlobalEnv)
    out <- as.list(environment())
    return(out)
}
