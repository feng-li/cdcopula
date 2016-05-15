##' The Main file for MCMC for the copula model.
##'
##' For details of individual variables, see the setup file.
##' @param MdlConfigFile "character".
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
##' @note Created: Thu Feb 02 13:33:06 CET 2012; Current: Fri Mar 27 12:08:32 CST 2015.
CplMain <- function(Mdl.Idx.training, MdlConfigFile)
{

###----------------------------------------------------------------------------
### DEBUGGING
###----------------------------------------------------------------------------

    DEBUGGING <- FALSE
    if(DEBUGGING)
    {## Turn warnings into error
        options(warn = 100)
    }
###----------------------------------------------------------------------------
### LOAD USER SETUP FILE
###----------------------------------------------------------------------------
    ## Load dependences
    require("mvtnorm", quietly = TRUE)
    require("numDeriv", quietly = TRUE)
    require("Rmpfr", quietly = TRUE)

    ## Load the sourceDir tool
    R_CPL_LIB_ROOT_DIR <- Sys.getenv("R_CPL_LIB_ROOT_DIR")
    if(nchar(R_CPL_LIB_ROOT_DIR) == 0L)
    {
        stop("R_CPL_LIB_ROOT_DIR is not set properly!")
    }

    R_CPL_NPARALLEL <- Sys.getenv("R_CPL_NPARALLEL")
    if(nchar(R_CPL_NPARALLEL) == 0L)
    {
        stop("R_CPL_NPARALLEL is not set properly!")
    }


    sys.source(file.path(R_CPL_LIB_ROOT_DIR, "R/flutils/sourceDir.R"),
               envir = .GlobalEnv)

    ## Load the whole library
    Cpl.source <- sourceDir(file.path(R_CPL_LIB_ROOT_DIR, "R"),
                            byte.compile = 0,
                            recursive = TRUE,
                            ignore.error = TRUE)

    ## Source the configuration file for the model.
    nObs <- NA
    nCross <- NA
    Mdl.X <- NA
    Mdl.Y <- NA
    Mdl.MargisType <- NA
    Mdl.MargisNM <- NA
    Mdl.parLink <- NA
    Mdl.priArgs <- NA
    MCMC.varSelArgs <- NA
    Mdl.betaInit <- NA
    Mdl.u <- NA

    MCMC.nIter <- NA
    MCMC.Update <- NA
    MCMC.track <- NA
    MCMC.burninProp <- NA
    MCMC.UpdateStrategy <- NA
    MCMC.UpdateOrder <- NA
    MCMC.propArgs <- NA
    MCMC.optimInit <- NA

    source(MdlConfigFile, local = TRUE)

    ## Parallel Setting up
    if(as.numeric(R_CPL_NPARALLEL)>1)
    {
        require("parallel", quietly = TRUE)
        require("snow", quietly = TRUE)
        cl4MCMC <- makeCluster(as.numeric(R_CPL_NPARALLEL), type = "MPI")

        setDefaultCluster(cl4MCMC)
        ce4MCMC <- clusterEvalQ(cl4MCMC,
        {sourceDir(file.path(Sys.getenv("R_CPL_LIB_ROOT_DIR"), "R"),
                   byte.compile = 0,
                   recursive = TRUE,
                   ignore.error = TRUE)
        })
    }

    ## Generating simple model information
    Starting.time <- Sys.time()
    ModelDescription <- paste(c(Mdl.MargisNM[-length(Mdl.MargisNM)],"+",
                                Mdl.MargisType, "+" , "nObs",
                                nObs, "nCross", nCross,  "+",
                                format(Starting.time, "%Y%m%d@%H.%M")),
                              collapse = "")
###----------------------------------------------------------------------------
### INITIALIZE THE DATA STRUCTURE AND INITIAL VALUES
###----------------------------------------------------------------------------
    ## Extract the training and testing data according to cross-validation
    subsetFun <- function(x, idx)x[idx, , drop = FALSE]

    if(exists("Mdl.Y"))
    {   ## Marginal Y supplies
        Mdl.Y.training <- rapply(object=Mdl.Y, f = subsetFun,
                                 idx = Mdl.Idx.training, how = "replace")
    }
    else
    {
        ## No marginal information supplies. Update copula component only
        Mdl.u.training <- Mdl.u[Mdl.Idx.training, , drop = FALSE]
    }

    ## if(tolower(Mdl.MargisType[length(Mdl.MargisType)]) %in% c("gogarch", "dccgarch"))
    if(exists("ForeignModelSpec") &&  ForeignModelSpec  != NA)
    {## Special case when a foreign multivariate model is introduced. Fit the model and
        ## quit the MCMC directly.
        Mdl.ForeignFitted <-ModelForeignEval(model  = Mdl.MargisType[length(Mdl.MargisType)],
                                             spec = ForeignModelSpec,
                                             data = Mdl.Y.training)
        print(Mdl.ForeignFitted) # Model summary

        out <- as.list(environment())
        return(out)
    }

    if(any(rapply(Mdl.X, class) != "matrix"))
    {
        ## Special case when a foreign marginal model is introduced.  Evaluating
        ## Foreign marginal models.
        cat("Evaluating foreign marginal models...\n")

        Mdl.X.training.Fitted <- MargiModelForeignEval(Mdl.MargisNM = Mdl.MargisNM,
                                                       Mdl.MargisType = Mdl.MargisType,
                                                       MargisForeignConfig = Mdl.X,
                                                       Mdl.Y = Mdl.Y.training)

        Mdl.X.training <- c(Mdl.X.training.Fitted[["Mdl.X"]],
                            rapply(object=Mdl.X[Mdl.MargisNM[length(Mdl.MargisNM)]],
                                   f = subsetFun,
                                   idx = Mdl.Idx.training,
                                   how = "replace"))
        Mdl.MargisForeignFitted <- Mdl.X.training.Fitted[["Mdl.MargisForeignFitted"]]
    }
    else
    {## Native marginal model structure
        Mdl.X.training <- rapply(object=Mdl.X, f = subsetFun,
                                 idx = Mdl.Idx.training, how = "replace")
    }

    ## Assign the initial values
    initParOut <- initPar(MCMC.varSelArgs = MCMC.varSelArgs, Mdl.priArgs = Mdl.priArgs,
                          Mdl.betaInit = Mdl.betaInit, Mdl.MargisType = Mdl.MargisType,
                          Mdl.X = Mdl.X.training, Mdl.Y = Mdl.Y.training,
                          Mdl.parLink = Mdl.parLink, MCMC.Update = MCMC.Update,
                          MCMC.optimInit = MCMC.optimInit)

    Mdl.betaIdx <- initParOut[["Mdl.betaIdx"]]
    Mdl.beta <- initParOut[["Mdl.beta"]]
###----------------------------------------------------------------------------
###  ALLOCATE THE STORAGE
###----------------------------------------------------------------------------

    ## The final parameters in current fold
    MCMC.beta <- MCMC.Update
    MCMC.betaIdx <- MCMC.Update
    MCMC.par <- MCMC.Update
    MCMC.AccProb <- MCMC.Update

    if(!exists("MCMC.density"))
    {
        ## MCMC.density <- list()
        ## MCMC.density[["d"]] <- array(NA, c(nTraining, length(Mdl.MargisNM) +1,  MCMC.nIter))
        ## MCMC.density[["u"]] <- array(NA, c(nTraining, length(Mdl.MargisNM),  MCMC.nIter))
        ## FIXME: This is really big ~ 1G
    }

    nTraining <- length(Mdl.Idx.training)
    for(CompCaller in names(MCMC.Update))
    {
        for(parCaller in names(MCMC.Update[[CompCaller]]))
        {
            ncolX.ij <- ncol(Mdl.X.training[[CompCaller]][[parCaller]])
            nPar.ij <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]
            namesX.ij <- rep(colnames(Mdl.X.training[[CompCaller]][[parCaller]]), nPar.ij)

            if((CompCaller %in% Mdl.MargisNM[length(Mdl.MargisNM)]) & nPar.ij != 1)
            {
                nDim <- length(Mdl.MargisNM)-1
                namesParFull.ij <- matrix(paste(matrix(1:nDim, nDim, nDim),
                                                matrix(1:nDim, nDim, nDim, byrow = TRUE), sep = "."), nDim)
                namesPar.ij <- namesParFull.ij[lower.tri(namesParFull.ij, diag = FALSE)]
            }
            else
            {
                namesPar.ij <- "1.1"
            }


            nMCMC <- ifelse(MCMC.Update[[CompCaller]][[parCaller]], MCMC.nIter, 1)

            ## The MCMC storage
            MCMC.betaIdx[[CompCaller]][[parCaller]] <- matrix(Mdl.betaIdx[[CompCaller]][[parCaller]],
                                                              nMCMC, ncolX.ij*nPar.ij,
                                                              dimnames = list(NULL, namesX.ij),
                                                              byrow = TRUE)
            MCMC.beta[[CompCaller]][[parCaller]] <- matrix(Mdl.beta[[CompCaller]][[parCaller]],
                                                           nMCMC, ncolX.ij*nPar.ij,
                                                           dimnames = list(NULL, namesX.ij),
                                                           byrow = TRUE)
            ## if(CompCaller == "MVT" & parCaller  == "rho") browser()
            MCMC.par[[CompCaller]][[parCaller]] <- array(NA, c(nMCMC, nTraining, nPar.ij),
                                                         dimnames = list(NULL, NULL, namesPar.ij))

            ## The Metropolis-Hasting acceptance rate
            MCMC.AccProb[[CompCaller]][[parCaller]] <- matrix(NA, nMCMC, 1)
        }
    }

###----------------------------------------------------------------------------
### THE GIBBS SAMPLER (WITH METROPOLIS-HASTINGS)
###----------------------------------------------------------------------------
    ## Dry run to obtain staticcache for the initial values. Again this time all the
    ## parameters should be updated.

    Mdl.DryRun <- logPost(Mdl.MargisType = Mdl.MargisType,
                          Mdl.Y = Mdl.Y.training,
                          Mdl.X = Mdl.X.training,
                          Mdl.beta = Mdl.beta,
                          Mdl.betaIdx = Mdl.betaIdx,
                          Mdl.parLink = Mdl.parLink,
                          MCMC.varSelArgs = MCMC.varSelArgs,
                          Mdl.priArgs = Mdl.priArgs,
                          parUpdate = MCMC.Update,
                          MCMC.UpdateStrategy = MCMC.UpdateStrategy)
    staticCache <- Mdl.DryRun[["staticCache"]]
    ## browser()

    cat("\nMEAN INITIAL VALUES FOR MODEL PARAMETERS:\n")

    Mdl.par <- staticCache[["Mdl.par"]]
    Mdl.par.mean <- rapply(Mdl.par, mean, how = "replace")
    print(unlist(Mdl.par.mean, recursive = TRUE))

    cat("\nPosterior sampling using Metropolis-Hastings within Gibbs\n")

    ## Clear all warnings during initial value optimization. NOTE: not working
    ## warningsClear(envir = environment())

    ## The updating matrix
    UpdateMat <- parCplRepCaller(parUpdate = MCMC.Update, parUpdateOrder = MCMC.UpdateOrder)
    nInner <- nrow(UpdateMat)

    for(iUpdate in 1:(nInner*MCMC.nIter))
    {
        iInner <- ifelse((iUpdate%%nInner) == 0, nInner, iUpdate%%nInner)
        iIter <- floor((iUpdate-1)/nInner)+1

        ## print(c(iUpdate, iInner, iIter))
        chainCaller <- UpdateMat[iInner, ]
        CompCaller <- UpdateMat[iInner, 1] ## SP100, SP600,  BB7
        parCaller <- UpdateMat[iInner, 2] ## mean,  df...

        ## Switch all the updating indicators OFF and switch current updating parameter
        ## indicator ON.
        parUpdate <- rapply(MCMC.Update, function(x) FALSE, how = "replace")
        parUpdate[[CompCaller]][[parCaller]] <- TRUE

        ## Call the Metropolis-Hastings algorithm
        MHOut <- MetropolisHastings(Mdl.MargisType = Mdl.MargisType,
                                    MCMC.propArgs = MCMC.propArgs,
                                    MCMC.varSelArgs = MCMC.varSelArgs,
                                    Mdl.priArgs = Mdl.priArgs,
                                    parUpdate = parUpdate,
                                    Mdl.Y = Mdl.Y.training,
                                    Mdl.X = Mdl.X.training,
                                    Mdl.beta = Mdl.beta,
                                    Mdl.betaIdx = Mdl.betaIdx,
                                    Mdl.parLink = Mdl.parLink,
                                    staticCache = staticCache,
                                    MCMC.UpdateStrategy = MCMC.UpdateStrategy)

        if(MHOut$errorFlag == FALSE)
        { ## Update the MH results to the current parameter structure
            staticCache <- MHOut[["staticCache"]]
            Mdl.beta[[CompCaller]][[parCaller]] <- MHOut[["beta"]]
            Mdl.betaIdx[[CompCaller]][[parCaller]] <- MHOut[["betaIdx"]]

            AccProbCurr <- MHOut[["accept.prob"]]
        }
        else
        { ## Keep everything unchanged. Only set acceptance probability to zero if this
            ## iteration fails.
            AccProbCurr <- 0
        }

        ## Export the parameters in each cross-validation fold
        MCMC.beta[[CompCaller]][[parCaller]][iIter, ] <- Mdl.beta[[CompCaller]][[parCaller]]
        MCMC.betaIdx[[CompCaller]][[parCaller]][iIter, ] <- Mdl.betaIdx[[CompCaller]][[parCaller]]
        MCMC.par[[CompCaller]][[parCaller]][iIter, ,] <- staticCache[["Mdl.par"]][[CompCaller]][[parCaller]]

        MCMC.AccProb[[CompCaller]][[parCaller]][iIter,] <- AccProbCurr

        ## Save the marginal densities.  FIXME: This is not updated for every iteration.
        ## MCMC.density[["u"]][, , iIter] <- staticCache[["Mdl.u"]]
        ## MCMC.density[["d"]][, , iIter] <- staticCache[["Mdl.d"]]

        ## MCMC trajectory
        if(MCMC.track == TRUE && iInner == nInner)
        {
            ## ProgressBar only available in interactive mode
            if(interactive()) progressbar(iIter = iIter, nIter = MCMC.nIter)

            CplMCMC.summary(iIter = iIter, MCMC.nIter = MCMC.nIter,
                            interval = 1, MCMC.burninProp = MCMC.burninProp,
                            OUT.MCMC = list(MCMC.beta = MCMC.beta,
                                            MCMC.betaIdx = MCMC.betaIdx,
                                            MCMC.par = MCMC.par,
                                            MCMC.AccProb = MCMC.AccProb,
                                            MCMC.Update = MCMC.Update,
                                            Starting.time = Starting.time)
                            )
        }

    }


    if(R_CPL_NPARALLEL>1)
    {
        stopCluster(cl4MCMC)
    }
    ## Fetch everything at current environment to a list
    ## list2env(out, envir = .GlobalEnv)
    ## GENERATING SHORT MODEL DESCRIPTION
    ModelDescription <- paste(c(Mdl.MargisNM[-length(Mdl.MargisNM)],"+",  Mdl.MargisType, "+" ,
                                "nObs", nObs, "nCross", nCross, "MCMC.nIter", MCMC.nIter, "+",
                                format(Starting.time, "%Y%m%d@%H.%M")),
                              collapse = "")

    gc()

    out <- as.list(environment())
    return(out)
}
