#' This function provide end-user function to run a covariate-dependent copula model.
#'
#' An equivalent bash script is available with the same name.
#' @title CDCopula end-user running function.
#' @param Mdl.ConfigFile Model configuration file
#' @param CDCOPULA_LIB_ROOT_DIR cdcopula package path
#' @param nCores Number of cores used. If nCores>1, parallel mode is applied.
#' @return List
#' @author Feng Li
#' @export
CplRun = function(Mdl.ConfigFile, CDCOPULA_LIB_ROOT_DIR = system.file(package = "cdcopula"), nCores = 1)
{
###----------------------------------------------------------------------------
### LOAD THE CODE LIBRARY AND INITIALIZE R ENVIRONMENT
###----------------------------------------------------------------------------
    if(CDCOPULA_LIB_ROOT_DIR == system.file(package = "cdcopula"))
    {
        DEV = FALSE
    }
    else
    {
        DEV = TRUE
    }

    message("Load cdcopula library from ", CDCOPULA_LIB_ROOT_DIR)
    if(DEV == FALSE)
    {
        require("flutils")
        require("cdcopula")
    }
    else
    {
        ## Load the sourceDir tool
        sys.source(file.path(CDCOPULA_LIB_ROOT_DIR, "../flutils", "R/systools/sourceDir.R"),
                   envir = .GlobalEnv)

        ## Load the whole Copula library
        Cpl.source <- sourceDir(file.path(CDCOPULA_LIB_ROOT_DIR, "R"),
                                file.path(CDCOPULA_LIB_ROOT_DIR, "../flutils/R"),
                                byte.compile = 0,
                                recursive = TRUE,
                                ignore.error = TRUE)
    }
###----------------------------------------------------------------------------
### PRINT THE NATIVE MODEL DESCRIPTION
###----------------------------------------------------------------------------
    ## Extract variables from configure files
    Mdl.ConfigEnv <- new.env()
    sys.source(Mdl.ConfigFile, envir = Mdl.ConfigEnv)
    attach(Mdl.ConfigEnv)

    if(!exists("ForeignModelSpec"))
    {
        cat(rep("-", getOption("width")-1), "\n", sep = "")
        cat("MODEL DESCRIPTION\n")
        cat(rep("-", getOption("width")-1), "\n", sep = "")


        cat("nObsRaw:", nObsRaw, "\n")
        cat("nObsUsed:", length(Mdl.dataUsedIdx), "\n")

        cat("MCMC.nIter: ", MCMC.nIter,   "\n")
        cat("MCMC.burninProp:", MCMC.burninProp, "\n")
        cat("MCMC.UpateStrategy: ", MCMC.UpdateStrategy, "\n")

        ## cat("nCross: ", nCross,  "\n")
        cat("Mdl.crossValidArgs:\n")
        print(unlist(Mdl.crossValidArgs))

        ## cat("Mdl.MargisNM: ", Mdl.MargisNM,   "\n")
        ## cat("Vine.CopulaType:\n")
        ## print(Vine.CopulaType)

        ## cat("Vine.TreeStruc:\n")
        ## print(Vine.TreeStruc)

        cat("Mdl.MargisType:\n")
        print(Mdl.MargisType)

        if(all(rapply(Mdl.X, class)  == "matrix"))
        {
            cat("No. of covariates used (including intercept):\n")
            print(rapply(Mdl.X, ncol))
        }
        cat("Mdl.varSelArgs:\n")
        print(unlist(Mdl.varSelArgs))
        cat("MCMC.Update:\n")
        print(unlist(MCMC.Update))

        cat("MCMC.optimInit: ", MCMC.optimInit, "\n", sep = "")
        cat(rep("-", getOption("width")-1), "\n", sep = "")

        ## Remove unneeded variables
        ## rm(list = setdiff(setdiff(ls(), obj),
        ##                   c("nCross", "Mdl.crossValidIdx", "LPDS", "CplMCMC.summary")))
    }
###----------------------------------------------------------------------------
### SETUP PARALLEL COMPUTING ENVIRONMENT
###----------------------------------------------------------------------------
    if(nCross > 1 & nCores >1)
    {
        CDCOPULA_NPARALLEL <- floor(nCores/nCross)

        if(CDCOPULA_NPARALLEL<1)
        {
            stop("No sufficient many cores (",  nCores,  ") for running ",  nCross,
                 " cross-validations")
        }
    } else
    {
        CDCOPULA_NPARALLEL <- nCores
    }

    if(nCross>1 & nCores>1)
    {
        ## We have enough cores to do parallelized cross-validation
        require("parallel", quietly = TRUE)
        ## require("snow", quietly = TRUE)
        ## require("Rmpi", quietly = TRUE) # mpi.quit("no")
        ## cl4CV <- makeCluster(nCross, type = "MPI")
        if(nCross>detectCores())
        {
            stop("Parallel with SOCK method does not support this.")
        }

        cl4CV <- makeCluster(nCross)
        clusterExport(cl = cl4CV, envir = environment(),
                      c("DEV", "CDCOPULA_LIB_ROOT_DIR", "CDCOPULA_NPARALLEL", "sourceDir"))
        ce4CV <- clusterEvalQ(cl4CV,{
            Sys.setenv(DEV = DEV,
                       CDCOPULA_LIB_ROOT_DIR = CDCOPULA_LIB_ROOT_DIR,
                       CDCOPULA_NPARALLEL = CDCOPULA_NPARALLEL)

            if (DEV == TRUE)
            {
                sourceDir(file.path(CDCOPULA_LIB_ROOT_DIR, "R"),
                          file.path(CDCOPULA_LIB_ROOT_DIR, "../flutils/R"),
                          byte.compile = 0,
                          recursive = TRUE,
                          ignore.error = TRUE)
            } else
            {
                require("flutils")
                require("cdcopula")
            }
        })
    } else
    {
        Sys.setenv(DEV = DEV,
                   CDCOPULA_LIB_ROOT_DIR = CDCOPULA_LIB_ROOT_DIR,
                   CDCOPULA_NPARALLEL = CDCOPULA_NPARALLEL)
    }
###----------------------------------------------------------------------------
### RUN MCMC
###----------------------------------------------------------------------------
    ## Recording starting time
    Starting.time <- Sys.time()
    if(nCross>1 & nCores>1)
    {
        sink.parallel(cl = cl4CV)
        OUT.FITTED <- parLapply(cl4CV,
                                Mdl.crossValidIdx[["training"]],
                                CplMain,
                                Mdl.ConfigEnv = Mdl.ConfigEnv)
        sink.parallel(cl = cl4CV, file = NULL)
    } else
    {
        OUT.FITTED <- lapply(X = Mdl.crossValidIdx[["training"]],
                             FUN = CplMain,
                             Mdl.ConfigEnv = Mdl.ConfigEnv)
    }

###----------------------------------------------------------------------------
### POSTERIOR INFERENCE, PREDICTION ETC.
###----------------------------------------------------------------------------

    ## Final Summary of MCMC when paralleled.
    if(nCross>1 & nCores >1)
    {
        summary.Cpl <- list()
        for(iCross in 1:nCross)
        {   ## TODO: change the loop
            cat("iCross: ", iCross, "\n")
            summary.Cpl[[iCross]] <- CplMCMC.summary(OUT.FITTED[[iCross]])
        }
    }


    ## Calculate the LPDS
    cat("Calculating predictions...")
    OUT.PRED.LPDS <- list()
    LPDS <- OUT.FITTED[[1]][["Mdl.ConfigEnv"]][["LPDS"]]
    {
        if(nCross>1 & nCores >1)
        {
            Mdl.PredLst <- clusterMap(cl = cl4CV,
                                      fun = logDensPred,
                                      CplFitted = OUT.FITTED,
                                      Mdl.Idx.testing = Mdl.crossValidIdx[["testing"]])
        } else
        {
            ## browser()
            Mdl.PredLst <- mapply(FUN = logDensPred,
                                  CplFitted = OUT.FITTED,
                                  Mdl.Idx.testing = Mdl.crossValidIdx[["testing"]],
                                  SIMPLIFY = FALSE)
        }
        if(length(LPDS)>1 || !is.null(LPDS))
        {
            for(iLPDS in LPDS)
            {   ## TODO: parallel
                OUT.PRED.LPDS[[iLPDS]] <- logPredDensScore(
                    lapply(Mdl.PredLst, function(x) x[["Mdl.logPredDens"]][, iLPDS, drop = FALSE]))
            }
        }

        OUT.PRED.MVSK <- lapply(Mdl.PredLst, function(x) x[["Mdl.PredMVSK"]])
        OUT.PRED.RESID <- lapply(Mdl.PredLst, function(x) x[["Mdl.PredRESID"]])
        cat("done\n")
        print(OUT.PRED.LPDS)
###----------------------------------------------------------------------------
### FINISH THE WORK AND SAVE THE OUTPUT
###----------------------------------------------------------------------------
        ## Stop the parallel environment
        if(nCross>1 & nCores >1)
        {
            stopCluster(cl4CV)
        }

        save.all(save.output = Mdl.ConfigEnv[["save.output"]],
                 ModelDescription = OUT.FITTED[[1]][["ModelDescription"]])
    }
}
