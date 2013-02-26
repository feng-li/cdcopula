#! /usr/bin/env Rscript
require("methods")

###----------------------------------------------------------------------------
### LOAD THE CODE LIBRARY AND INITIALIZE R ENVIRONMENT
###----------------------------------------------------------------------------

rm(list = ls()); gc()

## PATH TO THE ROOT DIRECTORY OF THE MODEL LIBRARY
pathLibRoot <- "~/workspace/copula/code/"

## Load the sourceDir tool
sys.source(file.path(pathLibRoot, "R/flutils/sourceDir.R"),
           envir = .GlobalEnv)

## Load the whole library
Cpl.source <- sourceDir(file.path(pathLibRoot, "R"),
                        byte.compile = 0,
                        recursive = TRUE,
                        ignore.error = TRUE)

## Load dependences
require("mvtnorm")

###----------------------------------------------------------------------------
### MCMC
###----------------------------------------------------------------------------

## PATH TO THE MODEL CONFIGURATION FILE
CplConfigFile <- file.path(pathLibRoot, "inst/config/config.main.sp100-600.R")
source(CplConfigFile, local = TRUE)

## Temporally Disabled for debugging mode
parallel <- FALSE
if(parallel == TRUE)
  {
    require(parallel)
    ## Use mcmlapply function
  }else
  {
    ## Sequential loops over folds.
    ## mapply(FUN = CplMain)
    CplMain.CrossOut <- list()
    for(iCross in 1:length(crossValidIdx[["training"]]))
      {
        CplMain.CrossOut[[iCross]] <- CplMain(
            CplConfigFile = CplConfigFile,
            Training.Idx = crossValidIdx[["training"]][[iCross]])
      }

  }

###----------------------------------------------------------------------------
### POSTERIOR INFERENCE, PREDICTION ETC
###----------------------------------------------------------------------------

## Temporally Disabled for debugging mode
parallel <- FALSE
if(parallel == TRUE)
  {
    require(parallel)
    ## Use mcmlapply function
  }else
  {
    ## Sequential loops over folds.
    ## mapply(FUN = CplMain)
    logPredLst <- list()
    for(iCross in 1:length(crossValidIdx[["testing"]]))
      {
        logPredLst[[iCross]] <- logPredDens(
            CplMain.out = CplMain.CrossOut[[iCross]],
            Testing.Idx = crossValidIdx[["testing"]][[iCross]])
      }
  }

LPDS <- logPredDensScore(logPredLst = logPredLst)

###----------------------------------------------------------------------------
### SAVE THE OUTPUT AND QUIT
###----------------------------------------------------------------------------
