#! /usr/bin/env Rscript
require("methods")

###----------------------------------------------------------------------------
### LOAD THE CODE LIBRARY AND INITIALIZE R ENVIRONMENT
###----------------------------------------------------------------------------

rm(list = ls()); gc()

## PATH TO THE ROOT DIRECTORY OF THE MODEL LIBRARY
pathLibRoot <- "~/workspace/copula/code/"

## PATH TO THE MODEL CONFIGURATION FILE
CplConfigFile <- file.path(pathLibRoot, "inst/config/config.main.sp100-600.R")

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
CplMain(configfile = CplConfigFile)

###----------------------------------------------------------------------------
### POSTERIOR INFERENCE, PREDICTION ETC
###----------------------------------------------------------------------------

## LPDS
logPredMatrix <-

LogPredDensScore(logPredMatrix = logPredMatrix)

###----------------------------------------------------------------------------
### SAVE THE OUTPUT AND QUIT
###----------------------------------------------------------------------------
