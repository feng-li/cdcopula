#! /usr/bin/env Rscript
require("methods")

###----------------------------------------------------------------------------
### LOAD THE CODE LIBRARY AND INITIALIZE R ENVIRONMENT
###----------------------------------------------------------------------------

rm(list = ls()); gc()

## PATH TO THE ROOT DIRECTORY OF THE MODEL LIBRARY
pathLibRoot <- "~/workspace/copulas/CovarDepenCopula/R/"

## PATH TO THE MODEL CONFIGURATION FILE
configfile <- file.path(pathLibRoot, "config/config.main.sp100-600-n100.R")

## Load the sourceDir tool
sys.source(file.path(pathLibRoot, "flutils/stable/sourceDir.R"),
           envir = .GlobalEnv)

## Load the whole library
sourceDir(file.path(pathLibRoot, c("flutils/stable", "mcmc", "models", "simul")),
          byte.compile = 0, recursive = TRUE, ignore.error = TRUE)

## LOAD DEPENDENCES
require("mvtnorm")

###----------------------------------------------------------------------------
### MCMC
###----------------------------------------------------------------------------
CplMain(configfile)

###----------------------------------------------------------------------------
###
###----------------------------------------------------------------------------
