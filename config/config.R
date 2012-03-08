###############################################################################
### TITLE
###        BIVARIATE COPULA MODEL WITH COVARIATE--DEPENDENT
###
### SYSTEM REQUIREMENTS
###        R > = 2.14.0 with packages ``mvtnorm'' 
###        BLAS (optional)
###
### INPUT  VARIABLES
###        The variables with comments in capital letters are user input
###        variables.
###
### OUTPUT VARIABLES
###        The output variables are always started with ``OUT.''
###
### GENERALIZATION 
###        If you are interested in developing new copula models based on the
###        existing code. Look into the ``model'' folder.     
###
### DATE
###        CREATED: Mon Jan 09 17:12:04 CET 2012
###        CURRENT: Mon Jan 23 15:20:01 CET 2012
###############################################################################

###----------------------------------------------------------------------------
### LOAD THE CODE LIBRARY AND INITIALIZE R ENVIRONMENT
###----------------------------------------------------------------------------

## rm(list = ls()); gc()

## PATH TO THE ROOT DIRECTORY OF THE MODEL LIBRARY
pathLibRoot <- "~/workspace/copulas/CovarDepenCopula/R/"

## Load the sourceDir tool
sys.source(file.path(pathLibRoot, "rutils/stable/sourceDir.R"),
           envir = attach(NULL, name = "sourceDir"))

## Load the whole library
sourceDir(file.path(pathLibRoot, c("rutils/stable", "mcmc", "models", "simul")),
          byte.compile = FALSE, ignore.error = TRUE)        

## LOAD DEPENDENCES
require("mvtnorm")

###----------------------------------------------------------------------------
### SPECIFY THE MODEL
###----------------------------------------------------------------------------

## SHORT MODEL DESCRIPTION
ModelDescription <- "bb7_copula_no_variable_selection"

## COPULA DENSITY NAME AND PARAMETERS
CplNM <- "BB7"
CplParNM <- list(c("tau", "lambdaL"))

## MARGINAL MODELS NAME, TYPE AND PARAMETERS
MargisNM <- c("SP500", "NASDAQ100")
MargisTypes <- c("GAUSSIAN", "GAUSSIAN")
MargisParNM <- list(c("mu", "sigma"), 
                    c("mu", "sigma"))

## Attribute name on the arguments 
names(CplParNM) <- CplNM
names(MargisTypes) <- MargisNM
names(MargisParNM) <- MargisNM

## The object structure for the model components
MdlDataStruc <- initDataStruc(CplParNM, MargisParNM)

## Generating the numerical tabular for the inverse Kendall's tau
tauTabular <- kendalltauTabular(CplNM = CplNM, tol = 0.005)

###----------------------------------------------------------------------------
### THE DATA AND MODEL
###----------------------------------------------------------------------------

## THE DATASET
DGPCpl(configfile = "config/config.DGPCpl.R", export = "list")

## COVARIATES USED FOR THE MARGINAL AND COPULA PARAMETERS
Mdl.X <- MdlDataStruc
Mdl.X[[1]][[1]] <- cbind(1, X[[1]])
Mdl.X[[1]][[2]] <- cbind(1, X[[1]])
Mdl.X[[2]][[1]] <- cbind(1, X[[2]])
Mdl.X[[2]][[2]] <- cbind(1, X[[2]])
Mdl.X[[3]][[1]] <- cbind(1, X[[1]], X[[2]])
Mdl.X[[3]][[2]] <- cbind(1, X[[1]], X[[2]])


## THE LINK FUNCTION USED IN THE MODEL
Mdl.parLink <- MdlDataStruc
Mdl.parLink[[1]][[1]] <- "identity"
Mdl.parLink[[1]][[2]] <- "log"
Mdl.parLink[[2]][[1]] <- "identity"
Mdl.parLink[[2]][[2]] <- "log"
Mdl.parLink[[3]][[1]] <- "logit"
Mdl.parLink[[3]][[2]] <- "logit"


## THE VARIABLE SELECTION SETTINGS AND STARTING POINT
## Variable selection candidates, NULL: no variable selection use full
## covariates. ("all-in", "all-out", "random", or user-input)

varSelArgs <- MdlDataStruc
varSelArgs[[1]][[1]] <- list(cand = c(2, 3),
                             init = "all-in") 
varSelArgs[[1]][[2]] <- list(cand = c(2, 3),
                             init = "all-out")
varSelArgs[[2]][[1]] <- list(cand = c(2, 4),
                             init = "random")
varSelArgs[[2]][[2]] <- list(cand = c(2, 4),
                             init = "all-out")
varSelArgs[[3]][[1]] <- list(cand = c(2, 3, 5, 6),
                             init = c(2, 3))
varSelArgs[[3]][[2]] <- list(cand = c(3, 5, 6),
                             init = "random")

###----------------------------------------------------------------------------
### THE MCMC CONFIGURATION 
###----------------------------------------------------------------------------

## NUMBER OF MCMC ITERATIONS 
nIter <- 20

## BURN-IN RATIO
burnin <- 0.1 # zero indicates no burn-in

## SAVE OUTPUT PATH
save.output <- FALSE # "save.output = FALSE" will not save anything

## MCMC TRAJECTORY
track.MCMC = TRUE

## WHICH VARIABLE SHOULD BE UPDATED?
MCMCUpdate <- MdlDataStruc
MCMCUpdate[[1]][[1]] <- TRUE
MCMCUpdate[[1]][[2]] <- TRUE
MCMCUpdate[[2]][[1]] <- TRUE
MCMCUpdate[[2]][[2]] <- TRUE
MCMCUpdate[[3]][[1]] <- TRUE
MCMCUpdate[[3]][[2]] <- TRUE

## THE METROPOLIS-HASTINGS ALGORITHM PROPOSAL ARGUMENTS 
propArgs <- MdlDataStruc
propArgs[[1]][[1]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"), 
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.2)) 
propArgs[[1]][[2]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"), 
       "beta" = list(type = "mvt", df = 6), 
       "indicators" = list(type = "binom", prob = 0.2))
propArgs[[2]][[1]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"), 
       "beta" = list(type = "mvt", df = 6), 
       "indicators" = list(type = "binom", prob = 0.2))
propArgs[[2]][[2]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"), 
       "beta" = list(type = "mvt", df = 6), 
       "indicators" = list(type = "binom", prob = 0.2))
propArgs[[3]][[1]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"), 
       "beta" = list(type = "mvt", df = 6), 
       "indicators" = list(type = "binom", prob = 0.2))
propArgs[[3]][[2]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
       "beta" = list(type = "mvt", df = 6), 
       "indicators" = list(type = "binom", prob = 0.2))

## CROSS-VALIDATION
## If N.subsets  =  0, no cross-validation 
crossValidArgs <- list(N.subsets = 1,    # Folds for cross-validation. 
                       partiMethod = "time-series", # How to partition the data
                       testRatio = 0.2)   # Testing percent if "time-series"

## Indices for training and testing sample according to cross-validation
## settings 
crossValidIdx <- set.crossvalid(nObs,crossValidArgs)
nCrossFold <- length(crossValidIdx[["training"]])

## SAMPLER PROPORTION FOR LPDS
LPDS.sampleProp = 0.05          

###----------------------------------------------------------------------------
### PRIOR SETTINGS
###----------------------------------------------------------------------------

## PRIOR FOR THE COPULA PARAMETERS
## NOTE: The variable are recycled if needed. For example
## indicators$prob can be a scaler or a vector with same length of variable
## section candidates.

priArgs <- MdlDataStruc
priArgs[[1]][[1]] <- 
  list("beta" = list(
         "intercept" = list(type = "custom", 
           input = list(type = "norm",  mean = 0.5, variance = 1), 
           output = list(type = "norm", shrinkage = 1)), 
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = 1)), 
       "indicators" = list(type = "bern", 
         prob = 0.5))
priArgs[[1]][[2]] <- 
  list("beta" = list(
         "intercept" = list(type = "custom", 
           input = list(type = "lognorm",  mean = 0.5, variance = 1), 
           output = list(type = "norm", shrinkage = 1)), 
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = 1)), 
       "indicators" = list(type = "bern", 
         prob = 0.5))
priArgs[[2]][[1]] <- 
  list("beta" = list(
         "intercept" = list(type = "custom", 
           input = list(type = "norm",  mean = 0.5, variance = 1), 
           output = list(type = "norm", shrinkage = 1)), 
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = 1)), 
       "indicators" = list(type = "bern", 
         prob = 0.5))
priArgs[[2]][[2]] <- 
  list("beta" = list(
         "intercept" = list(type = "custom", 
           input = list(type = "lognorm",  mean = 0.5, variance = 1), 
           output = list(type = "norm", shrinkage = 1)), 
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = 1)), 
       "indicators" = list(type = "bern", 
         prob = 0.5))
priArgs[[3]][[1]] <- 
  list("beta" = list(
         "intercept" = list(type = "custom", 
           input = list(type = "beta",  mean = 0.5, variance = 1), 
           output = list(type = "norm", shrinkage = 1)), 
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = 1)), 
       "indicators" = list(type = "bern", 
         prob = 0.5))
priArgs[[3]][[2]] <- 
  list("beta" = list(
         "intercept" = list(type = "custom", 
           input = list(type = "beta",  mean = 0.5, variance = 1), 
           output = list(type = "norm", shrinkage = 1)), 
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = 1)), 
       "indicators" = list(type = "bern", 
         prob = 0.5))

###----------------------------------------------------------------------------
### THE PARAMETERS FOR INITIAL AND CURRENT MCMC ITERATION
### The parameters in the current MCMC iteration. For the first iteration, it
### is set as the initial values
###----------------------------------------------------------------------------

## THE PARAMETER COEFFICIENTS STARTING POINT
## ("random", or user-input)

betaInit <- MdlDataStruc
betaInit[[1]][[1]] <- "random"
betaInit[[1]][[2]] <- "random"
betaInit[[2]][[1]] <- "random"
betaInit[[2]][[2]] <- "random"
betaInit[[3]][[1]] <- "random"
betaInit[[3]][[2]] <- "random"

################################################################################
###                                  THE END
################################################################################
