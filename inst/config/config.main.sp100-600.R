###############################################################################
### TITLE
###        COPULA MODEL WITH COVARIATE--DEPENDENT
###
### SYSTEM REQUIREMENTS
###        R > = 2.15.3 with packages ``mvtnorm'',  ``parallel''
###        BLAS (optional)
###
### INPUT  VARIABLES
###        The variables with comments in capital letters are user input
###        variables.
###
### OUTPUT VARIABLES
###        The output variables are always started with ``MCMC.''
###
### GENERALIZATION
###        If you are interested in developing new copula models based on the
###        existing code. Look into the ``model'' folder.
###
### DATE
###        CREATED: Mon Jan 09 17:12:04 CET 2012
###        CURRENT: Tue Mar 05 17:22:35 CET 2013
###############################################################################

###----------------------------------------------------------------------------
### SPECIFY THE MODEL
###----------------------------------------------------------------------------

## SHORT MODEL DESCRIPTION
ModelDescription <- "bb7_copula_with_vs"

## COPULA DENSITY NAME AND PARAMETERS
CplNM <- "BB7"
CplParNM <- list(c("tau", "lambdaL"))

## MARGINAL MODELS NAME, TYPE AND PARAMETERS
MargisNM <- c("SP600", "SP100")
MargisTypes <- c("SPLITT", "SPLITT")
MargisParNM <- list(c("mu", "phi", "df", "lmd"),
                    c("mu", "phi", "df", "lmd"))

## Attribute name on the arguments
names(CplParNM) <- CplNM
names(MargisTypes) <- MargisNM
names(MargisParNM) <- MargisNM

## The object structure for the model components
MdlDataStruc <- initDataStruc(CplParNM, MargisParNM)

###----------------------------------------------------------------------------
### THE DATA AND MODEL
###----------------------------------------------------------------------------

## THE DATASET
##-----------------------------------------------------------------------------
## The dataset should either provided via DGP or the real dataset. The dataset
## should be in the following structure:
## Mdl.X: "list" each list contains the covariates in each margin or copula.
## Mdl.Y: "list" each list contains the response variable of that margin.

load(file.path(R_CPL_LIB_ROOT_DIR, "data/SP100-SP600-20130116.Rdata"))

## No. of Total Observations
nObsRaw <- length(Y[[1]])

## Data subset used
nObsIdx <- (0 + nObsRaw-nObsRaw):nObsRaw

## No. of used Observations
nObs <- length(nObsIdx)

## COVARIATES USED FOR THE MARGINAL AND COPULA PARAMETERS
Mdl.X <- MdlDataStruc
Mdl.X[[1]][["mu"]] <- cbind(1, X[[1]][, 1:9])[nObsIdx, 1:10, drop = FALSE]
Mdl.X[[1]][["phi"]] <- cbind(1, X[[1]][, 1:9])[nObsIdx, 1:10, drop = FALSE]
Mdl.X[[1]][["df"]] <- cbind(1, X[[1]][, 1:9])[nObsIdx, 1:10, drop = FALSE]
Mdl.X[[1]][["lmd"]] <- cbind(1, X[[1]][, 1:9])[nObsIdx, 1:10, drop = FALSE]

Mdl.X[[2]][["mu"]] <- cbind(1, X[[2]][, 1:9])[nObsIdx, 1:10, drop = FALSE]
Mdl.X[[2]][["phi"]] <- cbind(1, X[[2]][, 1:9])[nObsIdx, 1:10, drop = FALSE]
Mdl.X[[2]][["df"]] <- cbind(1, X[[2]][, 1:9])[nObsIdx, 1:10, drop = FALSE]
Mdl.X[[2]][["lmd"]] <- cbind(1, X[[2]][, 1:9])[nObsIdx, 1:10, drop = FALSE]

Mdl.X[[3]][["tau"]] <- cbind(1, X[[1]][, 1:9], X[[2]][, 1:9])[nObsIdx, 1:19, drop = FALSE]
Mdl.X[[3]][["lambdaL"]] <- cbind(1, X[[1]][, 1:9], X[[2]][, 1:9])[nObsIdx, 1:19, drop = FALSE]

## THE RESPONSE VARIABLES
Mdl.Y <- lapply(Y, function(x, idx)x[idx, ,drop = FALSE], nObsIdx)
names(Mdl.Y) <- MargisNM

## THE LINK FUNCTION USED IN THE MODEL
Mdl.parLink <- MdlDataStruc
Mdl.parLink[[1]][["mu"]] <- list(type = "identity")
Mdl.parLink[[1]][["phi"]] <- list(type = "log")
Mdl.parLink[[1]][["df"]] <- list(type = "glog", a = 4)
Mdl.parLink[[1]][["lmd"]] <- list(type = "log")


Mdl.parLink[[2]][["mu"]] <- list(type = "identity")
Mdl.parLink[[2]][["phi"]] <- list(type = "log")
Mdl.parLink[[2]][["df"]] <- list(type = "glog", a = 4)
Mdl.parLink[[2]][["lmd"]] <- list(type = "log")

Mdl.parLink[[3]][["tau"]] <- list(type = "glogit", a = 0.01, b = 0.99)
Mdl.parLink[[3]][["lambdaL"]] <- list(type = "glogit", a = 0.01, b = 0.99)

## THE VARIABLE SELECTION SETTINGS AND STARTING POINT
## Variable selection candidates, NULL: no variable selection use full
## covariates. ("all-in", "all-out", "random", or user-input)

varSelArgs <- MdlDataStruc
varSelArgs[[1]][["mu"]] <- list(cand = 2:10,
                                init = "all-in")
varSelArgs[[1]][["phi"]] <- list(cand = 2:10,
                                 init = "all-in")
varSelArgs[[1]][["df"]] <- list(cand = 2:10,
                                init = "all-in")
varSelArgs[[1]][["lmd"]] <- list(cand = 2:10,
                                 init = "all-in")

varSelArgs[[2]][["mu"]] <- list(cand = 2:10,
                                init = "all-in")
varSelArgs[[2]][["phi"]] <- list(cand = 2:10,
                                 init = "all-in")
varSelArgs[[2]][["df"]] <- list(cand = 2:10,
                                init = "all-in")
varSelArgs[[2]][["lmd"]] <- list(cand = 2:10,
                                 init = "all-in")

varSelArgs[[3]][["tau"]] <- list(cand = 2:19,
                                 init = "all-in")
varSelArgs[[3]][["lambdaL"]] <- list(cand = 2:19,
                                     init = "all-in")

###----------------------------------------------------------------------------
### THE MCMC CONFIGURATION
###----------------------------------------------------------------------------

## NUMBER OF MCMC ITERATIONS
nIter <- 100

## SAVE OUTPUT PATH
##-----------------------------------------------------------------------------
## "save.output = FALSE" it will not save anything.
## "save.output = "path-to-directory"" it will save the working directory in
## the given directory.
save.output <- FALSE

## MCMC TRAJECTORY
##-----------------------------------------------------------------------------
## If TRUE,  the MCMC should be tracked during the evaluation.
track.MCMC = TRUE

## WHICH VARIABLE SHOULD BE UPDATED?
MCMCUpdate <- MdlDataStruc
MCMCUpdate[[1]][[1]] <- T
MCMCUpdate[[1]][[2]] <- T
MCMCUpdate[[1]][[3]] <- T
MCMCUpdate[[1]][[4]] <- T

MCMCUpdate[[2]][[1]] <- T
MCMCUpdate[[2]][[2]] <- T
MCMCUpdate[[2]][[3]] <- T
MCMCUpdate[[2]][[4]] <- T

MCMCUpdate[[3]][[1]] <- T
MCMCUpdate[[3]][[2]] <- T

MCMCUpdateOrder <- MdlDataStruc
MCMCUpdateOrder[[1]][[1]] <- 1
MCMCUpdateOrder[[1]][[2]] <- 2
MCMCUpdateOrder[[1]][[3]] <- 3
MCMCUpdateOrder[[1]][[4]] <- 4

MCMCUpdateOrder[[2]][[1]] <- 5
MCMCUpdateOrder[[2]][[2]] <- 6
MCMCUpdateOrder[[2]][[3]] <- 7
MCMCUpdateOrder[[2]][[4]] <- 8

MCMCUpdateOrder[[3]][[1]] <- 9
MCMCUpdateOrder[[3]][[2]] <- 10

## THE METROPOLIS-HASTINGS ALGORITHM PROPOSAL ARGUMENTS
propArgs <- MdlDataStruc
propArgs[[1]][[1]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))
propArgs[[1]][[2]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))
propArgs[[1]][[3]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))
propArgs[[1]][[4]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))

propArgs[[2]][[1]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))
propArgs[[2]][[2]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))
propArgs[[2]][[3]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))
propArgs[[2]][[4]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))

propArgs[[3]][[1]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))
propArgs[[3]][[2]] <-
  list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
       "beta" = list(type = "mvt", df = 6),
       "indicators" = list(type = "binom", prob = 0.5))


## POSTERIOR INFERENCE OPTIONS
##-----------------------------------------------------------------------------

## CROSS VALIDATION
## "N.subsets" is no. of folds for cross-validation. If N.subsets  =  0, no
## cross-validation. And "partiMethod" tells how to partition the data. Testing
## percent is used if partiMethod is "time-series". (use the old data to
## predict the new interval)

nCross <- 1
crossValidArgs <- list(N.subsets = nCross,
                       partiMethod = "time-series",
                       testRatio = 0.001)

## Indices for training and testing sample according to cross-validation
crossValidIdx <- set.crossvalid(nObs,crossValidArgs)
## nCrossFold <- length(crossValidIdx[["training"]])

## SAMPLER PROPORTION FOR POSTERIOR INFERENCE,
sampleProp <- 0.8

## BURN-IN RATIO
burnin <- 0.1 # zero indicates no burn-in

###----------------------------------------------------------------------------
### PRIOR SETTINGS
###----------------------------------------------------------------------------

## PRIOR FOR THE COPULA PARAMETERS
##-----------------------------------------------------------------------------
## NOTE: The variable are recycled if needed. For example indicators$prob can
## be a scaler or a vector with same length of variable section
## candidates. There might be connections between parameters in the models but
## is will not affect the prior settings on the coefficients as long as we use
## a dynamic link function.

priArgs <- MdlDataStruc
priArgs[[1]][["mu"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "norm",  mean = 0, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[1]][["phi"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "lognorm",  mean = 1, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[1]][["df"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "glognorm",  mean = 5, variance = 10, a = 4),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[1]][["lmd"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "lognorm",  mean = 1, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))

priArgs[[2]][["mu"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "norm",  mean = 0, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][["phi"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "lognorm",  mean = 1, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][["df"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "glognorm",  mean = 5, variance = 10, a = 4),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][["lmd"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "lognorm",  mean = 1, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))

priArgs[[3]][["tau"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "gbeta",  mean = 0.2, variance = 0.05, a = 0.01, b = 0.95),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[3]][["lambdaL"]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "gbeta",  mean = 0.2, variance = 0.05, a = 0.01, b = 0.95),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "identity", shrinkage = 1)),
       "indicators" = list(type = "bern", prob = 0.5))

###----------------------------------------------------------------------------
### THE PARAMETERS FOR INITIAL AND CURRENT MCMC ITERATION
### The parameters in the current MCMC iteration. For the first iteration, it
### is set as the initial values
###----------------------------------------------------------------------------

## THE PARAMETER COEFFICIENTS STARTING POINT
## The possible inputs are ("random", "ols"  or user-input).
betaInit <- MdlDataStruc
betaInit[[1]][[1]] <- "random"
betaInit[[1]][[2]] <- "random"
betaInit[[1]][[3]] <- "random"
betaInit[[1]][[4]] <- "random"

betaInit[[2]][[1]] <- "random"
betaInit[[2]][[2]] <- "random"
betaInit[[2]][[3]] <- "random"
betaInit[[2]][[4]] <- "random"

betaInit[[3]][[1]] <- "random"
betaInit[[3]][[2]] <- "random"

## betaInit[[1]][[1]] <- "random"
## betaInit[[1]][[2]] <- "random"
## betaInit[[1]][[3]] <- "random"
## betaInit[[1]][[4]] <- -0.11706760

## betaInit[[2]][[1]] <- "random"
## betaInit[[2]][[2]] <- "random"
## betaInit[[2]][[3]] <- "random"
## betaInit[[2]][[4]] <- -0.11707829

## betaInit[[3]][[1]] <- -3.13519286
## betaInit[[3]][[2]] <- 0.60447399


 ## [1]  0.11535526 -0.07386206 -3.60276562 -0.11706760  0.14084234  0.02443366
 ## [7] -2.53630083 -0.11707829 -3.13519286  0.60447399

################################################################################
###                                  THE END
################################################################################
