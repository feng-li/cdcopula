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
###        CURRENT: Sun Jan 04 10:21:56 CST 2015
###############################################################################


###----------------------------------------------------------------------------
### SPECIFY THE MODEL
###----------------------------------------------------------------------------
## MARGINAL MODELS NAME, TYPE AND PARAMETERS
MargisType <- c("SPLITT", "SPLITT", "BB7")
MargisNM <- c("M1", "M2", "BB7")

MCMCUpdate <- list(list("mu" = T, "phi" = T, "df" = T, "lmd" = T),
                   list("mu" = T, "phi" = T, "df" = T, "lmd" = T),
                   list("lambdaL" = T, "lambdaU" = T))

names(MCMCUpdate) <- MargisNM

## THE MODEL EVALUATION CRITERION
## Set this to NULL to turn of evaluation.
LPDS <- c("joint", MargisNM)

## The object structure for the model components
names(MargisType) <-  MargisNM

###----------------------------------------------------------------------------
### THE DATA AND MODEL
###----------------------------------------------------------------------------

## THE DATASET
##-----------------------------------------------------------------------------
## The dataset should either provided via DGP or the real dataset. The dataset
## should be in the following structure:
## Mdl.X: "list" each list contains the covariates in each margin or copula.
## Mdl.Y: "list" each list contains the response variable of that margin.
## Mdl.u: "list" each list contains the CDF for response variable of that margin.
DGPCpl(DGPconfigfile = file.path(R_CPL_LIB_ROOT_DIR,
                                 "inst/config/dgp/config.Cpl.DGP.R"),
       export = "parent.env")

## No. of Total Observations
nObsRaw <- nrow(MdlDGP.u[["u"]])

## Data subset used
nObsIdx <- (1 + nObsRaw-nObsRaw):nObsRaw

## No. of used Observations
nObs <- length(nObsIdx)

## THE RESPONSE VARIABLES (LOADED VIA DGP)

## COVARIATES USED FOR THE MARGINAL AND COPULA PARAMETERS ((LOADED VIA DGP))
## ------------------------------------------------------------------------------

## A trick to include foreign marginal models in the estimation which are hard to directly
## put into the "MargiModel()" is do the following settings: (1) Let "MCMCUpdate" be FALSE
## in all marginal densities.  (2) Estimate the density features in foreign models and set
## the features in "Mdl.X" directly.  (3) Set MCMCUpdateStrategy be "two-stage". (4) Set
## "betaInit" be one in all marginal features.

## THE LINK FUNCTION USED IN THE MODEL (LOADED VIA DGP)

## THE VARIABLE SELECTION SETTINGS AND STARTING POINT
## Variable selection candidates, NULL: no variable selection use full
## covariates. ("all-in", "all-out", "random", or user-input)

varSelArgs <- MCMCUpdate
varSelArgs[[1]][["mu"]] <- list(cand = "2:end", init = "all-in")
varSelArgs[[1]][["phi"]] <- list(cand = "2:end", init = "all-in")
varSelArgs[[1]][["df"]] <- list(cand = "2:end", init = "all-in")
varSelArgs[[1]][["lmd"]] <- list(cand = "2:end", init = "all-in")

varSelArgs[[2]][["mu"]] <- list(cand = "2:end", init = "all-in")
varSelArgs[[2]][["phi"]] <- list(cand = "2:end", init = "all-in")
varSelArgs[[2]][["df"]] <- list(cand = "2:end", init = "all-in")
varSelArgs[[2]][["lmd"]] <- list(cand = "2:end", init = "all-in")

varSelArgs[[3]][["lambdaL"]] <- list(cand = "2:end", init = "all-in")
varSelArgs[[3]][["lambdaU"]] <- list(cand = "2:end", init = "all-in")

###----------------------------------------------------------------------------
### THE MCMC CONFIGURATION
###----------------------------------------------------------------------------

## NUMBER OF MCMC ITERATIONS
MCMC.nIter <- 1000

## SAVE OUTPUT PATH
##-----------------------------------------------------------------------------
## "save.output = FALSE" it will not save anything.
## "save.output = "path-to-directory"" it will save the working directory in
## the given directory.
save.output <- "~/running"

## MCMC TRAJECTORY
##-----------------------------------------------------------------------------
## If TRUE,  the MCMC should be tracked during the evaluation.
MCMC.track <- TRUE

MCMCUpdateOrder <- MCMCUpdate
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

## MCMC UPDATING STRATEGY
##-----------------------------------------------------------------------------
## "joint"    : Update the joint posterior w.r.t. MCMCUpdate and MCMCUpdateOrder
## "margin"   : the marginal posterior.
## "twostage" : Update the joint posterior but using a two stage approach.

## NOTE: If one want to use "margin" or "two-stage" options just to to estimate the copula
## density. A variable "MCMC.density[["u"]]" must provide. "MCMC.density" consists of CDF of
## margins (i.e. u1,  u2, ...)

MCMCUpdateStrategy <- "joint"

## THE METROPOLIS-HASTINGS ALGORITHM PROPOSAL ARGUMENTS
propArgs <- MCMCUpdate
propArgs[[1]][[1]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))
propArgs[[1]][[2]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))
propArgs[[1]][[3]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))
propArgs[[1]][[4]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))

propArgs[[2]][[1]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))
propArgs[[2]][[2]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))
propArgs[[2]][[3]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))
propArgs[[2]][[4]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))

propArgs[[3]][[1]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
                           "beta" = list(type = "mvt", df = 6),
                           "indicators" = list(type = "binom", prob = 0.5))
propArgs[[3]][[2]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 3, hess = "outer"),
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
                       testRatio = 0.2)

## Indices for training and testing sample according to cross-validation
crossValidIdx <- set.crossvalid(nObs,crossValidArgs)
## nCrossFold <- length(crossValidIdx[["training"]])

## SAMPLER PROPORTION FOR POSTERIOR INFERENCE,
MCMC.sampleProp <- 0.8

## BURN-IN RATIO
MCMC.burninProp <- 0.1 # zero indicates no burn-in

###----------------------------------------------------------------------------
### PRIOR SETTINGS
###----------------------------------------------------------------------------

## PRIOR FOR THE COPULA PARAMETERS
## -----------------------------------------------------------------------------
## NOTE: The variable are recycled if needed. For example indicators$prob can be a scaler
## or a vector with same length of variable section candidates. There might be connections
## between parameters in the models but is will not affect the prior settings on the
## coefficients as long as we use a dynamic link function.

priArgs <- MCMCUpdate
priArgs[[1]][["mu"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                              input = list(type = "norm",  mean = 0, variance = 1),
                                                              output = list(type = "norm", shrinkage = 1)),
                                           "slopes" = list(type = "cond-mvnorm",
                                                           mean = 0, covariance = "identity", shrinkage = 1)),
                             "indicators" = list(type = "bern", prob = 0.5))
priArgs[[1]][["phi"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                               input = list(type = "lognorm",  mean = 1, variance = 1),
                                                               output = list(type = "norm", shrinkage = 1)),
                                            "slopes" = list(type = "cond-mvnorm",
                                                            mean = 0, covariance = "identity", shrinkage = 1)),
                              "indicators" = list(type = "bern", prob = 0.5))
priArgs[[1]][["df"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                              input = list(type = "glognorm",  mean = 5, variance = 10, a = 4),
                                                              output = list(type = "norm", shrinkage = 1)),
                                           "slopes" = list(type = "cond-mvnorm",
                                                           mean = 0, covariance = "identity", shrinkage = 1)),
                             "indicators" = list(type = "bern", prob = 0.5))
priArgs[[1]][["lmd"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                               input = list(type = "lognorm",  mean = 1, variance = 1),
                                                               output = list(type = "norm", shrinkage = 1)),
                                            "slopes" = list(type = "cond-mvnorm",
                                                            mean = 0, covariance = "identity", shrinkage = 1)),
                              "indicators" = list(type = "bern", prob = 0.5))

priArgs[[2]][["mu"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                              input = list(type = "norm",  mean = 0, variance = 1),
                                                              output = list(type = "norm", shrinkage = 1)),
                                           "slopes" = list(type = "cond-mvnorm",
                                                           mean = 0, covariance = "identity", shrinkage = 1)),
                             "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][["phi"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                               input = list(type = "lognorm",  mean = 1, variance = 1),
                                                               output = list(type = "norm", shrinkage = 1)),
                                            "slopes" = list(type = "cond-mvnorm",
                                                            mean = 0, covariance = "identity", shrinkage = 1)),
                              "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][["df"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                              input = list(type = "glognorm",  mean = 5, variance = 10, a = 4),
                                                              output = list(type = "norm", shrinkage = 1)),
                                           "slopes" = list(type = "cond-mvnorm",
                                                           mean = 0, covariance = "identity", shrinkage = 1)),
                             "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][["lmd"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                               input = list(type = "lognorm",  mean = 1, variance = 1),
                                                               output = list(type = "norm", shrinkage = 1)),
                                            "slopes" = list(type = "cond-mvnorm",
                                                            mean = 0, covariance = "identity", shrinkage = 1)),
                              "indicators" = list(type = "bern", prob = 0.5))

priArgs[[3]][["lambdaL"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                                   input = list(type = "gbeta",  mean = 0.2, variance = 0.05,
                                                                                a = 0.01, b = 0.99),
                                                                   output = list(type = "norm", shrinkage = 1)),
                                                "slopes" = list(type = "cond-mvnorm",
                                                                mean = 0, covariance = "identity", shrinkage = 1)),
                                  "indicators" = list(type = "bern", prob = 0.5))
priArgs[[3]][["lambdaU"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                                   input = list(type = "gbeta",  mean = 0.2, variance = 0.05,
                                                                                a = 0.01, b = 0.99),
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
betaInit <- MCMCUpdate
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

optimInit <- TRUE
################################################################################
###                                  THE END
################################################################################
