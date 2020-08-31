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
Mdl.MargisType <- c("STOCHVOL", "STOCHVOL", "MVT")
Mdl.MargisNM <- c("^SML", "^OEX", "MVT")

MCMC.Update <- list(list("mu" = F, "phi" = F),
                   list("mu" = F, "phi" = F),
                   list("df" = T, "rho" = T))


names(MCMC.Update) <- Mdl.MargisNM

## THE MODEL EVALUATION CRITERION
## Set this to NULL to turn of evaluation.
LPDS <- c("joint", Mdl.MargisNM)

## The object structure for the model components
names(Mdl.MargisType) <-  Mdl.MargisNM

###----------------------------------------------------------------------------
### THE DATA AND MODEL
###----------------------------------------------------------------------------

## THE DATASET
##-----------------------------------------------------------------------------
## The dataset should either provided via DGP or the real dataset. The dataset
## should be in the following structure:
## Mdl.X: "list" each list contains the covariates in each margin or copula.
## Mdl.Y: "list" each list contains the response variable of that margin.

load(file.path(CDCOPULA_LIB_ROOT_DIR, "data/SP100-SP400-SP600-20150206.Rdata"))

## No. of Total Observations
nObsRaw <- length(Y[[1]])

## Data subset used
Mdl.dataUsedIdx <- (1 + nObsRaw-nObsRaw):nObsRaw
Mdl.algorithm <- "full"

## THE RESPONSE VARIABLES
Mdl.Y <- lapply(Y[Mdl.MargisNM[-length(Mdl.MargisNM)]], function(x, idx)x[idx, ,drop = FALSE], Mdl.dataUsedIdx)

## The name of respond variables
names(Mdl.Y) <- Mdl.MargisNM[-length(Mdl.MargisNM)]

## COVARIATES USED FOR THE MARGINAL AND COPULA PARAMETERS
## ------------------------------------------------------------------------------

## A trick to include foreign marginal models in the estimation which are hard to directly
## put into the "MargiModel()" is do the following settings: (1) Let "MCMC.Update" be FALSE
## in all marginal densities.  (2) Estimate the density features in foreign models and set
## the features in "Mdl.X" directly.  (3) Set MCMC.UpdateStrategy be "twostage". (4) Set
## "Mdl.betaInit" be one in all marginal features.
Mdl.X <- MCMC.Update


Mdl.X[[1]] <- list(draws  =  1000,
                   burnin  =  100,
                   designmatrix  =  NA,
                   priormu  =  c(0,  100),
                   priorphi  =  c(5,  1.5),
                   priorsigma  =  1,
                   priornu  =  NA,
                   priorbeta  =  c(0,  10000),
                   thinpara  =  1,
                   thinlatent  =  1,
                   thintime  =  1,
                   keeptau  =  FALSE,
                   quiet  =  FALSE)

Mdl.X[[2]] <- list(draws  =  1000,
                   burnin  =  100,
                   designmatrix  =  NA,
                   priormu  =  c(0,  100),
                   priorphi  =  c(5,  1.5),
                   priorsigma  =  1,
                   priornu  =  NA,
                   priorbeta  =  c(0,  10000),
                   thinpara  =  1,
                   thinlatent  =  1,
                   thintime  =  1,
                   keeptau  =  FALSE,
                   quiet  =  FALSE)

Mdl.X[[3]][["df"]] <- cbind(1, X[[1]][Mdl.dataUsedIdx, 1:9], X[[2]][Mdl.dataUsedIdx, 1:9])#[, 1, drop = FALSE]
Mdl.X[[3]][["rho"]] <- cbind(1, X[[1]][Mdl.dataUsedIdx, 1:9], X[[2]][Mdl.dataUsedIdx, 1:9])#[, 1, drop = FALSE]

## THE LINK FUNCTION USED IN THE MODEL
Mdl.parLink <- MCMC.Update
Mdl.parLink[[1]][["mu"]] <- list(type = "identity", nPar = 1)
Mdl.parLink[[1]][["phi"]] <- list(type = "identity", nPar = 1)

Mdl.parLink[[2]][["mu"]] <- list(type = "identity", nPar = 1)
Mdl.parLink[[2]][["phi"]] <- list(type = "identity", nPar = 1)

Mdl.parLink[[3]][["df"]] <- list(type = "glog", nPar = 1, a = 2, b = 30)
Mdl.parLink[[3]][["rho"]] <- list(type = "glogit", nPar = 1, a = 0.01, b = 0.99)

## THE VARIABLE SELECTION SETTINGS AND STARTING POINT
## Variable selection candidates, NULL: no variable selection use full
## covariates. ("all-in", "all-out", "random", or user-input)

Mdl.varSelArgs <- MCMC.Update
Mdl.varSelArgs[[1]][["mu"]] <- list(cand = NULL, init = "all-in")
Mdl.varSelArgs[[1]][["phi"]] <- list(cand = NULL, init = "all-in")

Mdl.varSelArgs[[2]][["mu"]] <- list(cand = NULL, init = "all-in")
Mdl.varSelArgs[[2]][["phi"]] <- list(cand = NULL, init = "all-in")

Mdl.varSelArgs[[3]][["df"]] <- list(cand = NULL, init = "all-in")
Mdl.varSelArgs[[3]][["rho"]] <- list(cand = NULL, init = "all-in")

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

MCMC.UpdateOrder <- MCMC.Update
MCMC.UpdateOrder[[1]][[1]] <- 1
MCMC.UpdateOrder[[1]][[2]] <- 2

MCMC.UpdateOrder[[2]][[1]] <- 3
MCMC.UpdateOrder[[2]][[1]] <- 4

MCMC.UpdateOrder[[3]][[1]] <- 5
MCMC.UpdateOrder[[3]][[2]] <- 6

## MCMC UPDATING STRATEGY
##-----------------------------------------------------------------------------
## "joint"    : Update the joint posterior w.r.t. MCMC.Update and MCMC.UpdateOrder
## "margin"   : the marginal posterior.
## "twostage" : Update the joint posterior but using a two stage approach.

## NOTE: If one want to use "margin" or "twostage" options just to to estimate the copula
## density. A variable "MCMC.density[["u"]]" must provide. "MCMC.density" consists of CDF of
## margins (i.e. u1,  u2, ...)

MCMC.UpdateStrategy <- "twostage"

## THE METROPOLIS-HASTINGS ALGORITHM PROPOSAL ARGUMENTS
MCMC.propArgs <- MCMC.Update
MCMC.propArgs[[1]][[1]] <- NA
MCMC.propArgs[[1]][[2]] <- NA

MCMC.propArgs[[2]][[1]] <- NA
MCMC.propArgs[[2]][[2]] <- NA

MCMC.propArgs[[3]][[1]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
                                "beta" = list(type = "mvt", df = 6),
                                "indicators" = list(type = "binom", prob = 0.5))
MCMC.propArgs[[3]][[2]] <- list("algorithm" = list(type = "GNewtonMove", ksteps = 1, hess = "outer"),
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
Mdl.crossValidArgs <- list(N.subsets = nCross, partiMethod = "time-series", testRatio = 0.2)

## Indices for training and testing sample according to cross-validation
Mdl.crossValidIdx <- set.crossvalid(length(Mdl.dataUsedIdx),Mdl.crossValidArgs)
## nCrossFold <- length(Mdl.crossValidIdx[["training"]])

## SAMPLER PROPORTION FOR POSTERIOR INFERENCE,
MCMC.sampleProp <- 1

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

Mdl.priArgs <- MCMC.Update

Mdl.priArgs[[1]][["mu"]] <- NA
Mdl.priArgs[[1]][["phi"]] <- NA

Mdl.priArgs[[2]][["mu"]] <- NA
Mdl.priArgs[[2]][["phi"]] <- NA


Mdl.priArgs[[3]][["df"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                                  input = list(type = "glognorm",  mean = 5, variance = 10, a = 2),
                                                                  output = list(type = "norm", shrinkage = 1)),
                                               "slopes" = list(type = "cond-mvnorm",
                                                               mean = 0, covariance = "identity", shrinkage = 1)),
                                 "indicators" = list(type = "bern", prob = 0.5))
Mdl.priArgs[[3]][["rho"]] <- list("beta" = list("intercept" = list(type = "custom",
                                                                   input = list(type = "gbeta",  mean = 0.5, variance = 0.05,
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
Mdl.betaInit <- MCMC.Update
Mdl.betaInit[[1]][[1]] <- 1
Mdl.betaInit[[1]][[2]] <- 1

Mdl.betaInit[[2]][[1]] <- 1
Mdl.betaInit[[2]][[2]] <- 1

Mdl.betaInit[[3]][[1]] <- "random"
Mdl.betaInit[[3]][[2]] <- "random"

MCMC.optimInit <- 0.1
################################################################################
###                                  THE END
################################################################################
