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
###        CURRENT: Mon Apr 23 18:00:22 CEST 2012
###############################################################################

###----------------------------------------------------------------------------
### SPECIFY THE MODEL
###----------------------------------------------------------------------------

## SHORT MODEL DESCRIPTION
ModelDescription <- "bb7_copula_with_variable_selection"

## COPULA DENSITY NAME AND PARAMETERS
CplNM <- "BB7"
CplParNM <- list(c("tau", "lambdaL"))

## MARGINAL MODELS NAME, TYPE AND PARAMETERS
MargisNM <- c("SP600", "SP100")
MargisTypes <- c("GAUSSIAN", "GAUSSIAN")
MargisParNM <- list(c("mu", "sigma"),
                    c("mu", "sigma"))

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

load(file.path(pathLibRoot, "data/SP100-SP600-n100.Rdata"))

## COVARIATES USED FOR THE MARGINAL AND COPULA PARAMETERS
Mdl.X <- MdlDataStruc
Mdl.X[[1]][[1]] <- cbind(1, X[[1]][, 1:3])
Mdl.X[[1]][[2]] <- cbind(1, X[[1]][, 1:3])
Mdl.X[[2]][[1]] <- cbind(1, X[[2]][, 1:3])
Mdl.X[[2]][[2]] <- cbind(1, X[[2]][, 1:3])
Mdl.X[[3]][[1]] <- cbind(1, X[[1]][, 1:3], X[[2]][, 1:3])
Mdl.X[[3]][[2]] <- cbind(1, X[[1]][, 1:3], X[[2]][, 1:3])

## THE RESPONSE VARIABLES
Mdl.Y <- Y

## THE LINK FUNCTION USED IN THE MODEL
Mdl.parLink <- MdlDataStruc
Mdl.parLink[[1]][[1]] <- "identity"
Mdl.parLink[[1]][[2]] <- "log"
Mdl.parLink[[2]][[1]] <- "identity"
Mdl.parLink[[2]][[2]] <- "log"
Mdl.parLink[[3]][[1]] <- "glogit"
Mdl.parLink[[3]][[2]] <- "logit"

## THE VARIABLE SELECTION SETTINGS AND STARTING POINT
## Variable selection candidates, NULL: no variable selection use full
## covariates. ("all-in", "all-out", "random", or user-input)

varSelArgs <- MdlDataStruc
varSelArgs[[1]][[1]] <- list(cand = c(2, 3),
                             init = "all-in")
varSelArgs[[1]][[2]] <- list(cand = c(2, 3),
                             init = "all-out")
varSelArgs[[2]][[1]] <- list(cand = c(2, 3),
                             init = "random")
varSelArgs[[2]][[2]] <- list(cand = c(2, 3),
                             init = "all-out")
varSelArgs[[3]][[1]] <- list(cand = c(2, 3, 4),
                             init = c(2, 3))
varSelArgs[[3]][[2]] <- list(cand = c(2, 4),
                             init = "random")

###----------------------------------------------------------------------------
### THE MCMC CONFIGURATION
###----------------------------------------------------------------------------

## NUMBER OF MCMC ITERATIONS
nIter <- 20

## BURN-IN RATIO
burnin <- 0.1 # zero indicates no burn-in

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
MCMCUpdate[[1]][[1]] <- F
MCMCUpdate[[1]][[2]] <- F
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
##-----------------------------------------------------------------------------
## "N.subsets" is no. of folds for cross-validation. If N.subsets  =  0, no
## cross-validation. And "partiMethod" tells how to partition the data. Testing
## percent is used if partiMethod is "time-series". (use the old data to
## predict the new interval)
crossValidArgs <- list(N.subsets = 1,
                       partiMethod = "time-series",
                       testRatio = 0.2)

## SAMPLER PROPORTION FOR LPDS
LPDS.sampleProp = 0.05

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

nObs <- length(Mdl.Y[[1]])

priArgs <- MdlDataStruc
priArgs[[1]][[1]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "norm",  mean = 0, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = nObs)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[1]][[2]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "lognorm",  mean = 1, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = nObs)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][[1]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "norm",  mean = 0, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = nObs)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[2]][[2]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "lognorm",  mean = 1, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = nObs)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[3]][[1]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "gbeta",  mean = 0.5, variance = 1, a = 0.1, b = 0.3),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = nObs)),
       "indicators" = list(type = "bern", prob = 0.5))
priArgs[[3]][[2]] <-
  list("beta" = list(
         "intercept" = list(type = "custom",
           input = list(type = "beta",  mean = 0.5, variance = 1),
           output = list(type = "norm", shrinkage = 1)),
         "slopes" = list(type = "cond-mvnorm",
           mean = 0, covariance = "g-prior", shrinkage = nObs)),
       "indicators" = list(type = "bern", prob = 0.5))

###----------------------------------------------------------------------------
### THE PARAMETERS FOR INITIAL AND CURRENT MCMC ITERATION
### The parameters in the current MCMC iteration. For the first iteration, it
### is set as the initial values
###----------------------------------------------------------------------------

## THE PARAMETER COEFFICIENTS STARTING POINT
## The possible inputs are ("random", or user-input).

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
