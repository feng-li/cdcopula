###############################################################################
### Configuration file for the copula data generating process
###############################################################################

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

## NO. OF OBSERVATIONS
nObs <- 15

## THE LINK FUNCTION USED IN THE MODEL
MdlDGP.parLink <- MdlDataStruc
MdlDGP.parLink[[1]][[1]] <- "identity"
MdlDGP.parLink[[1]][[2]] <- "log"
MdlDGP.parLink[[2]][[1]] <- "identity"
MdlDGP.parLink[[2]][[2]] <- "log"
MdlDGP.parLink[[3]][[1]] <- "logit"
MdlDGP.parLink[[3]][[2]] <- "logit"

##-----------------------------------------------------------------------------
## THE TRUE PARAMETER VALUES IN THE DGP
## -----------------------------------------------------------------------------
## The parameters are the features of the model, e.g. mean, variance, ...  The
## parameters are observation specified (each observation has its own
## features.) This is the desired feature, the TRUE DGP might be not exactly
## the same as it. See "Mdl.par" for the comparison.
## -----------------------------------------------------------------------------
MdlDGP.par <- MdlDataStruc

## The first margin
MdlDGP.par[[1]][[1]] <- matrix(rnorm(n = nObs, mean = 0, sd = 1))
MdlDGP.par[[1]][[2]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1))

## The second margin
MdlDGP.par[[2]][[1]] <- matrix(rnorm(n = nObs, mean = 0, sd = 1))
MdlDGP.par[[2]][[2]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1))

## The copula
MdlDGP.par[[3]][[1]] <- matrix(rbeta2(n = nObs, mean = 0.7, sd = 0.1))
MdlDGP.par[[3]][[2]] <- matrix(rbeta2(n = nObs, mean = 0.3, sd = 0.1))

##------------------------------------------------------------------------------
## THE TRUE COVARIATE-DEPENDENT PARAMETER VALUES IN THE DGP
##------------------------------------------------------------------------------
## Note that (a) The intercept is always included.  (b) Variable selection is
## represented if the beta coefficient is zero.  (c) If there are no covariates
## dependent, only intercept should be set.
##------------------------------------------------------------------------------

## INTERCEPT INDICATOR
## If "TRUE",  the intercept should be in the covariate-dependent parameter
## structure.
MdlDGP.intercept <- MdlDataStruc

MdlDGP.intercept[[1]][[1]] <- TRUE
MdlDGP.intercept[[1]][[2]] <- TRUE
MdlDGP.intercept[[2]][[1]] <- TRUE
MdlDGP.intercept[[2]][[2]] <- TRUE
MdlDGP.intercept[[3]][[1]] <- TRUE
MdlDGP.intercept[[3]][[2]] <- TRUE

## THE COEFFICIENTS
## When the coefficient are not NA, the parameter are fixed. Otherwise it was
## determined by the spline covariates. If the intercept is included, the first
## entry should always be "NA".
## MdlDGP.beta <- MdlDataStruc

## ## The first margin
## MdlDGP.beta[[1]][[1]] <- matrix(c(NA, 1,  -1,  NA, NA, NA))
## MdlDGP.beta[[1]][[2]] <- matrix(c(NA, 1,  -1,  NA, NA, NA))

## ## The second margin
## MdlDGP.beta[[2]][[1]] <- matrix(c(NA, NA, NA,  NA, NA))
## MdlDGP.beta[[2]][[2]] <- matrix(c(NA, NA, NA,  NA, NA))

## ## The copula
## MdlDGP.beta[[3]][[1]] <- matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
## MdlDGP.beta[[3]][[2]] <- matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))

## NUMBER OF COVARIATES (EXCLUDING INTERCEPT)
MdlDGP.nCovs <- MdlDataStruc

## The first margin
MdlDGP.nCovs[[1]][[1]] <- list(total = 4, fixed = 2)
MdlDGP.nCovs[[1]][[2]] <- list(total = 4, fixed = 2)

## The second margin
MdlDGP.nCovs[[2]][[1]] <- list(total = 5, fixed = 2)
MdlDGP.nCovs[[2]][[2]] <- list(total = 5, fixed = 2)

## The copula
MdlDGP.nCovs[[3]][[1]] <- list(total = 10, fixed = 4)
MdlDGP.nCovs[[3]][[2]] <- list(total = 10, fixed = 4)

## Generating the numerical tabular for the inverse Kendall's tau
tauTabular <- kendalltauTabular(CplNM = CplNM, tol = 0.005)
