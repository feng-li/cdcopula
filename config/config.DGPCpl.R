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
Mdl.parLink <- MdlDataStruc
Mdl.parLink[[1]][[1]] <- "identity"
Mdl.parLink[[1]][[2]] <- "log"
Mdl.parLink[[2]][[1]] <- "identity"
Mdl.parLink[[2]][[2]] <- "log"
Mdl.parLink[[3]][[1]] <- "logit"
Mdl.parLink[[3]][[2]] <- "logit"

## THE TRUE PARAMETER VALUES IN THE DGP
##-----------------------------------------------------------------------------
## The parameters are the features of the model, e.g. mean,  variance, ...
## The parameters are observation specified (each observation has its own
## feature.)
##-----------------------------------------------------------------------------
MdlDGP.par <- MdlDataStruc

## The first margin
MdlDGP.par[[1]][[1]] <- matrix(rnorm(n = nObs, mean = 0, sd = 1))
MdlDGP.par[[1]][[2]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1))

## The second margin
MdlDGP.par[[2]][[1]] <- matrix(rnorm(n = nObs, mean = 0, sd = 1))
MdlDGP.par[[2]][[2]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1))

## The copula
MdlDGP.par[[3]][[1]] <- matrix(rbeta2(n = nObs, mean = 0.3, sd = 0.4))
MdlDGP.par[[3]][[2]] <- matrix(rbeta2(n = nObs, mean = 0.3, sd = 0.4))

## THE TRUE COVARIATE-DEPENDENT PARAMETER VALUES IN THE DGP
## -----------------------------------------------------------------------------
## Note that (a) The intercept is always included.  (b) Variable selection is
## represented if the beta coefficient is zero.  (c) If there are no covariates
## dependent, only intercept should be set.
## -----------------------------------------------------------------------------
MdlDGP.beta <- MdlDataStruc

## The first margin
MdlDGP.beta[[1]][[1]] <- matrix(c(1, -1, 1, -1, 0))
MdlDGP.beta[[1]][[2]] <- matrix(c(1, -1, 1, -1, 0))

## The second margin
MdlDGP.beta[[2]][[1]] <- matrix(c(1, -1, 1, -1, 0))
MdlDGP.beta[[2]][[2]] <- matrix(c(1, -1, 1, -1, 0))

## The copula
MdlDGP.beta[[3]][[1]] <- matrix(c(1, -1, 1, -1, 0, 1, -1, 1, -1, 0))
MdlDGP.beta[[3]][[2]] <- matrix(c(1, -1, 1, -1, 0, 1, -1, 1, -1, 0))

## Generating the numerical tabular for the inverse Kendall's tau
tauTabular <- kendalltauTabular(CplNM = CplNM, tol = 0.005)
