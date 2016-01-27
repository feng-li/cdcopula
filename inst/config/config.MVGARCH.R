###############################################################################
### TITLE
###        GARCH type models
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
###        CREATED: Mon Sep 21 15:22:40 CST 2015
###        CURRENT: Mon Sep 21 15:22:44 CST 2015
###############################################################################


###----------------------------------------------------------------------------
### SPECIFY THE MODEL
###----------------------------------------------------------------------------
## MARGINAL MODELS NAME, TYPE AND PARAMETERS
MargisType <- c("NA", "NA", "DCCGARCH")
MargisNM <- c("^SML", "^OEX", "DCCGARCH")

## THE MODEL EVALUATION CRITERION
## Set this to NULL to turn of evaluation.
LPDS <- c("joint")

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

load(file.path(R_CPL_LIB_ROOT_DIR, "data/SP100-SP400-SP600-20150206.Rdata"))

## No. of Total Observations
nObsRaw <- length(Y[[1]])

## Data subset used
nObsIdx <- (1 + nObsRaw-nObsRaw):nObsRaw

## No. of used Observations
nObs <- length(nObsIdx)

## THE RESPONSE VARIABLES
Mdl.Y <- lapply(Y[MargisNM[-length(MargisNM)]], function(x, idx)x[idx, ,drop = FALSE], nObsIdx)

## The name of respond variables
names(Mdl.Y) <- MargisNM[-length(MargisNM)]

### THE MODEL
###----------------------------------------------------------------------------
{if(tolower(MargisType[length(MargisType)]) == "gogarch")
 {
   ## GoGARCH
   require("rmgarch")

   mean.model <- list(model  =  "constant",  robust  = FALSE,
                      lag  =  1,  lag.max  =  NULL,
                      lag.criterion  = "AIC",
                      external.regressors  =  NULL,
                      robust.control  =  list("gamma"  =  0.25,
                                              "delta"  =  0.01,
                                              "nc" =  10,
                                              "ns"  =  500))
   variance.model <- list(model  =  "sGARCH",
                           garchOrder  =  c(1, 1),
                           submodel  =  NULL,
                           variance.targeting  =  FALSE)
   ForeignModelSpec <- gogarchspec(mean.model = mean.model,
                                   variance.model = variance.model,
                                   distribution.model  =  "mvnorm")
 }
 else if(tolower(MargisType[length(MargisType)]) == "dccgarch")
 {
   ## DCC-GARCH
   require("rmgarch")
   uspec <- ugarchspec()
   mspec <- multispec(replicate(2, uspec))

   ForeignModelSpec <- dccspec(uspec = mspec, VAR = FALSE,
                               robust = FALSE, lag = 1, lag.max = NULL,
                               lag.criterion = "AIC",
                               external.regressors = NULL,
                               robust.control = list("gamma" = 0.25, "delta" = 0.01,
                                                     "nc" = 10, "ns" = 500),
                               dccOrder = c(1,1), model = "DCC",
                               groups = rep(1, length(uspec@spec)),
                               distribution = c("mvnorm"),
                               start.pars = list(), fixed.pars = list())
 }
 else
 {
   stop("No such foreign models!")
 }
}


## SAVE OUTPUT PATH
##-----------------------------------------------------------------------------
## "save.output = FALSE" it will not save anything.
## "save.output = "path-to-directory"" it will save the working directory in
## the given directory.
save.output <- "~/running"


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

optimInit <- TRUE
################################################################################
###                                  THE END
################################################################################
