###############################################################################
### TITLE
###        Foreign Multivariate Models
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
## Mdl.MargisType <- c("", "", "DCCGARCH")
## Mdl.MargisType <- c("", "", "GOGARCH")
Mdl.MargisType <- c("", "", "VAR")

## Mdl.MargisNM <- c("000021", "000066", "VAR")
## Mdl.MargisNM <- c("000021", "600755", "VAR")
## Mdl.MargisNM <- c("000066", "600755", "VAR")


## Mdl.MargisNM <- c("000021", "000032", "VAR")
## Mdl.MargisNM <- c("000021", "000727", "VAR")
## Mdl.MargisNM <- c("000021", "000748", "VAR")
## Mdl.MargisNM <- c("000021", "600171", "VAR")
## Mdl.MargisNM <- c("000021", "600536", "VAR")
## Mdl.MargisNM <- c("000021", "600764", "VAR")
## Mdl.MargisNM <- c("000021", "600775", "VAR")

## Mdl.MargisNM <- c("000032", "000066", "VAR")
## Mdl.MargisNM <- c("000032", "000727", "VAR")
## Mdl.MargisNM <- c("000032", "000748", "VAR")
## Mdl.MargisNM <- c("000032", "600171", "VAR")
## Mdl.MargisNM <- c("000032", "600536", "VAR")
## Mdl.MargisNM <- c("000032", "600755", "VAR")
## Mdl.MargisNM <- c("000032", "600764", "VAR")
## Mdl.MargisNM <- c("000032", "600775", "VAR")

## Mdl.MargisNM <- c("000066", "000727", "VAR")
## Mdl.MargisNM <- c("000066", "000748", "VAR")
## Mdl.MargisNM <- c("000066", "600171", "VAR")
## Mdl.MargisNM <- c("000066", "600536", "VAR")
## Mdl.MargisNM <- c("000066", "600764", "VAR")
## Mdl.MargisNM <- c("000066", "600775", "VAR")

## Mdl.MargisNM <- c("000727", "000748", "VAR")
## Mdl.MargisNM <- c("000727", "600755", "VAR")
## Mdl.MargisNM <- c("000727", "600764", "VAR")

## Mdl.MargisNM <- c("000748", "600171", "VAR")
## Mdl.MargisNM <- c("000748", "600536", "VAR")
## Mdl.MargisNM <- c("000748", "600755", "VAR")
## Mdl.MargisNM <- c("000748", "600764", "VAR")
## Mdl.MargisNM <- c("000748", "600775", "VAR")

## Mdl.MargisNM <- c("600171", "600755", "VAR")
## Mdl.MargisNM <- c("600171", "600764", "VAR")

## Mdl.MargisNM <- c("600536", "600755", "VAR")
## Mdl.MargisNM <- c("600536", "600764", "VAR")

## Mdl.MargisNM <- c("600755", "600764", "VAR")
## Mdl.MargisNM <- c("600755", "600775", "VAR")

Mdl.MargisNM <- c("600764", "600775", "VAR")


## THE MODEL EVALUATION CRITERION
## Set this to NULL to turn of evaluation.
LPDS <- c("joint")

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

source(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/scripts/Get-RISK-CONTAGION.R"))

## No. of Total Observations
nObsRaw <- length(Y[[1]])

## Data subset used
Mdl.dataUsedIdx <- (1 + nObsRaw-nObsRaw):nObsRaw
Mdl.algorithm <- "full"

## THE RESPONSE VARIABLES
Mdl.Y <- lapply(Y[Mdl.MargisNM[-length(Mdl.MargisNM)]], function(x, idx)x[idx, ,drop = FALSE], Mdl.dataUsedIdx)

## The name of respond variables
names(Mdl.Y) <- Mdl.MargisNM[-length(Mdl.MargisNM)]

### THE MODEL
###----------------------------------------------------------------------------
{if(tolower(Mdl.MargisType[length(Mdl.MargisType)]) == "gogarch")
 {
     ## GoGARCH
     require("rmgarch", quietly = TRUE)

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
 else if(tolower(Mdl.MargisType[length(Mdl.MargisType)]) == "dccgarch")
    {
        ## DCC-GARCH
        require("rugarch", quietly = TRUE)
        require("rmgarch", quietly = TRUE)

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
 else if(tolower(Mdl.MargisType[length(Mdl.MargisType)]) == "msbvar")
    {
        require("MSBVAR")

        h = 1 # Number of regimes/states, an integer

        ForeignModelSpec <- list(
            Y = NA,  # A T x m multiple time series object created with ts().

            p = 1,  # Lag length, an integer

            h = h,  # Number of regimes/states, an integer

            lambda0 = 0.8, # Value in [0, 1], Overall tightness of the prior (discounting
                           # of prior scale)

            lambda1 = 0.15, # Value in [0, 1], Standard deviation or tightness of the
                            # prior around the AR(1) parameters.

            lambda3 = 1, # Lag decay (>0, with 1 = harmonic)

            lambda4 = 0.25, # Standard deviation or tightness around the intercept >0

            lambda5 = 1, # Standard deviation or tightness around the exogneous variable
                         # coefficients >0

            mu5 = 0, # Sum of coefficients prior weight ≥0. Larger values imply difference
                     # stationarity.

            mu6 = 0, # Dummy initial observations or drift prior ≥0. Larger values allow
                     # for common trends.

            qm = 12, # Frequency of the data for lag decay equivalence. Default is 4, and
                     # a value of 12 will match the lag decay of monthly to quarterly
                     # data. Other values have the same effect as "4"

            alpha.prior = (100*diag(h) + matrix(2,  h,  h)), # Prior for the Dirichlet
                                                             # process for the MS
                                                             # process. Default is 100 *
                                                             # diag(h) + matrix(2, h, h),
                                                             # but the model will be
                                                             # sensitive to this.))
            prior = 0,# One of three values: 0 = Normal-Wishart prior, 1 = Normal-flat
                      # prior, 2 = flat-flat prior (i.e., akin to MLE). The conjugate
                      # prior is the first one, which is the default.

            max.iter = 30, # One of three values: 0 = Normal-Wishart prior, 1 =
                           # Normal-flat prior, 2 = flat-flat prior (i.e., akin to
                           # MLE). The conjugate prior is the first one, which is the
                           # default.

            initialize.opt = NULL # Initial values for the block optimization
                                  # algorithm. If default = NULL initialize.msbvar is
                                  # called to provide values. User can specify values as
                                  # long as they conform to the structure produced by
                                  # initialize.msbvar.
                     )

    }
 else if(tolower(Mdl.MargisType[length(Mdl.MargisType)]) == "var")
    {
        require("MTS")

        ForeignModelSpec <- list(
            p  =  4, # Order of VAR model.
            output  =  TRUE,
            include.mean  =  TRUE,
            fixed  =  NULL # A logical matrix used in constrained estimation. It is used
                           # mainly in model simplifcation, e.g., removing insignificant
                           # estimates.
        )
    }

 else
    {
        stop("No such foreign marginal models!")
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
Mdl.crossValidArgs <- list(N.subsets = nCross,
                           partiMethod = "time-series",
                           testRatio = 0.2)

## Indices for training and testing sample according to cross-validation
Mdl.crossValidIdx <- set.crossvalid(length(Mdl.dataUsedIdx), Mdl.crossValidArgs)

MCMC.optimInit <- TRUE
################################################################################
###                                  THE END
################################################################################
