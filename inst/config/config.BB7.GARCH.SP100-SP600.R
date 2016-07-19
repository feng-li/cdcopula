###############################################################################
### TITLE
###        COPULA MODEL WITH COVARIATE--DEPENDENT
###
### SYSTEM REQUIREMENTS
###        R > = 2.15.3 WITH PACKAGES ``MVTNORM'',  ``PARALLEL''
###        BLAS (OPTIONAL)
###
### INPUT  VARIABLES
###        THE VARIABLES WITH COMMENTS IN CAPITAL LETTERS ARE USER INPUT
###        VARIABLES.
###
### OUTPUT VARIABLES
###        THE OUTPUT VARIABLES ARE ALWAYS STARTED WITH ``MCMC.''
###
### GENERALIZATION
###        IF YOU ARE INTERESTED IN DEVELOPING NEW COPULA MODELS BASED ON THE
###        EXISTING CODE. LOOK INTO THE ``MODEL'' FOLDER.
###
### DATE
###        CREATED: MON JAN 09 17:12:04 CET 2012
###        CURRENT: SUN JAN 04 10:21:56 CST 2015
###############################################################################


###----------------------------------------------------------------------------
### SPECIFY THE MODEL
###----------------------------------------------------------------------------
## MARGINAL MODELS NAME, TYPE AND PARAMETERS
MDL.MARGISTYPE <- C("GARCH", "GARCH", "BB7")
MDL.MARGISNM <- C("^SML", "^OEX", "BB7")

MCMC.UPDATE <- LIST(LIST("MU" = F, "PHI" = F),
                    LIST("MU" = F, "PHI" = F),
                    LIST("LAMBDAL" = T, "TAU" = T))

NAMES(MCMC.UPDATE) <- MDL.MARGISNM

## THE MODEL EVALUATION CRITERION
## SET THIS TO NULL TO TURN OF EVALUATION.
LPDS <- C("JOINT", MDL.MARGISNM)

## THE OBJECT STRUCTURE FOR THE MODEL COMPONENTS
NAMES(MDL.MARGISTYPE) <-  MDL.MARGISNM

###----------------------------------------------------------------------------
### THE DATA AND MODEL
###----------------------------------------------------------------------------

## THE DATASET
##-----------------------------------------------------------------------------
## THE DATASET SHOULD EITHER PROVIDED VIA DGP OR THE REAL DATASET. THE DATASET
## SHOULD BE IN THE FOLLOWING STRUCTURE:
## MDL.X: "LIST" EACH LIST CONTAINS THE COVARIATES IN EACH MARGIN OR COPULA.
## MDL.Y: "LIST" EACH LIST CONTAINS THE RESPONSE VARIABLE OF THAT MARGIN.

LOAD(FILE.PATH(R_CPL_LIB_ROOT_DIR, "DATA/SP100-SP400-SP600-20150206.RDATA"))

## NO. OF TOTAL OBSERVATIONS
NOBSRAW <- LENGTH(Y[[1]])

## DATA SUBSET USED
MDL.DATAUSEDIDX <- (1 + NOBSRAW-NOBSRAW):NOBSRAW


## THE RESPONSE VARIABLES
MDL.Y <- LAPPLY(Y[MDL.MARGISNM[-LENGTH(MDL.MARGISNM)]], FUNCTION(X, IDX)X[IDX, ,DROP = FALSE], MDL.DATAUSEDIDX)

## THE NAME OF RESPOND VARIABLES
NAMES(MDL.Y) <- MDL.MARGISNM[-LENGTH(MDL.MARGISNM)]

## COVARIATES USED FOR THE MARGINAL AND COPULA PARAMETERS
## ------------------------------------------------------------------------------

## A TRICK TO INCLUDE FOREIGN MARGINAL MODELS IN THE ESTIMATION WHICH ARE HARD TO DIRECTLY
## PUT INTO THE "MARGIMODEL()" IS DO THE FOLLOWING SETTINGS: (1) LET "MCMC.UPDATE" BE FALSE
## IN ALL MARGINAL DENSITIES.  (2) ESTIMATE THE DENSITY FEATURES IN FOREIGN MODELS AND SET
## THE FEATURES IN "MDL.X" DIRECTLY.  (3) SET MCMC.UPDATESTRATEGY BE "TWOSTAGE". (4) SET
## "MDL.BETAINIT" BE ONE IN ALL MARGINAL FEATURES.
MDL.X <- MCMC.UPDATE

MDL.X[[1]] <- LIST(INCLUDE.MEAN = TRUE, COND.DIST = "NORM", TRACE = TRUE)
MDL.X[[2]] <- LIST(INCLUDE.MEAN = TRUE, COND.DIST = "NORM", TRACE = TRUE)

MDL.X[[3]][["LAMBDAL"]] <- CBIND(1, X[[1]][MDL.DATAUSEDIDX, 1:9], X[[2]][MDL.DATAUSEDIDX, 1:9])
MDL.X[[3]][["TAU"]] <- CBIND(1, X[[1]][MDL.DATAUSEDIDX, 1:9], X[[2]][MDL.DATAUSEDIDX, 1:9])

## THE LINK FUNCTION USED IN THE MODEL
MDL.PARLINK <- MCMC.UPDATE
MDL.PARLINK[[1]][["MU"]] <- LIST(TYPE = "IDENTITY", NPAR = 1)
MDL.PARLINK[[1]][["PHI"]] <- LIST(TYPE = "IDENTITY", NPAR = 1)

MDL.PARLINK[[2]][["MU"]] <- LIST(TYPE = "IDENTITY", NPAR = 1)
MDL.PARLINK[[2]][["PHI"]] <- LIST(TYPE = "IDENTITY", NPAR = 1)

MDL.PARLINK[[3]][["LAMBDAL"]] <- LIST(TYPE = "GLOGIT", NPAR = 1, A = 0.01, B = 0.99)
MDL.PARLINK[[3]][["TAU"]] <- LIST(TYPE = "GLOGIT", NPAR = 1, A = 0.01, B = 0.99)

## THE VARIABLE SELECTION SETTINGS AND STARTING POINT
## VARIABLE SELECTION CANDIDATES, NULL: NO VARIABLE SELECTION USE FULL
## COVARIATES. ("ALL-IN", "ALL-OUT", "RANDOM", OR USER-INPUT)

MDL.VARSELARGS <- MCMC.UPDATE
MDL.VARSELARGS[[1]][["MU"]] <- LIST(CAND = NULL, INIT = "ALL-IN")
MDL.VARSELARGS[[1]][["PHI"]] <- LIST(CAND = NULL, INIT = "ALL-IN")

MDL.VARSELARGS[[2]][["MU"]] <- LIST(CAND = NULL, INIT = "ALL-IN")
MDL.VARSELARGS[[2]][["PHI"]] <- LIST(CAND = NULL, INIT = "ALL-IN")

MDL.VARSELARGS[[3]][["LAMBDAL"]] <- LIST(CAND = "2:END", INIT = "ALL-IN")
MDL.VARSELARGS[[3]][["TAU"]] <- LIST(CAND = "2:END", INIT = "ALL-IN")

###----------------------------------------------------------------------------
### THE MCMC CONFIGURATION
###----------------------------------------------------------------------------

## NUMBER OF MCMC ITERATIONS
MCMC.NITER <- 1000

## SAVE OUTPUT PATH
##-----------------------------------------------------------------------------
## "SAVE.OUTPUT = FALSE" IT WILL NOT SAVE ANYTHING.
## "SAVE.OUTPUT = "PATH-TO-DIRECTORY"" IT WILL SAVE THE WORKING DIRECTORY IN
## THE GIVEN DIRECTORY.
SAVE.OUTPUT <- "~/RUNNING"

## MCMC TRAJECTORY
##-----------------------------------------------------------------------------
## IF TRUE,  THE MCMC SHOULD BE TRACKED DURING THE EVALUATION.
MCMC.TRACK <- TRUE

MCMC.UPDATEORDER <- MCMC.UPDATE
MCMC.UPDATEORDER[[1]][[1]] <- 1

MCMC.UPDATEORDER[[2]][[1]] <- 2

MCMC.UPDATEORDER[[3]][[1]] <- 3
MCMC.UPDATEORDER[[3]][[2]] <- 4

## MCMC UPDATING STRATEGY
##-----------------------------------------------------------------------------
## "JOINT"    : UPDATE THE JOINT POSTERIOR W.R.T. MCMC.UPDATE AND MCMC.UPDATEORDER
## "MARGIN"   : THE MARGINAL POSTERIOR.
## "TWOSTAGE" : UPDATE THE JOINT POSTERIOR BUT USING A TWO STAGE APPROACH.

## NOTE: IF ONE WANT TO USE "MARGIN" OR "TWOSTAGE" OPTIONS JUST TO TO ESTIMATE THE COPULA
## DENSITY. A VARIABLE "MCMC.DENSITY[["U"]]" MUST PROVIDE. "MCMC.DENSITY" CONSISTS OF CDF OF
## MARGINS (I.E. U1,  U2, ...)

MCMC.UPDATESTRATEGY <- "TWOSTAGE"

## THE METROPOLIS-HASTINGS ALGORITHM PROPOSAL ARGUMENTS
MCMC.PROPARGS <- MCMC.UPDATE
MCMC.PROPARGS[[1]][[1]] <- NA
MCMC.PROPARGS[[1]][[2]] <- NA

MCMC.PROPARGS[[2]][[1]] <- NA
MCMC.PROPARGS[[2]][[2]] <- NA

MCMC.PROPARGS[[3]][[1]] <- LIST("ALGORITHM" = LIST(TYPE = "GNEWTONMOVE", KSTEPS = 3, HESS = "OUTER"),
                           "BETA" = LIST(TYPE = "MVT", DF = 6),
                           "INDICATORS" = LIST(TYPE = "BINOM", PROB = 0.5))
MCMC.PROPARGS[[3]][[2]] <- LIST("ALGORITHM" = LIST(TYPE = "GNEWTONMOVE", KSTEPS = 3, HESS = "OUTER"),
                           "BETA" = LIST(TYPE = "MVT", DF = 6),
                           "INDICATORS" = LIST(TYPE = "BINOM", PROB = 0.5))

## POSTERIOR INFERENCE OPTIONS
##-----------------------------------------------------------------------------

## CROSS VALIDATION
## "N.SUBSETS" IS NO. OF FOLDS FOR CROSS-VALIDATION. IF N.SUBSETS  =  0, NO
## CROSS-VALIDATION. AND "PARTIMETHOD" TELLS HOW TO PARTITION THE DATA. TESTING
## PERCENT IS USED IF PARTIMETHOD IS "TIME-SERIES". (USE THE OLD DATA TO
## PREDICT THE NEW INTERVAL)

NCROSS <- 1
MDL.CROSSVALIDARGS <- LIST(N.SUBSETS = NCROSS,
                       PARTIMETHOD = "TIME-SERIES",
                       TESTRATIO = 0.2)

## INDICES FOR TRAINING AND TESTING SAMPLE ACCORDING TO CROSS-VALIDATION
MDL.CROSSVALIDIDX <- SET.CROSSVALID(LENGTH(MDL.DATAUSEDIDX),MDL.CROSSVALIDARGS)
## NCROSSFOLD <- LENGTH(MDL.CROSSVALIDIDX[["TRAINING"]])

## SAMPLER PROPORTION FOR POSTERIOR INFERENCE,
MCMC.SAMPLEPROP <- 1

## BURN-IN RATIO
MCMC.BURNINPROP <- 0.1 # ZERO INDICATES NO BURN-IN

###----------------------------------------------------------------------------
### PRIOR SETTINGS
###----------------------------------------------------------------------------

## PRIOR FOR THE COPULA PARAMETERS
## -----------------------------------------------------------------------------
## NOTE: THE VARIABLE ARE RECYCLED IF NEEDED. FOR EXAMPLE INDICATORS$PROB CAN BE A SCALER
## OR A VECTOR WITH SAME LENGTH OF VARIABLE SECTION CANDIDATES. THERE MIGHT BE CONNECTIONS
## BETWEEN PARAMETERS IN THE MODELS BUT IS WILL NOT AFFECT THE PRIOR SETTINGS ON THE
## COEFFICIENTS AS LONG AS WE USE A DYNAMIC LINK FUNCTION.

MDL.PRIARGS <- MCMC.UPDATE

MDL.PRIARGS[[1]][["MU"]] <- NA
MDL.PRIARGS[[1]][["PHI"]] <- NA

MDL.PRIARGS[[2]][["MU"]] <- NA
MDL.PRIARGS[[2]][["PHI"]] <- NA


MDL.PRIARGS[[3]][["LAMBDAL"]] <- LIST("BETA" = LIST("INTERCEPT" = LIST(TYPE = "CUSTOM",
                                                                   INPUT = LIST(TYPE = "GBETA",  MEAN = 0.2, VARIANCE = 0.05,
                                                                                A = 0.01, B = 0.99),
                                                                   OUTPUT = LIST(TYPE = "NORM", SHRINKAGE = 1)),
                                                "SLOPES" = LIST(TYPE = "COND-MVNORM",
                                                                MEAN = 0, COVARIANCE = "IDENTITY", SHRINKAGE = 1)),
                                  "INDICATORS" = LIST(TYPE = "BERN", PROB = 0.5))
MDL.PRIARGS[[3]][["TAU"]] <- LIST("BETA" = LIST("INTERCEPT" = LIST(TYPE = "CUSTOM",
                                                                   INPUT = LIST(TYPE = "GBETA",  MEAN = 0.2, VARIANCE = 0.05,
                                                                                A = 0.01, B = 0.99),
                                                                   OUTPUT = LIST(TYPE = "NORM", SHRINKAGE = 1)),
                                                "SLOPES" = LIST(TYPE = "COND-MVNORM",
                                                                MEAN = 0, COVARIANCE = "IDENTITY", SHRINKAGE = 1)),
                                  "INDICATORS" = LIST(TYPE = "BERN", PROB = 0.5))

###----------------------------------------------------------------------------
### THE PARAMETERS FOR INITIAL AND CURRENT MCMC ITERATION
### THE PARAMETERS IN THE CURRENT MCMC ITERATION. FOR THE FIRST ITERATION, IT
### IS SET AS THE INITIAL VALUES
###----------------------------------------------------------------------------

## THE PARAMETER COEFFICIENTS STARTING POINT
## THE POSSIBLE INPUTS ARE ("RANDOM", "OLS"  OR USER-INPUT).
MDL.BETAINIT <- MCMC.UPDATE
MDL.BETAINIT[[1]][[1]] <- 1
MDL.BETAINIT[[1]][[2]] <- 1

MDL.BETAINIT[[2]][[1]] <- 1
MDL.BETAINIT[[2]][[2]] <- 1

MDL.BETAINIT[[3]][[1]] <- "RANDOM"
MDL.BETAINIT[[3]][[2]] <- "RANDOM"

MCMC.OPTIMINIT <- TRUE
################################################################################
###                                  THE END
################################################################################
