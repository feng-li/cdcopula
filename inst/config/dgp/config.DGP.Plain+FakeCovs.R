###############################################################################
### Configuration file for the copula data generating process
###############################################################################

###----------------------------------------------------------------------------
### SPECIFY THE MODEL
###----------------------------------------------------------------------------
## MARGINAL MODELS NAME, TYPE AND PARAMETERS
Mdl.MargisType <- c("SPLITT", "SPLITT", "BB7")
Mdl.MargisNM <- c("M1", "M2", "BB7")

MdlDGP.Update <- list(list("mu" = T, "phi" = T, "df" = T, "lmd" = T),
                      list("mu" = T, "phi" = T, "df" = T, "lmd" = T),
                      list("lambdaL" = T, "lambdaU" = T))

names(MdlDGP.Update) <- Mdl.MargisNM

## The object structure for the model components
names(Mdl.MargisType) <-  Mdl.MargisNM


## NO. OF OBSERVATIONS
nObs <- 1000

## THE LINK FUNCTION USED IN THE MODEL
Mdl.parLink <- MdlDGP.Update
Mdl.parLink[[1]][["mu"]] <- list(type = "identity", nPar = 1)
Mdl.parLink[[1]][["phi"]] <- list(type = "glog", a = 0.01, b = 100, nPar = 1)
Mdl.parLink[[1]][["df"]] <- list(type = "glog", nPar = 1,  a = 2, b = 30)
Mdl.parLink[[1]][["lmd"]] <- list(type = "glog", a = 0.01, b = 100, nPar = 1)

Mdl.parLink[[2]][["mu"]] <- list(type = "identity", nPar = 1)
Mdl.parLink[[2]][["phi"]] <- list(type = "glog", a = 0.01, b = 100, nPar = 1)
Mdl.parLink[[2]][["df"]] <- list(type = "glog", nPar = 1,  a = 2, b = 30)
Mdl.parLink[[2]][["lmd"]] <- list(type = "glog", a = 0.01, b = 100, nPar = 1)

Mdl.parLink[[3]][["lambdaL"]] <- list(type = "glogit", nPar = 1, a = 0.01, b = 0.99)
Mdl.parLink[[3]][["lambdaU"]] <- list(type = "glogit", nPar = 1, a = 0.01, b = 0.99)

## -----------------------------------------------------------------------------
## TRUE PARAMETER VALUES IN DGP
## -----------------------------------------------------------------------------
## The parameters are the features of the model, e.g. mean, variance, ...  The
## parameters are observation specified (each observation has its own
## features.) This is the desired feature, the TRUE DGP might be not exactly
## the same as it. See "Mdl.par" for the comparison.
## -----------------------------------------------------------------------------
MdlDGP.par <- MdlDGP.Update

## ## The first margin
## MdlDGP.par[[1]][["mu"]] <- matrix(0, nObs, 1)
## MdlDGP.par[[1]][["phi"]] <- matrix(0.5, nObs, 1)
## MdlDGP.par[[1]][["df"]] <- matrix(6, nObs, 1)
## MdlDGP.par[[1]][["lmd"]] <- matrix(1, nObs, 1)

## ## The second margin
## MdlDGP.par[[2]][["mu"]] <- matrix(0, nObs, 1)
## MdlDGP.par[[2]][["phi"]] <- matrix(0.5, nObs, 1)
## MdlDGP.par[[2]][["df"]] <- matrix(6, nObs, 1)
## MdlDGP.par[[2]][["lmd"]] <- matrix(1,  nObs, 1)

## ## The copula component
## MdlDGP.par[[3]][["lambdaL"]] <- matrix(0.1, nObs, 1)
## MdlDGP.par[[3]][["lambdaU"]] <- matrix(0.1, nObs, 1)

## The first margin
MdlDGP.par[[1]][["mu"]] <- matrix(rnorm(n = nObs, mean = 0, sd = 1), nObs, 1)
MdlDGP.par[[1]][["phi"]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1), nObs, 1)
MdlDGP.par[[1]][["df"]] <- matrix(rlnorm2(n = nObs, mean = 6, sd = 1), nObs, 1)
MdlDGP.par[[1]][["lmd"]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1), nObs, 1)

## The second margin
MdlDGP.par[[2]][["mu"]] <- matrix(rnorm(n = nObs, mean = 0, sd = 1), nObs, 1)
MdlDGP.par[[2]][["phi"]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1), nObs, 1)
MdlDGP.par[[2]][["df"]] <- matrix(rlnorm2(n = nObs, mean = 6, sd = 1), nObs, 1)
MdlDGP.par[[2]][["lmd"]] <- matrix(rlnorm2(n = nObs, mean = 1, sd = 1), nObs, 1)

## The copula component
MdlDGP.par[[3]][["lambdaL"]] <- matrix(rbeta2(n = nObs, mean = 0.3, sd = 0.1), nObs, 1)
MdlDGP.par[[3]][["lambdaU"]] <- matrix(rbeta2(n = nObs, mean = 0.5, sd = 0.1), nObs, 1)

## Special treatments for parameters under/over Mdl.parLink boundary
MdlDGP.par <- DGPCplRestrictPar(MdlDGP.par = MdlDGP.par, Mdl.parLink = Mdl.parLink)


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
MdlDGP.intercept <- MdlDGP.Update

MdlDGP.intercept[[1]][[1]] <- TRUE
MdlDGP.intercept[[1]][[2]] <- TRUE
MdlDGP.intercept[[1]][[3]] <- TRUE
MdlDGP.intercept[[1]][[4]] <- TRUE

MdlDGP.intercept[[2]][[1]] <- TRUE
MdlDGP.intercept[[2]][[2]] <- TRUE
MdlDGP.intercept[[2]][[3]] <- TRUE
MdlDGP.intercept[[2]][[4]] <- TRUE

MdlDGP.intercept[[3]][[1]] <- TRUE
MdlDGP.intercept[[3]][[2]] <- TRUE

## NUMBER OF COVARIATES (INCLUDING INTERCEPT)
## If MdlDGP.intercept is TRUE, the first element in MdlDGP.beta is the intercept.

beta0 <- MdlDGP.Update
for(i in names(MdlDGP.Update))
{
    for(j in names(MdlDGP.Update[[i]]))
    {
        ## beta0 is the only TRUE intercepts.
        beta0[[i]][[j]] = parLinkFun(mean(MdlDGP.par[[i]][[j]]),
                                     linkArgs = Mdl.parLink[[i]][[j]])
    }
}

MdlDGP.beta <- MdlDGP.Update

## ## The first margin
## MdlDGP.beta[[1]][[1]] <- matrix(c(beta0[[1]][[1]], 0, 0, 0, 0))
## MdlDGP.beta[[1]][[2]] <- matrix(c(beta0[[1]][[2]], 0, 0, 0, 0))
## MdlDGP.beta[[1]][[3]] <- matrix(c(beta0[[1]][[3]], 0, 0, 0, 0))
## MdlDGP.beta[[1]][[4]] <- matrix(c(beta0[[1]][[4]], 0, 0, 0, 0))

## ## The second margin
## MdlDGP.beta[[2]][[1]] <- matrix(c(beta0[[2]][[1]], 0, 0, 0, 0))
## MdlDGP.beta[[2]][[2]] <- matrix(c(beta0[[2]][[2]], 0, 0, 0, 0))
## MdlDGP.beta[[2]][[3]] <- matrix(c(beta0[[2]][[3]], 0, 0, 0, 0))
## MdlDGP.beta[[2]][[4]] <- matrix(c(beta0[[2]][[4]], 0, 0, 0, 0))

## ## The copula
## MdlDGP.beta[[3]][[1]] <- matrix(c(beta0[[3]][[1]], 1, -1, 1, -1, 0, 0, 0, 0))
## MdlDGP.beta[[3]][[2]] <- matrix(c(beta0[[3]][[2]], 1, -1, 1, -1, 0, 0, 0, 0))


MdlDGP.beta[[1]][[1]] <- matrix(c(beta0[[1]][[1]]))
MdlDGP.beta[[1]][[2]] <- matrix(c(beta0[[1]][[2]]))
MdlDGP.beta[[1]][[3]] <- matrix(c(beta0[[1]][[3]]))
MdlDGP.beta[[1]][[4]] <- matrix(c(beta0[[1]][[4]]))

##  The second margin
MdlDGP.beta[[2]][[1]] <- matrix(c(beta0[[2]][[1]]))
MdlDGP.beta[[2]][[2]] <- matrix(c(beta0[[2]][[2]]))
MdlDGP.beta[[2]][[3]] <- matrix(c(beta0[[2]][[3]]))
MdlDGP.beta[[2]][[4]] <- matrix(c(beta0[[2]][[4]]))

## The copula
## MdlDGP.beta[[3]][[1]] <- matrix(c(beta0[[3]][[1]]))
## MdlDGP.beta[[3]][[2]] <- matrix(c(beta0[[3]][[2]]))

## MdlDGP.beta[[3]][[1]] <- matrix(c(beta0[[3]][[1]], 0, 0, 0, 0, 0, 0, 0, 0))
## MdlDGP.beta[[3]][[2]] <- matrix(c(beta0[[3]][[2]], 0, 0, 0, 0, 0, 0, 0, 0))

MdlDGP.beta[[3]][[1]] <- matrix(c(beta0[[3]][[1]], 1, -1, 1, -1, 0, 0, 0, 0))
MdlDGP.beta[[3]][[2]] <- matrix(c(beta0[[3]][[2]], 1, -1, 1, -1, 0, 0, 0, 0))

################################################################################
###                                  THE END
################################################################################
