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

## THE TRUE PARAMETER VALUES IN THE DGP (INCLUDING INTERCEPTS)
DGP.betaTRUE <- MdlDataStruc
DGP.betaTRUE[[1]][[1]] <- c(0.2,  0.6, 0.7)
DGP.betaTRUE[[1]][[2]] <- c(0.2, -0.6, 0.7)
DGP.betaTRUE[[2]][[1]] <- c(0.2,  0.6, 0.7,  0.9)
DGP.betaTRUE[[2]][[2]] <- c(0.2, -0.6, 0.7, -0.9)
DGP.betaTRUE[[3]][[1]] <- c(0.2,  0.6, 0.7,  0.2, 0.6, 0.7)
DGP.betaTRUE[[3]][[2]] <- c(0.2,  0.6, 0.7,  0.2, 0.6, 0.7)

## Generating the numerical tabular for the inverse Kendall's tau
tauTabular <- kendalltauTabular(CplNM = CplNM, tol = 0.005)
