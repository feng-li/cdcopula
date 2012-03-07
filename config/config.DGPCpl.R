###############################################################################
### Configuration file for the copula data generating process 
###############################################################################

## NO. OF OBSERVATIONS
nObs <- 15

## THE VARIABLE SELECTION SETTINGS AND STARTING POINT
## Variable selection candidates, NULL: no variable selection use full
## covariates. ("all-in", "all-out", "random", or user-input)
varSelArgs <- MdlDataStruc
varSelArgs[[1]][[1]] <- list(cand = c(2, 3),
                             init = "all-in") 
varSelArgs[[1]][[2]] <- list(cand = c(2, 3),
                             init = "all-out")
varSelArgs[[2]][[1]] <- list(cand = c(2, 4),
                             init = "random")
varSelArgs[[2]][[2]] <- list(cand = c(2, 4),
                             init = "all-out")
varSelArgs[[3]][[1]] <- list(cand = c(2, 3, 5, 6),
                             init = c(2, 3))
varSelArgs[[3]][[2]] <- list(cand = c(3, 5, 6),
                             init = "random")

## THE LINK FUNCTION USED IN THE MODEL
Mdl.parLink <- MdlDataStruc
Mdl.parLink[[1]][[1]] <- "identity"
Mdl.parLink[[1]][[2]] <- "log"
Mdl.parLink[[2]][[1]] <- "identity"
Mdl.parLink[[2]][[2]] <- "log"
Mdl.parLink[[3]][[1]] <- "logit"
Mdl.parLink[[3]][[2]] <- "logit"

## THE TRUE PARAMETER VALUES IN THE DGP
DGP.betaTRUE <- MdlDataStruc
DGP.betaTRUE[[1]][[1]] <- c(0.2, 0.6, 0.7)
DGP.betaTRUE[[1]][[2]] <- c(0.2, -0.6, 0.7)
DGP.betaTRUE[[2]][[1]] <- c(0.2, 0.6, 0.7, 0.9)
DGP.betaTRUE[[2]][[2]] <- c(0.2, -0.6, 0.7, -0.9)
DGP.betaTRUE[[3]][[1]] <- c(0.2, 0.6, 0.7, 0.2, 0.6, 0.7)
DGP.betaTRUE[[3]][[2]] <- c(0.2, 0.6, 0.7, 0.2, 0.6, 0.7)
