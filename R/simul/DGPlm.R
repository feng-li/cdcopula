##' Simple data generating process for the deterministic linear model Y = XB given Y and
##' coefficients.
##'
##' Remember there is no error term in the model. This function is primarily used in the
##' covariate-dependent type of models to generate desired features.
##' @title Deterministic linear model DGP.
##' @param Y "scaler"
##' @param beta "scaler"
##' @param Xlim "vector"
##' @param intercept "logical" The intercept be considered in the DGP if the value for the
##'     argument is TRUE.
##' @return "matrix"
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note DEPENDS: flutils Created: Wed Mar 14 10:46:20 CET 2012; Current: Wed Mar 14
##'     10:46:25 CET 2012.
DGPlm <- function(Y, beta, Xlim, intercept = TRUE)
{
    n <- length(Y)
    q <- length(beta)
    B <- matrix(beta, q, 1) # Currently only allow univariate response

    ## Initialize the storage for the covariates. If intercept is in,  the first column is
    ## for intercept automatically.
    X <- matrix(1, n, q)

    ## Check if intercept is included
    if(intercept)
    {
        ## Indices for zero beta coefficients (except intercept) i.e. covariates not
        ## selected via variable selection.

        betaIdx.Zero <- which(abs(beta[-1]) ==0)+1
        betaIdx.NonZero <- which(abs(beta[-1]) !=0)+1
    }
    else
    {
        betaIdx.Zero <- which(abs(beta) ==0)
        betaIdx.NonZero <- which(abs(beta) !=0)
    }

    ## Indices Length for zero/non-zero beta coefficients (except intercept)
    betaIdxLen.Zero <- length(betaIdx.Zero)
    betaIdxLen.NonZero <- length(betaIdx.NonZero)

    if(length(betaIdxLen.Zero) > 0)
    {
        betaIdx.NoInt <- c(betaIdx.NonZero, betaIdx.Zero)
        ## variable selection used, fill those variables with garbage information
        X[, betaIdx.NoInt] <- runif(n*length(betaIdx.NoInt), min = Xlim[1], max = Xlim[2])
    }


    ## Variational Y with constant X is not possible when intercept in.
    if((betaIdxLen.NonZero == 0 & intercept & any(Y != beta[1])) ||
       (all(beta == 0) & any(Y != 0)))
    {
        stop("Variational Y with constant X is not possible!")
    }

    ## Algorithm that finds x_1, ...x_p for y = x_1*b1+...+x_p*b_p give y and b1, ...bp
    ## Randomly select some positions to be determined from all observations.
    if(length(betaIdx.NonZero) == 1)
    {
        IdxLastN <- rep(betaIdx.NonZero, n)
    }
    else
    {
        if(length(betaIdx.NonZero) == 0) browser()
        IdxLastN <- sample(x = betaIdx.NonZero, size = n, replace = TRUE)
    }

    IdxLastNLoc <- whichInd(arr.ind = cbind(1:n, IdxLastN), dims = c(n, NA))

    ## Make those positions be zero temporally.
    X[IdxLastNLoc] <- 0

    ## Calculate the appropriate values for those positions.
    Y4NoLastN <- X%*%B
    XLastN <- (Y-Y4NoLastN)/beta[IdxLastN]

    ## Fix the hole
    X[IdxLastNLoc] <- XLastN

    return(X)
}

###----------------------------------------------------------------------------
### TESTING
###----------------------------------------------------------------------------

## Y <- matrix(rnorm(20, 1))
## beta <- runif(6)
## Xlim <- c(-3, 3)

## DGPlm(Y, beta, Xlim, intercept = TRUE)
