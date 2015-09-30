##' Setup the initial values for the copula model
##'
##' This function is only once for the first iteration. The values are updated
##' via MCMC scheme.
##' @param varSelArgs "list".
##'
##'        Variable selection argument
##'
##' @param betaInit "character or numeric".
##'
##'        If is "character", the corresponding method are used to generate the
##' initial value.
##'
##'        If is "numeric", The initial value are taken as is.
##'
##' @param Mdl.X "list".
##'
##'        The covariates list
##'
##' @return "list"
##'
##'       Both the initial value for the beta coefficients ("Mdl.beta") and the
##' initial value for the variable selection indicator ("Mdl.betaIdx") are
##' returned.
##'
##' @references Li 2012
##' @author Feng Li, Central University of Finance and Economics.
##' @note Created: Thu Dec 22 15:57:14 CET 2011; Current: Tue Sep 29 23:50:07 CST 2015.
initPar <- function(varSelArgs, betaInit, Mdl.X, Mdl.Y, Mdl.parLink)
{
  ## The output structure.
  Mdl.betaIdx <- betaInit
  Mdl.beta <- betaInit
  ## Mdl.X <- Mdl.X
  ## Mdl.Y <- Mdl.Y


  ## Loop to assign the initial values
  CompNM <- names(betaInit)
  if(length(CompNM) == 0)
  {
    stop("Model component name must be specified first!")
  }


  for(i in CompNM)
  {
    CompParNM <- names(Mdl.parLink[[i]])
    for(j in CompParNM)
    {
      ## No. of col for covariates, including intercept if applicable.
      ncolX.ij <- ncol(Mdl.X[[i]][[j]])
      nPar.ij <- Mdl.parLink[[i]][[j]][["nPar"]]
###----------------------------------------------------------------------------
### Initialize the variable section indicator
###----------------------------------------------------------------------------
      ## Initial value for variable selection indicators which can be
      ## character or vector
      varSelCandConfigCurr <- varSelArgs[[i]][[j]][["cand"]]

      if(class(varSelCandConfigCurr) == "character" &&
         tolower(varSelCandConfigCurr) == "2:end")
      {
        varSelCandCurr <- 2:ncolX.ij
      }
      else
      {
        varSelCandCurr <- varSelCandConfigCurr
      }


      varSelInitCurr <- varSelArgs[[i]][[j]][["init"]]

      ## Check if variable selection candidates are all subsets of X.
      if(!all(varSelCandCurr %in% 1:ncolX.ij))
      {
        stop("Variable selection candidates are not subset of covariates in component: ",  i, "-", j, "!")
      }

      if(class(varSelInitCurr) == "character" &&
         tolower(varSelInitCurr) == "all-in")
      {
        Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)
      }
      else if(class(varSelInitCurr) == "character" &&
              tolower(varSelInitCurr) == "all-out")
      {
        Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)

        Mdl.betaIdx[[i]][[j]][varSelCandCurr, ] <- FALSE
      }
      else if(class(varSelInitCurr) == "character" &&
              tolower(varSelInitCurr) == "random")
      {
        Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)

        ## Randomly pick up half in
        betaIdxCurrSubOut <- sample(varSelCandCurr,
                                    round(length(varSelCandCurr)/2))
        Mdl.betaIdx[[i]][[j]][betaIdxCurrSubOut, ] <- FALSE
      }
      else # Do nothing, use user input
      {
        Mdl.betaIdx[[i]][[j]] <- matrix(TRUE, ncolX.ij, nPar.ij)
        Mdl.betaIdx[[i]][[j]][varSelCandCurr, ] <- FALSE
        Mdl.betaIdx[[i]][[j]][varSelInitCurr, ] <- TRUE
      }

###----------------------------------------------------------------------------
### Initial value for coefficients
###----------------------------------------------------------------------------
      betaInitCurr <- betaInit[[i]][[j]]
      if(class(betaInitCurr) == "character" &&
         tolower(betaInitCurr) == "random")
      {
        ## betaInitCurr = runif(ncolX.ij, -1, 1)
        ## Mdl.beta[[i]][[j]] <- array(betaInitCurr, c(ncolX.ij,
        ## 1)) Mdl.beta[[i]][[j]][!Mdl.betaIdx[[i]][[j]]] <- 0 let
        ## beta = 0 for non-selected variables NOTE: is this
        ## needed? -- YES.

        ## A simple version of random initial value. random for
        ## intercept and zero for the remaining values.
        Mdl.beta[[i]][[j]] <- array(0 , c(ncolX.ij, nPar.ij))
        Mdl.beta[[i]][[j]][1, ] <- runif(1, -1, 1)

      }
      else if (class(betaInitCurr) == "character" &&
               tolower(betaInitCurr) == "ols")
      {
        Y <- Mdl.Y[[i]]
        X <- Mdl.X[[i]][[j]]

        lmcoef <- lm(Y~0+X)$coef
        Mdl.beta[[i]][[j]] <- array(lmcoef, c(ncolX.ij, nPar.ij))
      }
      else # Do nothing and use user input
      {
        Mdl.beta[[i]][[j]] <- array(betaInitCurr, c(ncolX.ij, nPar.ij))
      }

    }
  }

  ## browser()
  out <- list(Mdl.beta = Mdl.beta,
              Mdl.betaIdx = Mdl.betaIdx)
  return(out)
}
