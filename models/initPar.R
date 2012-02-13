##' Setup the initial values for the copula model
##'
##' This function is only once for the first iteration. The values are updated
##' via MCMC scheme.
##' @title Setup copula model initial value.
##' @param varSelArgs 
##' @param betaInit 
##' @param Mdl.X 
##' @return "list"
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Dec 22 15:57:14 CET 2011;
##'       Current: Wed Jan 11 16:03:44 CET 2012.
initPar <- function(varSelArgs, betaInit, Mdl.X)
  {
    ## The output structure.
    Mdl.betaIdx <- betaInit
    Mdl.beta <- betaInit
    
    ## Loop to assign the initial values
    CompNM <- names(Mdl.X)
    for(i in CompNM)
      {
        CompParNM <- names(Mdl.X[[i]])
        for(j in CompParNM)
          {
            ## No. of col for covariates, including intercept if applicable.
            ncolX.ij <- ncol(Mdl.X[[i]][[j]])
###----------------------------------------------------------------------------
### Initialize the variable section indicator
###----------------------------------------------------------------------------            
            ## Initial value for variable selection indicators which can be
            ## character or vector
            varSelInitCurr <- varSelArgs[[i]][[j]][["init"]]
            
            if(class(varSelInitCurr) == "character" &&
               tolower(varSelInitCurr) == "all-in")
              {
                Mdl.betaIdx[[i]][[j]] <- array(1, c(ncolX.ij, 1))
              }
            else if(class(varSelInitCurr) == "character" &&
                    tolower(varSelInitCurr) == "all-out")
              {
                Mdl.betaIdx[[i]][[j]] <- array(1, c(ncolX.ij, 1))

                varSelCandCurr <- varSelArgs[[i]][[j]][["cand"]]
                Mdl.betaIdx[[i]][[j]][varSelCandCurr] <- 0
              }
            else if(class(varSelInitCurr) == "character" &&
                    tolower(varSelInitCurr) == "random")
              {
                Mdl.betaIdx[[i]][[j]] <- array(1, c(ncolX.ij, 1))
                varSelCandCurr <- varSelArgs[[i]][[j]][["cand"]]
                ## Randomly pick up half in
                betaIdxCurrSubOut <- sample(varSelCandCurr,
                                            round(length(varSelCandCurr)/2))
                Mdl.betaIdx[[i]][[j]][betaIdxCurrSubOut] <- 0
              }     
            else # Do nothing, use user input
              {
                Mdl.betaIdx[[i]][[j]] <- array(1, c(ncolX.ij, 1))
                Mdl.betaIdx[[i]][[j]][varSelInitCurr] <- 1

              }
            
###----------------------------------------------------------------------------
### Initial value for coefficients
###----------------------------------------------------------------------------
            betaInitCurr <- betaInit[[i]][[j]]
            if(class(betaInitCurr) == "character" &&
               tolower(betaInitCurr) == "random")
              {
                betaInitCurr = runif(ncolX.ij, -1, 1)
                Mdl.beta[[i]][[j]] <- array(betaInitCurr, c(ncolX.ij, 1))
                Mdl.beta[[i]][[j]][Mdl.betaIdx[[i]][[j]] == 0] <- 0
                ## let beta  =  0 for nu-selected variables NOTE: is this
                ## needed? -- YES.
              }
            else # Do nothing and use user input
              {
                Mdl.beta[[i]][[j]] <- array(betaInitCurr, c(ncolX.ij, 1))
              }  

          }
      }
    
    
    out <- list(Mdl.beta = Mdl.beta,
                Mdl.betaIdx = Mdl.betaIdx)
    return(out)
  }
