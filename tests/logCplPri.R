##' The log prior for the parameters in the copula density
##'
##' <details>
##' @title <short tile>
##' @param par "list"
##' @return 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Sun Oct 23 18:40:55 CEST 2011;
##'       Current: Sun Oct 23 18:41:03 CEST 2011.
logCplPri <- function(beta, priCplArgs, link)
  {
    ## The tail dependence parameters [0, 1]
    ## theta <- parCpl[["theta"]] #  the lower tail dependence 
    ## delta <- parCpl[["delta"]] # the upper tail dependence

###----------------------------------------------------------------------------
### PRIOR FOR KENDALL'S TAU
###----------------------------------------------------------------------------
    beta0 <- beta[1]
    priCplIntArgs <- priCplArgs[["intercept"]]
    
    ## Prior for intercept (via link function)
    logCplIntPri <- logPri(beta = beta0, priArgs=priCplIntArgs, link = link)
    
    ## Prior for slopes (via link function)
    logCplSlopePri <- logPri(...)

###----------------------------------------------------------------------------
### PRIOR FOR LOWER TAIL DEPENDENT PARAMETER lambda_L
###----------------------------------------------------------------------------

    ## Prior for intercept (via link function)
    beta0 <- beta[1] 
    logCplIntPri <- logPri(beta = beta0, priArgs=priIntArgs, link = link)
    
    ## Prior for slopes (via link function)
    logCplSlopePri <- logPri(...) 

###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------
    ## The output
    out <- logCplIntPri + logCplSlopePri
    return(out) 
  }
###----------------------------------------------------------------------------
### Testing 
###----------------------------------------------------------------------------
## par <- list(theta = 2, delta = 2)
## logCplPri(par)
