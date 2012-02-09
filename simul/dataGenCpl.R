##' <description>
##'
##' <details>
##' @title <short tile>
##' @param nObs 
##' @param MargisTypes 
##' @param MdlData.betaIdx 
##' @param MdlData.beta 
##' @return 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: ; Current: .
dataGenCpl <- function(nObs, MargisTypes,  MdlData.betaIdx,  MdlData.beta)
  {
    ## DATA COVARIATES for the marginal distrib 
    covarNoInt <- list(matrix(runif(nObs*4), nObs),
                       matrix(runif(nObs*5), nObs))
    
    ## Covariates with intercept
    covar <- lapply(covarNoInt, function(x) cbind(1, x))

    
  }
