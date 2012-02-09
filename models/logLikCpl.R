##' This function make the log-likelihood function of the copula component
##'
##' This function only update the copula density part. The whole log-likelihood
##' should consist of two parts: copula and marginal parts.
##' @title Compute the copula likelihood of copula function
##' @param u 
##' @param copula 
##' @param par 
##' @return "matrix"; 
##' @references Joe 1997, p. 153
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Oct 20 18:15:13 CEST 2011;
##'       Current: Thu Oct 20 18:19:34 CEST 2011. 
logLikCpl <- function(u, CplNM, parCpl, extArgs)
{
  ## The sum of log copula density
  if(tolower(CplNM) == "bb7")
    {
      ## The name of marginal model
      MargisNM <- names(u)
      
      ## Subtract the parameters list. 
      tau <- parCpl[["tau"]]
      lambdaL <- parCpl[["lambdaL"]]
      
      ## Transform the parameters into the standard form
      ## FIXME: Consider to speed it up if it is really slow
      parOut <- kendalltauInv(CplNM = CplNM, parRepCpl = parCpl,
                              tauTabular = extArgs[["tauTabular"]])
      delta <- parOut[["delta"]]
      theta <- parOut[["theta"]]

      ## delta <- -log(2)/log(lambda_L)
      ## theta <- log(2)/log(2-lambda_U)
      
      ## temporal data
      ## L1 <- 1-(1-u1)^theta
      ## L2 <- 1-(1-u2)^theta
      TC1 <- 1-(1-u)^theta
      ## L1 <- TC1[, 1]
      ## L2 <- TC1[, 2]
      
      ## L3 <- (1-u1)^(-1+theta)
      ## L4 <- (1-u2)^(-1+theta)
      TC2 <- (1-u)^(-1+theta)
      ## L3 <- TC2[, 1]
      ## L4 <- TC2[, 2]
      
      ## L5 <- -1+L1^(-delta)+L2^(-delta)
      L5 <- -1 + rowSums(TC1^(-delta))
      
      L6 <- 1-L5^(-1/delta) # FIXME: log(L6)->Inf when u->1,  v->1.
        
      logCpl <- (-1-delta)*rowSums(log(TC1))+
        rowSums(log(TC2))-
          2*(1+delta)/delta*log(L5)+
            (-2+1/theta)*log(L6)+
              log(-1+theta+L5^(1/delta)*L6*(1+delta)*theta)
      loglik <- sum(logCpl)
      
      out <- loglik
      ## if(is.infinite(out)) browser()
      
      }
  else if(tolower(CplNM) == "gaussian")
    {
      ## The Gaussian copula
    }
  else
    {
      stop("The copula is not implemented yet!")
    }  
  return(out)
  ## The sum of log marginal density   
}
