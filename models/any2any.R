##' <description>
##'
##' <details>
##' @title <short tile>
##' @param densArgs 
##' @param linkType 
##' @return 
##' @references 
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Jan 12 00:05:00 CET 2012;
##'       Current: Thu Jan 12 00:05:07 CET 2012.
any2any <- function(densArgs, linkType)
  {

    ## Input distribution type
    inType <- densArgs[["input"]][["type"]]

    ## Output distribution type
    outType <- densArgs[["output"]][["type"]]

###----------------------------------------------------------------------------
### The input distribution 
###----------------------------------------------------------------------------    

    if(tolower(inType) == "beta") 
      {
        ## Beta distribution
        ## Subtract the parameter input in a convenient way
        mean <- densArgs[["input"]][["mean"]] # mean of beta density
        variance <- densArgs[["input"]][["variance"]] # Covariates
        ## shrinkage <- densArgs[["shrinkage"]] # Shrinkage
                
        ## Transform to standard parametrization 
        alpha <- - mean*(mean^2 -mean + variance)/
          (variance)
        beta <- (mean-1)^2*mean/(variance) + mean -1

        if(tolower(linkType) == "logit")
          {
            ## Assume this normal, See Beta-rep.nb
            meanLinked<- digamma(alpha)-digamma(beta)
            varLinked <- trigamma(alpha)+trigamma(beta)
          }
        else
          {
            stop("Not implemented yet!")
          }
      }
    else if(tolower(inType) == "lognorm")
      {
        ## Lognormal distribution.
        ## Subtract the parameter input in a convenient way
        mean <- densArgs[["input"]][["mean"]] # mean of beta density
        variance <- densArgs[["input"]][["variance"]] # Covariates
                
        ## Transform to standard parametrization
        sigma2 <- log(1+variance/mean^2)
        mu <- log(mean) - sigma2/2
          

        if(tolower(linkType) == "log")
          {
            meanLinked<- mu
            varLinked <- sigma2
            
          }


      }
    else if(tolower(inType) == "norm")
      {
        ## Lognormal distribution.
        ## Subtract the parameter input in a convenient way
        mean <- densArgs[["input"]][["mean"]] # mean of beta density
        variance <- densArgs[["input"]][["variance"]] # Covariates
        ## shrinkage <- densArgs[["shrinkage"]] # Shrinkage
        
        ## Transform to standard parametrization
        mu <- mean
        sigma2 <- variance

        if(tolower(linkType) == "identity")
          {
            meanLinked<- mu
            varLinked <- sigma2
          }
      }

###----------------------------------------------------------------------------
### The output distribution
###----------------------------------------------------------------------------    

    if(tolower(outType) == "norm") 
      {
        out <- list(mean = meanLinked,
                    variance = varLinked)
      }
###----------------------------------------------------------------------------
### The final output
###----------------------------------------------------------------------------    
    return(out)
    
  }









        
  ##               else if(tolower(linkCurr) == "identity")
  ##                 {
  ##                   ## TODO:
                    
  ##                 }
  ##               else if(tolower(linkCurr) == "log")
  ##                 {
  ##                   ## TODO:
                    
  ##                 }

  ##             }            
  ##           else if(tolower(priArgsCurr[["type"]]) == "norm")
  ##             {
  ##               ## Normal approximation with beta input

  ##               ## Subtract the prior information 
  ##               mean <- priArgsCurr[["mean"]]
  ##               variance <- priArgsCurr[["variance"]]
  ##               shrinkage <- priArgsCurr[["shrinkage"]]
  ##                                       # density

  ##               ## Transform to standard parametrization 
  ##               alpha <- - mean*(mean^2 -mean + variance*shrinkage)/
  ##                 (variance*shrinkage)
  ##               beta <- (mean-1)^2*mean/(variance*shrinkage) + mean -1
                
  ##               if(tolower(linkCurr) == "logit")
  ##                 {
  ##                   ## Assume this normal, 
  ##                   ## See Beta-rep.nb
                    
  ##                   meanLinked<- digamma(alpha)-digamma(beta)
  ##                   varLinked <- trigamma(alpha)+trigamma(beta)

  ##                   ## skewnessLinked <- (psigamma(alpha, 2)-psigamma(beta, 2))/
  ##                   ##   (trigamma(alpha)+trigamma(beta))^(3/2)
  ##                   outCurr[["beta"]][["intercept"]] <-
  ##                     sum(dnorm(xCurr, meanLinked, varLinked, log = TRUE))

  ##                 }
  ##               else if(tolower(linkCurr) == "identity")
  ##                 {
  ##                   ## TODO:
                    
  ##                 }
  ##               else if(tolower(linkCurr) == "log")
  ##                 {
  ##                   ## TODO:
                    
  ##                 }

  ##             }

  ##           else if
    
    
  ## }
