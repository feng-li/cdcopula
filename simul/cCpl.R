##' Compute the copula density 
##'     
##' @title Calculate the density function for give bivariate copula.
##' 
##' @param u "matrix"
##'     Each column of the matrix is the value in the unit interval of the real
##' line for the corresponding margin of the copula function.
##' @param theta "list";
##'     The parameter list supplied to the copula function.
##' @param copula "character";
##'     This is the option for different types of copulas used in the function.
##' Currently "gaussian" for Gaussian copula.
##' @param par "list";
##'     Any further parameters need to pass to the copula. Which is equivalent
##' to the ... argument but less possible to make error.
##' @return "list" see below.
##'     \item {copula} {"matrix"; The copula function with the same length of
##' that in each column entries in "u".}
##'     \item {tao} {The theoretical Kendall's tao.}
##'     \item {rho} {The theoretical Spearman's rho.}
##' @references
##'     Trivedi and Zimmer (2007).
##' @author
##'     Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     Created: Mon Sep 12 13:36:10 CEST 2011;
##'     Current: Tue Sep 13 11:38:16 CEST 2011.
##' TODO: the first argument maybe should be x instead of u,
##'       Return log density instead of current form.
cCpl <- function(u, theta, copula, par = par)
  {
    if(tolower(copula) == "bb7")
      {
        u1 <- u[, 1]
        u2 <- u[, 2]
        theta1 <- theta[1]
        theta2 <- theta[2]
        
        ## Local temporal variables 
        L1 <- 1-(1-u1)^theta1
        L2 <- 1-(1-u2)^theta1
        L3 <- (1-u1)^(-1+theta1)
        L4 <- (1-u2)^(-1+theta1)
        L5 <- -1 + L1^(-theta2) + L2^(-theta2)
        L6 <- 1 - L5^(-1/theta2)
        
        density <- (L1*L2)^(-1-theta2)*L3*L4*L5^(-2*(1+theta2)/theta2)*
          L6^(-2+1/theta1)*(-1+theta1+L5^(1/theta2)*L6*(1 + theta2)*theta1)

        out <- matrix(density)
      }
    else if(tolower(copula) == "gaussian")
      {                
        ## The quantile for normal CDF 
        u.quantile <- qnorm(u)
        x1 <- u.quantile[, 1]
        x2 <- u.quantile[, 2]
        ## x1 <- X[, 1]
        ## x2 <- X[, 2]

        ## The copula density function C_12(u1, u2)
        ## TODO: verify this
        density <- 1/sqrt(1-theta^2)*
          exp(-1/2*1/(1-theta^2)*(x1^2+x2^2-2*theta*x1*x2))*
            exp(1/2*(x1^2+x2^2))
        
        ## The output
        out <- matrix(density)
      }
    else if(tolower(copula) == "mvt") # The multivariate t-copula
      {
        df <- par[["df"]]
        p <- ncol(u)
        ## The quantile for *univariate* t 
        u.quantile <- qt(u, df = df)
        ## u.quantile <- X
        
        P.mat <- matrix(1, p, p) # The correlation matrix
        P.mat[lower.tri(P.mat)] <- theta
        P.mat[upper.tri(P.mat)] <- t(P.mat)[upper.tri(P.mat)]
        corr <- P.mat # The correlation matrix with scale 1.

        ## The copula density function C_12(u1, u2)  
        dst0 <- dmvt(u.quantile, sigma = corr, df = df)
        dst1 <- matrix(apply(dt(u.quantile, df = df), 1, prod))
        density <- dst0/dst1
        
        ## The output
        out <- matrix(density)
      }
   else if(tolower(copula) == "fgm")
      {
        u1 <- u[, 1]
        u2 <- u[, 2]
        density <- 1+theta*(1-2*u1)*(1-2*u2)
        out <- density
      }
    else if(tolower(copula) == "gumbel")
      {
        pctl <- uCpl(u, theta, copula = "gumbel")
        u.tilde <- -log(u)
        
        u.tildeProd <- u.tilde[, 1]*u.tilde[, 2]
        u.tildeSumtheta <- rowSums(u.tilde^theta)
        
        density <- pctl/(u[, 1]*u[, 2])*u.tildeProd^(theta-1)/
          u.tildeSumtheta^(2-1/theta)*(u.tildeSumtheta^(1/theta)+theta-1)       
        out <- density
      }    
    else if(tolower(copula)  == "frank")
      {                
        ## The quantile for normal CDF 
        u1 <- u[, 1]
        u2 <- u[, 2]
        ## The copula function
        density <- -1/theta*log(1+(exp(-theta*u1)-1)*(exp(-theta*u2)-1)/
                                (exp(-theta)-1))
        ## The output
        out <- density
      }
      else stop("Given copula name is not implemented.")
    
    return(out)
  }
