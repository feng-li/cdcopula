##' Generate random variable from a positive stable distribution
##'
##' Positive stable distribution generator
##' @title Positive stable distribution generator
##' @param n 
##' @param alpha alpha = 1/theta
##' @return 
##' @references
##'     Chambers et al (1976)
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note
##'     Created: Mon Sep 26 19:05:37 CEST 2011;
##'     Current: Mon Sep 26 19:05:43 CEST 2011.
rps <- function(n, alpha)
  {
    eta <- runif(n, 0, pi)
    omega <- rexp(n, rate = 1)
    z <- (sin((1-alpha)*eta)*sin(alpha*eta)^(alpha/(1-alpha)))/
      sin(eta)^(1/(1-alpha))
    gamma <- (z/omega)^((1-alpha)/alpha)
    return(gamma)
  }
