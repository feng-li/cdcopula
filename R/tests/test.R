
## for(i in 1:10)
##   {

##     if(i == 5) next
##     print(i)

##   }

################################################################################
### TESTING
################################################################################

## Sigma <- cor(matrix(rnorm(100), 20))
## theta <- Sigma[lower.tri(Sigma)]

## copula = "mvt"
## par = list("df" = 5)
## n <- 10

## u <- ruCpl(n = n, theta = theta, copula = copula, par = par)
## c <- cCpl(u = u, theta = theta, copula = copula, par = par)
## p <- uCpl(u = u, theta = theta, copula = copula, par = par)

## ###----------------------------------------------------------------------------
## ### BB7 Copula
## ###----------------------------------------------------------------------------
## n <- 100
## theta <- c(1.2, 1)
## copula <- "BB7"

## U0 <- ruCpl(n = n, theta = theta, copula = copula)
## c0 <- cCpl(u = U0$u, theta = theta, copula = copula)

## x0 <- qnorm(U0$u)
## xlim <- c(-3, 3)
## par(mfrow = c(2, 2))
## plot(U0$u, main = round(c(U0$lambdaL, U0$lambdaU), 2))
## plot(x0, xlim = xlim, ylim = xlim)

## require("mvtnorm")
## rho <- 0.6
## plot(rmvnorm(n, sigma = matrix(c(1, rho, rho, 1), 2)), xlim = xlim, ylim =
##      xlim, main = rho)

## theta1 <- theta[1]
## theta2 <- theta[2]

## p <- 2 #Hard code for bivariate copula
## n <- 5000
## v <- matrix(runif(n*p, 0, 1), n, p)
## theta <- c(0.5, 0.8)
## u <- matrix(NA, n, p)
## u[, 1] <- v[, 1]

## for(i in 1:n)
## {
##   u[i, 2] <- uniroot(function(x, theta, v0)
##           {
##             theta1 = theta[1]
##             theta2 = theta[2]
##             u1 <- v0[1]
##             v2 <- v0[2]
##             u2 <- x

##             L1 <- 1-(1-u1)^theta1
##             L2 <- 1-(1-u2)^theta1
##             L5 <- -1 + L1^(-theta2) + L2^(-theta2)
##             L6 <- 1- L5^(-1/theta2)
##             out <- L1^(-1-theta2)*L5^(-1-1/theta2)*
##               L6^(-1+1/theta1)*(1-u1)^(-1+theta1)-v2
##             return(out)
##           },
##           interval = c(0, 1), theta = theta, v0 = v[i, ])[["root"]]
## }
## print(proc.time()-a)


## Apply type function will not save much time with self-defined functions.
## gc()
## a <- proc.time()
## out <- apply(matrix(1:n), 1,
##       function(x, theta, v)
##       {
##         uniroot(function(x, theta, v0)
##                 {
##                   theta1 = theta[1]
##                   theta2 = theta[2]
##                   u1 <- v0[1]
##                   v2 <- v0[2]
##                   u2 <- x

##                   L1 <- 1-(1-u1)^theta1
##                   L2 <- 1-(1-u2)^theta1
##                   L5 <- -1 + L1^(-theta2) + L2^(-theta2)
##                   L6 <- 1- L5^(-1/theta2)
##                   out <- L1^(-1-theta2)*L5^(-1-1/theta2)*
##                     L6^(-1+1/theta1)*(1-u1)^(-1+theta1)-v2
##                   return(out)
##                 },
##                 interval = c(0, 1), theta = theta, v0 = v[x, ])[["root"]]
##       },theta = theta, v = v)

## print(proc.time()-a)
