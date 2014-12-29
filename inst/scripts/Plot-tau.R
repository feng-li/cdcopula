## This script plots the contour plot for Kendall's tau with respect to lower
## and upper tail dependence.

CplNM <- "BB7"
xlim <- c(0, 1)
ylim <- c(0, 1)
xat <- seq(0, 1, 0.1)
yat <- xat

lambdaL <- seq(0.001, 0.999, 0.01)
lambdaU <- lambdaL

lambdaLU <- mesh.grid(lambdaL, lambdaU)

delta <- -log(2)/log(lambdaLU[, 1])
theta <- log(2)/log(2-lambdaLU[, 2]) # ff(delta)

parCpl <- list(delta = delta, theta = theta)

tau <- matrix(kendalltau(CplNM, parCpl), length(lambdaL))

###----------------------------------------------------------------------------
### Debuging code to check if inverse Kendall's tau is correctly calculated
###----------------------------------------------------------------------------
lambdaU.e <- kendalltauInv(CplNM = "BB7",
                         parRepCpl = list(lambdaL = lambdaLU[, 1], tau = tau),
                         method = "tabular")
theta.e <- log(2)/log(2-lambdaU.e)
## tau.e <- kendalltau(CplNM = "BB7",
##                     parCpl = list(theta = theta1, delta = delta))
idx <- which(abs(lambdaU.e-lambdaLU[, 2])>0.1)



###----------------------------------------------------------------------------
### Plot counter
###----------------------------------------------------------------------------
## Colored contour
filled.contour(x = lambdaL, y = lambdaU, z = tau,
               xlim = xlim, ylim = xlim,
               xlab = expression(lambda[L]),
               ylab = expression(lambda[U]),
               key.title = title(main = expression(tau)),
               key.axes  =  {axis(4,  at  =  seq(0,  1,  0.1))},
               plot.axes = {
                 axis(1,  at  =  xat)
                 axis(2,  at  =  yat)

                 contour(lambdaL,  lambdaU,  tau,  add  =  TRUE,
                         col  = "black",  lwd  =  0.8,  lty  = "solid",
                         xlim  = xlim,  ylim  =  ylim)

                 ## DEBUGGING CODE:
                 points(lambdaLU[idx, ], col = "blue", pch = 20)

               })


## Contour line only
contour(lambdaL,  lambdaU,  tau,  add  =  FALSE,
        col  = "black",  lwd  =  0.8,  lty  = "solid",
        xlim  = xlim,  ylim  =  ylim)
