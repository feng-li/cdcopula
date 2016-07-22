## This script plots the daily return of SP100 and SP600 data
iCross <- 1
MCMC.par <- parCplMCMC(MCMC.beta = OUT.FITTED[[iCross]][["MCMC.beta"]],
                       Mdl.X = OUT.FITTED[[iCross]][["Mdl.X.training"]],
                       Mdl.parLink = OUT.FITTED[[iCross]][["Mdl.parLink"]],
                       MCMC.Update = OUT.FITTED[[iCross]][["MCMC.Update"]],
                       MCMC.sampleIdx = MCMC.sampleIdx)


## Plot the empirical copula
nPoints <- 20
u1 <- seq(0, 1, length.out = nPoints)
u2 <- u1
u <- mesh.grid(u1)

## empDist <- matrix(hatCpl(u = u[, 1], v = u[, 2],
##                          x = Y[[1]], y = Y[[2]]),
##                   nPoints)

## tau <- MCMC.par[[3]][["tau"]]
lambdaL <- MCMC.par[[3]][["lambdaL"]]
lambdaU <- MCMC.par[[3]][["lambdaL"]]
## lambdaU <- as.vector(kendalltauInv(
##           CplNM = CplNM, parRepCpl = MCMC.par[[3]]))

      ## cat(lambdaL[1], lambdaU[1], tau[1], "\n")
      ## points( lambdaL[1], lambdaU[2], col = "red", pch = 20)
      ## Sys.sleep(0.5)

      ## if(tau[1]>0.9) browser()

      ## cat("tau", tau, "lambdaL", lambdaL, "lambdaU", lambdaU, "\n")
      ## The standard copula parameters (recycled if necessary, should not have
      ## dimension attributed).

delta <- -log(2)/log(lambdaL)
theta <- log(2)/log(2-lambdaU)

##########################

xlim <- c(0, 1)
ylim <- c(0, 1)

## frechet.l <- matrix(uCpl(u = u, CplNM = "frechet-lower"), nPoints)
## frechet.u <- matrix(uCpl(u = u, CplNM = "frechet-upper"), nPoints)

nPlot <- 5
par(mfrow = c(nPlot, 4), mar = c(4, 2, .5, .5))

##ID.cand <- c(seq(1, 3000, length.out = 5), seq(3001, nTraining, length.out = 5))
ad <- order(delta)
ID.cand <- cbind(ad[seq(1, 1000, length.out = nPlot)], ad[seq(length(ad)-200, length(ad), length.out = nPlot)])
normaldays <- c(1, 3, 5, 7, 9)
## ID.card <- 3000:3500
j <- 0
for(i in as.vector(t(ID.cand)))
  {
    j <- j+1
    theta1 <- theta[i]
    delta1 <- delta[i]

    empDist <- matrix(pCpl(u = u, CplNM = CplNM,
                           theta = c(theta1, delta1)), nPoints)
    empDens <- matrix(dCpl(u = u, CplNM = CplNM, theta = c(theta1, delta1)), nPoints)

    if(j %in% normaldays)
      {
        col = "blue"
      }
    else
      {
        col = "red"
      }

    contour(as.matrix(u1),  as.matrix(u2),  empDist,  add  =  FALSE,
            col  = col,  lwd  =  0.8,  lty  = "solid",
            xlim  = xlim,  ylim  =  ylim,
            xlab = "", ylab = "")
    contour(as.matrix(u1),  as.matrix(u2),  empDens,  add  =  FALSE,
            col  = col,  lwd  =  0.8,  lty  = "solid",
            xlim  = xlim,  ylim  =  ylim,
            xlab = "", ylab = "")


  }
