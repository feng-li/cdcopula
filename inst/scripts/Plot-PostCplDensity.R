## This script plots the daily return of SP100 and SP600 data
load(file.path("~/running/", "^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross1MCMC.nIter1000twostage+20160522@23.11.c168ab.Rdata"))

sourceDir("~/code/cdcopula/R", recursive = TRUE)

iCross <- 1


MCMC.Update = OUT.FITTED[[iCross]][["MCMC.Update"]]
MCMC.burninProp <- OUT.FITTED[[iCross]][["MCMC.burninProp"]]
MCMC.sampleProp <- OUT.FITTED[[iCross]][["MCMC.sampleProp"]]
MCMC.nIter <- OUT.FITTED[[iCross]][["MCMC.nIter"]]

CplNM <- names(MCMC.Update)[length(MCMC.Update)]


## n.burn <- round(MCMC.nIter*MCMC.burninProp)
## MCMC.sampleIdx <- round(seq(n.burn+1, MCMC.nIter,
##                             length.out = round((MCMC.nIter-n.burn)*MCMC.sampleProp)))
## MCMC.par <- parCplMCMC(MCMC.beta = OUT.FITTED[[iCross]][["MCMC.beta"]],
##                        Mdl.X = OUT.FITTED[[iCross]][["Mdl.X.training"]],
##                        Mdl.parLink = OUT.FITTED[[iCross]][["Mdl.parLink"]],
##                        MCMC.Update = OUT.FITTED[[iCross]][["MCMC.Update"]],
##                        MCMC.sampleIdx = MCMC.sampleIdx)


summary.Cpl <- CplMCMC.summary(OUT.MCMC = OUT.FITTED[[iCross]])

lambdaL <- summary.Cpl[["par.summary"]][["ts.mean"]][[CplNM]][["lambdaL"]]
lambdaU <- summary.Cpl[["par.summary"]][["ts.mean"]][[CplNM]][["lambdaU"]]

###----------------------------------------------------------------------------
###
###----------------------------------------------------------------------------

## The copula density scaling
nPoints <- 20
u1 <- seq(0, 1, length.out = nPoints)
u2 <- u1
u <- mesh.grid(u1)


delta <- -log(2)/log(lambdaL)
theta <- log(2)/log(2-lambdaU)

##########################

xlim <- c(0, 1)
ylim <- c(0, 1)


nPlot <- 5
par(mfrow = c(nPlot, 4), mar = c(4, 2, .5, .5))

ad <- order(lambdaL, decreasing = T)

## ID.cand <- cbind(ad[seq(1, 1000, length.out = nPlot)],
##                  ad[seq(length(ad)-200, length(ad), length.out = nPlot)])
normaldays <- c(1, 3, 5, 7, 9)
## ID.card <- 3000:3500
j <- 0
for(i in as.vector(t(ID.cand)))
  {
    j <- j+1
    theta1 <- theta[i]
    delta1 <- delta[i]

    empDist <- matrix(pCpl(u = u, CplNM = CplNM,
                           parCpl = c(theta = theta1, delta = delta1)), nPoints)
    empDens <- matrix(dCpl(u = u, CplNM = CplNM,
                           parCpl = c(theta = theta1, delta = delta1)), nPoints)

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
