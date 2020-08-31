## This script makes the posterior plot
## load(file.path("~/running/", "JOB835.^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross0MCMC.nIter10000joint+20160720@14.44.5c1da7.Rdata"))

load(file.path("~/running/","JOB835.^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross0MCMC.nIter10000joint+20160720@14.44.5c1da7-withsummary.Rdata"))

sourceDir("~/code/cdcopula/R", recursive = TRUE)

iCross <- 1

###----------------------------------------------------------------------------
### Plot of Time Series
###----------------------------------------------------------------------------

Mdl.Y <- OUT.FITTED[[iCross]][["Mdl.Y"]]
ID <- as.Date(OUT.FITTED[[iCross]][["ID"]])

ylim <- c(-5, 5)

## Plot of returns
par(mar = c(4, 2, 0, 0)+0.1, mfrow = c(1, 2), cex = 0.7)

plot(ID, Mdl.Y[[1]], type = "l", col = "blue", xlab = "S&P 600", ylab = "", ylim = ylim)
plot(ID, Mdl.Y[[2]], type = "l", col = "blue", xlab = "S&P 100", ylab = "", ylim = ylim)


## Plot of fitted values and HPD
par(mar = c(2, 4, 0, 0)+0.1, mfrow = c(2, 1))

plot(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
       OUT.PRED.MVSK[[iCross]][["mean"]][[1]][, 1],
       type = "l", col = "red", ylim = c(-1, 2))

points(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
       OUT.PRED.MVSK[[iCross]][["HPDL"]][[1]][, 1],
       type = "l", col = "grey")
points(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
       OUT.PRED.MVSK[[iCross]][["HPDU"]][[1]][, 1],
       type = "l", col = "grey")


plot(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
       OUT.PRED.MVSK[[iCross]][["mean"]][[2]][, 1],
       type = "l", col = "red", ylim = c(-1, 2))

points(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
       OUT.PRED.MVSK[[iCross]][["HPDL"]][[2]][, 1],
       type = "l", col = "grey")
points(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
       OUT.PRED.MVSK[[iCross]][["HPDU"]][[2]][, 1],
       type = "l", col = "grey")


## Residual plot
par(mar = c(2, 4, 0, 0)+0.1, mfrow = c(2, 1))
plot(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
     OUT.PRED.RESID[[iCross]][["mean"]][[1]],
     type = "l", col = "blue")
plot(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
     OUT.PRED.RESID[[iCross]][["mean"]][[2]],
     type = "l", col = "blue")


## dev.copy2eps(file = "~/workspace/cdcopula/paper/SP100-SP600.eps")
###----------------------------------------------------------------------------
### Plot the Posterior Summary
###----------------------------------------------------------------------------

## The basic model plot summary
## summary.Cpl <- CplMCMC.summary(OUT.MCMC = OUT.FITTED[[iCross]])


## Plot time series of Kendall's tau (when lambdaL, lambdaU are parameterized)
MCMC.Update = OUT.FITTED[[iCross]][["MCMC.Update"]]
MCMC.burninProp <- OUT.FITTED[[iCross]][["MCMC.burninProp"]]
MCMC.sampleProp <- OUT.FITTED[[iCross]][["MCMC.sampleProp"]]
MCMC.nIter <- OUT.FITTED[[iCross]][["MCMC.nIter"]]
CplNM <- names(MCMC.Update)[length(MCMC.Update)]


n.burn <- round(MCMC.nIter*MCMC.burninProp)
MCMC.sampleIdx <- round(seq(n.burn+1, MCMC.nIter,
                            length.out = round((MCMC.nIter-n.burn)*MCMC.sampleProp)))

plotCplParTS(MCMC.Update = MCMC.Update,
             MCMC.parSummary = summary.Cpl[["par.summary"]],
             MdlDGP.par = NULL, ObsIdx4Plot = NA)


MCMC.par <- parCplMCMC(MCMC.beta = OUT.FITTED[[iCross]][["MCMC.beta"]],
                       Mdl.X = OUT.FITTED[[iCross]][["Mdl.X.training"]],
                       Mdl.parLink = OUT.FITTED[[iCross]][["Mdl.parLink"]],
                       MCMC.Update = MCMC.Update,
                       MCMC.sampleIdx = MCMC.sampleIdx)

summary.tau <- parCplMCMCSummary4Tau(MCMC.par, MCMC.sampleIdx = 1:length(MCMC.sampleIdx))
plotCplParTS(MCMC.Update = list("BB7" = list("tau" = TRUE)),
             MCMC.parSummary = summary.tau,
             MdlDGP.par = NULL, ObsIdx4Plot = NA)


###----------------------------------------------------------------------------
### Print the Posterior summary with LaxTex marked
###----------------------------------------------------------------------------

options(width = 120)
rapply(summary.Cpl[["beta.mean"]], function(x) round(x, 3), how = "replace")
rapply(summary.Cpl[["betaIdx.mean"]], function(x) round(x, 2), how = "replace")


mean(unlist(summary.Cpl[["beta.ineff"]]))
mean(unlist(summary.Cpl[["accept.prob.mean"]]))
