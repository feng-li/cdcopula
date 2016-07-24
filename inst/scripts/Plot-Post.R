## This script make the posterior plot
sourceDir("~/code/cdcopula/R", recursive = TRUE)
load(file.path("~/running/", "^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross1MCMC.nIter1000twostage+20160522@23.11.c168ab.Rdata"))
iCross <- 1

###----------------------------------------------------------------------------
### Plot of Time Series
###----------------------------------------------------------------------------
par(mar = c(2, 4, 0, 0)+0.1, mfrow = c(2, 1))

Mdl.Y <- OUT.FITTED[[iCross]][["Mdl.Y"]]
ID <- as.Date(OUT.FITTED[[iCross]][["ID"]])

ylim <- c(-5, 5)
plot(ID, Mdl.Y[[1]], type = "l", col = "blue", ylab = "S&P 600", ylim = ylim)

points(ID[Mdl.crossValidIdx[["testing"]][[iCross]]],
       OUT.PRED.MVSK[[iCross]][["mean"]][[1]], col = "red")

plot(ID, Mdl.Y[[2]], type = "l", col = "blue", ylab = "S&P 100", ylim = ylim)

dev.copy2eps(file = "~/workspace/cdcopula/paper/SP100-SP600.eps")
###----------------------------------------------------------------------------
### Plot the Posterior Summary
###----------------------------------------------------------------------------

## The basic model plot summary
summary.Cpl <- CplMCMC.summary(OUT.MCMC = OUT.FITTED[[iCross]])


## Plot time series of Kendall's tau (when lambdaL, lambdaU are parameterized)
MCMC.Update = OUT.FITTED[[iCross]][["MCMC.Update"]]
MCMC.burninProp <- OUT.FITTED[[iCross]][["MCMC.burninProp"]]
MCMC.sampleProp <- OUT.FITTED[[iCross]][["MCMC.sampleProp"]]
MCMC.nIter <- OUT.FITTED[[iCross]][["MCMC.nIter"]]
CplNM <- names(MCMC.Update)[length(MCMC.Update)]


n.burn <- round(MCMC.nIter*MCMC.burninProp)
MCMC.sampleIdx <- round(seq(n.burn+1, MCMC.nIter,
                            length.out = round((MCMC.nIter-n.burn)*MCMC.sampleProp)))

MCMC.par <- parCplMCMC(MCMC.beta = OUT.FITTED[[iCross]][["MCMC.beta"]],
                       Mdl.X = OUT.FITTED[[iCross]][["Mdl.X.training"]],
                       Mdl.parLink = OUT.FITTED[[iCross]][["Mdl.parLink"]],
                       MCMC.Update = MCMC.Update,
                       MCMC.sampleIdx = MCMC.sampleIdx)

summary.tau <- parCplMCMCSummary4Tau(MCMC.par)

plotCplParTS(MCMC.Update = list("BB7" = list("tau" = TRUE)),
             MCMC.parSummary = summary.tau,
             MdlDGP.par = NULL, ObsIdx4Plot = NA)

###----------------------------------------------------------------------------
### Plot the predictive summary
###----------------------------------------------------------------------------

plot(Mdl.Y[[1]], type = "l")
