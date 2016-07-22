## This script make the posterior plot
sourceDir("~/code/cdcopula/R", recursive = TRUE)
load(file.path("~/running/", "^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross1MCMC.nIter1000twostage+20160522@23.11.c168ab.Rdata"))
iCross <- 1

###----------------------------------------------------------------------------
### Plot of Time Series
###----------------------------------------------------------------------------
par(mar = c(3, 4, 0, 0)+0.1, mfrow = c(2, 1))

Mdl.Y <- OUT.FITTED[[iCross]][["Mdl.Y"]]
ID <- as.Date(OUT.FITTED[[iCross]][["ID"]])

ylim <- c(-5, 5)
plot(ID, Mdl.Y[[1]], type = "l", col = "blue", ylab = "S&P 600", ylim = ylim)
plot(ID, Mdl.Y[[2]], type = "l", col = "blue", ylab = "S&P 100", ylim = ylim)

dev.copy2eps(file = "~/workspace/cdcopula/paper/SP100-SP600.eps")
###----------------------------------------------------------------------------
### Plot the Posterior Summary
###----------------------------------------------------------------------------

summary.Cpl <- CplMCMC.summary(OUT.MCMC = OUT.FITTED[[iCross]])
