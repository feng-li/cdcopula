## This script make the posterior plot
load(file.path("~/running/", "^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross1MCMC.nIter1000twostage+20160522@23.11.c168ab.Rdata"))
iCross <- 1
summary.Cpl <- CplMCMC.summary(OUT.MCMC = OUT.FITTED[[iCross]])
