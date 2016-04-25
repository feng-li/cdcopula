## This script make the posterior plot
load(file.path("~/running/",
               "M1M2+SPLITTSPLITTBB7+nObs1000nCross1MCMC.nIter1000+20160425@20.57.34b75b.Rdata"))

iCross <- 1
summary.Cpl <- CplMCMC.summary(iIter = MCMC.nIter, MCMC.nIter = MCMC.nIter,
                               interval = 0.1, MCMC.burninProp = MCMC.burninProp,
                               OUT.MCMC = OUT.CplCross[[iCross]])
