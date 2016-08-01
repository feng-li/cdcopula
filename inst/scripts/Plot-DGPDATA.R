## This script makes the posterior plot for DGP data

sourceDir("~/code/cdcopula/R", recursive = TRUE)
load(file.path("~/running/",
               ## "DGP.Rdata" ## 0.3, 0.7 variational
               "M1M2+SPLITTSPLITTBB7+nObsNAnCross1MCMC.nIter1000+20160515@15.13.9ba81f.Rdata"

               ))


iCross <- 1

###----------------------------------------------------------------------------
### Plot the predictive summary for DGP data
###----------------------------------------------------------------------------


## The basic model plot summary
summary.Cpl <- CplMCMC.summary(OUT.MCMC = OUT.FITTED[[iCross]],
                               ObsIdx4Plot = 1:100)

plotCplParTS(MCMC.Update = MCMC.Update,
             MCMC.parSummary = summary.Cpl,
             MdlDGP.par = NULL, ObsIdx4Plot = 1:100)


## Hard coded for static model posterior plugin
## abline(h = 0.682, col = "black", lwd = 1)
## abline(h = 0.682+0.05, col = "black", lty = "dashed", lwd = 1)
## abline(h = 0.682-0.05, col = "black", lty = "dashed", lwd = 1)
## legend("bottomright", ncol = 2, legend = c("Posterior mean (No covariate-dependent)", "95% HPD"), lty = c("solid", "dashed"), col = "black", lwd = 2)

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
             MdlDGP.par = NULL, ObsIdx4Plot = 1:100)
