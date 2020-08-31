## This script makes the posterior plot for DGP data
rm(list = ls())
load(file.path("~/running/",
               ## "DGP2.Rdata" ## 0.3, 0.7 variational
               ## "M1M2+SPLITTSPLITTBB7+nObsNAnCross1MCMC.nIter1000+20160515@15.13.9ba81f.Rdata"
               "JOB899.M1M2+SPLITTSPLITTBB7+nObs1000nCross4MCMC.nIter1000joint+20160730@12.57.283761.Rdata"
               ))

source("~/code/cdcopula/R/flutils/R/systools/sourceDir.R")
sourceDir("~/code/cdcopula/R", recursive = TRUE)

iCross <- 4

###----------------------------------------------------------------------------
### Plot the predictive summary for DGP data
###----------------------------------------------------------------------------

## Convert the old data to adapt new model structure
if(length(OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]]) == 0)
{
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]] <- list()
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.burninProp"]] <- OUT.FITTED[[iCross]][["MCMC.burninProp"]]
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.sampleProp"]] <- OUT.FITTED[[iCross]][["MCMC.sampleProp"]]
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.Update"]] <- OUT.FITTED[[iCross]][["MCMC.Update"]]
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.nIter"]] <- OUT.FITTED[[iCross]][["MCMC.nIter"]]
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.parLink"]] <- OUT.FITTED[[iCross]][["Mdl.parLink"]]
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.crossValidIdx"]] <- OUT.FITTED[[iCross]][["Mdl.crossValidIdx"]]
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.Y"]] <- OUT.FITTED[[iCross]][["Mdl.Y"]]

}

MCMC.Update <- OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.Update"]]

## The basic model plot summary
summary.Cpl <- CplMCMC.summary(OUT.MCMC = OUT.FITTED[[iCross]],
                               ObsIdx4Plot = 1:100)

plotCplParTS(MCMC.Update = MCMC.Update,
             MCMC.parSummary = summary.Cpl[["par.summary"]],
             MdlDGP.par = NULL, ObsIdx4Plot = 1:100)


## Hard coded for static model posterior plugin
## abline(h = 0.682, col = "black", lwd = 1)
## abline(h = 0.682+0.05, col = "black", lty = "dashed", lwd = 1)
## abline(h = 0.682-0.05, col = "black", lty = "dashed", lwd = 1)
## legend("bottomright", ncol = 2, legend = c("Posterior mean (No covariate-dependent)", "95% HPD"), lty = c("solid", "dashed"), col = "black", lwd = 2)

MCMC.Update = OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.Update"]]
MCMC.burninProp <- OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.burninProp"]]
MCMC.sampleProp <- OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.sampleProp"]]
MCMC.nIter <- OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.nIter"]]
CplNM <- names(MCMC.Update)[length(MCMC.Update)]


n.burn <- round(MCMC.nIter*MCMC.burninProp)
MCMC.sampleIdx <- round(seq(n.burn+1, MCMC.nIter,
                            length.out = round((MCMC.nIter-n.burn)*MCMC.sampleProp)))

MCMC.par <- parCplMCMC(MCMC.beta = OUT.FITTED[[iCross]][["MCMC.beta"]],
                       Mdl.X = OUT.FITTED[[iCross]][["Mdl.X.training"]],
                       Mdl.parLink = OUT.FITTED[[iCross]][["Mdl.parLink"]],
                       MCMC.Update = MCMC.Update,
                       MCMC.sampleIdx = MCMC.sampleIdx)

summary.tau <- parCplMCMCSummary4Tau(MCMC.par,
                                     MCMC.sampleIdx = 1:length(MCMC.sampleIdx))

plotCplParTS(MCMC.Update = list("BB7" = list("tau" = TRUE)),
             MCMC.parSummary = summary.tau,
             MdlDGP.par = NULL, ObsIdx4Plot = 1:100)
