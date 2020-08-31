#! /usr/bin/env Rscript

## This script makes the posterior plot
rm(list = ls())
require("methods", quietly = TRUE)
require("parallel", quietly = TRUE)

Models <- list()
## Models[["Joe-Clayton"]] <- "~/running/^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross1MCMC.nIter1000joint+20160523@22.58.f299f8.Rdata"
## Models[["Clayton"]] <- "~/running/^SML^OEX+SPLITTSPLITTCLAYTON+nObs6557nCross1MCMC.nIter5000joint+20160523@13.43.d13c0b.Rdata"
## Models[["Gumbel"]] <- "~/running/^SML^OEX+SPLITTSPLITTGUMBEL+nObs6557nCross1MCMC.nIter5000joint+20160523@23.04.fbf7d2.Rdata"
## Models[["bivariate DCC-GARCH"]] <- "~/running/^SML^OEX+NANADCCGARCH+nObsNAnCross1+20160519@20.58.58aebe.Rdata"

Models[["BB7+MixtureT"]] = "~/downloads/JOB16640.^SML^MID+TEIGENTEIGENBB7+nObs1490nCross1MCMC.nIterNA+20190402@15.17.44b61f.Rdata"

source("~/code/flutils/R/systools/sourceDir.R")
sourceDir("~/code/cdcopula/R", "~/code/flutils/R", recursive = TRUE)
iCross <- 1


## Make a parallel cluster
cl <- makeCluster(detectCores())
setDefaultCluster(cl)
ce <- clusterEvalQ(cl,{
    source("~/code/flutils/R/systools/sourceDir.R")
    sourceDir("~/code/cdcopula/R", "~/code/flutils/R", recursive = TRUE)
})

weights <- c(0.5, 0.5)
probs <- c(0.05, 0.01)

## Extract VARS
VARS <- list()
for(iModel in 1:length(Models))
{
    load(Models[[iModel]])
    OUT.FITTED[[iCross]][["MCMC.sampleProp"]] <- 0.1
    ## rm(rCpl)
    sourceDir("~/code/cdcopula/R", "~/code/flutils/R", recursive = TRUE)

    ## Convert the old data to adapt new model structure
    if(length(OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]]) == 0)
    {
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]] <- list()
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.burninProp"]] <- OUT.FITTED[[iCross]][["MCMC.burninProp"]]
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.sampleProp"]] <- 0.1
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.Update"]] <- OUT.FITTED[[iCross]][["MCMC.Update"]]
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.nIter"]] <- OUT.FITTED[[iCross]][["MCMC.nIter"]]
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.parLink"]] <- OUT.FITTED[[iCross]][["Mdl.parLink"]]
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.crossValidIdx"]] <- OUT.FITTED[[iCross]][["Mdl.crossValidIdx"]]
        OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.Y"]] <- OUT.FITTED[[iCross]][["Mdl.Y"]]
    }


    Mdl.crossValidIdx <- OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.crossValidIdx"]]
    Mdl.Idx.testing <- Mdl.crossValidIdx[["testing"]][[iCross]]

    ## Setting the MCMC sample to reduce computing time.
    OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["MCMC.sampleProp"]] <- 0.1


    ## Special case when a foreign multivariate model is introduced. Use the model's
    ## predict method.
    ForeignModelSpec <- OUT.FITTED[[iCross]][["ForeignModelSpec"]]
    if(length(ForeignModelSpec)>0 &&  class(ForeignModelSpec) != "logical")
    {
        replicate = 1
    } else
    {
        replicate = 10
    }

    VARS[[iModel]] <- CplMCMC.VaR(CplFitted = OUT.FITTED[[iCross]], probs = probs,
                                  Mdl.Idx.testing = Mdl.Idx.testing,
                                  args = list(method = "linear",
                                              weights = weights,
                                              plot = FALSE,
                                              replicate = replicate))

    cat("Calculating VaR for model ", names(Models)[iModel], ", ",
        iModel, " of ", length(Models), "model done!\n")
}

nModels <- length(VARS)

Portfolio.Real <- VARS[[1]][["Portfolio.Real"]]

Portfolio.LeftOut <- sapply(as.list(1:nModels),
                            function(i, x) rowSums(VARS[[i]][["VaR"]]>x)/length(x), x = Portfolio.Real)

## Plot VARS
Mdl.crossValidIdx <- OUT.FITTED[[iCross]][["Mdl.ConfigEnv"]][["Mdl.crossValidIdx"]]
if(length(OUT.FITTED[[iCross]][["ID"]]) == 0)
{
    ## ID = 1:length(Portfolio.Real)
    ID4VARS <- Mdl.crossValidIdx[["testing"]][[iCross]]

} else
{
    ID <- as.Date(OUT.FITTED[[iCross]][["ID"]])
    ID4VARS <- ID[Mdl.crossValidIdx[["testing"]][[iCross]]]

}


## ylim = c(quantile(Portfolio.Real, 0.001),
##          quantile(Portfolio.Real, 0.999))
ylim = c(-10, 6)
plotIdx = 1:length(ID4VARS)

## plotIdx = (length(ID4VARS)-500+1):length(ID4VARS)

col = c("red",  "black", "purple", "blue")


pdf(file = "~/VaR001005.pdf", width = 10, height = 8)

nProbs = length(probs)
par(mfrow = c(nProbs, 1), mar = c(4, 4, 0.1, 0.1), cex = 0.8, ps = 12)
for(i in 1:nProbs)
{
    plot(ID4VARS[plotIdx],
         VARS[[1]][["VaR"]][i, plotIdx],
         type = "l", col = "red", ylim = ylim,
         ylab = "Portofolio and VaR",
         xlab = "Date")

    if(nModels>1)
    {
        for(iModel in 2:nModels)
        {
            points(ID4VARS[plotIdx], VARS[[iModel]][["VaR"]][i, plotIdx],
                   type = "l", col = col[iModel], lty = 1 )
        }
    }
    points(ID4VARS[plotIdx], Portfolio.Real[plotIdx], type = "p", pch = 20, col = "green")
    legend("topright",
           legend = c(paste("Portfolio with weights: (",
                            weights[1], ", " , weights[2], ")", sep = ""),
                      paste(probs[i]*100,  "% VaR for ", names(Models), " (PER: ",
                            round(Portfolio.LeftOut[i, ]*100), "%)", sep = "")
                      ),
           col = c("green", col[1:nModels]),
           ## lty = c(3, 1:nModels),
           lty = c(3, rep(1, nModels)),
           lwd = c(3, rep(1.5, nModels)),
           ncol = 1)
}
dev.off()
stopCluster(cl)
