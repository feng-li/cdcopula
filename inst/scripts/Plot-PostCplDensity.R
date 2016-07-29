## This script plots the daily return of SP100 and SP600 data
load(file.path("~/running/", "JOB835.^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross0MCMC.nIter10000joint+20160720@14.44.5c1da7.Rdata"))
## load(file.path("~/running/", "^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross1MCMC.nIter1000twostage+20160522@23.11.c168ab.Rdata"))

sourceDir("~/code/cdcopula/R", recursive = TRUE)

iCross <- 1


MCMC.Update = OUT.FITTED[[iCross]][["MCMC.Update"]]
MCMC.burninProp <- OUT.FITTED[[iCross]][["MCMC.burninProp"]]
MCMC.sampleProp <- OUT.FITTED[[iCross]][["MCMC.sampleProp"]]
MCMC.nIter <- OUT.FITTED[[iCross]][["MCMC.nIter"]]
Mdl.Y <- OUT.FITTED[[iCross]][["Mdl.Y"]]
CplNM <- names(MCMC.Update)[length(MCMC.Update)]



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

##########################


## Plot The Returns
nObs <- length(lambdaL)
nDataWindow <- 6
nPlot <- 5
nPoints <- 100


par(mfrow = c(nDataWindow, nPlot), mar = c(2, 2, 0, 0)+0.1)
ylim <- c(-4, 4)

dataWindowIdx <- data.partition(nObs = nObs,
                                args = list(partiMethod = "ordered",
                                            N.subsets = nDataWindow))
for(i in 1:nDataWindow)
{

    idx <- dataWindowIdx[[i]]
    idx.order <- order(lambdaL[idx, ])

    ID.cand <- idx[idx.order[seq(1, length(idx), length.out = nPlot)]]

    ## normaldays <- c(1, 3, 5, 7, 9)
    ## ID.card <- 3000:3500

    ## print(lambdaL[ID.cand])

    whichCol <- 0
    for(j in ID.cand)
    {

        whichCol <- whichCol+1

        if(whichCol == nPlot & i == 5)
        {
            ylim = c(-30, 20)
            col = "red"
        }
        else if(whichCol == nPlot)
        {
            ylim = c(-4, 4)
            col = "orange"

        }
        else
        {
            ylim = c(-2, 2)
            if(i == 5)
            {
                col = "red"
            }
            else
            {
                col = "blue"
            }
        }



        y1 <- matrix(seq(ylim[1], ylim[2], length.out = nPoints))
        y2 <- y1
        y <- mesh.grid(y1)


        Mdl.Y.grid <- list(y1 = y[, 1, drop = FALSE],
                           y2 = y[, 2, drop = FALSE])
        names(Mdl.Y.grid) <- names(MCMC.Update[1:2])




        Mdl.par.curr <- rapply(summary.Cpl[["par.summary"]][["ts.mean"]],
                               function(x) x[j, drop = FALSE], how = "replace")
        Mdl.ud <- logDens(Mdl.MargisType = OUT.FITTED[[iCross]][["Mdl.MargisType"]],
                          Mdl.Y = Mdl.Y.grid,
                          Mdl.par = Mdl.par.curr,
                          parUpdate = OUT.FITTED[[iCross]][["MCMC.Update"]],
                          MCMC.UpdateStrategy = OUT.FITTED[[iCross]][["MCMC.UpdateStrategy"]])

        dens.grid <- matrix(exp(rowSums(Mdl.ud[["Mdl.d"]])), nPoints)


        contour(as.matrix(y1),  as.matrix(y2),  dens.grid,  add  =  FALSE,
                col  = col,  lwd  =  0.8,  lty  = "solid",
                xlim  = ylim,  ylim  =  ylim,
                xlab = "", ylab = "", drawlabels = TRUE)
        legend("topleft", ncol = 1, bty = "n",
               legend = c(paste("(", round(Mdl.par.curr[[CplNM]][["lambdaL"]], 2),", ",
                                round(Mdl.par.curr[[CplNM]][["lambdaU"]], 2), ")", sep = "")))

        legend("bottomright", ncol = 1, bty = "n", legend = paste(ID[j], "  ", sep = ""))
    }
}
