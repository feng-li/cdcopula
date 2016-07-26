plotCplParTS <- function(MCMC.Update, MCMC.parSummary, MdlDGP.par = NULL, ObsIdx4Plot = NA)
{

    has.Display <- (nchar(Sys.getenv("DISPLAY"))>0)
    par.ts.mean <- MCMC.parSummary[["ts.mean"]]
    par.ts.hpd95 <- MCMC.parSummary[["ts.hpd95"]]
    if(has.Display)
    {
        nPlot <- sum(sapply(MCMC.Update, function(x) any(unlist(x) == TRUE)))
        nDev <- length(dev.list())
        if( nDev < nPlot)
        {
            replicate(nPlot-nDev, dev.new())

        }
        jDev <- 1

        if(any(is.na(ObsIdx4Plot)))
        {
            nObs <- nrow(par.ts.mean[[1]][[1]])
            ObsIdx4Plot <- 1:nObs
        }

    }


    for(i in names(MCMC.Update))
    {
        npar <- sum(unlist(MCMC.Update[[i]]))
        if(has.Display & npar>0)
        {
            dev.set(dev.list()[jDev])
            if(npar <= 2)
            {
                par(mfrow = c(npar, 1))
            }
            else
            {
                ## par(mfrow = c(ceiling(npar/2), 2))
                par(mfrow = c(npar, 1))

            }
            jDev <- jDev+1
        }
        for(j in names(MCMC.Update[[i]]))
        {
            if(MCMC.Update[[i]][[j]])
            {
                if(has.Display)
                {
                    if(ncol(par.ts.mean[[i]][[j]]) == 1)
                    {
                        ## TODO: Add smoothing plot for polygon.

                        par(mar = c(2, 4, 0, 0)+0.1)

                        hpd95 <- par.ts.hpd95[[i]][[j]][, ObsIdx4Plot, 1] #2-by-nobs
                        ylim <- c(quantile(hpd95[1, ], 0.001), quantile(hpd95[2, ], 0.999))

                        ## Initial plot to draw the plot window
                        if(length(colnames(hpd95)) == 0)
                        {
                            ## No time stamp attributed,  use 1, 2, ...
                            x <- 1:ncol(hpd95)
                        }
                        else
                        {
                            x <- as.Date(colnames(hpd95))
                        }


                        y <- par.ts.mean[[i]][[j]][ObsIdx4Plot, 1]


                        plot(x = x, y = y, type = "l", lty = "solid", col = "white",
                             ylim = ylim, ylab = j, xlab = "",
                             log = ifelse(j == "phi", "y", ""))

                        ## HPD Polygon
                        hpd95.smoothL <- spline(x, hpd95[1, ],
                                                n = length(ObsIdx4Plot)*10)
                        hpd95.smoothU <- spline(x, hpd95[2, ],
                                                n = length(ObsIdx4Plot)*10)

                        polygon(x = c(hpd95.smoothL$x, rev(hpd95.smoothU$x)),
                                y = c(hpd95.smoothL$y, rev(hpd95.smoothU$y)),
                                border = "grey", col = "grey")

                        ## polygon(x = c(hpd95.smoothL$x, rev(hpd95.smoothU$x)),
                        ##         y = c(hpd95.smoothL$y, rev(hpd95.smoothU$y)),
                        ##         border = "grey", col = "grey")


                        ## Posterior Mean
                        points(x = x, y = par.ts.mean[[i]][[j]][ObsIdx4Plot, 1],
                               type = "l", lty = "solid", col = "blue", lwd = 2)

                        ## DGP (Only for DGP data)
                        ## MdlDGP.par <- OUT.MCMC[["MdlDGP.par"]]
                        if(!(length(MdlDGP.par)  == 0 & is.null(MdlDGP.par)))
                        {
                            ## Mdl.Idx.training <- OUT.MCMC[["Mdl.Idx.training"]]
                            points(x = x, y = MdlDGP.par[[i]][[j]][ObsIdx4Plot],
                                   type = "l", lty = "dashed", col = "red", lwd = 2)
                            legend.idx <- 1:3
                        } else
                        {
                            legend.idx <- 1:2
                        }

                        ## Legend
                        legend("topleft",ncol = 1, # length(legend.idx),
                               lty = c("solid", "dotted", "dashed")[legend.idx],
                               lwd = c(3, 20, 3)[legend.idx],
                               col = c("blue", "grey", "red")[legend.idx],
                               legend = c("Posterior mean", "95% HPD", "DGP values")[legend.idx])
                    }
                    else
                    {
                        plot(0, main = i, ylab = j,
                             xlab = "(multivariate parameters such as covariance matrix are not plotted)")
                    }
                }
            }
        }
    }
}
