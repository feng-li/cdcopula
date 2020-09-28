#' @export
plotCplParTS <- function(MCMC.Update, MCMC.parSummary, MdlDGP.par = NULL,
                         ObsIdx4Plot = NA, ModelDescription, save.output, smoothHPD = FALSE)
{

    has.Display <- (nchar(Sys.getenv("DISPLAY"))>0)
    ## has.Display <- FALSE
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

        ## if(any(is.na(ObsIdx4Plot)))
        ## {
        ##     nObs <- nrow(par.ts.mean[[1]][[1]])
        ##     ObsIdx4Plot <- 1:nObs
        ## }

    }
    else
    {
        ## Batch mode,  save figures into a pdf file
        JOB_ID <- Sys.getenv("SLURM_JOB_ID")
        if(nchar(JOB_ID)>0)
        {
            JOB.str <- paste("JOB", JOB_ID, ".", sep = "")
        }
        else
        {
            JOB.str <- ""
        }

        outfile <- file.path(save.output,
                             paste(JOB.str, ModelDescription, ".pdf" ,   sep  =  ""))

        pdf(file = outfile, width = 8, height = 4)

    }

    for(i in names(MCMC.Update))
    {
        npar <- sum(unlist(MCMC.Update[[i]]))
        if(has.Display & npar>0)
        {
            ## Interactive mode
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
                ## if(has.Display)
                ## {
                if(ncol(par.ts.mean[[i]][[j]]) == 1)
                {

                    if(any(is.na(ObsIdx4Plot)))
                    {
                        nObs <- nrow(par.ts.mean[[1]][[1]])
                        ObsIdx4Plot <- 1:ncol(par.ts.hpd95[[i]][[j]])
                    }

                    par(mar = c(2, 4, 0.5, 0.5), las = 1)
                    hpd95 <- par.ts.hpd95[[i]][[j]][, ObsIdx4Plot, 1] #2-by-nobs

                    hardlim = FALSE
                    if(hardlim == TRUE)
                        {
                            if(j  == "mu") {
                                ylim = c(-0.5, 1)
                                ylab = expression(mu)
                            }
                            else if(j  == "phi") {
                                ylim = c(0.3, 15)
                                ylab = expression(phi)
                            }
                            else if(j  == "df") {
                                ylim = c(2, 30)
                                ylab = expression(nu)
                            }
                            else if(j  == "lmd") {
                                ylim = c(0.5, 1.5)
                                ylab = expression(kappa)
                            }
                            else  if(j  == "lambdaL") {
                                ylim = c(0.3, 0.9)
                                ylab = expression(lambda[L])
                            }
                            else if(j  == "lambdaU") {
                                ylim = c(0.3, 0.9)
                                ylab = expression(lambda[U])
                            }
                            else if(j  == "tau") {
                                ylim = c(0.3, 0.9)
                                ylab = expression(tau)
                            }
                            else
                            {
                            }
                        }
                     else
                    {
                        dgp <- MdlDGP.par[[i]][[j]][ObsIdx4Plot]
                        ylim <- c(quantile(c(hpd95[1, ], dgp), 0.001),
                                  quantile(c(hpd95[2, ], dgp), 0.999))
                        ylab = j
                    }
                    ## Initial plot to draw the plot window
                    if(length(colnames(hpd95)) == 0)
                    {
                        ## No time stamp attributed,  use 1, 2, ...
                        x <- ObsIdx4Plot
                    }
                    else
                    {
                        x <- as.Date(colnames(hpd95))
                    }


                    y <- par.ts.mean[[i]][[j]][ObsIdx4Plot, 1]


                    plot(x = x, y = y, type = "l", lty = "solid", col = "white",
                         ylim = ylim, ylab = ylab, xlab = "",
                         log = ifelse(j %in% c("phi", "df"), "y", ""))

                    ## HPD Polygon
                    if(smoothHPD == "TRUE")
                    {
                        nSpline = length(ObsIdx4Plot)*10
                        polyIdx = round(seq(1, nSpline, by = 10))

                        hpd95.smoothL <- spline(x, hpd95[1, ], nSpline)
                        hpd95.smoothU <- spline(x, hpd95[2, ], nSpline)
                        if(j == "df") browser()
                        polygon(x = c(hpd95.smoothL$x[polyIdx], rev(hpd95.smoothU$x[polyIdx])),
                                y = c(hpd95.smoothL$y[polyIdx], rev(hpd95.smoothU$y[polyIdx])),
                                border = "grey", col = "grey")
                    }
                    else
                    {
                        polygon(x = c(x, rev(x)),
                                y = c(hpd95[1, ], rev(hpd95[2, ])),
                                border = "grey", col = "grey")
                    }

                    ## Posterior Mean
                    points(x = x, y = par.ts.mean[[i]][[j]][ObsIdx4Plot, 1],
                           type = "l", lty = "solid", col = "blue", lwd = 2)

                    ## Hard code
                    ## browser()
                    ## abline(v = 0.6556)
                    ## abline(v = 0.682)

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
                    legend("topleft",ncol = 3, # length(legend.idx),
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
                ## }
            }
        }

    }
    if(!has.Display)
    {
        ## Save plots to PDF file
        dev.off()
    }
    else
    {
        ## plot is not saved automatically. Use dev.copy2pdf()
    }

}
