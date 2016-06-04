##' Trajectory MCMC.
##'
##' This function can also be used for summarizing the posterior results.
##' @param MCMC.nIter "integer"
##' @param iIter "integer"
##' @param interval e.g. 10%
##' @param ... Other arguments used
##' @return A summary object
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Fri Feb 01 14:49:15 CET 2013; Current: Mon Mar 30 16:32:00 CST 2015.
CplMCMC.summary <- function(MCMC.nIter, iIter, interval = 0.1, MCMC.burninProp, OUT.MCMC,
                            maxcovprint = 20, ObsIdx4Plot = NA)
{
    dev.width <- getOption("width")
    has.Display <- (nchar(Sys.getenv("DISPLAY"))>0)


    ## Set default values for missing arguments
    if(missing(MCMC.nIter))
    {
        MCMC.burninProp <- OUT.MCMC[["MCMC.burninProp"]]

        MCMC.nIter <- OUT.MCMC[["MCMC.nIter"]]

        if(length(MCMC.burninProp) == 0 | length(MCMC.nIter) == 0)
        {
            stop("MCMC.burninProp & MCMC.nIter are not int OUT.MCMC. Supply them manually.")
        }
    }

    if(missing(iIter))
    {
        iIter <- MCMC.nIter
    }


    ## Progress statistics
    Starting.time <- OUT.MCMC[["Starting.time"]]

    TimeUsed <- difftime(Sys.time(), Starting.time, units = "hours")
    TimeToGo <-  round(TimeUsed/(iIter-1)*(MCMC.nIter-1-iIter), 2) # take away first
                                        # waring-up iteration.

    TimeUsed.units <- ifelse(abs(TimeUsed) < 1, "mins", "hours")
    TimeToGo.units <- ifelse(abs(TimeToGo) < 1, "mins", "hours")


    if(TimeToGo<0)
    {
        TimeToGo <- 0
        TimeUsed <- NA
    }


    donePercent <- round(iIter/MCMC.nIter*100)

    progress <- paste(date(), ": ", round(as.numeric(TimeUsed, units = TimeUsed.units), 1),
                      " ", TimeUsed.units, " elapsed, ", donePercent, "% done, ",
                      round(as.numeric(TimeToGo, units = TimeToGo.units), 1),
                      " ", TimeToGo.units, " to go...\n", sep = "")


    printInterval2 <- ifelse(MCMC.nIter*interval < 1,  1,
                             floor(MCMC.nIter*0.1))

    printIter1 <- 4
    printIter2 <- seq(from = printInterval2, to = MCMC.nIter, by = printInterval2)
    if(printIter2[1] >= printIter1)
    {
        printIter2 <- c(printIter1, printIter2)
    }

    if(iIter  %in% printIter2)
    {
        cat(progress)
    }

    ## The burning
    n.burn.default <- round(MCMC.nIter*MCMC.burninProp)
    n.burn <- ifelse(iIter>n.burn.default, n.burn.default, 0)

    ## Set the print interval and consider burnin. If print interval is to narrow, set to
    ## one.
    printInterval <- ifelse(floor(MCMC.nIter*interval) == 0,  1,
                            floor(MCMC.nIter*interval))
    printIter <- seq(from = printInterval, to = MCMC.nIter, by = printInterval)
    if(printIter[length(printIter)] != MCMC.nIter)
    {# Always print summary at last iteration.
        printIter <- c(printIter, MCMC.nIter)
    }

    if(iIter %in% printIter)
    {
        ## par <- list(...)
        MCMC.betaIdx <- OUT.MCMC[["MCMC.betaIdx"]]
        MCMC.beta <- OUT.MCMC[["MCMC.beta"]]
        MCMC.par <- OUT.MCMC[["MCMC.par"]]
        MCMC.AccProb <- OUT.MCMC[["MCMC.AccProb"]]
        MCMC.Update <- OUT.MCMC[["MCMC.Update"]]
        MCMC.sampleProp <- OUT.MCMC[["MCMC.sampleProp"]]

        welcome <- paste("MCMC SUMMARY: ", donePercent, "% (",
                         round(n.burn/MCMC.nIter*100), "% MCMC.burninProp) "
                       , TimeToGo, " hours to go.\n", sep = "")
        cat("\n", rep("=", dev.width), "\n", sep = "")
        cat(welcome)
        cat(rep("=", dev.width), "\n", sep = "")

        MCMC.sampleIdx <- round(seq(n.burn+1, iIter,
                                    length.out = round((iIter-n.burn)*MCMC.sampleProp)))
        nMCMCSample <- length(MCMC.sampleIdx)

        ## Regenerate ``MCMC.par`` if not supplied
        MCMC.par <- NULL
        if(length(MCMC.par) == 0)
        {
            MCMC.par <- list()
            Mdl.X.training <- OUT.MCMC[["Mdl.X.training"]]
            Mdl.parLink <- OUT.MCMC[["Mdl.parLink"]]


            nTraining <- nrow(Mdl.X.training[[1]][[1]])

            Mdl.MargisNM <- names(MCMC.Update)
            for(CompCaller in names(MCMC.Update))
            {
                for(parCaller in names(MCMC.Update[[CompCaller]]))
                {
                    ncolX.ij <- ncol(Mdl.X.training[[CompCaller]][[parCaller]])
                    nPar.ij <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]
                    namesX.ij <- rep(colnames(Mdl.X.training[[CompCaller]][[parCaller]]), nPar.ij)


                    if((CompCaller %in% Mdl.MargisNM[length(Mdl.MargisNM)]) & nPar.ij != 1)
                    {
                        nDim <- length(Mdl.MargisNM)-1
                        namesParFull.ij <- matrix(paste(matrix(1:nDim, nDim, nDim),
                                                        matrix(1:nDim, nDim, nDim,
                                                               byrow = TRUE), sep = "."), nDim)
                        namesPar.ij <- namesParFull.ij[lower.tri(namesParFull.ij,
                                                                 diag = FALSE)]
                    }
                    else
                    {
                        namesPar.ij <- "1.1"
                    }


                    MCMC.par[[CompCaller]][[parCaller]] <- array(NA, c(nMCMCSample, nTraining, nPar.ij),
                                                                 dimnames = list(NULL, NULL, namesPar.ij))
                }
            }

            subsetFun4beta <- function(x, idx)
            {
                if((dim(x)[1] == 1 && length(idx)>1) ||
                   (dim(x)[1] == 1 && length(idx) == 1 && idx != 1))
                {# check whether some parameters are not updated
                    out <- x
                }
                else
                {
                    out <- x[idx, , drop = FALSE]
                }
                return(out)
            }

            for(iMCMC.sampleIdx in 1:nMCMCSample)
            {
                Mdl.beta.curr <- rapply(object=MCMC.beta, f = subsetFun4beta,
                                        idx = MCMC.sampleIdx[iMCMC.sampleIdx],
                                        how = "replace")

                Mdl.par.curr <- parCplMeanFun(Mdl.X = Mdl.X.training,
                                              Mdl.parLink = Mdl.parLink,
                                              Mdl.beta = Mdl.beta.curr,
                                              parUpdate = MCMC.Update)

                for(CompCaller in names(MCMC.Update))
                {
                    for(parCaller in names(MCMC.Update[[CompCaller]]))
                    {
                        MCMC.par[[CompCaller]][[parCaller]][iMCMC.sampleIdx, ,] <- Mdl.par.curr[[CompCaller]][[parCaller]]
                    }
                }
            }
        }


        subFun3 <- function(obj, fun, dim, ...)
        {
            if(any(is.na(obj)))
            {
                out <- NA
            }
            else
            {
                out <- apply(obj, dim, fun, ...)
            }
            return(out)
        }
        par.mean <- rapply(MCMC.par, subFun3, how = "replace", fun = mean, dim = 3)
        par.median <- rapply(MCMC.par, subFun3, how = "replace", fun = median, dim = 3,)
        par.sd <- rapply(MCMC.par, subFun3, how = "replace",   fun = sd, dim = 3)
        par.ts.mean <- rapply(MCMC.par, subFun3, how = "replace", fun = mean, dim = c(2, 3))
        par.ts.median <- rapply(MCMC.par, subFun3, how = "replace", fun = median,
                                dim = c(2, 3))
        par.ts.sd <- rapply(MCMC.par, subFun3, how = "replace", fun = sd, dim = c(2, 3))
        par.ts.hpd95 <- rapply(MCMC.par, subFun3, how = "replace", fun = quantile,
                               dim = c(2, 3), probs = c(0.025, 0.975))


        subFun <- function(obj, fun)
        {
            if(any(is.na(obj)))
            {
                out <- NA
            }
            else
            {
                out <- apply(obj, 2, fun)
            }
            return(out)
        }
        accept.prob.mean <- rapply(MCMC.AccProb, subFun, how = "replace", fun = mean)
        betaIdx.mean <- rapply(MCMC.betaIdx, subFun, how = "replace", fun = mean)
        betaIdx.median <- rapply(MCMC.betaIdx, subFun, how = "replace", fun = median)
        beta.mean <- rapply(MCMC.beta, subFun, how = "replace", fun = mean)
        beta.median <- rapply(MCMC.beta, subFun, how = "replace", fun = median)
        beta.sd <- rapply(MCMC.beta, subFun, how = "replace",   fun = sd)
        beta.ineff <- rapply(MCMC.beta, subFun, how = "replace", fun = ineff)


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
                ObsIdx4Plot <- 1:nTraining
            }

        }

        for(i in names(MCMC.beta))
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
                    par(mfrow = c(ceiling(npar/2), 2))
                }
                jDev <- jDev+1
            }
            for(j in names(MCMC.beta[[i]]))
            {
                if(MCMC.Update[[i]][[j]])
                {
                    if(has.Display)
                    {
                        if(ncol(par.ts.mean[[i]][[j]]) == 1)
                        {
                            ## TODO: Add smoothing plot for polygon.

                            par(mar = c(3, 4, 0, 0)+0.1)

                            hpd95 <- par.ts.hpd95[[i]][[j]][, ObsIdx4Plot, 1]
                            ylim <- c(min(hpd95[1, ]), max(hpd95[2, ]))
                            ## browser()
                            ## Initial plot to draw the plot window
                            plot(par.ts.mean[[i]][[j]][ObsIdx4Plot, 1], type = "l",
                                 lty = "solid", col = "blue",
                                 ylim = ylim, ylab = j, xlab = "")

                            ## HPD Polygon
                            hpd95.smoothL <- spline(1:length(ObsIdx4Plot), hpd95[1, ],
                                                    n = length(ObsIdx4Plot)*10)
                            hpd95.smoothU <- spline(1:length(ObsIdx4Plot), hpd95[2, ],
                                                    n = length(ObsIdx4Plot)*10)

                            polygon(x = c(hpd95.smoothL$x, rev(hpd95.smoothU$x)),
                                    y = c(hpd95.smoothL$y, rev(hpd95.smoothU$y)),
                                    border = "grey", col = "grey")

                            ## Posterior Mean
                            points(par.ts.mean[[i]][[j]][ObsIdx4Plot, 1],
                                   type = "l", lty = "solid", col = "blue", lwd = 2)

                            ## DGP (Only for DGP data)
                            MdlDGP.par <- OUT.MCMC[["MdlDGP.par"]]
                            if(!(length(MdlDGP.par)  == 0 & is.null(MdlDGP.par)))
                            {
                                ## Mdl.Idx.training <- OUT.MCMC[["Mdl.Idx.training"]]
                                points(MdlDGP.par[[i]][[j]][ObsIdx4Plot],
                                       type = "l", lty = "dashed", col = "red", lwd = 2)
                                legend.idx <- 1:3
                            } else
                            {
                                legend.idx <- 1:2
                            }

                            ## Legend
                            legend("topright",ncol = length(legend.idx),
                                   lty = c("dotted", "solid", "dashed")[legend.idx],
                                   lwd = c(30, 1, 1)[legend.idx],
                                   col = c("grey", "blue", "red")[legend.idx],
                                   legend = c("95% HPD", "Posterior mean", "DGP values")[legend.idx])
                        }
                        else
                        {
                            plot(0, main = i, ylab = j,
                                 xlab = "(multivariate parameters such as covariance matrix are not plotted)")
                        }
                    }

                    obj.par <- rbind(par.mean[[i]][[j]],
                                     par.median[[i]][[j]],
                                     par.sd[[i]][[j]],
                                     round(accept.prob.mean[[i]][[j]], 2))
                    rownames(obj.par) <- c("par.mean", "par.median", "par.sd", "acc.prob")
                    colnames(obj.par) <- names(par.mean[[i]][[j]])

                    obj <- rbind(beta.mean[[i]][[j]],
                                 beta.median[[i]][[j]],
                                 beta.sd[[i]][[j]],
                                 betaIdx.mean[[i]][[j]],
                                 beta.ineff[[i]][[j]])

                    rownames(obj) <- c("beta.mean", "beta.median", "beta.sd",
                                       "betaIdx.mean", "beta.ineff")
                    colnames(obj) <- paste(rep(colnames(obj.par),
                                               each = ncol(obj)/ncol(obj.par)),
                                           colnames(obj), sep = "|")
                    ## colnames(obj) <- paste(j, rep(1:ncol(obj.par),
                    ##                               each = ncol(obj)/ncol(obj.par)), colnames(obj),
                    ##                        sep = ".")

                    cat("\n", rep("-", dev.width-1), "\n", sep = "")
                    cat(i, j, "(", donePercent, "% )\n")

                    covprint <- min(maxcovprint, ncol(obj))
                    covrintCol <- order(obj["betaIdx.mean", ], decreasing = TRUE)[0:covprint]
                    print(obj.par)
                    cat("\n")
                    print(obj[, covrintCol, drop = FALSE])
                    if(maxcovprint< ncol(obj))
                    {
                        cat("(Too many covariates used. The first 20 largest `betaIdx.mean` printed.)\n" )
                    }
                }
            }
        }

        if(iIter == MCMC.nIter)
        {
            cat(rep("-", dev.width), "\n\n",  sep = "")
        }
        out <- list(par.mean = par.mean,
                    par.median = par.median,
                    par.sd = par.sd,
                    par.ts.mean = par.ts.mean,
                    par.ts.median = par.ts.median,
                    par.ts.sd = par.ts.sd,
                    par.ts.hpd95 = par.ts.hpd95,
                    beta.mean = beta.mean,
                    beta.median = beta.median,
                    beta.sd = beta.sd,
                    betaIdx.mean = betaIdx.mean,
                    beta.ineff = beta.ineff)

        invisible(out)
    }

}
