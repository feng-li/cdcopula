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
##' TODO: write this function as a summary
CplMCMC.summary <- function(MCMC.nIter, iIter = MCMC.nIter, interval = 0.1, MCMC.burninProp, OUT.MCMC, maxcovprint = 20)
{
    ## Set the print interval and consider burnin. If print interval is to narrow, set to
    ## one.
    printInterval <- ifelse(floor(MCMC.nIter*interval) == 0,  1, floor(MCMC.nIter*interval))
    printIter <- seq(from = printInterval,
                     to = MCMC.nIter,
                     by = printInterval)
    if(printIter[length(printIter)] != MCMC.nIter)
    {# Always print summary at last iteration.
        printIter <- c(printIter, MCMC.nIter)
    }

    dev.width <- getOption("width")
    has.Display <- (nchar(Sys.getenv("DISPLAY"))>0)

    ## The burning
    n.burn.default <- round(MCMC.nIter*MCMC.burninProp)
    n.burn <- ifelse(iIter>n.burn.default, n.burn.default, 0)

    subFun <- function(x, iIter, fun){
        candIdx <- (n.burn+1):iIter
        if((dim(x)[1] == 1 && length(candIdx)>1) ||
           (dim(x)[1] == 1 && length(candIdx) == 1 && candIdx != 1))
        {
            out <- NA
        }
        else
        {
            obj <- x[candIdx, ,drop = FALSE]
            if(any(is.na(obj)))
            {
                out <- NA
            }
            else
            {
                out <- apply(obj, 2, fun)
            }
        }

        return(out)
    }

    subFun3 <- function(x, iIter, fun, dim, ...)
    {
        candIdx <- (n.burn+1):iIter
        if((dim(x)[1] == 1 && length(candIdx)>1) ||
           (dim(x)[1] == 1 && length(candIdx) == 1 && candIdx != 1))
        {
            out <- NA
        }
        else
        {
            obj <- x[candIdx, , ,drop = FALSE]
            if(any(is.na(obj)))
            {
                out <- NA
            }
            else
            {
                out <- apply(obj, dim, fun, ...)
            }
        }
        return(out)
    }



    if(iIter %in% printIter)
    {
        ## par <- list(...)
        MCMC.betaIdx <- OUT.MCMC$MCMC.betaIdx
        MCMC.beta <- OUT.MCMC$MCMC.beta
        MCMC.par <- OUT.MCMC$MCMC.par
        MCMC.AccProb <- OUT.MCMC$MCMC.AccProb
        MCMCUpdate <- OUT.MCMC$MCMCUpdate
        Starting.time <- OUT.MCMC$Starting.time

        TimeToGo <-  round(difftime(Sys.time(), Starting.time,
                                    units = "hours")/(iIter-1)*(MCMC.nIter-iIter), 2)

        if(iIter  == printIter[1] && iIter != MCMC.nIter)
        {
            cat("\nabout", TimeToGo, " hours to go.\n", sep = " ")
            return()
        }

        donePercent <- round(iIter/MCMC.nIter*100)
        welcome <- paste("MCMC SUMMARY: ", donePercent, "% (",
                         round(n.burn/MCMC.nIter*100), "% MCMC.burninProp) "
                       , TimeToGo, " hours to go.\n", sep = "")
        cat("\n", rep("=", dev.width), "\n", sep = "")
        cat(welcome)
        cat(rep("=", dev.width), "\n", sep = "")

        accept.prob.mean <- rapply(MCMC.AccProb, subFun, how = "replace",
                                   iIter = iIter,  fun = mean)

        par.mean <- rapply(MCMC.par, subFun3, how = "replace",
                           iIter = iIter,  fun = mean, dim = 3)
        par.median <- rapply(MCMC.par, subFun3, how = "replace",
                             iIter = iIter,  fun = median, dim = 3)
        par.sd <- rapply(MCMC.par, subFun3, how = "replace",
                         iIter = iIter,  fun = sd, dim = 3)

        betaIdx.mean <- rapply(MCMC.betaIdx, subFun, how = "replace",
                               iIter = iIter, fun = mean)
        betaIdx.median <- rapply(MCMC.betaIdx, subFun, how = "replace",
                                 iIter = iIter, fun = median)
        beta.mean <- rapply(MCMC.beta, subFun, how = "replace",
                            iIter = iIter,  fun = mean)
        beta.median <- rapply(MCMC.beta, subFun, how = "replace",
                              iIter = iIter,  fun = median)
        beta.sd <- rapply(MCMC.beta, subFun, how = "replace",
                          iIter = iIter,  fun = sd)

        par.ts.mean <- rapply(MCMC.par, subFun3, how = "replace", iIter = iIter,
                              fun = mean, dim = c(2, 3))
        par.ts.median <- rapply(MCMC.par, subFun3, how = "replace", iIter = iIter,
                                fun = median, dim = c(2, 3))
        par.ts.sd <- rapply(MCMC.par, subFun3, how = "replace", iIter = iIter,
                            fun = sd, dim = c(2, 3))
        par.ts.hpd95 <- rapply(MCMC.par, subFun3, how = "replace", iIter = iIter,
                               fun = quantile, dim = c(2, 3), probs = c(0.025, 0.975))

        ## Efficiency factor of MCMC
        beta.ineff <- rapply(MCMC.beta, subFun, how = "replace", iIter = iIter,
                             fun = ineff)


        if(has.Display)
        {
            nPlot <- sum(sapply(MCMCUpdate, function(x) any(unlist(x) == TRUE)))
            nDev <- length(dev.list())
            if( nDev < nPlot)
            {
                replicate(nPlot-nDev, dev.new())

            }
            jDev <- 1
        }

        for(i in names(MCMC.beta))
        {
            npar <- sum(unlist(MCMCUpdate[[i]]))
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
                if(MCMCUpdate[[i]][[j]])
                {

                    if(has.Display)
                    {
                        if(ncol(par.ts.mean[[i]][[j]]) == 1)
                        {
                            hpd95 <- par.ts.hpd95[[i]][[j]][, , 1]
                            ylim <- c(min(hpd95[1, ]), max(hpd95[2, ]))

                            ## Initial plot to draw the plot window
                            plot(par.ts.mean[[i]][[j]][, 1], type = "l", lty = "solid",
                                 col = "blue", ylim = ylim, ylab = j)


                            ## HPD Polygon
                            polygon(x = c(1:ncol(hpd95), ncol(hpd95):1),
                                    y = c(hpd95[1, ], rev(hpd95[2, ])),
                                    border = "grey", col = "grey")

                            points(par.ts.mean[[i]][[j]][, 1],
                                   type = "l", lty = "solid", col = "blue", lwd = 2)

                            MdlDGP.par <- OUT.MCMC[["MdlDGP.par"]]
                            if(!(length(MdlDGP.par) = 1 & is.null(MdlDGP.par)))
                            {
                                Mdl.Idx.training <- OUT.MCMC[["Mdl.Idx.training"]]
                                points(MdlDGP.par[[i]][[j]][Mdl.Idx.training],
                                       type = "l", lty = "solid", col = "red", lwd = 2)
                                legend.idx <- 1:3
                            }
                            else
                            {
                                legend.idx <- 1:2
                            }

                            legend("topright",ncol = length(legend.idx),
                                   lty = c("dotted", "solid", "dashed")[legend.idx],
                                   lwd = c(10, 1, 1)[legend.idx],
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
