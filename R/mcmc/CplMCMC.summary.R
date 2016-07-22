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
CplMCMC.summary <- function(MCMC.nIter, iIter, interval = 0.1, MCMC.burninProp,
                            MCMC.sampleProp, OUT.MCMC,
                            maxcovprint = 20, ObsIdx4Plot = NA)
{
    dev.width <- getOption("width")

    ## Set default values for missing arguments
    if(missing(MCMC.nIter))
    {
        MCMC.burninProp <- OUT.MCMC[["MCMC.burninProp"]]

        MCMC.nIter <- OUT.MCMC[["MCMC.nIter"]]

        if(length(MCMC.burninProp) == 0 | length(MCMC.nIter) == 0)
        {
            stop("MCMC.burninProp & MCMC.nIter are not available in OUT.MCMC. Supply them manually.")
        }
    }

    if(missing(iIter))
    {
        iIter <- MCMC.nIter
    }

    if(missing(MCMC.sampleProp))
    {
        MCMC.sampleProp <- OUT.MCMC[["MCMC.sampleProp"]]

        if(length(MCMC.sampleProp) == 0)
        {
            stop("MCMC.sampleProp is not available in OUT.MCMC. Supply it manually.")
        }

    }


    ## Progress statistics
    Starting.time <- OUT.MCMC[["Starting.time"]]

    TimeUsed <- difftime(Sys.time(), Starting.time, units = "hours")
    TimeToGo <- TimeUsed/iIter*(MCMC.nIter-iIter)
    TimeUsed.units <- ifelse(abs(TimeUsed) < 1, "mins", "hours")
    TimeToGo.units <- ifelse(abs(TimeToGo) < 1, "mins", "hours")


    if(TimeToGo<0)
    {
        TimeToGo <- 0
        TimeUsed <- NA
    }

    ## if(TimeUsed <- 0.01) browser()

    donePercent <- round(iIter/MCMC.nIter*100)

    progress <- paste(format(Sys.time(), "%Y-%b-%d@%H:%M"), ": ",
                      round(as.numeric(TimeUsed, units = TimeUsed.units), 1),
                      " ", TimeUsed.units, " elapsed, ", donePercent, "% done, ",
                      round(as.numeric(TimeToGo, units = TimeToGo.units), 1),
                      " ", TimeToGo.units, " to go...\n", sep = "")


    printInterval2 <- ifelse(MCMC.nIter*interval < 1,  1, floor(MCMC.nIter*0.1))

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
        ## MCMC.par <- OUT.MCMC[["MCMC.par"]]
        MCMC.AccProb <- OUT.MCMC[["MCMC.AccProb"]]
        MCMC.Update <- OUT.MCMC[["MCMC.Update"]]
        ## MCMC.sampleProp <- OUT.MCMC[["MCMC.sampleProp"]]

        welcome <- paste("MCMC SUMMARY: ", donePercent, "% (",
                         round(n.burn/MCMC.nIter*100), "% MCMC.burninProp) "
                       , TimeToGo, " hours to go.\n", sep = "")
        cat("\n", rep("=", dev.width), "\n", sep = "")
        cat(welcome)
        cat(rep("=", dev.width), "\n", sep = "")

        MCMC.sampleIdx <- round(seq(n.burn+1, iIter,
                                    length.out = round((iIter-n.burn)*MCMC.sampleProp)))

        ## Generate ``MCMC.par summary
        MCMC.parSummary <- parCplMCMC(MCMC.beta = MCMC.beta,
                                      Mdl.X = OUT.MCMC[["Mdl.X.training"]],
                                      Mdl.parLink = OUT.MCMC[["Mdl.parLink"]],
                                      MCMC.Update = MCMC.Update,
                                      MCMC.sampleIdx = MCMC.sampleIdx)
        ## Plot Post Summary
        plotCplParTS(MCMC.Update = MCMC.Update,
                     MCMC.parSummary = MCMC.parSummary,
                     ObsIdx4Plot = ObsIdx4Plot,
                     MdlDGP.par = OUT.MCMC[["MdlDGP.par"]])


        ## Print Post Summary
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

        for(i in names(MCMC.Update))
        {
            for(j in names(MCMC.Update[[i]]))
            {
                if(MCMC.Update[[i]][[j]])
                {
                    obj.par <- rbind(MCMC.parSummary[["mean"]][[i]][[j]],
                                     MCMC.parSummary[["median"]][[i]][[j]],
                                     MCMC.parSummary[["sd"]][[i]][[j]],
                                     round(accept.prob.mean[[i]][[j]], 2))
                    rownames(obj.par) <- c("par.mean", "par.median", "par.sd", "acc.prob")
                    ##  colnames(obj.par) <- names(par.mean[[i]][[j]])

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

                    if(maxcovprint < ncol(obj))
                    {
                        cat("(Too many covariates in the model. The first ", maxcovprint,
                            " largest `betaIdx.mean` printed.)\n" )
                    }

                }
            }
        }


        if(iIter == MCMC.nIter)
        {
            cat(rep("-", dev.width), "\n\n",  sep = "")
        }


        out <- list(par.summary = MCMC.parSummary,
                    beta.mean = beta.mean,
                    ## beta.median = beta.median,
                    beta.sd = beta.sd,
                    betaIdx.mean = betaIdx.mean,
                    beta.ineff = beta.ineff,
                    accept.prob.mean = accept.prob.mean)

        invisible(out)
    }

}
