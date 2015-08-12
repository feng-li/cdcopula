##' Trajectory MCMC.
##'
##' This function can also be used for summarizing the posterior results.
##' @param nIter "integer"
##' @param iIter "integer"
##' @param interval e.g. 10%
##' @param ... Other arguments used
##' @return A summary object
##' @references NA
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Fri Feb 01 14:49:15 CET 2013; Current: Mon Mar 30 16:32:00 CST 2015.
##' TODO: write this function as a summary
CplMCMC.summary <- function(nIter, iIter = nIter, interval = 0.1, burnin, OUT.MCMC)
{
  ## Set the print interval and consider burnin
  printIter <- c(5, seq(from = floor(nIter*interval),
                        to = nIter,
                        by = floor(nIter*interval)))
  if(printIter[length(printIter)] != nIter)
  {
    printIter <- c(printIter, nIter)
  }

  dev.width <- getOption("width")
  has.Display <- (nchar(Sys.getenv("DISPLAY"))>0)

  ## The burning
  n.burn.default <- round(nIter*burnin)
  n.burn <- ifelse(iIter>n.burn.default, n.burn.default, 0)


  subFun <- function(x, iIter, fun){
    obj <- x[(n.burn+1):iIter, ,drop = FALSE]
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

  subFun3 <- function(x, iIter, fun, dim, ...){
    obj <- x[(n.burn+1):iIter, , ,drop = FALSE]
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
                                units = "hours")/(iIter-1)*(nIter-iIter), 2)

    if(iIter  == printIter[1])
    {
      cat(TimeToGo, " hours to go.\n", sep = "")
      return()
    }

    donePercent <- round(iIter/nIter*100)
    welcome <- paste("MCMC SUMMARY: ", donePercent, "% (",
                     round(n.burn/nIter*100), "% burnin) "
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
        par(mfrow = c(npar, 1))
        jDev <- jDev+1
      }
      for(j in names(MCMC.beta[[i]]))
      {
        if(MCMCUpdate[[i]][[j]])
        {

          if(has.Display && ncol(par.ts.mean[[i]][[j]]) == 1)
          {
            hpd95 <- par.ts.hpd95[[i]][[j]][, , 1]
            ylim <- c(min(hpd95[1, ]), max(hpd95[2, ]))
            plot(hpd95[1, ], type = "l", lty = "dotted", col = "red",
                 ylim = ylim , ylab = j, main = i)
            points(hpd95[2, ], type = "l", lty = "dotted", col = "red")
            points(par.ts.mean[[i]][[j]][, 1], type = "l", lty = "solid", col = "blue")
            points(par.ts.median[[i]][[j]][,1], type = "l", lty = "dashed", col = "black")

            legend("topright",ncol = 3,bg = "grey",
                   lty = c("dotted", "solid", "dashed"),
                   col = c("red", "blue", "black"),
                   legend = c("95% HPD", "Posterior mean", "Posterior median"))
          }

          obj.par <- rbind(round(accept.prob.mean[[i]][[j]], 2),
                           par.mean[[i]][[j]],
                           par.median[[i]][[j]],
                           par.sd[[i]][[j]])
          rownames(obj.par) <- c("acc.prob", "par.mean",
                                 "par.median", "par.sd")
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
          print(obj.par)
          cat("\n")
          print(obj)
        }
      }
    }

    if(iIter == nIter)
    {
      cat(rep("-", dev.width), "\n\n",  sep = "")
    }
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
