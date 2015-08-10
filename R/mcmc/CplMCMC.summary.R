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

  ## The burning
  n.burn.default <- round(nIter*burnin)
  n.burn <- ifelse(iIter>n.burn.default, n.burn.default, 0)


  subFun <- function(x, iIter, fun){
    apply(x[(n.burn+1):iIter, ,drop = FALSE], 2, fun)
  }

  subFun3 <- function(x, iIter, fun){
    apply(x[(n.burn+1):iIter, , ,drop = FALSE], 3, fun)
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

    par.mean <- rapply(MCMC.par, subFun3, how = "replace", iIter = iIter,  fun = mean)
    par.sd <- rapply(MCMC.par, subFun3, how = "replace", iIter = iIter,  fun = sd)

    betaIdx.mean <- rapply(MCMC.betaIdx, subFun, how = "replace", iIter = iIter,
                           fun = mean)
    beta.mean <- rapply(MCMC.beta, subFun, how = "replace", iIter = iIter,  fun = mean)
    beta.sd <- rapply(MCMC.beta, subFun, how = "replace", iIter = iIter,  fun = sd)


    ## Efficiency factor of MCMC
    beta.ineff <- rapply(MCMC.beta, subFun, how = "replace", iIter = iIter,
                         fun = ineff)

    for(i in names(MCMC.beta))
    {
      for(j in names(MCMC.beta[[i]]))
      {
        if(MCMCUpdate[[i]][[j]])
        {
          obj.par <- rbind(round(accept.prob.mean[[i]][[j]], 2),
                           par.mean[[i]][[j]],
                           par.sd[[i]][[j]])
          rownames(obj.par) <- c("acc.prob", "par.mean",
                                 "par.sd")
          colnames(obj.par) <- names(par.mean[[i]][[j]])

          obj <- rbind(beta.mean[[i]][[j]],
                       beta.sd[[i]][[j]],
                       betaIdx.mean[[i]][[j]],
                       beta.ineff[[i]][[j]])

          rownames(obj) <- c("beta.mean", "beta.sd", "betaIdx.mean", "beta.ineff")
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

}
