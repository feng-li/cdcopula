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
##' @note Initial: Fri Feb 01 14:49:15 CET 2013;
##'       Current: Fri Feb 01 14:49:33 CET 2013.
##' TODO: write this function as a summary
CplMCMC.summary <- function(nIter, iIter = nIter, interval = 0.1, burnin, ...)
{
    ## Set the print interval and consider burnin
    printIter <- seq(from = floor(nIter*interval),
                     to = nIter,
                     by = floor(nIter*interval))
    if(printIter[length(printIter)] != nIter)
        {
            printIter <- c(printIter, nIter)
        }

    dev.width <- getOption("width")

    ## The burning
    n.burn.default <- round(nIter*burnin)
    n.burn <- ifelse(iIter>n.burn.default, n.burn.default, 0)

    if(iIter %in% printIter)
        {
            donePercent <- round(iIter/nIter*100)

            welcome <- paste("MCMC SUMMARY: ", donePercent, "% (",
                             round(n.burn/nIter*100), "% burnin)\n", sep = "")
            cat("\n", rep("=", dev.width), "\n", sep = "")
            cat(welcome)
            cat(rep("=", dev.width), "\n", sep = "")

            ## format.name0 <- c("Done(%)|", parsList)
            ## format.name <- format(c("Done(%)|", parsList), width = 13, justify = "right")
            ## cat(format.name, "\n")
            ## }

            par <- list(...)
            MCMC.betaIdx <- par$MCMC.betaIdx
            MCMC.beta <- par$MCMC.beta
            MCMC.par <- par$MCMC.par
            MCMC.AccProb <- par$MCMC.AccProb
            MCMCUpdate <- par$MCMCUpdate

            accept.prob.mean <- rapply(MCMC.AccProb,
                                       function(x, iIter) mean(x[(n.burn+1):iIter]),
                                       how = "replace", iIter = iIter)
            par.mean <- rapply(MCMC.par,
                               function(x, iIter) mean(x[(n.burn+1):iIter]),
                               how = "replace", iIter = iIter)
            par.sd <- rapply(MCMC.par,
                             function(x, iIter) sd(x[(n.burn+1):iIter]),
                             how = "replace", iIter = iIter)

            betaIdx.mean <- rapply(MCMC.betaIdx,
                                   function(x, iIter){
                                       colMeans(x[(n.burn+1):iIter, , drop = FALSE])},
                                   how = "replace", iIter = iIter)

            beta.mean <- rapply(MCMC.beta,
                                function(x, iIter){
                                    colMeans(x[(n.burn+1):iIter, , drop = FALSE])},
                                how = "replace", iIter = iIter)

            beta.sd <- rapply(MCMC.beta,
                              function(x, iIter){
                                  colSds(x[(n.burn+1):iIter, , drop = FALSE])},
                              how = "replace", iIter = iIter)

            ## Efficiency factor of MCMC
            colIneffs <- function(x) {apply(x, 2, ineff)}
            beta.ineff <- rapply(MCMC.beta,
                                 function(x, iIter){
                                     colIneffs(x[(n.burn+1):iIter, , drop = FALSE])},
                                 how = "replace", iIter = iIter)

            for(i in names(MCMC.beta))
                {
                    for(j in names(MCMC.beta[[i]]))
                        {
                            if(MCMCUpdate[[i]][[j]])
                                {
                                    obj.par <- rbind(round(accept.prob.mean[[i]][[j]], 2),
                                                     par.mean[[i]][[j]],
                                                     par.sd[[i]][[j]])
                                    rownames(obj.par) <- c("acc.prob", "par.mean", "par.sd")
                                    colnames(obj.par) <- ""

                                    obj <- rbind(beta.mean[[i]][[j]],
                                                 beta.sd[[i]][[j]],
                                                 betaIdx.mean[[i]][[j]],
                                                 beta.ineff[[i]][[j]])
                                    rownames(obj) <- c("beta.mean", "beta.sd", "betaIdx.mean", "beta.ineff")

                                    cat("\n", i, j, "(", donePercent, "% )\n")
                                    cat(rep("-", dev.width), sep = "")
                                    print(obj.par)
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
