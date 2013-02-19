##' Trajectory MCMC.
##'
##' This function can also be used for summarizing the posterior results.
##' @param 0Iter "integer"
##' @param nIter "integer"
##' @param iIter "integer"
##' @param interval e.g. 10%
##' @param ...
##' @return A summary object
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Initial: Fri Feb 01 14:49:15 CET 2013;
##'       Current: Fri Feb 01 14:49:33 CET 2013.
##' TODO: write this function as a summary
MCMC.trajectory <- function(nIter, iIter = nIter, interval = 0.1, ...)
{
  ## Set the print interval
  printIter <- seq(from = floor(nIter*interval),
                   to = nIter,
                   by = floor(nIter*interval))
  if(printIter[length(printIter)] != nIter)
    {
      printIter <- c(printIter, nIter)
    }

  if(iIter %in% printIter)
    {
      donePercent <- round(iIter/nIter*100)

      welcome <- paste("MCMC trajectory summary:", donePercent, "%\n", sep = " ")
      cat("\n=====================================================================\n")
      cat(welcome)
      cat(  "=====================================================================\n")

      ## format.name0 <- c("Done(%)|", parsList)
      ## format.name <- format(c("Done(%)|", parsList), width = 13, justify = "right")
      ## cat(format.name, "\n")
      ## }
      par <- list(...)
      MCMC.betaIdx <- par$MCMC.betaIdx
      MCMC.beta <- par$MCMC.beta
      MCMC.par <- par$MCMC.par
      MCMC.AccProb <- par$MCMC.AccProb

      accept.prob.mean <- rapply(MCMC.AccProb,
                                 function(x, iIter) mean(x[1:iIter]),
                                 how = "replace", iIter = iIter)
      par.mean <- rapply(MCMC.par,
                         function(x, iIter) mean(x[1:iIter]),
                         how = "replace", iIter = iIter)
      par.sd <- rapply(MCMC.par,
                   function(x, iIter) sd(x[1:iIter]),
                   how = "replace", iIter = iIter)

      betaIdx.mean <- rapply(MCMC.betaIdx,
                             function(x, iIter){
                               colMeans(x[1:iIter, , drop = FALSE])},
                             how = "replace", iIter = iIter)

      beta.mean <- rapply(MCMC.beta,
                          function(x, iIter){
                            colMeans(x[1:iIter, , drop = FALSE])},
                          how = "replace", iIter = iIter)

      beta.sd <- rapply(MCMC.beta,
                          function(x, iIter){
                            colSds(x[1:iIter, , drop = FALSE])},
                          how = "replace", iIter = iIter)

      for(i in names(MCMC.beta))
        {
          for(j in names(MCMC.beta[[i]]))
          {
            obj.par <- rbind(round(accept.prob.mean[[i]][[j]], 2),
                             par.mean[[i]][[j]],
                             par.sd[[i]][[j]])
            rownames(obj.par) <- c("acc.prob", "par.mean", "par.sd")
            colnames(obj.par) <- ""

            obj <- rbind(beta.mean[[i]][[j]], beta.sd[[i]][[j]],
                         betaIdx.mean[[i]][[j]])
            rownames(obj) <- c("beta.mean", "beta.sd", "betaIdx.mean")

            cat("\n", i, j, "\n")
            cat(".....................................................................")
            print(obj.par)
            print(obj)
          }
        }


      if(iIter == nIter)
            {
              cat("----------------------------------------------------------------------\n\n")
            }
    }

}
