#' Value at Risk based on MCMC draws
#'
#' Calculate the Value at Risk and Expected Shortfall
#' @param CplFitted Fitted Results
#' @param Mdl.Idx.testing Testing index
#' @param probs Probability
#' @param args Other arguments
#' @return NA
#' @references NA
#' @author Feng Li, Central University of Finance and Economics.
#' @export
CplMCMC.VaR <- function(CplFitted, Mdl.Idx.testing, probs, args)
{
    method <- args[["method"]]

    if(tolower(method) ==  "linear")
    {
        ## Huang et al (2009) Insurance: Mathematics and Economics
        weights <- args[["weights"]]
        if(sum(weights) != 1)
        {
            stop("Weights for portfolio should sum to one.")
        }
        Mdl.ConfigEnv <- CplFitted[["Mdl.ConfigEnv"]]
        if(length(Mdl.ConfigEnv) == 0)
        {
            stop("Mld.ConfigEnv should be supplied via Cplfitted")
        }

        ## attach(CplFitted[["Mdl.ConfigEnv"]])
        Mdl.Y <- Mdl.ConfigEnv[["Mdl.Y"]]


        Portfolio.Real <- (Mdl.Y[[1]][Mdl.Idx.testing]*weights[1]+
                           Mdl.Y[[2]][Mdl.Idx.testing]*weights[2])



        replicate = args[["replicate"]]

        nVaR = length(probs)
        VaR.ary = array(NA, c(nVaR, length(Mdl.Idx.testing), replicate))
        for(i in 1:replicate)
        {
            YPredSiml <- rCplPred(CplFitted, Mdl.Idx.testing)
            Portfolio.Sim <- YPredSiml[, 1, ]*weights[1] + YPredSiml[, 2, ]*weights[2]
            VaR.ary[, , i] = apply(Portfolio.Sim, 2, quantile, probs = probs)

        }
        VaR = apply(VaR.ary, c(1, 2), mean)

        plot = args[["plot"]]
        if(length(plot)>0 && plot == TRUE)
        {
            ylim = c(quantile(c(Portfolio.Real, Portfolio.Sim), 0.01),
                     quantile(Portfolio.Real, 0.99))
            par(cex = 0.8)

            col = c("red", "purple", "black", "green", colors())

            plot(VaR[1, ], type = "l", col = "red", ylim = ylim, ylab = "VaR")

            if(nVaR>1)
            {
                for(iVaR in 2:nVaR)
                {
                    points(VaR[iVaR, ], type = "l", col = col[iVaR], lty = i )
                }
            }
            points(Portfolio.Real, type = "p", pch = 20, col = "blue")
            legend("topleft",
                   legend = c(paste("Portfolio with weights: (",
                                    weights[1], ", " , weights[2], ")", sep = ""),
                              paste("VaR ", probs*100,  "%", sep = "")),
                   col = c("blue", col[1:nVaR]),
                   lty = c(3, 1:nVaR),
                   lwd = c(2, rep(1.5, nVaR)),
                   ncol = 1,
                   bg = "grey")
        }


    }
    else if(tolower(method) ==  "multivariate")
    {
        ## Cousin Bernardino, 2013 JMVA
        stop("Not implemented yet.")
    }
    else
    {
        stop("No such method!")
    }

    out <- list(Portfolio.Real = Portfolio.Real,
                VaR = VaR)
    return(out)
}


## TESTING
## rm(list = ls())
## load(file.path("~/running/", "JOB835.^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross0MCMC.nIter10000joint+20160720@14.44.5c1da7.Rdata"))
## load(file.path("~/running/", "^SML^OEX+SPLITTSPLITTBB7+nObs6557nCross1MCMC.nIter1000twostage+20160522@23.11.c168ab.Rdata"))

## rm(list = ls())
## load(file.path("~/running/", "^SML^OEX+NANADCCGARCH+nObsNAnCross1+20160519@20.58.58aebe.Rdata"))
## sourceDir("~/code/cdcopula/R", recursive = TRUE)

## CplMCMC.VaR(OUT.FITTED[[1]], probs = c(0.01, 0.05, 0.1, 0.2),
##         args = list(method = "linear", weights = c(0.5, 0.5), plot = TRUE))
