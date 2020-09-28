#' Log predictive density score
#'
#' @param logPredDensLst "list"
#'
#'        The log posterior predictive sampler from MCMC for n-fold.
#'
#' @return "matrix"
#'
#'        The LPDS and the numerical standard error.
#'
#' @references Li 2010
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @note First: Tue Nov 30 23:18:18 CET 2010; Current: Tue Sep 15 17:49:03 CST 2015
#' @export
logPredDensScore <- function(logPredDensLst)
{
###----------------------------------------------------------------------------
### Calculate the LPDS and its numerical standard error
###----------------------------------------------------------------------------
    ## NOTE: LPDS in a copula model: LPDS_full = sum(LPDS_marginals) + LPDS_copula

    ## Check if the log PredLst is from time series.  One fold and within each fold
    ## multi-columns
    nFold <- length(logPredDensLst)
    if(!is.matrix(logPredDensLst[[1]]))
    {
        stop("PredDensLst must be a list containing matrices.")
    }
    dim0 <- logPredDensLst[[1]]

    nMCMCSample <- dim(logPredDensLst[[1]])[1]

    logPredMatrix <- matrix(unlist(logPredDensLst), nMCMCSample, nFold)

    scaleFactors <- apply(logPredMatrix,2, median)
    scaleMatrix <- matrix(scaleFactors, nMCMCSample, nFold)

    expPredMatrix <- exp(logPredMatrix-scaleMatrix)
    expPredMatrix[!is.finite(expPredMatrix)] <- NA # TODO: Think about it carefully
    expPredMean <- colMeans(expPredMatrix, na.rm = TRUE)

    ## The LPDS and numerical standard error
    LPDS <- mean(scaleFactors + log(expPredMean))

    if (sum(!is.finite(expPredMatrix)) > nMCMCSample*.05)
    {
        warning("Too many NA/NaN or Inf predict densities produced. The LPDS may not be accurate.")
    }

    ## The numerical standard error for LPDS.
    ## See Box Jenkins (2007)  p.31
    ACFxx <- matrix(0, 1, nFold)
    for(k in 1:nFold)
    {
        acfSample.tmp0 <- expPredMatrix[, k]
        acfSample.tmp <- acfSample.tmp0[acfSample.tmp0<1e100] # Just let it numerically stable.

        predACF <- acf(acfSample.tmp, plot = FALSE,  na.action = na.pass)$acf
        nlagk <- length(predACF)
        for(kk in 0:(nlagk-1))
        {
            ACFxx[k] <- ACFxx[k] + (1 - kk/nlagk)*predACF[kk +1]
        }
    }

    expPredMatrix.tmp <- expPredMatrix
    expPredMatrix.tmp[expPredMatrix.tmp>1e100] <- NA # numerically stable
    predVar <- apply(expPredMatrix.tmp, 2, var, na.rm = TRUE)

    ## var.MeanexpPredMatrix <- predVar/nUsed*(1+2*ACFxx)
    var.MeanexpPredMatrix <- predVar/nMCMCSample*(1+2*ACFxx)
    nvarLPDS <- 1/nFold^2*sum(1/(expPredMean)^2*var.MeanexpPredMatrix)
    LPDS.nse <- sqrt(nvarLPDS)

    out <- cbind(LPDS = LPDS, LPDS.nse = LPDS.nse)
    if(is.nan(LPDS.nse))
    {
        warning("NaN in LPDS.nse for this component probably caused by ill-behaved MCMC samples.")
    }

    return(out)
}
