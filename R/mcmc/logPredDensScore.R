##' Log predictive density score
##'
##' @param logPredLst "list"
##'
##'        The log posterior predictive sampler from MCMC.
##'
##' @return "list"
##'
##'        The LPDS and the numerical standard error.
##'
##' @references Li 2010
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Nov 30 23:18:18 CET 2010;
##'       Current:       Tue Nov 30 23:18:25 CET 2010.
##' TODO: Don't use the full draws when nIter is large.
##'       Still doubt the nseLPDS
logPredDensScore <- function(logPredLst)
{

###----------------------------------------------------------------------------
### Calculate the LPDS and its numerical standard error
###----------------------------------------------------------------------------

  ## Check if the log PredLst is from time series.
  ## One fold and within each fold multi-columns
  nFold <- length(logPredLst)
  dim0 <- logPredLst[[1]]
  ## LikLst.len <- dim0[2]
  nMCMCSample <- dim(logPredLst[[1]])[1]

  ## LPDS.vec <- matrix(NA, LikLst.len, 1)
  ## LPDS.nse.vec <- matrix(NA, LikLst.len, 1)
  ## for(i in 1:LikLst.len)
  ##   {

  ## if(partiMethod  == "time-series")
  ##   {
  ##     logPredMatrix <- logPredLst[[1]][, i]
  ##   }
  ## else
  ##   {
  ##     ## matrix NOT from time series settings or only one observation in
  ##     ## the prediction
  logPredMatrix <- matrix(unlist(logPredLst), nMCMCSample, nFold)
  ## }

  ## The LPDS
  scaleFactors <- apply(logPredMatrix,2, median)
  scaleMatrix <- matrix(scaleFactors, nMCMCSample, nFold)

  expPredMatrix <- exp(logPredMatrix-scaleMatrix)
  expPredMatrix[is.infinite(expPredMatrix)] <- NA # TODO: Think about it carefully
  expPredMean <- colMeans(expPredMatrix, na.rm = TRUE)
  LPDS <- mean(scaleFactors + log(expPredMean))

  if (sum(is.infinite(expPredMatrix)) > nMCMCSample*.05)
    {
      warning("Too many infinite predict densities produced. The LPDS may not be accurate.")
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

  ## LPDS.nse.vec[i] <- nseLPDS
  ## }

  ## LPDS <- sum(LPDS.vec)
  ## LPDS.nse <- sum(LPDs.nse.vec)

  out <- list(LPDS = LPDS, LPDS.nse = LPDS.nse)
  return(out)
}
