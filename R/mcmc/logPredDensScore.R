##' Log predictive density score
##'
##' @param logPredMatrix
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
logPredDensScore <- function(logPredMatrix)
{

###----------------------------------------------------------------------------
### Calculate the LPDS and its numerical standard error
###----------------------------------------------------------------------------
  mat.dim <- dim(logPredMatrix)

  nSample <- mat.dim[1]
  nFold <- mat.dim[2]
  ## nUsed <- xx # No. of MCMC used taking out burnin.

  ## The LPDS
  scaleFactors <- apply(logPredMatrix, 2, median)
  scaleMatrix <- matrix(scaleFactors, nSample, nFold)
  expPredMatrix <- exp(logPredMatrix-scaleMatrix)
  expPredMatrix[is.infinite(expPredMatrix)] <- NA # TODO: Think about it carefully
  expPredMean <- colMeans(expPredMatrix, na.rm = TRUE)
  LPDS <- mean(scaleFactors + log(expPredMean))

  ## If at least 5% are infinity, maybe due to unsuccessful MCMC.
  if (sum(is.infinite(expPredMatrix)) > nSample*.05)
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
  var.MeanexpPredMatrix <- predVar/nSample*(1+2*ACFxx)
  nvarLPDS <- 1/nFold^2*sum(1/(expPredMean)^2*var.MeanexpPredMatrix)
  nseLPDS <- sqrt(nvarLPDS)

  out <- list(LPDS = LPDS, nseLPDS = nseLPDS)
  return(out)
}
