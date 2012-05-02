##' Calculate the CDF and PDF of the marginal distribution.
##'
##' The detailed description can be found in the main setting file for each
##' input variable.
##' @title CDF of marginal model.
##' @param Mdl.Y "list".
##'        The response variables for the marginal model.
##' @param MargisTypes "vector" with character input.
##'        The types of the marginal model.
##' @param parMargis "list".
##'        The parameters input for the marginal model.
##' @param whichMargis
##'        Which marginals are going to update?
##' @return "matrix" with marginal names attributed.
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Jan 17 19:27:25 CET 2012;
##'       Current: Tue Jan 17 19:27:30 CET 2012.
MargisModels <- function(Mdl.Y, MargisTypes, parMargis, whichMargis =
                         names(Mdl.Y),  staticArgs)
  {
    ## The marginal name to be updated.
    ## MargisUpNM <- names(parMargis)[whichMargis]

    ## No. of observations
    ## nObs <- length(Mdl.Y[[1]])

    ## Fetch previous values for u and d.
    Mdl.ud <- staticArgs[c("Mdl.u", "Mdl.d")]

    ## Loop over all possible marginal models
    for(i in whichMargis)
      {
        if(tolower(MargisTypes[i]) == "gaussian")
          {
            ## The observed data in the marginal model
            y <- Mdl.Y[[i]]

            ## The mean and standard deviation for Gaussian density
            mu <- parMargis[[i]][["mu"]] # scaler
            sigma <- parMargis[[i]][["sigma"]]   # scaler

            ## The percentile representation
            Mdl.ud[["Mdl.u"]][, i] <- pnorm(y, mean = mu, sd = sigma, log = FALSE)

            ## The quantile representation
            Mdl.ud[["Mdl.d"]][, i] <- dnorm(y, mean = mu, sd = sigma, log = TRUE)
          }
        else if (tolower(margiType) == "student-t")
          {
            ## the student t version
          }
        else
          {
            stop("This type of margin is not implemented.")
          }
      }

    ## The output
    return(Mdl.ud)
  }
