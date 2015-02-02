qCpl <- function(u, parMargis, MargisTypes)
  {
    p <- length(MargisTypes)
    xOut <- vector("list", p) # The output matrix for latent variable
    names(xOut) <- names(parMargis)
    for(i in 1:p)
      {
        margiType <- MargisTypes[i]
        if(tolower(margiType) == "gaussian") # continues of course
          {
            mu <- parMargis[[i]][["mu"]] # scaler
            sigma <- parMargis[[i]][["sigma"]]   # scaler
            xOut[[i]] <- matrix(qnorm(u[, i], mean = mu, sd = sigma))
          }
        else if(tolower(margiType) == "poisson") # the discrete case
          {
            ## use latent variable instead.
          }
        else
          {
            stop("The quantile for the marginal density is not implemented yet.")
          }
      }

    ## The output
    return(xOut)
  }
