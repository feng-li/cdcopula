ModelForeignPred <- function(model, fitted.model, data.pred, ...)
{
  if(tolower(model)  == "gogarch")
  {
    dataMat <- do.call("cbind", data.pred)
    n.ahead <- nrow(dataMat)
    PredDens <- matrix(NA, n.ahead, 1)


    pred <- gogarchforecast(fitted.model, n.ahead = n.ahead, n.roll = 0,
                            external.forecasts = list(mregfor = NULL),
                            cluster = NULL, ...)

    ## "mean", "variance", "skewness",  "kurtosis" for p(Y_p|Y_old)
    MVSK <- list()
    MVSK[["mean"]] <- pred@mforecast[["mu"]][, , 1]
    MVSK[["var"]] <- pred@mforecast[["factor.sigmas"]][, , 1]

    if(fitted.model@model[["modeldesc"]][["distribution"]]  == "mvnorm")
    {
      for(i in 1:n.ahead)
      {

        PredDens[i, ] <- dmvnorm(x = dataMat[i, ],
                                 mean = MVSK[["mean"]][i, ],
                                 sigma = diag(MVSK[["var"]][i ,], nrow = ncol(dataMat)),
                                 log = TRUE)
      }
    }
    else
    {
      stop("Not implemented yet!")
    }
    Mdl.logPredDens <- matrix(sum(PredDens))


    RESID <- list()
    RESID[["mean"]] <- dataMat-MVSK[["mean"]]
  }
  else if(tolower(model) == "dccgarch")
  {
    dataMat <- do.call("cbind", data.pred)
    n.ahead <- nrow(dataMat)

    PredDens <- matrix(NA, n.ahead, 1)

    pred <- dccforecast(fitted.model, n.ahead = n.ahead, n.roll = 0,
                        external.forecasts = list(mregfor = NULL,
                                                  vregfor = NULL),
                        cluster = NULL, ...)


    ## "mean", "variance", "skewness",  "kurtosis" for p(Y_p|Y_old)
    MVSK <- list()
    MVSK[["mean"]] <- pred@mforecast[["mu"]][, , 1]
    MVSK[["var"]] <- pred@mforecast[["H"]][[1]]


    if(fitted.model@model[["modeldesc"]][["distribution"]]  == "mvnorm")
    {
      for(i in 1:n.ahead)
      {
        PredDens[i, ] <- dmvnorm(x = dataMat[i, ],
                                 mean = MVSK[["mean"]][i, ],
                                 sigma = MVSK[["var"]][, , i],
                                 log = TRUE)
      }
    }
    else
    {
      stop("Not implemented yet!")
    }

    Mdl.logPredDens <- matrix(sum(PredDens))

    RESID <- list()
    RESID[["mean"]] <- dataMat-MVSK[["mean"]]

  }
  else
  {
    stop("No such model implemented!")
  }


  out <- list(MVSK = MVSK,
              RESID = RESID,
              Mdl.logPredDens = Mdl.logPredDens)
  return(out)
}
