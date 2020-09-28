#' @export
ModelForeignPred <- function(model, fitted.model, data.pred, ...)
{
    MVSK <- list()
    MVSK.fitted <- list()
    RESID <- list()
    Mdl.logPredDens <- NA

    if(tolower(model)  == "gogarch")
    {
        if(length(data.pred)>0)
        {
            dataMat <- do.call("cbind", data.pred)
            n.ahead <- nrow(dataMat)
            PredDens <- matrix(NA, n.ahead, 1)

            pred <- gogarchforecast(fitted.model, n.ahead = n.ahead, n.roll = 0,
                                    external.forecasts = list(mregfor = NULL),
                                    cluster = NULL, ...)

            ## "mean", "variance", "skewness",  "kurtosis" for p(Y_p|Y_old)
            MVSK[["mean"]] <- pred@mforecast[["mu"]][, , 1]
            MVSK[["var"]] <- pred@mforecast[["factor.sigmas"]][, , 1]

            if(fitted.model@model[["modeldesc"]][["distribution"]]  == "mvnorm")
            {
                for(i in 1:n.ahead)
                {
                    PredDens[i, ] <- dmvnorm(x = dataMat[i, ],
                                             mean = MVSK[["mean"]][i, ],
                                             sigma = diag(MVSK[["var"]][i ,],
                                                          nrow = ncol(dataMat)),
                                             log = TRUE)
                }
            }
            else
            {
                stop("Not implemented yet!")
            }
            Mdl.logPredDens <- matrix(sum(PredDens), 1, 1, dimnames = list(NULL, "joint"))
            RESID[["mean"]] <- dataMat-MVSK[["mean"]]
        }

        ## TODO: in-sample prediction

    }
    else if(tolower(model) == "dccgarch")
    {
        ## out-of-sample prediction
        if(length(data.pred)>0)
        {
            dataMat <- do.call("cbind", data.pred)
            n.ahead <- nrow(dataMat)

            PredDens <- matrix(NA, n.ahead, 1)

            pred <- dccforecast(fitted.model, n.ahead = n.ahead, n.roll = 0,
                                external.forecasts = list(mregfor = NULL, vregfor = NULL),
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

            Mdl.logPredDens <- matrix(sum(PredDens), 1, 1, dimnames = list(NULL, "joint"))

            RESID[["mean"]] <- dataMat-MVSK[["mean"]]
        }

        ## In-sample mean and variance
        MVSK.fitted <- list()
        MVSK.fitted[["mean"]] <- fitted.model@model[["mu"]]
        MVSK.fitted[["var"]] <- fitted.model@mfit[["H"]][[1]]

    }
    else if(tolower(model) == "msbvar")
    {
        ## browser()

        dataMat <- do.call("cbind", data.pred)
        n.ahead <- nrow(dataMat)
        PredDens <- matrix(NA, n.ahead, 1)

        n.burnin <- 1000
        nDraws <- 1000

        pred <- forecast(fitted.model, nsteps = n.ahead, N1 = n.burnin, N2 = nDraws)

        stop("Not done yet!")

    }
    else if(tolower(model) == "var")
    {
        ## out-of-sample prediction
        if(length(data.pred)>0)
        {

            dataMat <- do.call("cbind", data.pred)
            n.ahead <- nrow(dataMat)
            nDim <- ncol(dataMat)

            PredDens <- matrix(NA, n.ahead, 1)

            pred <- VARpred(fitted.model,  h  =  n.ahead,  orig  =  0,
                            Out.level  =  FALSE)

            MVSK <- list()
            MVSK[["mean"]] <- pred[["pred"]]
            MVSK[["var"]] <-  array(apply(pred$se.err, 1, diag), c(nDim, nDim, n.ahead))


            for(i in 1:n.ahead)
            {
                PredDens[i, ] <- dmvnorm(x = dataMat[i, ],
                                         mean = MVSK[["mean"]][i, ],
                                         sigma = MVSK[["var"]][, , i],
                                         log = TRUE)
            }

            Mdl.logPredDens <- matrix(sum(PredDens), 1, 1, dimnames = list(NULL, "joint"))
            RESID[["mean"]] <- dataMat-MVSK[["mean"]]

        }

        ## In-sample mean and variance
        MVSK.fitted <- list()
        MVSK.fitted[["mean"]] <- NA
        MVSK.fitted[["var"]] <- NA

    }

    else
    {
        stop("No such model implemented!")
    }


    out <- list(MVSK = MVSK,
                MVSK.fitted = MVSK.fitted,
                RESID = RESID,
                Mdl.logPredDens = Mdl.logPredDens)
    return(out)
}
