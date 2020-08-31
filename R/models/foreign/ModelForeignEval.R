#' @export
ModelForeignEval <- function(model, spec, data, ...)
{
    if(tolower(model)  == "gogarch")
    {
        dataMat <- do.call("cbind", data)
        out <- gogarchfit(spec = spec, data = dataMat, ...)
    }
    else if(tolower(model) == "dccgarch")
    {
        dataMat <- do.call("cbind", data)
        out <- dccfit(spec = spec,
                      data = dataMat,
                      out.sample = 0,
                      solver = "solnp",
                      solver.control = list(),
                      fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE),
                      cluster = NULL, fit = NULL, VAR.fit = NULL, realizedVol = NULL)
    }
    else if(tolower(model) == "msbvar")
    {
        dataMat <- do.call("cbind", data)

        spec[["Y"]] <- as.ts(dataMat)
        model <- do.call(msbvar, spec)
        out <- gibbs.msbvar(model)
    }
    else if(tolower(model) == "var")
    {
        dataMat <- do.call("cbind", data)
        spec[["x"]] <- dataMat
        out <- do.call(VAR, spec)
    }
    else
    {
        stop("No such foreign multivariate model implemented!")
    }
    return(out)
}
