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
  else
  {
    stop("No such model implemented!")
  }
  return(out)
}
