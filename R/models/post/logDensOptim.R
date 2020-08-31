## This is the optimization version. Note that this should only work with one observation.
#' @export
logDensOptim <- function(x, jPar, Mdl.MargisType, Mdl.Y,
                         Mdl.par, Mdl.u, Mdl.d, parUpdate,
                         MCMC.UpdateStrategy)
{
  ## There should be only one chain in parUpdate

  chainCaller <- parCplRepCaller(parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  Mdl.par[[CompCaller]][[parCaller]][, jPar] <- x

  Mdl.ud <- logDens(Mdl.MargisType = Mdl.MargisType,
                    Mdl.Y = Mdl.Y,
                    Mdl.par = Mdl.par,
                    Mdl.u = Mdl.u,
                    Mdl.d = Mdl.d,
                    parUpdate = parUpdate,
                    MCMC.UpdateStrategy = MCMC.UpdateStrategy)

  Mdl.d <- Mdl.ud[["Mdl.d"]]
  Mdl.PostComp <- Mdl.ud[["Mdl.PostComp"]]
  out <- sum(Mdl.d[, unlist(Mdl.PostComp)])

  return(out)
}
