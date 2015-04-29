## This is the optimization version. Note that this should only work with one observation.
logDensOptim <- function(x, jPar, chainCaller, MargisType, Mdl.Y,
                         Mdl.par, Mdl.u, Mdl.d, parUpdate,
                         MCMCUpdateStrategy)
{
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  Mdl.par[[CompCaller]][[parCaller]][, jPar] <- x

  Mdl.ud <- logDens(MargisType = MargisType,
                    Mdl.Y = Mdl.Y,
                    Mdl.par = Mdl.par,
                    Mdl.u = Mdl.u,
                    Mdl.d = Mdl.d,
                    parUpdate = parUpdate,
                    MCMCUpdateStrategy = MCMCUpdateStrategy)
  Mdl.d <- Mdl.ud[["Mdl.d"]]
  Mdl.PostComp <- Mdl.ud[["Mdl.PostComp"]]

  out <- sum(Mdl.d[, unlist(Mdl.PostComp)])
  return(out)
}
