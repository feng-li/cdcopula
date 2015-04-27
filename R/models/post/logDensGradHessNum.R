## APPROACH TWO: This version calculates the numerical gradient (for log
## posterior) via the full log likelihood function
logDensGradHessNum <- function(par,jPar,iObs,...)
  {
    ## List and import all necessary variables. This order is very
    ## important. The data in local environment should always be nested inside
    ## global environment.
    ## list2env(data.global.env, envir = environment())
    ## list2env(data.parent.env, envir = environment())

    require("numDeriv")

    Mdl.u <- staticCache[["Mdl.u"]][iObs, , drop = FALSE]
    Mdl.d <- staticCache[["Mdl.d"]][iObs, , drop = FALSE]

    gradTry <-  try(grad(func = logDensitiesOptim,
                         x = ipar,
                         jPar = jPar,
                         CplNM = CplNM,
                         Mdl.Y = lapply(Mdl.Y, function(x, i)x[i, , drop = FALSE], i = iObs),
                         Mdl.par = lapply(Mdl.par, function(x, i)x[i, , drop = FALSE], i = iObs),
                         Mdl.u = lapply(Mdl.u, function(x, i)x[i, , drop = FALSE], i = iObs),
                         Mdl.d = lapply(Mdl.d, function(x, i)x[i, , drop = FALSE], i = iObs),
                         parUpdate = parUpdate,
                         chainCaller = chainCaller,
                         MCMCUpdateStrategy = MCMCUpdateStrategy)
                   ,silent = TRUE)

    if(is(gradTry, "try-error"))
      {
        out <- NA
      }
    else
      {
        out <- gradTry
      }


    logLikelihoodGrad.num <- mapply(FUN = logDensitiesGradNumFun,
                                    par = Mdl.par[[CompCaller]][[parCaller]],
                                    jPar = 1:ncol(Mdl.par[[CompCaller]][[parCaller]]),
                                    MoreArgs = list(
                                        CplNM = CplNM,
                                        Mdl.Y = Mdl.Y,
                                        Mdl.par = Mdl.par,
                                        staticCache = staticCache.curr,
                                        parUpdate = parUpdate,
                                        chainCaller = chainCaller,
                                        MCMCUpdateStrategy = MCMCUpdateStrategy))
    return(out)
  }
