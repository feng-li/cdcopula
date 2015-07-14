## A simple wrapper that calculates the numerical gradient for given parameters
## NOTE: The numeric gradient depends on numDeriv package
MargiModelGradNum <- function(y, par, type, parCaller, densCaller)
{
  MargiModelGradNumFun.subtask <- function(subtask, parCaller, parCurr, yCurr, typeCurr, densCaller)
  {
    nSubTask <- length(subtask)
    out <- matrix(NA, nSubTask, length(densCaller), dimnames = list(NULL, densCaller))

    require("numDeriv")
    MargiModelGradNumFun <- function(x, parCaller, parCurr, yCurr, typeCurr, densCallerCurr)
    {## browser()
      parCurr[[parCaller]] <- x

      ## We only have one output, so [[1]] did the trick.
      MargiLogLikObs <- MargiModel(par = parCurr, y = yCurr,
                                   type = typeCurr, densCaller = densCallerCurr)[[1]]
      out <- MargiLogLikObs
      return(out)
    }

    for(i in 1:nSubTask)
    {
      for(j in densCaller)
      {
        gradTry <-  try(grad(func = MargiModelGradNumFun,
                             x = par[[parCaller]][i],
                             parCaller = parCaller,
                             parCurr = lapply(par, function(x, i)x[i], i = i),
                             yCurr = y[subtask][i],
                             typeCurr = type,
                             densCallerCurr = j), silent = TRUE)

        if(is(gradTry, "try-error"))
        {
          out[i, j] <- NA
        }
        else
        {
          out[i, j] <- gradTry
        }
      }

    }
    return(out)
  }

  ## This code is used for checking if the analytical gradients are correctly coded, not
  ## designed for parallel computing. To parallel the numerical gradients, use
  ## "logDensGradHessNum()"

  out <- MargiModelGradNumFun.subtask(subtask = 1:length(y), parCaller = parCaller,
                                      yCurr = y, typeCurr = type, densCaller = densCaller)
  return(out)
}
