parCplCaller <- function(parUpdate, parUpdateOrder)
{
  ## check which parameter need to update
  CompNM <- names(parUpdate)
  out.init <- NULL
  for(i in CompNM)
    {
      parNM <- names(parUpdate[[i]])
      parUpdateIdx <- parNM[parUpdate[[i]] == TRUE]
      for(j in parUpdateIdx)
        {

          out.init <- rbind(out.init , c(i, j))
        }
    }


  ## Order the parameters if asked
    if(!missing(parUpdateOrder))
    {
      parOrder <- unlist(parUpdateOrder)
      parUpdate <- unlist(parUpdate)
      parUpdateOrderVec <- order(parOrder[parUpdate])
      out <- out.init[parUpdateOrderVec, , drop = FALSE]
    }
    else
      {
        out <- matrix(out.init, , 2)
      }

  return(out)
}
