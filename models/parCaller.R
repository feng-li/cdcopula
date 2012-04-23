parCaller <- function(parUpdate)
{
  ## check which parameter need to update
  CompNM <- names(parUpdate)
  for(i in CompNM)
    {
      parNM <- names(parUpdate[[i]])
      parUpdateIdx <- parNM[parUpdate[[i]] == TRUE]
      for(j in parUpdateIdx)
        {
          out <- c(i, j)
        }
    }
  return(out)
}
