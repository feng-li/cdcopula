parCaller <- function(parCurrUpdate)
{
  ## check which parameter need to update
  CompNM <- names(parCurrUpdate)
  for(i in CompNM)
    {
      parNM <- names(parCurrUpdate[[i]])
      parUpdateIdx <- parNM[parCurrUpdate[[i]] == TRUE] 
      for(j in parUpdateIdx)
        {
          out <- c(i, j) 
        }
    }
  return(out)
}
