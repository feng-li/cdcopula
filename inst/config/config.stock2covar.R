###----------------------------------------------------------------------------
### Options of configures
###----------------------------------------------------------------------------

## THE SOURCE OF DATA
RawDataPath <- list(SP600 = "data/raw/SP600-SML.csv",
                    SP100 = "data/raw/SP100-OEX.csv")

## EXTRACT SUBSETS TO USE
nObs <- 100

## STANDARDIZE THE DATA
StandardizeData <- "norm-0-1"

## SAVE TO RData FORMAT
save2disk <- "data/SP100-SP600-n100.Rdata"

###----------------------------------------------------------------------------
### The script
###----------------------------------------------------------------------------

## Convert the raw data into covariates
Scovar <- lapply(X = RawDataPath, FUN = stock2covariates)

## Subtract the subset of the file
nTotal <- nrow(Scovar[[1]])
SubObs <- (nTotal-nObs+1):nTotal
XRaw <- list()
YRaw <- list()
ID <- Scovar[[1]][SubObs, "Date"]

for(i in names(Scovar))
{
  CovCurr <- Scovar[[i]]
  nCol <- ncol(CovCurr)
  XRaw[[i]] <- as.matrix(CovCurr[SubObs, 3:nCol, drop = FALSE])
  YRaw[[i]] <- as.matrix(CovCurr[SubObs, "Returns"])
}

## Standardize the data
if(!(StandardizeData == FALSE))
{
  XNew <- lapply(XRaw, StdData, method = StandardizeData)
  YNew <- lapply(YRaw, StdData, method = StandardizeData)

  X <- lapply(XNew, function(x) x$data)
  Y <- lapply(YNew, function(x) x$data)
  X.config <- lapply(XNew, function(x) x$config)
  Y.config <- lapply(YNew, function(x) x$config)
}else
{
  X <- XRaw
  Y <- YRaw
  X.config <- NA
  Y.config <- NA
}

## Save to file
if(is.character(save2disk))
  {
    save(ID, X, Y,  X.config, Y.config, file = save2disk)
  }
###----------------------------------------------------------------------------
### End of the script
###----------------------------------------------------------------------------
