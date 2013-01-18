## This script is used to construct the Financial data from Yahoo Finance for
## the copula model
################################################################################
## In Yahoo Finance,  the first month is zero
## The starting day: a=07&b=2&c=1982  reads: 1982 August 02
## The ending dat:   d=00&e=17&f=2013 reads: 2013 January 17

StartingDay <- "a=07&b=16&c=1995" #The earliest day is "a=07&b=16&c=1995"
EndingDay <- "d=00&e=17&f=2013"

## The Data ticker used in Yahoo Finance
DataTicker <- c("OEX", "SML") # SP100 and SP600

## Standardized data method
StdData.method <- "norm-0-1"

## Save File
SavePath <- "SP100-SP600-20130116.Rdata"

## Download the data
Data <- list()
X <- list()
Y <- list()
X.config <- list()
X.date <- list()
for(i in DataTicker)
  {
    DataPath.i <- paste("http://ichart.finance.yahoo.com/table.csv?s=%5E",
                        i, "&", StartingDay, "&", EndingDay,
                        "&g=d&ignore=.csv", sep ="")
    Data[[i]] <- stock2covariates(file = DataPath.i)
  }

## Match the date exactly
## There might be situations that stock is not available in some days
ID <- intersect(Data[[1]]$ID, Data[[2]]$ID)

for(i in DataTicker)
  {
    Idx <- Data[[i]]$ID%in%ID
    X.iraw <- Data[[i]]$X[Idx,  , drop = FALSE]
    Data.i <- StdData(X.iraw, method = StdData.method)
    X[[i]] <- Data.i$data
    X.config[[i]] <- Data.i$config
    Y[[i]] <- as.matrix(Data[[i]]$Y[Idx])
  }

## Save the ID in date format
X.ID <- as.Date(ID, "%Y-%m-%d")

save(X, Y, X.config, X.ID = X.ID, file = SavePath)
## Construct Copula model style dataset
