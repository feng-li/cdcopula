
## Save path
save2diskPath <- "~/code/cdcopula/data/Stocks-Texts.Rdata"


## Load Discrete Data
load("~/code/stockandtexts/data/日期情感单词矩阵.Rdata")
rownames(DateSen) <- DateSen[, 1]

## Load Continuous Stock Data
load("~/code/cdcopula/data/BABA-Texts.Rdata")
Data.Stocks <- MVStocks(from = "20140801", to = "20150922", stocks = "BABA")

## Merge Data
ID.discrete <-  DateSen[, 1]
ID.continuous <- Data.Stocks[["ID"]]
ID <- intersect(ID.discrete, ID.continuous)


X.discrete <- list("BABA" = as.matrix(DateSen[ID, 5:ncol(DateSen)]))
Y.discrete <- list("BABA" = as.matrix(DateSen[ID, 2:4, drop = FALSE]))
X.discrete.config <- NA
Y.discrete.config = NA


X.continuous <- list("BABA" = Data.Stocks[["X"]][["BABA"]][ID, , drop = FALSE])
Y.continuous <- list("BABA" = Data.Stocks[["Y"]][["BABA"]][ID, , drop = FALSE])
X.continuous.config <- Data.Stocks[["X.config"]]
Y.continuous.config <- Data.Stocks[["Y.config"]]


save(ID, X.discrete, Y.discrete,
     X.continuous, Y.continuous,
     X.discrete.config , Y.discrete.config ,
     X.continuous.config, Y.continuous.config,
     file = save2diskPath, compress = "xz")
