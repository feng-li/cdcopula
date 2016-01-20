load("~/code/stockandtexts/data/日期情感单词矩阵.Rdata")
save2diskPath <- "~/code/cdcopula/data/BABA-Texts.Rdata"

rownames(DateSen) <- DateSen[, 1]

Data.Stocks <- MVStocks(from = "20140801", to = "20150922", stocks = "BABA")

ID.Texts <-  DateSen[, 1]
ID.Stocks <- Data.Stocks[["ID"]]

ID <- intersect(ID.Texts, ID.Stocks)

X <- list("BABA" = Data.Stocks[["X"]][["BABA"]][ID, , drop = FALSE],
          "TEXTS"  =  as.matrix(DateSen[ID, 5:ncol(DateSen)]))
Y <- list("BABA"  = Data.Stocks[["Y"]][["BABA"]][ID, , drop = FALSE],
          "TEXTS" = as.matrix(DateSen[ID, 2, drop = FALSE])) # Positive
X.config <- Data.Stocks[["X.config"]]
Y.config <- Data.Stocks[["Y.config"]]

save(ID, X, Y,  X.config, Y.config, file = save2diskPath, compress = "xz")
