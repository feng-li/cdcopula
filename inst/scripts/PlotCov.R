
load("data/BABA-Texts.Rdata")
par(mfrow = c(2, 3))

plot(Y[["BABA"]], X[["BABA"]][, "RMA1"],
     xlab = "Last Day's Return", ylab = "Return", main = "",
     col = "blue", pch = 20)

plot(Y[["BABA"]], X[["BABA"]][, "RMA5"],
     xlab = "Last Five Day's Return", ylab = "Return", main = "",
     col = "blue", pch = 20)

plot(Y[["BABA"]], X[["BABA"]][, "RMA20"],
     xlab = "Last Month's Return", ylab = "Return", main = "",
     col = "blue", pch = 20)

plot(Y[["BABA"]], X[["BABA"]][, "MaxMin95"],
     xlab = "MaxMin95", ylab = "Return", main = "",
     col = "blue", pch = 20)

plot(Y[["BABA"]], X[["BABA"]][, "CloseAbs95"],
     xlab = "CloseAbs95", ylab = "Return", main = "",
     col = "blue", pch = 20)

plot(Y[["BABA"]], X[["BABA"]][, "CloseSqrt95"],
     xlab = "CloseSqrt95", ylab = "Return", main = "",
     col = "blue", pch = 20)

dev.copy2eps(file = "BABACovPlot.eps")

par(mfrow = c(3, 3), mar = c(2, 1, 2, 0))
for(i in 1:9)
  {
    plot(X[[1]][, i], X[[2]][, i],
         main = colnames(X[[1]])[i],
         xlab = "", ylab = "",
         col = "blue", pch = 20)
  }
