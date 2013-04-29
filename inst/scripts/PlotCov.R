par(mfrow = c(1, 2))
plot(X[[1]][, "CloseAbs80"], X[[2]][, "CloseAbs80"],
     xlab = "SP100", ylab = "SP600", main = "CloseAbs80",
     col = "blue", pch = 20)
plot(X[[1]][, "MaxMin95"], X[[2]][, "MaxMin95"],
     xlab = "SP100", ylab = "", main = "MaxMin95",
     col = "blue", pch = 20)

dev.copy2eps(file = "MarginCorrPlot.eps")

par(mfrow = c(3, 3), mar = c(2, 1, 2, 0))
for(i in 1:9)
  {
    plot(X[[1]][, i], X[[2]][, i],
         main = colnames(X[[1]])[i],
         xlab = "", ylab = "",
         col = "blue", pch = 20)
  }
