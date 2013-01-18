## This script plots the daily return of SP100 and SP600 data

nObs <- length(X.ID)

ylim <- c(-10, 10)

par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
plot(X.ID, Y[[1]], col = "blue", ylim = ylim, type = "l",
     xlab = "", ylab = "SP100", axes = TRUE)

plot(X.ID, Y[[2]], col = "blue", ylim = ylim, type = "l",
     xlab = "Time", ylab = "SP600", axes = TRUE)
