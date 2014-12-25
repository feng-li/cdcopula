## This script plots the daily return of SP100 and SP600 data

## Plot The Returns
nObs <- length(ID)

ylim <- c(-10, 10)

par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
plot(as.Date(ID), Y[[1]], col = "blue", ylim = ylim, type = "l",
     xlab = "", ylab = "SP100", axes = TRUE)

plot(as.Date(ID), Y[[2]], col = "blue", ylim = ylim, type = "l",
     xlab = "Time", ylab = "SP600", axes = TRUE)

## Plot the empirical copula
nPoints <- 100
u1 <- seq(0, 1, length.out = nPoints)
u2 <- u1
u <- mesh.grid(u1)

empDist <- matrix(hatCpl(u = u[, 1], v = u[, 2],
                         x = Y[[1]], y = Y[[2]]),
                  nPoints)
xlim <- c(0, 1)
ylim <- c(0, 1)

frechet.l <- matrix(uCpl(u = u, CplNM = "frechet-lower"), nPoints)
frechet.u <- matrix(uCpl(u = u, CplNM = "frechet-upper"), nPoints)


par(mfrow = c(1, 3), mar = c(4, 2, .5, .5))
contour(u1,  u2,  frechet.l,  add  =  FALSE,
        col  = "black",  lwd  =  0.8,  lty  = "solid",
        xlim  = xlim,  ylim  =  ylim,
        xlab = "Fréchet-Hoeffding lower bound copula", ylab = "")

contour(u1,  u2,  empDist,  add  =  FALSE,
        col  = "black",  lwd  =  0.8,  lty  = "solid",
        xlim  = xlim,  ylim  =  ylim,
        xlab = "empirical copula", ylab = "")

contour(u1,  u2,  frechet.u,  add  =  FALSE,
        col  = "black",  lwd  =  0.8,  lty  = "solid",
        xlim  = xlim,  ylim  =  ylim,
        xlab = "Fréchet-oeffding upper bound copula", ylab = "")

               ## })
