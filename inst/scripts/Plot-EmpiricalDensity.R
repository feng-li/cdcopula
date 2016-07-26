## This script plots the daily return of SP100 and SP600 data
sourceDir("~/code/cdcopula/R", recursive = TRUE)

load("data/SP100-SP400-SP600-20150206.Rdata")

## Plot The Returns
nObs <- length(ID)


## The empirical copula
nPoints <- 150
u1 <- seq(0, 1, length.out = nPoints)
u2 <- u1
u <- mesh.grid(u1)




nDataWindow <- 6

xlim <- c(0, 1)
ylim <- c(0, 1)


par(mfrow = c(2, nDataWindow/2), mar = c(5, 2, .5, .5))

dataWindowIdx <- data.partition(nObs = nObs,
                                args = list(partiMethod = "ordered",
                                            N.subsets = nDataWindow))
for(i in 1:nDataWindow)
{
    ## The empirical density
    require("MASS")
    idx <- dataWindowIdx[[i]]
    bivn.kde <- kde2d(Y[[1]][idx, ],  Y[[2]][idx, ],  n  =  nPoints)

    contour(bivn.kde, xlim = c(-2.2, 2.2), ylim = c(-2.2, 2.2),
            xlab = paste("from ", ID[idx[1]]," to ",  ID[idx[length(idx)]], sep = ""),
            col = ifelse(i == 5, "red", "blue"))


    ## op <- par(fig  =  c(.06, .46, .53, .99), new  =  TRUE, cex = 0.6)
    ## empDist <- matrix(hatCpl(u = u[, 1], v = u[, 2],
    ##                      x = Y[[1]][idx, ], y = Y[[2]][idx, ]),
    ##                   nPoints)
    ## contour(u1,  u2,  empDist,  add  =  FALSE,
    ##         col  = "black",  lwd  =  0.4,  lty  = "solid",
    ##         xlim  = xlim,  ylim  =  ylim,
    ##         xlab = "", ylab = "")
    ## par(op)
}



## Plot the Lower and Uper bounds
frechet.l <- matrix(pCpl(u = u, CplNM = "frechet-lower"), nPoints)
frechet.u <- matrix(pCpl(u = u, CplNM = "frechet-upper"), nPoints)

contour(u1,  u2,  frechet.l,  add  =  FALSE,
        col  = "black",  lwd  =  0.8,  lty  = "solid",
        xlim  = xlim,  ylim  =  ylim,
        xlab = "Fréchet-Hoeffding lower bound copula", ylab = "")


contour(u1,  u2,  frechet.u,  add  =  FALSE,
        col  = "black",  lwd  =  0.8,  lty  = "solid",
        xlim  = xlim,  ylim  =  ylim,
        xlab = "Fréchet-Hoeffding upper bound copula", ylab = "")

               ## })
