#' @export
lambdaTabular <- function(CplNM, tol = 1e-4)
{          ## tau, df -> lambda
    CplNM0 <- strsplit(CplNM, split = "_")[[1]][1]

    if(tolower(CplNM0) == "mvt")
    {

        x1lim <- c(0, 1)
        x2lim <- c(2, 30)

        x1Grid <- seq(x1lim[1]+tol, x1lim[2]-tol, tol)
        x2Grid <- exp(seq(log(x2lim[1])+tol, log(x2lim[2])-tol, tol))

        nGrid1 <- length(x1Grid)
        nGrid2 <- length(x2Grid)

        Mat <- matrix(NA, nGrid1, nGrid2)

        ## Big table takes huge amount of memory. We split the calculation if we
        ## require a very precise table.
        ## Split the calculations
        MaxLenCurr <- round(min(nGrid1*nGrid2, 1e6)/nGrid1)
        LoopIdx <- c(seq(1, nGrid2, MaxLenCurr), nGrid2)
        LoopIdx[1] <- 0

        tauIdxCurr0 <- 0

        nLoops <- length(LoopIdx)-1


        fun <- function(x1, x2)
        {
            tau <- x1
            df <- x2

            c0 <- sqrt(1-sin(pi*tau/2))/sqrt(1+sin(pi*tau/2))
            out <- (1/(1+c0^2))^((1+df)/2)/sqrt(df)/beta(df/2, 1/2)
            return(out)
        }

        for(j in 1:nLoops)
        {
            IdxCurr0 <- LoopIdx[j]+1
            IdxCurr1 <- LoopIdx[j+1]

            x1 <- rep(x1Grid, times = IdxCurr1-IdxCurr0+1)
            x2 <- rep(x2Grid[IdxCurr0:IdxCurr1], each = nGrid1)

            ## delta <- -log(2)/log(x1)
            ## theta <- log(2)/log(2-x2)
            ## parCpl <- list(theta = theta, delta = delta)

            Mat[, IdxCurr0:IdxCurr1] <- fun(x1, x2)
        }

        out <- list(Mat = Mat,
                    nGrid1 = nGrid1,
                    nGrid2 = nGrid2,
                    tol = tol,
                    x2Grid = x2Grid
                    )
    }
    return(out)
}
