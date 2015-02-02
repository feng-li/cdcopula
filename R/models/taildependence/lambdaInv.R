lambdaInv <- function(
    CplNM,
    parRepCpl,
    method = c("tabular", "iterative")[1])
    {
        if(tolower(method) == "tabular")
            {
                ## If the tauTabular not exist, create it
                if(!exists("tauTabular", envir = .GlobalEnv))
                    {
                        lambdaTabular <<- lambdaTabular(CplNM = CplNM, tol = 1e-3)
                    }
                out <- lambdaInv.tab(CplNM = CplNM,
                                     parRepCpl = parRepCpl,
                                     lambdaTabular = lambdaTabular)
            }
        else if(tolower(method) == "iterative")
            {
                out <- lambdaInv.iter(CplNM = CplNM,
                                      parRepCpl = parRepCpl,
                                      parCaller = "theta")
            }
        return(out)
    }

lambdaInv.tab <- function(CplNM, parRepCpl, lambdaTabular)
    {
        if(tolower(CplNM) == "mvt")
            {
                ## out <- vector("list", length(parRepCpl))
                lambda <- parRepCpl[["lambda"]]
                tau <- parRepCpl[["tau"]]

                if(length(lambdaL) !=length(tau))
                    {
                        stop("The input parameters should be of the same length.")
                    }

                ## The dictionary look up method for the upper tail dependence given
                ## lower tail dependence and Kendall's tau.
                tol <- lambdaTabular$tol
                nGrid1 <- lambdaTabular$nGrid1
                nGrid2 <- lambdaTabular$nGrid2

                Mat <- lambdaTabular$Mat

                x2Grid <- lambdaTabular$x2Grid

                ## The lower tail dependence indices.
                lambdaIdxRaw <- lambda/tol
                lambdaIdxFloor <- round(lambdaIdxRaw)

                ## Extra work to avoid under and over flow
                lambdaIdxFloor1 <- (lambdaIdxFloor < 1)
                lambdaIdxFloor2 <- (lambdaIdxFloor > nGridL)
                if(any(lambdaIdxFloor1))
                    {
                        lambdaIdxFloor[lambdaIdxFloor1] <- 1
                    }
                if(any(lambdaIdxFloor2))
                    {
                        lambdaIdxFloor[lambdaIdxFloor2] <- nGridL
                    }
                lambdaMatTabFloor <- Mat[lambdaIdxFloor, ,drop = FALSE]

                ## Find the indices of the closed values close to tau's left and right side
                nObs <- length(lambda)
                x2Test <- matrix(lambda, nObs, nGrid2)

                tauFloorDev0 <- - abs(x2Test-lambdaMatTabFloor)


                ## The indices of lambdaU
                ## This is the bottom neck of speed.
                x2FloorIdx0 <- max.col(tauFloorDev0)
                x2Floor0 <- x2Grid[x2FloorIdx0]

                ## Make sure the output format is same as the input
                out <- lambda
                out[1:length(out)] <- x2Floor0
            }
        return(out)
    }
