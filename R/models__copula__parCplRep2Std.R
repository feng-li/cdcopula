#' Convert the parameters to standard notation
#'
#' See Li and Kang 2018 paper.
#' @param CplNM NA
#' @param parCplRep NA
#' @return NA
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @export
parCplRep2Std <- function(CplNM, parCplRep)
{
    CplNM0 <- strsplit(CplNM, split = "_")[[1]][1]

    ## Initialize output Storage and name it
    out <- list()
    parNM <- names(parCplRep)

    if(tolower(CplNM0) %in% c("bb7", "sjc"))
    {
        ## The reparameterized parameters input
        if("lambdaL" %in% parNM)
        {
            lambdaL <- parCplRep[["lambdaL"]]
            delta <- (-log(2)/log(lambdaL))
        }

        if("lambdaU" %in% parNM)
        {
            lambdaU <- parCplRep[["lambdaU"]]
            theta <- log(2)/log(2-lambdaU)
        }

        if("tau" %in% parNM)
        {
            tau <- parCplRep[["tau"]]

            ## The new approach with funinv2d()
            if("lambdaL" %in% parNM)
            {
                ## The old approach.
                ## lambdaU <- as.vector(kendalltauInv(CplNM0 = CplNM0, parCplRep = parCplRep))
                ## theta <- log(2)/log(2-lambdaU)
                FUN <- function(x1, x2)
                {
                    parCpl <- list(delta = x1, theta = x2)
                    kendalltau(CplNM = "bb7", parCpl = parCpl)
                }
                theta <- funinv2d(FUN = FUN,
                                  x1 = delta,
                                  y = tau,
                                  x1lim = c(0, 30),
                                  x2lim = c(1, 30),
                                  method = "tabular",
                                  tol = 1e-3) # n-by-lq
            }
            else if("lambdaU" %in% parNM)
            {
                FUN <- function(x1, x2)
                {
                    parCpl <- list(delta = x2, theta = x1)
                    kendalltau(CplNM = "bb7", parCpl = parCpl)
                }

                delta <- funinv2d(FUN = FUN,
                                  x1 = theta,
                                  y = tau,
                                  x1lim = c(1, 30),
                                  x2lim = c(0, 30),
                                  method = "tabular",
                                  tol = 1e-3) # n-by-lq
            }

        }

        ## The first parameter
        out[["delta"]] <- delta

        ## The second parameter
        out[["theta"]] <- theta
    }
    else if(tolower(CplNM0) == "mvt")
    {
        ## Convert lower tail dependence and Kendall's tau to degrees of freedom and
        ## scale-matrix.
        if("tau" %in% parNM)
        {
            ## Kendall's Tau
            tau <- parCplRep[["tau"]]

            rho <- sin(tau*pi/2) # n-by-lq correlation
            ## TODO: Think about a better way to make correlation matrix always positive
            ## definite. Ping Ma: Simga  =  Sigma_0 + l*I
        }

        if("lambdaL" %in% parNM)
        {
            ## The lower tail dependence
            lambdaL <- parCplRep[["lambdaL"]]

            FUN <- function(x1, x2)
            {
                ## x1 = rho; x2  =  df
                parCpl <- list(rho = x1, df = x2)
                lambda(CplNM = "mvt", parCpl = parCpl)[["lambdaL"]]
            }

            df0 <- funinv2d(FUN = FUN,
                            method = "tabular",
                            x1 = rho,
                            y = lambdaL,
                            x1lim = c(0, 1),
                            x2lim = c(1, 30),
                            tol = 1e-3) # n-by-lq

            ## In multivariate case, lambda_ij = f(rho_ij, df) where df is independent of i
            ## and j. To make the calculate stable. let df_ij = df and make the mean of
            ## them. See Demarta & McNeil (2005) and Hult & Lindskog (2002)

            df <- as.matrix(rowMeans(df0)) # n-by-1
        }


        ## No re-parameterization
        if("rho" %in% parNM)
        {
            rho  =  parCplRep[["rho"]]
        }
        if("df" %in% parNM)
        {
            df  =  parCplRep[["df"]]
        }

        out[["rho"]] <- rho
        out[["df"]] <- df
    }
    else if(tolower(CplNM0) == "gumbel")
    {
        if("lambdaU" %in% parNM)
        {
            lambdaU <- parCplRep[["lambdaU"]]

            delta <- log(2)/log(2-lambdaU)
            out[["delta"]] <- delta
        }
        else if("tau" %in% parNM)
        {
            tau <- parCplRep[["tau"]]

            delta <- 1/(1-tau)
            out[["delta"]] <- delta
        }

    }
    else if(tolower(CplNM0) == "clayton")
    {

        if("tau" %in% parNM)
        {
            tau <- parCplRep[["tau"]]
            delta <- 2*tau/(1-tau)
            out[["delta"]] <- delta
        }
        else if("lambdaL" %in% parNM)
        {
            lambdaL <- parCplRep[["lambdaL"]]

            delta <- (-log(2)/log(lambdaL))

            out[["delta"]] <- delta
        }



    }
    else
    {
        stop("This copula ", CplNM,  " is not implemented!")
    }
    return(out)
}
