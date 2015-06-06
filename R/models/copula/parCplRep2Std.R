##' <description>
##'
##' <details>
##' @title <short tile>
##' @param copula
##' @param parCplRep
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Sun Feb 08 17:22:15 CST 2015; Current: Sun Feb 08 17:22:21 CST 2015.
parCplRep2Std <- function(CplNM, parCplRep)
  {
    ## Initialize output Storage and name it
    out <- list()

    if(tolower(CplNM) == "bb7")
      {
        ## The reparameterized parameters input
        lambdaL <- parCplRep[["lambdaL"]]
        tau <- parCplRep[["tau"]]

        lambdaU <- as.vector(kendalltauInv(
          CplNM = CplNM, parCplRep = parCplRep))

        ## The standard parametrization
        delta <- -log(2)/log(lambdaL)
        theta <- log(2)/log(2-lambdaU)

        ## The first parameter
        out[["delta"]] <- delta

        ## The second parameter
        out[["theta"]] <- theta
      }
    else if(tolower(CplNM) == "mvt")
      {
        ## Convert lower tail dependence and Kendall's tau to degrees of freedom and
        ## scale-matrix.

        ## The lower tail dependence
        lambdaL <- parCplRep[["lambdaL"]]

        ## Kendall's Tau
        tau <- parCplRep[["tau"]]

        rho <- sin(tau*pi/2) # n-by-lq correlation

        FUN <- function(x1, x2)
          {
            ## x1 = rho; x2  =  df
            parCpl <- list(rho = x1, df = x2)
            lambda(CplNM = "mvt", parCpl = parCpl)[["lambdaL"]]
          }

        ## browser()

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

        out[["rho"]] <- rho
        out[["df"]] <- df
      }
    else
      {
        stop("This copula is not implemented!")
      }
    return(out)
  }
